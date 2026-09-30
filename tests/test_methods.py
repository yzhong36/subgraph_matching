#!/usr/bin/env python3
"""Every method against plain whole-graph VF2 on the synthetic example (examples/make_example.py), node set by node set:
seeded exact and hop_filter must equal VF2.  If the published tools are built (build/bin, or BIN), each tool's mapping
count must equal hop_filter's and RI's distinct node sets must equal VF2's.

  python tests/test_methods.py            (or pytest tests/)
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
BIN = Path(os.environ.get('BIN', ROOT / 'build' / 'bin'))
_WORK: dict = {}


def work() -> Path:
    """The example and every method's output, made once."""
    if 'dir' not in _WORK:
        d = Path(tempfile.mkdtemp(prefix='sgm_test_'))
        subprocess.run([PY, ROOT / 'examples' / 'make_example.py', d], check=True)
        t, p = d / 'target.gml', d / 'patterns.txt'
        subprocess.run([PY, ROOT / 'methods' / 'vf2.py', t, p, '--workers', '2', '--output', d / 'vf2.jsonl'], check=True)
        for name, args in (('seeded', ['seeded_exact.py']), ('hop', ['hop_filter.py'])):
            subprocess.run([PY, ROOT / 'methods' / args[0], *args[1:], t, p, '--workers', '2', '--sets',
                            '--output', d / f'{name}.jsonl'], check=True)
        _WORK['dir'] = d
    return _WORK['dir']


def sets(path: Path) -> dict:
    return {r['pattern']: {tuple(s) for s in r['sets']} for r in map(json.loads, path.open())}


def test_seeded_exact_equals_vf2():
    d = work()
    assert sets(d / 'seeded.jsonl') == sets(d / 'vf2.jsonl')


def test_hop_filter_equals_vf2():
    d = work()
    assert sets(d / 'hop.jsonl') == sets(d / 'vf2.jsonl')


def relabelled(kind: str) -> Path:
    """The example with its labels rewritten (strings, or large sparse integers) and three patterns added: single nodes
    with a label the target has and with one it lacks, and a triangle using the lacking label; patterns as a folder."""
    import igraph as ig
    src, d = work(), work() / kind
    if not d.exists():
        (d / 'patterns').mkdir(parents=True)
        new = (lambda l: f'type_{int(l)}') if kind == 'strings' else (lambda l: 1_000_003 * int(l) + 7)
        for rel in ['target.gml'] + [f'patterns/{p.name}' for p in sorted((src / 'patterns').glob('*.gml'))]:
            g = ig.Graph.Read_GML(str(src / rel))
            g.vs['label'] = [new(l) for l in g.vs['label']]
            g.write_gml(str(d / rel))
        for name, labs, edges in (('single_present', [0], []), ('single_absent', [99], []),
                                  ('triangle_absent', [0, 99, 0], [(0, 1), (1, 2), (0, 2)])):
            g = ig.Graph(n=len(labs), edges=edges)
            g.vs['id'], g.vs['label'] = list(range(len(labs))), [new(l) for l in labs]
            g.write_gml(str(d / 'patterns' / f'{name}.gml'))
        for m in ('seeded_exact.py', 'hop_filter.py', 'vf2.py'):
            subprocess.run([PY, ROOT / 'methods' / m, d / 'target.gml', d / 'patterns', '--workers', '2',
                            '--output', d / f'{m[:-3]}.jsonl'], check=True, stdout=subprocess.DEVNULL)
    return d


def test_labels_of_any_kind():
    ref = {p: len(s) for p, s in sets(work() / 'vf2.jsonl').items()}
    ref['single_present.gml'] = sum(1 for l in __import__('igraph').Graph.Read_GML(str(work() / 'target.gml')).vs['label']
                                    if l == 0)
    ref['single_absent.gml'] = ref['triangle_absent.gml'] = 0
    for kind in ('strings', 'sparse'):
        d = relabelled(kind)
        for m in ('seeded_exact', 'hop_filter', 'vf2'):
            got = {r['pattern']: r['occurrences'] for r in map(json.loads, (d / f'{m}.jsonl').open())}
            assert got == ref, (kind, m)


def test_published_tools():
    """On the example, and on its string-labelled variant with single-node patterns (patterns as a folder)."""
    if not (BIN / 'ri_batch').exists():
        print(f'  published tools not built in {BIN}: skipped')
        return
    for d, plist, hop_out in ((work(), work() / 'patterns.txt', 'hop.jsonl'),
                              (relabelled('strings'), relabelled('strings') / 'patterns', 'hop_filter.jsonl')):
        subprocess.run([PY, ROOT / 'benchmark' / 'convert.py', d / 'target.gml', plist, d / 'input'], check=True,
                       stdout=subprocess.DEVNULL)
        names = (d / 'input' / 'patterns.txt').read_text().split()
        hop = {r['pattern']: r for r in map(json.loads, (d / hop_out).open())}
        for exe, fmt in (('ri_batch', 'gfu'), ('vf3_batch_VF3', 'vft'), ('vf3_batch_VF3L', 'vft'), ('daf_batch', 'daf'),
                         ('fastiso_batch', 'gfu')):
            out = d / f'{exe}.tsv'
            subprocess.run([BIN / exe, d / 'input' / f'target.{fmt}', d / 'input' / f'queries_{fmt}.txt', '2', out, '60'],
                           check=True, stdout=subprocess.DEVNULL)
            rows = {int(r[0]): r[1:] for r in (l.split('\t') for l in out.read_text().splitlines() if l.strip())}
            assert len(rows) == len(names), exe
            for j, p in enumerate(names):
                assert int(rows[j][1]) == hop[p]['mappings'], (d.name, exe, p)
                if exe == 'ri_batch':
                    assert int(rows[j][2]) == hop[p]['occurrences'], (d.name, exe, p)

if __name__ == '__main__':
    failed = 0
    for name, fn in [(k, v) for k, v in globals().items() if k.startswith('test_')]:
        try:
            fn()
            print(f'ok      {name}')
        except AssertionError as e:
            failed += 1
            print(f'FAILED  {name} {e}')
    d = work()
    total = sum(len(s) for s in sets(d / 'vf2.jsonl').values())
    print(f'{total} occurrences of {len(sets(d / "vf2.jsonl"))} patterns in the example; work dir {d}')
    sys.exit(1 if failed else 0)
