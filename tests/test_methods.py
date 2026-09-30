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


def test_published_tools():
    if not (BIN / 'ri_batch').exists():
        print(f'  published tools not built in {BIN}: skipped')
        return
    d = work()
    subprocess.run([PY, ROOT / 'benchmark' / 'convert.py', d / 'target.gml', d / 'patterns.txt', d / 'input'], check=True)
    names = (d / 'input' / 'patterns.txt').read_text().split()
    hop = {r['pattern']: r for r in map(json.loads, (d / 'hop.jsonl').open())}
    ref = sets(d / 'vf2.jsonl')
    for exe, fmt in (('ri_batch', 'gfu'), ('vf3_batch_VF3', 'vft'), ('vf3_batch_VF3L', 'vft'), ('daf_batch', 'daf'),
                     ('fastiso_batch', 'gfu')):
        out = d / f'{exe}.tsv'
        subprocess.run([BIN / exe, d / 'input' / f'target.{fmt}', d / 'input' / f'queries_{fmt}.txt', '2', out, '60'],
                       check=True, stdout=subprocess.DEVNULL)
        rows = {int(r[0]): r[1:] for r in (l.split('\t') for l in out.read_text().splitlines() if l.strip())}
        assert len(rows) == len(names), exe
        for j, p in enumerate(names):
            assert int(rows[j][1]) == hop[p]['mappings'], (exe, p)
            if exe == 'ri_batch':
                assert int(rows[j][2]) == len(ref[p]), (exe, p)


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
