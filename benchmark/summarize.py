#!/usr/bin/env python3
"""Summarise run_all.sh outputs: exactness and wall time, per result directory and worker count.

Exactness, per pattern:
  reference   plain VF2's node sets where it ran and finished, else hop_filter's (itself checked against VF2)
  ours        node-set counts of seeded exact and hop_filter against the reference; exact_match filter's shortfall
              (occurrences that stick out of their owner window)
  tools       RI's distinct node sets against the reference; every tool's mapping count against hop_filter's (hop counts
              every mapping onto each node set it owns)
Time: wall for every query of the target.  The tools read and prepare the target before their clock; so are
hop_filter and seeded exact shown as their search (hop's setup in brackets), exact_match filter as setup + search, VF2
as its wall.  '*': some queries hit the per-query limit and count at the limit, so the wall is a lower bound.

  summarize.py <out_dir> [<out_dir> ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

TOOLS = [('ri', 'RI'), ('vf3', 'VF3'), ('vf3l', 'VF3L'), ('daf', 'DAF'), ('fastiso', 'FASTiso')]


def jl(p: Path):
    return {r['pattern']: r for r in map(json.loads, p.open())} if p.exists() and p.stat().st_size else None


def js(p: Path):
    return json.loads(p.read_text()) if p.exists() and p.stat().st_size else None


def tsv(p: Path):
    if not p.exists():
        return None
    return {int(r[0]): r[1:] for r in (l.split('\t') for l in p.read_text().splitlines() if l.strip())}


def report(R: Path, W: int) -> list[str]:
    names = (R / 'input' / 'patterns.txt').read_text().split()
    n = len(names)
    hop, vf2, seeded, filt = (jl(R / f'w{W}_{m}.jsonl') for m in ('hop', 'vf2', 'seeded', 'filter'))
    if hop is None:
        return [f'  {W} workers: no hop_filter run']
    ref = {p: hop[p]['occurrences'] for p in names}
    src = 'hop_filter'
    if vf2:
        src = 'VF2'
        for p in names:
            if not vf2[p]['capped']:
                ref[p] = vf2[p]['occurrences']
    same = lambda rows: sum(rows[p]['occurrences'] == ref[p] for p in names)            # noqa: E731
    lines = [f'  {W} workers, reference {src}: {sum(ref.values()):,} occurrences']
    checks = [f'hop {same(hop)}/{n}']
    if vf2:
        checks.insert(0, f'VF2 finished {sum(not r["capped"] for r in vf2.values())}/{n}')
    if seeded:
        checks.append(f'seeded {same(seeded)}/{n}')
    if filt:
        checks.append(f'filter {same(filt)}/{n} (short {sum(ref[p] - filt[p]["occurrences"] for p in names):,})')
    hmaps = [hop[p]['mappings'] for p in names]
    for key, name in TOOLS:
        t = tsv(R / f'w{W}_{key}.tsv')
        if t is None or len(t) < n:
            checks.append(f'{name} missing' if t is None else f'{name} returned {len(t)}/{n}')
            continue
        c = f'{name} mappings {sum(int(t[j][1]) == hmaps[j] for j in range(n))}/{n}'
        if key == 'ri':
            c += f', node sets {sum(int(t[j][2]) == ref[names[j]] for j in range(n))}/{n}'
        to = sum(int(t[j][-1]) for j in range(n))
        checks.append(c + (f', {to} timed out' if to else ''))
    lines.append('    exact: ' + '; '.join(checks))
    times = []
    for key, name in TOOLS:
        w = R / f'w{W}_{key}.wall'
        if w.exists() and w.stat().st_size:
            t = tsv(R / f'w{W}_{key}.tsv')
            star = '*' if sum(int(r[-1]) for r in t.values()) else ''
            times.append(f'{name} {float(w.read_text()):.2f}{star}')
    h = js(R / f'w{W}_hop.json')
    times.append(f"hop {h['search_s']:.2f} (+{h['prepare_s']:.2f})")
    for key, name, field in (('filter', 'filter', 'total_s'), ('seeded', 'seeded exact', 'search_s')):
        d = js(R / f'w{W}_{key}.json')
        if d:
            times.append(f'{name} {d[field]:.2f}')
    v = js(R / f'w{W}_vf2.json')
    if v:
        times.append(f"VF2 {v['wall_s']:.2f}{'*' if v['capped'] else ''}")
    lines.append('    wall s: ' + ', '.join(times))
    return lines


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for d in map(Path, sys.argv[1:]):
        ws = sorted({int(p.name[1:].split('_')[0]) for p in d.glob('w*_host.txt')}, reverse=True)
        hosts = {w: (d / f'w{w}_host.txt').read_text().strip() for w in ws}
        print(f'== {d.name}')
        for w in ws:
            print(f'  host: {hosts[w]}')
            print("\n".join(report(d, w)))


if __name__ == '__main__':
    main()
