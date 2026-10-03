#!/usr/bin/env python3
"""Markdown wall-time tables (one per worker count) from run_all.sh result directories, fastest per target in bold,
'*' where some queries hit the per-query limit (a lower bound).  hop_filter and seeded exact are their search time,
the tools and VF2 their wall.

  tables.py results/small results/large ... [--names "Small" "Large" ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROWS = [('vf3', 'VF3'), ('vf3l', 'VF3L'), ('fastiso', 'FASTiso'), ('daf', 'DAF'), ('ri', 'RI 3.6'),
        ('hop', 'hop_filter (ours)'), ('seeded', 'seeded_exact (ours)'),
        ('vf2', 'plain VF2')]
FIELD = {'hop': 'search_s', 'seeded': 'search_s', 'vf2': 'wall_s'}


def value(d: Path, w: int, key: str):
    if key in FIELD:
        p = d / f'w{w}_{key}.json'
        if not (p.exists() and p.stat().st_size):
            return None
        x = json.loads(p.read_text())
        return x[FIELD[key]], bool(x.get('capped'))
    p = d / f'w{w}_{key}.wall'
    if not (p.exists() and p.stat().st_size):
        return None
    rows = [l.split('\t') for l in (d / f'w{w}_{key}.tsv').read_text().splitlines() if l.strip()]
    return float(p.read_text()), any(int(r[-1]) for r in rows)


def fmt(v: float) -> str:
    return f'{v:,.0f}' if v >= 100 else f'{v:.1f}' if v >= 1 else f'{v:.2f}'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+', type=Path)
    ap.add_argument('--names', nargs='*', default=None)
    args = ap.parse_args()
    names = args.names or [d.name for d in args.dirs]
    workers = sorted({int(p.name[1:].split('_')[0]) for d in args.dirs for p in d.glob('w*_host.txt')})
    for w in workers:
        vals = {k: [value(d, w, k) for d in args.dirs] for k, _n in ROWS}
        best = [min((v[i][0] for v in vals.values() if v[i] is not None), default=None) for i in range(len(args.dirs))]
        print(f'\n**{w} worker{"s" if w > 1 else ""}**: seconds for every query of the target\n')
        print('| Method | ' + ' | '.join(names) + ' |')
        print('| --- |' + ' --- |' * len(names))
        for k, label in ROWS:
            if all(v is None for v in vals[k]):
                continue
            nodes = []
            for i, v in enumerate(vals[k]):
                if v is None:
                    nodes.append('–')
                    continue
                s = fmt(v[0]) + ('\\*' if v[1] else '')
                nodes.append(f'**{s}**' if v[0] == best[i] else s)
            print(f'| {label} | ' + ' | '.join(nodes) + ' |')
    hosts = sorted({(d / f'w{w}_host.txt').read_text().split(',')[1].strip() for d in args.dirs for w in workers
                    if (d / f'w{w}_host.txt').exists()})
    print(f'\nCPU: {", ".join(hosts)}; each target\'s methods on one machine.')


if __name__ == '__main__':
    main()
