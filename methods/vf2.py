#!/usr/bin/env python3
"""The reference: plain whole-graph igraph VF2 per pattern, the target's labels against the pattern's, occurrences
deduplicated by node set (the ground-truth definition).  Shares no matching code with seeded_exact.py or
hop_filter.py (only common.py's input reading).  One pattern per worker (with --recorded, an earlier run's
times order them slowest first); a pattern still running after --limit seconds is recorded as capped, with no count.

  vf2.py target.gml patterns.txt --workers 16 --limit 3600 --output vf2.jsonl   (lines: pattern, occurrences,
  seconds, capped, sets as sorted target node ids)"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import signal
import time
from pathlib import Path

import igraph as ig

import common

_T: dict = {}


class Capped(Exception):
    pass


def _alarm(_signum, _frame):
    raise Capped()


def _one(j: int):
    P = _T['patterns'][j]
    signal.signal(signal.SIGALRM, _alarm)
    start = time.perf_counter()
    signal.alarm(_T['limit'])
    try:
        maps = _T['G'].get_subisomorphisms_vf2(P, color1=_T['G'].vs['label'], color2=P.vs['label'])
        signal.alarm(0)
    except Capped:
        return j, time.perf_counter() - start, None
    ids = _T['ids']
    return j, time.perf_counter() - start, sorted({tuple(sorted(ids[v] for v in m)) for m in maps})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('target', type=Path)
    ap.add_argument('patterns', nargs='+', help='pattern GMLs, directories of them, or text files listing them')
    ap.add_argument('--workers', type=int, default=len(os.sched_getaffinity(0)))
    ap.add_argument('--limit', type=int, default=3600, help='seconds per pattern')
    ap.add_argument('--recorded', type=Path, default=None, help='an earlier vf2.py output, to start the slowest patterns first')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    G, paths, patterns = common.load(args.target, args.patterns)
    _T.update(G=G, ids=[int(v) for v in G.vs['id']], patterns=patterns, limit=args.limit)
    order = list(range(len(paths)))
    if args.recorded and args.recorded.exists():
        took = {r['pattern']: r['seconds'] for r in map(json.loads, args.recorded.open())}
        order.sort(key=lambda j: -took.get(paths[j].name, 0.))
    rows = {}
    with mp.get_context('fork').Pool(args.workers) as pool:
        start = time.perf_counter()
        for j, sec, sets in pool.imap_unordered(_one, order, chunksize=1):
            rows[j] = {'pattern': paths[j].name, 'occurrences': None if sets is None else len(sets), 'seconds': sec,
                       'capped': sets is None, 'sets': sets or []}
        wall = time.perf_counter() - start
    with args.output.open('w', encoding='utf-8') as fh:
        for j in range(len(paths)):
            fh.write(json.dumps(rows[j]) + '\n')
    print(json.dumps({'method': 'vf2', 'target': str(args.target), 'patterns': len(paths), 'workers': args.workers,
                      'limit_s': args.limit, 'wall_s': wall, 'capped': sum(r['capped'] for r in rows.values()),
                      'occurrences': sum(r['occurrences'] or 0 for r in rows.values())}), flush=True)


if __name__ == '__main__':
    main()
