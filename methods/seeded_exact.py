#!/usr/bin/env python3
"""Seeded exact matching: every occurrence of each pattern in a labelled target graph.

An occurrence is a set of target nodes onto which the pattern maps one-to-one, each node onto a node of its label and
each pattern edge onto a target edge (a monomorphism: extra target edges are allowed), counted once per node set.
Patterns must be connected.

Every occurrence holds a node of the pattern's rarest label (the anchor) and lies within the anchor's eccentricity of
it, since an embedding never lengthens a path.  So VF2 runs in that ball around each such node, the anchor pinned by a
label nothing else has.  One pattern per worker.  Exact: checked node set by node set against whole-graph VF2 (vf2.py).

  seeded_exact.py target.gml patterns.txt --workers 16 --output hits.jsonl [--sets]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import igraph as ig
import numpy as np

import common

PIN = 1 << 20          # a label no node uses: pins the anchor
_W: dict = {}          # set before the workers fork, so they share it: target and patterns


def seeded(G: ig.Graph, labels: np.ndarray, P: ig.Graph) -> set:
    """Every occurrence of P in G, as sorted tuples of G's vertex indices."""
    pl = P.vs['label']
    counts = np.bincount(labels, minlength=max(pl) + 1)
    anchor = min(range(P.vcount()), key=lambda u: (counts[pl[u]], u))
    radius = int(P.eccentricity(anchor))
    colour2 = list(pl)
    colour2[anchor] = PIN
    found = set()
    for node in np.flatnonzero(labels == pl[anchor]).tolist():
        ball = sorted(G.neighborhood(node, order=radius))      # induced_subgraph numbers vertices in id order
        if len(ball) < P.vcount():
            continue
        colour1 = labels[ball].tolist()
        colour1[ball.index(node)] = PIN
        sub = G.induced_subgraph(ball, implementation='create_from_scratch')
        for m in sub.get_subisomorphisms_vf2(P, color1=colour1, color2=colour2):
            found.add(tuple(sorted(ball[i] for i in m)))
    return found


def _exact_worker(j: int):
    start = time.perf_counter()
    found = seeded(_W['G'], _W['labels'], _W['patterns'][j])
    return j, time.perf_counter() - start, list(found) if _W['sets'] else len(found)


def run_exact(G: ig.Graph, patterns: list, workers: int, sets: bool) -> tuple[list, dict]:
    t0 = time.perf_counter()
    bare = ig.Graph(n=G.vcount(), edges=G.get_edgelist())       # induced_subgraph would copy every GML attribute
    _W.update(G=bare, labels=np.asarray(G.vs['label'], np.int64), patterns=patterns, sets=sets)
    prepare = time.perf_counter() - t0
    with mp.get_context('fork').Pool(workers) as pool:
        t1 = time.perf_counter()
        rows = sorted(pool.imap_unordered(_exact_worker, range(len(patterns)), chunksize=1))
        search = time.perf_counter() - t1
    return ([{'occurrences': f if isinstance(f, int) else len(f), 'seconds': s, **({'sets': f} if sets else {})}
             for _j, s, f in rows], {'prepare_s': prepare, 'search_s': search})



def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('target', type=Path, help='target GML: integer `label` per node')
    ap.add_argument('patterns', nargs='+', help='pattern GMLs, directories of them, or text files listing them')
    ap.add_argument('--workers', type=int, default=len(os.sched_getaffinity(0)))
    ap.add_argument('--sets', action='store_true', help='also write each occurrence as its sorted target node ids')
    ap.add_argument('--output', type=Path, default=None,
                    help='JSON lines, one per pattern: occurrences, seconds (its search), and sets with --sets')
    args = ap.parse_args()

    t0 = time.perf_counter()
    G, paths, patterns = common.load(args.target, args.patterns)
    info = {'method': 'seeded_exact', 'target': str(args.target), 'nodes': G.vcount(), 'edges': G.ecount(),
            'patterns': len(patterns), 'workers': args.workers, 'read_s': time.perf_counter() - t0}
    rows, timing = run_exact(G, patterns, args.workers, args.sets)
    info.update(timing, total_s=timing['prepare_s'] + timing['search_s'],
                occurrences=sum(r['occurrences'] for r in rows))
    if args.output:
        ids = np.asarray(G.vs['id'], np.int64)
        with args.output.open('w', encoding='utf-8') as fh:
            for p, row in zip(paths, rows):
                if args.sets:
                    row['sets'] = sorted(sorted(ids[list(s)].tolist()) for s in row['sets'])
                fh.write(json.dumps({'pattern': p.name, **row}) + '\n')
    print(json.dumps(info), flush=True)


if __name__ == '__main__':
    main()
