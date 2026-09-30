#!/usr/bin/env python3
"""A toy target graph and patterns, for trying the methods and for the tests.

Target: nodes at the points of a jittered triangular lattice, each joined to its lattice neighbours (a planar
triangulation, like the Delaunay graph of a point set), each with one of a few integer labels, some common and some
rare, and x, y coordinates.  Patterns: connected pieces grown from random target nodes (so each occurs at least once),
3-8 nodes with every target edge among them, and every eighth with its labels shuffled (often absent).

  make_example.py out_dir [--side 40] [--labels 5] [--patterns 40] [--seed 0]
  -> out_dir/target.gml (id, label, x, y), out_dir/patterns/*.gml (id, label), out_dir/patterns.txt
"""
from __future__ import annotations

import argparse
from pathlib import Path

import igraph as ig
import numpy as np


def lattice(side: int, rng: np.random.Generator) -> tuple[np.ndarray, list[tuple[int, int]]]:
    idx = lambda r, c: r * side + c                                                   # noqa: E731
    xy = np.array([[c + .5 * (r % 2), r * np.sqrt(3) / 2] for r in range(side) for c in range(side)])
    xy += rng.normal(0, .08, xy.shape)
    edges = []
    for r in range(side):
        for c in range(side):
            if c + 1 < side:
                edges.append((idx(r, c), idx(r, c + 1)))
            if r + 1 < side:                                   # the two nodes above, offset by the row's parity
                for dc in ((0, 1) if r % 2 else (-1, 0)):
                    if 0 <= c + dc < side:
                        edges.append((idx(r, c), idx(r + 1, c + dc)))
    return xy, edges


def grow(g: ig.Graph, size: int, rng: np.random.Generator) -> list[int]:
    nodes = [int(rng.integers(g.vcount()))]
    while len(nodes) < size:
        frontier = sorted({u for v in nodes for u in g.neighbors(v)} - set(nodes))
        nodes.append(int(rng.choice(frontier)))
    return sorted(nodes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out', type=Path)
    ap.add_argument('--side', type=int, default=40)
    ap.add_argument('--labels', type=int, default=5)
    ap.add_argument('--patterns', type=int, default=40)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    xy, edges = lattice(args.side, rng)
    p = np.array([2.0 ** -k for k in range(args.labels)])                                # common to rare labels
    labels = rng.choice(args.labels, size=len(xy), p=p / p.sum())
    g = ig.Graph(n=len(xy), edges=edges)
    g.vs['id'] = list(range(len(xy)))
    g.vs['label'] = labels.tolist()
    g.vs['x'], g.vs['y'] = xy[:, 0].tolist(), xy[:, 1].tolist()
    (args.out / 'patterns').mkdir(parents=True, exist_ok=True)
    g.write_gml(str(args.out / 'target.gml'))
    names = []
    for j in range(args.patterns):
        nodes = grow(g, int(rng.integers(3, 9)), rng)
        P = g.induced_subgraph(nodes)
        if j % 8 == 7:
            P.vs['label'] = rng.permutation(P.vs['label']).tolist() if P.vcount() > 3 else rng.integers(args.labels, size=3).tolist()
        P.vs['id'] = list(range(P.vcount()))
        name = f'pattern_{j:03d}_n{P.vcount()}.gml'
        P.write_gml(str(args.out / 'patterns' / name))
        names.append(f'patterns/{name}')
    (args.out / 'patterns.txt').write_text(''.join(n + '\n' for n in names))
    print(f'{args.out}: {g.vcount()} nodes, {g.ecount()} edges, {args.labels} labels, {len(names)} patterns')


if __name__ == '__main__':
    main()
