#!/usr/bin/env python3
"""A target GML and its patterns -> the published tools' input files.

  target.gfu, q/NNNN.gfu    RI's gfu (RI, FASTiso): vertex count, one label per vertex, edge count, one edge per line
  target.vft, q/NNNN.vft    vf3lib's text format (VF3, VF3L, read undirected): vertex count, "id label" per vertex, then
                            per vertex its edge count and each undirected edge once, from its lower end
  target.daf, q/NNNN.daf    DAF's format: "t n m", "v id label" per vertex, "e a b" per edge
  queries_{gfu,vft,daf}.txt the query files in pattern order;  patterns.txt  the pattern file names in that order

The .vft and .daf targets carry one isolated placeholder vertex for every label below the largest label (of the target
and the patterns) that the target lacks (as the benchmarked inputs did, for label-indexed tools).  Placeholders
have no edges and every query node has degree >= 1, so no mapping can use one; counts are unchanged.

  convert.py target.gml patterns.txt out_dir
"""
from __future__ import annotations

import sys
from pathlib import Path

import igraph as ig


def read(path) -> tuple[list[int], list[tuple[int, int]]]:
    g = ig.Graph.Read_GML(str(path))
    return [int(x) for x in g.vs['label']], g.get_edgelist()


def write(stem: Path, labels: list[int], edges: list[tuple[int, int]], placeholders: list[int] = ()) -> None:
    n = len(labels)
    stem.with_suffix('.gfu').write_text(f'#g\n{n}\n' + ''.join(f'{l}\n' for l in labels) + f'{len(edges)}\n'
                                        + ''.join(f'{a} {b}\n' for a, b in edges))
    lab = list(labels) + list(placeholders)
    below = {}
    for a, b in edges:
        below.setdefault(min(a, b), []).append(max(a, b))
    with stem.with_suffix('.vft').open('w') as f:
        f.write(f'{len(lab)}\n' + ''.join(f'{i} {l}\n' for i, l in enumerate(lab)))
        for i in range(len(lab)):
            nb = below.get(i, [])
            f.write(f'{len(nb)}\n' + ''.join(f'{i} {j}\n' for j in nb))
    stem.with_suffix('.daf').write_text(f't {len(lab)} {len(edges)}\n' + ''.join(f'v {i} {l}\n' for i, l in enumerate(lab))
                                        + ''.join(f'e {a} {b}\n' for a, b in edges))


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    target, plist, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]).resolve()
    (out / 'q').mkdir(parents=True, exist_ok=True)
    paths = [plist.parent / line.strip() for line in plist.read_text().splitlines() if line.strip()]
    tl, te = read(target)
    queries = [read(p) for p in paths]
    top = max(set(tl) | {l for ql, _qe in queries for l in ql})
    write(out / 'target', tl, te, sorted(set(range(top + 1)) - set(tl)))
    for j, (ql, qe) in enumerate(queries):
        write(out / 'q' / f'{j:04d}', ql, qe)
    for ext in ('gfu', 'vft', 'daf'):
        (out / f'queries_{ext}.txt').write_text(''.join(f"{out / 'q' / f'{j:04d}.{ext}'}\n" for j in range(len(paths))))
    (out / 'patterns.txt').write_text(''.join(p.name + '\n' for p in paths))
    print(f'{target.name}: {len(tl):,} nodes, {len(te):,} edges, {len(paths)} patterns, '
          f'{top + 1 - len(set(tl))} placeholder labels -> {out}', flush=True)


if __name__ == '__main__':
    main()
