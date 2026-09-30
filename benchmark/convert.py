#!/usr/bin/env python3
"""A target GML and its patterns -> the published tools' input files.

  target.gfu, q/NNNN.gfu    RI's gfu (RI, FASTiso): vertex count, one label per vertex, edge count, one edge per line
  target.vft, q/NNNN.vft    vf3lib's text format (VF3, VF3L, read undirected): vertex count, "id label" per vertex, then
                            per vertex its edge count and each undirected edge once, from its lower end
  target.daf, q/NNNN.daf    DAF's format: "t n m", "v id label" per vertex, "e a b" per edge
  queries_{gfu,vft,daf}.txt the query files in pattern order;  patterns.txt  the pattern file names in that order

Inputs are read as the methods read them (methods/common.py load): patterns given as GML files, folders or lists, and
labels (integers of any size, or strings) recoded as 0..L-1 over the target and the patterns together, so every tool
sees the same labels as our methods.

  convert.py target.gml patterns out_dir        (patterns: a list file, a folder of .gml files, or .gml files)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'methods'))
import common  # noqa: E402


def write(stem: Path, labels: list[int], edges: list[tuple[int, int]]) -> None:
    n = len(labels)
    stem.with_suffix('.gfu').write_text(f'#g\n{n}\n' + ''.join(f'{l}\n' for l in labels) + f'{len(edges)}\n'
                                        + ''.join(f'{a} {b}\n' for a, b in edges))
    below = {}
    for a, b in edges:
        below.setdefault(min(a, b), []).append(max(a, b))
    with stem.with_suffix('.vft').open('w') as f:
        f.write(f'{n}\n' + ''.join(f'{i} {l}\n' for i, l in enumerate(labels)))
        for i in range(n):
            nb = below.get(i, [])
            f.write(f'{len(nb)}\n' + ''.join(f'{i} {j}\n' for j in nb))
    stem.with_suffix('.daf').write_text(f't {n} {len(edges)}\n' + ''.join(f'v {i} {l}\n' for i, l in enumerate(labels))
                                        + ''.join(f'e {a} {b}\n' for a, b in edges))


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    target, out = Path(sys.argv[1]), Path(sys.argv[-1]).resolve()
    G, paths, patterns = common.load(target, sys.argv[2:-1])
    (out / 'q').mkdir(parents=True, exist_ok=True)
    write(out / 'target', G.vs['label'], G.get_edgelist())
    for j, P in enumerate(patterns):
        write(out / 'q' / f'{j:04d}', P.vs['label'], P.get_edgelist())
    for ext in ('gfu', 'vft', 'daf'):
        (out / f'queries_{ext}.txt').write_text(''.join(f"{out / 'q' / f'{j:04d}.{ext}'}\n" for j in range(len(paths))))
    (out / 'patterns.txt').write_text(''.join(p.name + '\n' for p in paths))
    print(f'{target.name}: {G.vcount():,} nodes, {G.ecount():,} edges, {len(paths)} patterns, '
          f'{len(set(G.vs["label"]))} labels -> {out}', flush=True)


if __name__ == '__main__':
    main()
