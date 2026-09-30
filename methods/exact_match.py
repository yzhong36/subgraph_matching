#!/usr/bin/env python3
"""Exact labelled subgraph matching on a spatial node graph, two ways.

An occurrence is a set of target nodes onto which the pattern maps one-to-one, each node onto a node of its
label and each pattern edge onto a target edge (a monomorphism: extra target edges are allowed), counted once
per node set.  Patterns must be connected.

  exact    seeded exact matching on the whole graph.  Every occurrence holds a node of the pattern's rarest label
           (the anchor) and lies within the anchor's eccentricity of it, since an embedding never lengthens a
           path.  So VF2 runs in that ball around each such node, the anchor pinned by a label nothing else
           has.  One pattern per worker.
  filter   refined filter + exact enumeration.  The plane is cut into square cores, each read through a window
           twice as wide: the core plus half a core of apron on every side.  For each (pattern, window) it
           checks label counts, one-hop candidates, candidate refinement, and a connected candidate group
           reaching the core.  Each test is necessary, so no owned occurrence is skipped.  Passing pairs are
           enumerated with LAD, each pattern node limited to window nodes of its label.  An occurrence is kept
           by the one window whose core holds its bounding-box centre.  One window per worker.

Both are exact, except that `filter` misses an occurrence that sticks out of its owner's window (one with a long
edge; rare on geometric graphs such as Delaunay triangulations); it never reports a false one.  hop_filter.py is the
exact, coordinate-free version of `filter`.  All were checked node set by node set against whole-graph VF2 (vf2.py).

  exact_match.py exact  target.gml patterns/ --workers 16 --output hits.jsonl
  exact_match.py filter target.gml list.txt  --workers 16 --output hits.jsonl --sets
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import igraph as ig
import numba as nb
import numpy as np

PIN = 1 << 20          # a label no node uses: pins the anchor
_W: dict = {}          # set before the workers fork, so they share it: target, patterns, tiling


def read_gml(path) -> ig.Graph:
    """A simple undirected GML graph with an integer `label` per node (and `x`, `y` for the target)."""
    g = ig.Graph.Read_GML(str(path))
    if g.is_directed() or not g.is_simple():
        raise ValueError(f'{path}: expected a simple undirected graph')
    g.vs['label'] = [int(v) for v in g.vs['label']]
    return g


def pattern_paths(args: list[str]) -> list[Path]:
    """.gml files, directories of them, or text files listing one path per line (relative to the list)."""
    out = []
    for a in map(Path, args):
        if a.is_dir():
            out += sorted(a.glob('*.gml'))
        elif a.suffix == '.gml':
            out.append(a)
        else:
            out += [a.parent / line.strip() for line in a.read_text().splitlines() if line.strip()]
    return out


# --------------------------------------------------------------------------- exact: seeded matching

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


# --------------------------------------------------------------------------- filter: graphs and tiling

def csr(n: int, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Adjacency of an undirected edge list: the neighbours of v are idx[ptr[v]:ptr[v + 1]]."""
    both = np.concatenate([edges, edges[:, ::-1]])
    both = both[np.argsort(both[:, 0], kind='stable')]
    return np.searchsorted(both[:, 0], np.arange(n + 1)), both[:, 1].copy()


def nbr_counts(labels: np.ndarray, ptr: np.ndarray, idx: np.ndarray, nl: int) -> np.ndarray:
    """out[v, l] = how many neighbours of v carry label l."""
    out = np.zeros((len(labels), nl), np.int32)
    np.add.at(out, (np.repeat(np.arange(len(labels)), np.diff(ptr)), labels[idx]), 1)
    return out


def pack(patterns: list, nl: int) -> tuple:
    """The patterns as one disjoint union for numba: node offsets, labels, neighbour-label counts, CSR, label counts."""
    off = np.concatenate([[0], np.cumsum([P.vcount() for P in patterns])]).astype(np.int64)
    lab = np.concatenate([P.vs['label'] for P in patterns]).astype(np.int64)
    edges = np.concatenate([np.asarray(P.get_edgelist(), np.int64).reshape(-1, 2) + o for P, o in zip(patterns, off)])
    ptr, idx = csr(len(lab), edges)
    cnt = np.stack([np.bincount(lab[a:b], minlength=nl) for a, b in zip(off[:-1], off[1:])]).astype(np.int64)
    return off, lab, nbr_counts(lab, ptr, idx, nl), ptr, idx, cnt


def tiling(xy: np.ndarray, core: float) -> dict:
    """Square cores of side `core` over the nodes' bounding box; nodes sorted by core."""
    x0 = xy.min(0)
    ncx, ncy = ((xy.max(0) - x0) // core).astype(np.int64) + 1
    b = ((xy - x0) // core).astype(np.int64)
    key = b[:, 1] * ncx + b[:, 0]
    order = np.argsort(key, kind='stable')
    return {'x0': x0, 'core': core, 'shape': (int(ncy), int(ncx)), 'order': order,
            'bptr': np.searchsorted(key[order], np.arange(ncx * ncy + 1))}


def window(t: dict, xy: np.ndarray, r: int, c: int) -> np.ndarray:
    """The nodes of core (r, c)'s window (the core and half a core around it, bounds included), sorted."""
    ncy, ncx = t['shape']
    bptr, lo_c, hi_c = t['bptr'], max(c - 1, 0), min(c + 1, ncx - 1)
    cand = np.concatenate([t['order'][bptr[rr * ncx + lo_c]:bptr[rr * ncx + hi_c + 1]]
                           for rr in range(max(r - 1, 0), min(r + 2, ncy))])
    lo = t['x0'] + (np.array([c, r]) - .5) * t['core']
    return np.sort(cand[((xy[cand] >= lo) & (xy[cand] <= lo + 2 * t['core'])).all(1)])


def calibrate(xy: np.ndarray, nodes: int) -> float:
    """The core side at which the median occupied window (at least nodes / 4) holds `nodes` nodes: the tiling
    the measurements used (make_tiles.py --target-nodes 128, 512 px windows at 256 px stride)."""
    core = float(np.ptp(xy, 0).max() * np.sqrt(nodes / len(xy)) / 2)       # from the bounding-box density
    for _ in range(5):
        t = tiling(xy, core)
        n = np.array([len(window(t, xy, r, c)) for r in range(t['shape'][0]) for c in range(t['shape'][1])])
        median = float(np.median(n[n >= nodes // 4]))
        if abs(median / nodes - 1) <= .02:
            break
        core *= np.sqrt(nodes / median)
    return core


# --------------------------------------------------------------------------- filter: necessary conditions

@nb.njit(cache=True)
def _cover(u, v, C, qptr, qidx, tptr, tidx, match_l, match_r, parent_r, visited, queue):
    """Can u's neighbours be matched to distinct neighbours of v, each u' to a w with C[u', w]?"""
    d = qptr[u + 1] - qptr[u]
    r0 = tptr[v]
    r = tptr[v + 1] - r0
    if d > r:
        return False
    for j in range(r):
        match_r[j] = -1
    for i in range(d):
        match_l[i] = -1
    for i in range(d):                     # augmenting paths, breadth first
        for j in range(r):
            visited[j] = False
        head, tail = 0, 1
        queue[0] = i
        found = -1
        while head < tail and found < 0:
            x = queue[head]
            head += 1
            ux = qidx[qptr[u] + x]
            for j in range(r):
                if not visited[j] and C[ux, tidx[r0 + j]]:
                    visited[j] = True
                    parent_r[j] = x
                    if match_r[j] == -1:
                        found = j
                        break
                    queue[tail] = match_r[j]
                    tail += 1
        if found < 0:
            return False
        y = found
        while y != -1:
            x = parent_r[y]
            prev = match_l[x]
            match_l[x] = y
            match_r[y] = x
            y = prev
    return True


@nb.njit(cache=True)
def _alive(C, ql, qcounts):
    """Every pattern node has a candidate, and each label has at least as many distinct candidate nodes as the
    pattern has nodes of that label."""
    m, n = C.shape
    for u in range(m):
        if not C[u].any():
            return False
    for lab in range(qcounts.shape[0]):
        if qcounts[lab] > 1:
            union = 0
            for v in range(n):
                for u in range(m):
                    if ql[u] == lab and C[u, v]:
                        union += 1
                        break
            if union < qcounts[lab]:
                return False
    return True


@nb.njit(cache=True)
def _refine(C, qptr, qidx, tptr, tidx, ql, qcounts, buf):
    """Drop C[u, v] where _cover fails, until nothing changes; False once a pattern node runs out."""
    m, n = C.shape
    match_l, match_r, parent_r, visited, queue = buf
    changed = True
    while changed:
        changed = False
        for u in range(m):
            for v in range(n):
                if C[u, v] and not _cover(u, v, C, qptr, qidx, tptr, tidx, match_l, match_r, parent_r, visited, queue):
                    C[u, v] = False
                    changed = True
        if not _alive(C, ql, qcounts):
            return False
    return True


@nb.njit(cache=True)
def _reaches_core(C, tptr, tidx, uv, lo, hi):
    """Is there a connected group of candidate nodes holding a candidate for every pattern node whose bounding box
    meets [lo, hi)^2?  Every occurrence lies in one such group, and its bounding-box centre lies in its box."""
    m, n = C.shape
    cand = np.zeros(n, np.bool_)
    for v in range(n):
        for u in range(m):
            if C[u, v]:
                cand[v] = True
                break
    comp = np.full(n, -1, np.int64)
    queue = np.empty(n, np.int64)
    has = np.zeros(m, np.bool_)
    for s0 in range(n):
        if not cand[s0] or comp[s0] >= 0:
            continue
        comp[s0] = s0
        head, tail = 0, 1
        queue[0] = s0
        x0 = x1 = uv[s0, 0]
        y0 = y1 = uv[s0, 1]
        has[:] = False
        while head < tail:
            v = queue[head]
            head += 1
            for u in range(m):
                if C[u, v]:
                    has[u] = True
            x0, x1 = min(x0, uv[v, 0]), max(x1, uv[v, 0])
            y0, y1 = min(y0, uv[v, 1]), max(y1, uv[v, 1])
            for k in range(tptr[v], tptr[v + 1]):
                w = tidx[k]
                if cand[w] and comp[w] < 0:
                    comp[w] = s0
                    queue[tail] = w
                    tail += 1
        if has.all() and x1 >= lo and x0 < hi and y1 >= lo and y0 < hi:
            return True
    return False


@nb.njit(cache=True)
def _passing(tl, tptr, tidx, tnbr, tcnt, lord, loff, uv, qoff, qlab, qnbr, qptr_all, qidx_all, qcnt):
    """Which patterns can have an occurrence owned by this window?  Label counts; one-hop candidates C[u, v]
    (v has u's label and at least u's count of neighbours of every label); refinement; a candidate group
    reaching the core.  uv: node positions in core sides from the core's corner, so the core is [0, 1)^2."""
    nq, nl = qcnt.shape
    n = tl.shape[0]
    deg, ml = 0, 0
    for v in range(n):
        deg = max(deg, tptr[v + 1] - tptr[v])
    for q in range(nq):
        ml = max(ml, qoff[q + 1] - qoff[q])
    buf = (np.empty(ml + 1, np.int64), np.empty(deg + 1, np.int64), np.empty(deg + 1, np.int64),
           np.empty(deg + 1, np.bool_), np.empty(ml + 1, np.int64))
    out = np.zeros(nq, np.bool_)
    for q in range(nq):
        ok = True
        for lab in range(nl):
            if qcnt[q, lab] > tcnt[lab]:
                ok = False
                break
        if not ok:
            continue
        o, m = qoff[q], qoff[q + 1] - qoff[q]
        ql = qlab[o:o + m]
        C = np.zeros((m, n), np.bool_)
        for u in range(m):
            found = False
            for s in range(loff[ql[u]], loff[ql[u] + 1]):
                v = lord[s]
                good = True
                for k in range(nl):
                    if tnbr[v, k] < qnbr[o + u, k]:
                        good = False
                        break
                if good:
                    C[u, v] = True
                    found = True
            if not found:
                ok = False
                break
        if not ok:
            continue
        qptr = qptr_all[o:o + m + 1] - qptr_all[o]
        qidx = qidx_all[qptr_all[o]:qptr_all[o + m]] - o
        if (_refine(C, qptr, qidx, tptr, tidx, ql, qcnt[q], buf)
                and _reaches_core(C, tptr, tidx, uv, -1e-6, 1 + 1e-6)):     # the margin only absorbs rounding
            out[q] = True
    return out


# --------------------------------------------------------------------------- filter: enumeration

def _filter_worker(k: int):
    """Cut window k, filter every pattern, enumerate the passing ones, keep what this core owns."""
    r, c = divmod(k, _W['tiling']['shape'][1])
    ptr, idx, nl = _W['ptr'], _W['idx'], _W['nl']
    members = window(_W['tiling'], _W['xy'], r, c)
    lens = ptr[members + 1] - ptr[members]
    src = np.repeat(np.arange(len(members)), lens)
    dst = idx[np.repeat(ptr[members] - np.cumsum(lens) + lens, lens) + np.arange(lens.sum())]
    pos = np.minimum(np.searchsorted(members, dst), len(members) - 1)
    inside = members[pos] == dst
    src, dst = src[inside], pos[inside]                 # in-window edges, both directions, sorted by src
    tl = _W['labels'][members]
    tptr = np.searchsorted(src, np.arange(len(members) + 1))
    lord = np.argsort(tl, kind='stable')
    loff = np.searchsorted(tl[lord], np.arange(nl + 1))
    guv = _W['uv'][members]
    ok = _passing(tl, tptr, dst, nbr_counts(tl, tptr, dst, nl), np.bincount(tl, minlength=nl), lord, loff,
                  guv - (c, r), *_W['packed'])
    out = []
    if not ok.any():
        return k, out
    G = ig.Graph(n=len(members), edges=np.column_stack([src, dst])[src < dst].tolist())
    by_label = [lord[loff[lab]:loff[lab + 1]].tolist() for lab in range(nl)]
    for j in np.flatnonzero(ok).tolist():
        P = _W['patterns'][j]
        start = time.perf_counter()
        maps = G.get_subisomorphisms_lad(P, domains=[by_label[lab] for lab in P.vs['label']], induced=False)
        owned = np.zeros((0, P.vcount()), np.int64)
        if maps:
            found = np.array(list({tuple(sorted(m)) for m in maps}), np.int64)       # one row per node set
            box = guv[found]                            # global core units, so every window computes the same owner
            owned = members[found[(np.floor((box.min(1) + box.max(1)) / 2) == (c, r)).all(1)]]
        out.append((j, time.perf_counter() - start, owned.tolist() if _W['sets'] else len(owned)))
    return k, out


def run_filter(G: ig.Graph, patterns: list, workers: int, sets: bool, core: float) -> tuple[list, dict]:
    t0 = time.perf_counter()
    labels = np.asarray(G.vs['label'], np.int64)
    nl = int(max(labels.max(), max(max(P.vs['label']) for P in patterns))) + 1
    xy = np.column_stack([G.vs['x'], G.vs['y']]).astype(np.float64)
    ptr, idx = csr(G.vcount(), np.asarray(G.get_edgelist(), np.int64).reshape(-1, 2))
    t = tiling(xy, core)
    ncy, ncx = t['shape']
    load = np.pad(np.diff(t['bptr']).reshape(ncy, ncx), 1)           # nodes in the 3x3 cores around each core
    load = sum(load[1 + dy:1 + dy + ncy, 1 + dx:1 + dx + ncx] for dy in (-1, 0, 1) for dx in (-1, 0, 1)).ravel()
    tasks = [int(k) for k in np.argsort(-load, kind='stable') if load[k]]          # the heaviest windows first
    _W.update(labels=labels, xy=xy, uv=(xy - t['x0']) / core, ptr=ptr, idx=idx, nl=nl, tiling=t,
              patterns=patterns, packed=pack(patterns, nl), sets=sets)
    prepare = time.perf_counter() - t0
    t1 = time.perf_counter()
    _filter_worker(tasks[0])                            # numba compile (or cache load) before the workers fork
    one_off = {'numba_s': time.perf_counter() - t1}
    rows = [{'occurrences': 0, 'seconds': 0., 'pairs': 0, **({'sets': []} if sets else {})} for _ in patterns]
    with mp.get_context('fork').Pool(workers) as pool:
        t2 = time.perf_counter()
        for _k, hits in pool.imap_unordered(_filter_worker, tasks, chunksize=1):
            for j, s, owned in hits:
                row = rows[j]
                row['pairs'] += 1
                row['seconds'] += s
                row['occurrences'] += len(owned) if sets else owned
                if sets:
                    row['sets'] += owned
        search = time.perf_counter() - t2
    return rows, {'prepare_s': prepare, 'search_s': search, 'core': core, 'grid': [ncy, ncx],
                  'windows': len(tasks), **one_off}


# --------------------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('method', choices=('exact', 'filter'))
    ap.add_argument('target', type=Path, help='target GML: integer `label`, and `x`, `y` for filter, per node')
    ap.add_argument('patterns', nargs='+', help='pattern GMLs, directories of them, or text files listing them')
    ap.add_argument('--workers', type=int, default=len(os.sched_getaffinity(0)))
    ap.add_argument('--core', type=float, default=None, help='filter: core side in target units (default: calibrated)')
    ap.add_argument('--window-nodes', type=int, default=128,
                    help='filter, without --core: the median nodes per occupied window to calibrate the core to')
    ap.add_argument('--sets', action='store_true', help='also write each occurrence as its sorted target node ids')
    ap.add_argument('--output', type=Path, default=None,
                    help='JSON lines, one per pattern: occurrences, seconds (exact: its search; filter: its enumeration '
                         'summed over its passing windows, with pairs = how many), and sets with --sets')
    args = ap.parse_args()

    t0 = time.perf_counter()
    G = read_gml(args.target)
    paths = pattern_paths(args.patterns)
    patterns = [read_gml(p) for p in paths]
    for p, P in zip(paths, patterns):
        if not P.is_connected():
            raise ValueError(f'{p}: pattern is not connected')
    info = {'method': args.method, 'target': str(args.target), 'nodes': G.vcount(), 'edges': G.ecount(),
            'patterns': len(patterns), 'workers': args.workers, 'read_s': time.perf_counter() - t0}
    if args.method == 'exact':
        rows, timing = run_exact(G, patterns, args.workers, args.sets)
    else:
        if not {'x', 'y'} <= set(G.vs.attributes()):
            raise SystemExit(f'{args.target}: filter tiles the plane, so every target node needs x and y')
        core = args.core
        if core is None:
            t1 = time.perf_counter()
            core = calibrate(np.column_stack([G.vs['x'], G.vs['y']]).astype(np.float64), args.window_nodes)
            info['calibrate_s'] = time.perf_counter() - t1
        rows, timing = run_filter(G, patterns, args.workers, args.sets, core)
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
