#!/usr/bin/env python3
"""Candidate filtering + exact enumeration on hop-distance windows: every occurrence of each pattern in a labelled
target graph, with no coordinates needed.  The matching definition is that of seeded_exact.py and vf2.py.

  cores       a k-hop net of seed nodes (taken in node order, each more than k hops from every earlier seed), and
              every node assigned to its nearest seed (graph Voronoi, ties to the earlier seed)
  windows     a pattern of radius r (the eccentricity of its most central nodes) is searched in the nodes within
              r hops of a core
  ownership   an occurrence S belongs to the core holding o(S): the smallest node index that any mapping onto S
              sends a central pattern node to.  Every node of S lies within r hops of o(S), since an embedding never
              lengthens a path.  So the owner's window holds all of S, and LAD finds every mapping onto S there.
              Every other window that holds S computes the same o(S), so each occurrence is counted exactly once.
  filter      per (pattern, window): label counts, one-hop candidates, refinement (common.py); then, on the
              graph, a connected group of candidates covering every pattern node and holding a core node that is a
              candidate for a central pattern node (o(S) is one)

  hop_filter.py target.gml patterns.txt --workers 16 --output hits.jsonl [--sets] [--core-hops 4]
"""
from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import igraph as ig
import numba as nb
import numpy as np

import common

_H: dict = {}          # set before the workers fork


def gather(ptr: np.ndarray, idx: np.ndarray, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(source, neighbour) for every edge leaving `nodes`."""
    lens = ptr[nodes + 1] - ptr[nodes]
    return np.repeat(nodes, lens), idx[np.repeat(ptr[nodes] - np.cumsum(lens) + lens, lens) + np.arange(lens.sum())]


@nb.njit(cache=True)
def _net(ptr, idx, k):
    """Seed nodes in node order, each more than k hops from every earlier seed."""
    n = ptr.shape[0] - 1
    covered = np.zeros(n, np.bool_)
    stamp = np.full(n, -1, np.int64)                # the seed whose search last reached the node
    depth = np.zeros(n, np.int64)
    queue = np.empty(n, np.int64)
    seeds = np.empty(n, np.int64)
    ns = 0
    for v in range(n):
        if covered[v]:
            continue
        seeds[ns] = v
        ns += 1
        head, tail = 0, 1
        queue[0] = v
        stamp[v], depth[v], covered[v] = v, 0, True
        while head < tail:
            x = queue[head]
            head += 1
            if depth[x] == k:
                continue
            for e in range(ptr[x], ptr[x + 1]):
                y = idx[e]
                if stamp[y] != v:
                    stamp[y], depth[y], covered[y] = v, depth[x] + 1, True
                    queue[tail] = y
                    tail += 1
    return seeds[:ns]


def make_cores(ptr: np.ndarray, idx: np.ndarray, k: int) -> np.ndarray:
    """Each node's core: its nearest seed of a k-hop net, as a core number."""
    n = len(ptr) - 1
    core = np.full(n, -1, np.int64)
    frontier = _net(ptr, idx, k)
    core[frontier] = np.arange(len(frontier))
    while len(frontier):
        src, dst = gather(ptr, idx, frontier)
        new = core[dst] < 0
        src, dst = src[new], dst[new]
        if not len(dst):
            break
        o = np.lexsort((core[src], dst))                     # per node, the earliest seed reaching it first
        src, dst = src[o], dst[o]
        first = np.flatnonzero(np.r_[True, dst[1:] != dst[:-1]])
        core[dst[first]] = core[src[first]]
        frontier = dst[first]
    return core


# --------------------------------------------------------------------------- filter

@nb.njit(cache=True)
def _holds_core(C, tptr, tidx, incore, central):
    """Is there a connected group of candidate nodes holding a candidate for every pattern node and a core node
    that is a candidate for a central pattern node?"""
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
        has[:] = False
        hit = False
        while head < tail:
            v = queue[head]
            head += 1
            for u in range(m):
                if C[u, v]:
                    has[u] = True
                    if incore[v] and central[u]:
                        hit = True
            for k in range(tptr[v], tptr[v + 1]):
                w = tidx[k]
                if cand[w] and comp[w] < 0:
                    comp[w] = s0
                    queue[tail] = w
                    tail += 1
        if hit and has.all():
            return True
    return False


@nb.njit(cache=True)
def _passing(tl, tptr, tidx, tnbr, tcnt, lord, loff, incore, qoff, qlab, qnbr, qptr_all, qidx_all, qcnt, qcentral):
    """Which patterns can have an occurrence owned by this window: label counts; one-hop candidates C[u, v] (v has
    u's label and at least u's count of neighbours of every label); refinement; a candidate group holding a core
    node that is a candidate for a central pattern node (_holds_core)."""
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
        if (common._refine(C, qptr, qidx, tptr, tidx, ql, qcnt[q], buf)
                and _holds_core(C, tptr, tidx, incore, qcentral[o:o + m])):
            out[q] = True
    return out


# --------------------------------------------------------------------------- windows

@nb.njit(cache=True)
def _grow(ptr, idx, core_nodes, R, dist, queue):
    """The nodes within R hops of the core, sorted, and each one's hop distance from it.  dist: scratch, all -1, left so."""
    tail = 0
    for v in core_nodes:
        dist[v] = 0
        queue[tail] = v
        tail += 1
    head = 0
    while head < tail:
        x = queue[head]
        head += 1
        if dist[x] == R:
            continue
        for e in range(ptr[x], ptr[x + 1]):
            y = idx[e]
            if dist[y] < 0:
                dist[y] = dist[x] + 1
                queue[tail] = y
                tail += 1
    nodes = np.sort(queue[:tail])
    d = np.empty(tail, np.int64)
    for i in range(tail):
        d[i] = dist[nodes[i]]
    for i in range(tail):
        dist[queue[i]] = -1
    return nodes, d


@nb.njit(cache=True)
def _level(ptr, idx, labels, nl, members, local):
    """The window graph on `members` (sorted nodes): labels, local CSR, neighbour-label counts, label counts, and the
    nodes by label (lord, loff).  local: scratch, all -1, left so."""
    n = members.shape[0]
    for i in range(n):
        local[members[i]] = i
    tl = np.empty(n, np.int64)
    tptr = np.zeros(n + 1, np.int64)
    for i in range(n):
        v = members[i]
        tl[i] = labels[v]
        c = 0
        for e in range(ptr[v], ptr[v + 1]):
            if local[idx[e]] >= 0:
                c += 1
        tptr[i + 1] = tptr[i] + c
    tidx = np.empty(tptr[n], np.int64)
    nbr = np.zeros((n, nl), np.int32)
    k = 0
    for i in range(n):
        v = members[i]
        for e in range(ptr[v], ptr[v + 1]):
            j = local[idx[e]]
            if j >= 0:
                tidx[k] = j
                k += 1
                nbr[i, tl[j]] += 1
    for i in range(n):
        local[members[i]] = -1
    tcnt = np.zeros(nl, np.int64)
    for i in range(n):
        tcnt[tl[i]] += 1
    lord = np.argsort(tl, kind='mergesort')
    loff = np.searchsorted(tl[lord], np.arange(nl + 1))
    return tl, tptr, tidx, nbr, tcnt, lord, loff


# --------------------------------------------------------------------------- enumeration

def _worker(c: int):
    """Grow core c's window; at each pattern radius, filter that radius's patterns, enumerate the passing ones, keep
    what this core owns."""
    ptr, idx, nl = _H['ptr'], _H['idx'], _H['nl']
    nodes, d = _grow(ptr, idx, _H['order'][_H['cptr'][c]:_H['cptr'][c + 1]], _H['levels'][-1][0], _H['dist'], _H['queue'])
    out, sizes = [], []
    for r, js, packed in _H['levels']:
        inner = d <= r
        members = nodes[inner]
        sizes.append(len(members))
        tl, tptr, tidx, nbr, tcnt, lord, loff = _level(ptr, idx, _H['labels'], nl, members, _H['local'])
        ok = _passing(tl, tptr, tidx, nbr, tcnt, lord, loff, d[inner] == 0, *packed)
        if not ok.any():
            continue
        src = np.repeat(np.arange(len(members)), np.diff(tptr))
        G = ig.Graph(n=len(members), edges=np.column_stack([src, tidx])[src < tidx].tolist())
        by_label = [lord[loff[lab]:loff[lab + 1]].tolist() for lab in range(nl)]
        for q in np.flatnonzero(ok).tolist():
            j = js[q]
            P, central = _H['patterns'][j], _H['central'][j]
            start = time.perf_counter()
            maps = G.get_subisomorphisms_lad(P, domains=[by_label[lab] for lab in P.vs['label']], induced=False)
            owned, nmaps = np.zeros((0, P.vcount()), np.int64), 0
            if maps:
                g = members[np.asarray(maps, np.int64)]          # global node indices, one row per mapping
                key = np.sort(g, axis=1)
                own = g[:, central].min(1)
                o = np.lexsort((own, *key.T[::-1]))              # by node set, then by o: the first row of a set has o(S)
                key, own = key[o], own[o]
                first = np.flatnonzero(np.r_[True, (key[1:] != key[:-1]).any(1)])
                mine = _H['core'][own[first]] == c
                owned = key[first][mine]
                nmaps = int(np.diff(np.r_[first, len(key)])[mine].sum())     # every mapping onto an owned set is here
            out.append((j, time.perf_counter() - start, owned.tolist() if _H['sets'] else len(owned), nmaps))
    return c, out, sizes


def run_hop(G: ig.Graph, patterns: list, workers: int, sets: bool, k: int) -> tuple[list, dict]:
    tw = time.perf_counter()
    _net(np.zeros(2, np.int64), np.zeros(0, np.int64), 1)       # numba compile (or cache load), not charged
    warm = time.perf_counter() - tw
    t0 = time.perf_counter()
    labels = np.asarray(G.vs['label'], np.int64)
    nl = int(max(labels.max(), max(max(P.vs['label']) for P in patterns))) + 1
    edges = np.fromiter(itertools.chain.from_iterable(G.get_edgelist()), np.int64, count=2 * G.ecount()).reshape(-1, 2)
    ptr, idx = common.csr(G.vcount(), edges)
    core = make_cores(ptr, idx, k)
    ncores = int(core.max()) + 1
    order = np.argsort(core, kind='stable')
    cptr = np.searchsorted(core[order], np.arange(ncores + 1))
    ecc = [np.asarray(P.eccentricity(), np.int64) for P in patterns]
    radius = [int(e.min()) for e in ecc]
    central = [np.flatnonzero(e == e.min()) for e in ecc]
    levels = []
    for r in sorted(set(radius)):
        js = [j for j in range(len(patterns)) if radius[j] == r]
        packed = common.pack([patterns[j] for j in js], nl)
        qcentral = np.concatenate([ecc[j] == radius[j] for j in js])
        levels.append((r, js, (*packed, qcentral)))
    tasks = [int(c) for c in np.argsort(-np.diff(cptr), kind='stable')]            # the largest cores first
    _H.update(labels=labels, ptr=ptr, idx=idx, nl=nl, core=core, order=order, cptr=cptr, levels=levels,
              patterns=patterns, central=central, sets=sets, dist=np.full(G.vcount(), -1, np.int64),
              local=np.full(G.vcount(), -1, np.int64), queue=np.empty(G.vcount(), np.int64))
    prepare = time.perf_counter() - t0
    t1 = time.perf_counter()
    _worker(tasks[0])                                   # numba compile (or cache load) before the workers fork
    numba_s = warm + time.perf_counter() - t1
    rows = [{'occurrences': 0, 'seconds': 0., 'pairs': 0, 'mappings': 0, **({'sets': []} if sets else {})} for _ in patterns]
    window = np.zeros(len(levels))
    with mp.get_context('fork').Pool(workers) as pool:
        t2 = time.perf_counter()
        # ~16 batches per worker: one pipe round trip per core costs ms when the pool sits on idle cores (tens of
        # thousands of light cores), while batches stay small enough to balance (the largest cores go first)
        chunk = max(1, len(tasks) // (16 * workers))
        for _c, hits, sizes in pool.imap_unordered(_worker, tasks, chunksize=chunk):
            window += sizes
            for j, s, owned, nmaps in hits:
                row = rows[j]
                row['pairs'] += 1
                row['mappings'] += nmaps
                row['seconds'] += s
                row['occurrences'] += len(owned) if sets else owned
                if sets:
                    row['sets'] += owned
        search = time.perf_counter() - t2
    return rows, {'prepare_s': prepare, 'search_s': search, 'numba_s': numba_s, 'core_hops': k, 'cores': ncores,
                  'mean_core_nodes': G.vcount() / ncores,
                  'mean_window_nodes': {str(r): float(w / ncores) for (r, _j, _p), w in zip(levels, window)}}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('target', type=Path, help='target GML: integer `label` per node (no coordinates needed)')
    ap.add_argument('patterns', nargs='+', help='pattern GMLs, directories of them, or text files listing them')
    ap.add_argument('--workers', type=int, default=len(os.sched_getaffinity(0)))
    ap.add_argument('--core-hops', type=int, default=4, help='seed spacing of the cores, in hops (about 34 nodes per core on a Delaunay graph)')
    ap.add_argument('--sets', action='store_true', help='also write each occurrence as its sorted target node ids')
    ap.add_argument('--output', type=Path, default=None, help='JSON lines, one per pattern: occurrences, mappings, pairs, seconds, and sets with --sets')
    args = ap.parse_args()
    t0 = time.perf_counter()
    G, paths, patterns = common.load(args.target, args.patterns)
    info = {'method': 'hop', 'target': str(args.target), 'nodes': G.vcount(), 'edges': G.ecount(),
            'patterns': len(patterns), 'workers': args.workers, 'read_s': time.perf_counter() - t0}
    rows, timing = run_hop(G, patterns, args.workers, args.sets, args.core_hops)
    info.update(timing, total_s=timing['prepare_s'] + timing['search_s'], occurrences=sum(r['occurrences'] for r in rows))
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
