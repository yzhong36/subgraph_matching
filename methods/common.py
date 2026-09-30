#!/usr/bin/env python3
"""Shared by the methods: GML input, and the candidate-filter kernels of hop_filter.py.

  load                      the target and the patterns, their labels recoded as 0..L-1 over all of them
  read_gml, pattern_paths   one simple undirected GML graph with a `label` per node; pattern lists
  csr, nbr_counts, pack     adjacency arrays, neighbour-label counts, and the patterns packed for numba
  _refine (_cover, _alive)  candidate refinement: drop a candidate v of pattern node u unless u's neighbours can take
                            distinct neighbours of v, each one of its own candidates (bipartite matching), and require
                            enough distinct candidates per label; repeat until nothing changes.  Every test is one an
                            occurrence must pass, so no occurrence is ever filtered out.
"""
from __future__ import annotations

from pathlib import Path

import igraph as ig
import numba as nb
import numpy as np


def read_gml(path) -> ig.Graph:
    """A simple undirected GML graph with a `label` on every node, labels as GML gives them (numbers or strings)."""
    g = ig.Graph.Read_GML(str(path))
    if g.is_directed() or not g.is_simple():
        raise ValueError(f'{path}: expected a simple undirected graph (no direction, loops or repeated edges)')
    if 'label' not in g.vs.attributes() or any(v is None or v != v or v == '' for v in g.vs['label']):
        raise ValueError(f'{path}: every node needs a label')
    return g


def label_key(v):
    """A GML label as a key: GML numbers arrive as floats, so whole ones become ints; strings stay strings."""
    return int(v) if isinstance(v, float) and v.is_integer() else v


def load(target, patterns: list[str]) -> tuple[ig.Graph, list[Path], list[ig.Graph]]:
    """The target, the pattern paths and the patterns, with every label recoded as 0..L-1 over all of them (the
    target's labels first), so labels may be integers of any size or strings.  Labels match by value.  Patterns must
    be connected."""
    G = read_gml(target)
    paths = pattern_paths(patterns)
    names = [p.name for p in paths]
    if len(set(names)) < len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f'pattern file names must be unique (results are keyed by them): {dup[:5]}')
    Ps = [read_gml(p) for p in paths]
    for p, P in zip(paths, Ps):
        if not P.is_connected():
            raise ValueError(f'{p}: pattern is not connected')
    code: dict = {}
    for g in (G, *Ps):
        for v in g.vs['label']:
            code.setdefault(label_key(v), len(code))
    for g in (G, *Ps):
        g.vs['label'] = [code[label_key(v)] for v in g.vs['label']]
    return G, paths, Ps


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



# --------------------------------------------------------------------------- adjacency and patterns

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



# --------------------------------------------------------------------------- candidate refinement

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
