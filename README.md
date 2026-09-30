# subgraph_matching

Exact labelled subgraph matching on large graphs, and a harness that runs published exact matchers on the same inputs so all methods can be checked against each other and timed.

**The problem.** You have a target graph whose nodes carry integer labels and a small, connected pattern graph with labelled nodes. The task is to find every occurrence of the pattern. An occurrence is a set of target nodes onto which the pattern maps one-to-one, with each node going to a node of the same label and each pattern edge to a target edge. Extra target edges are allowed, so this is a labelled monomorphism (non-induced subgraph matching). Each node set counts once.

## Methods in this repository

All are Python (numpy, python-igraph, numba) and run in parallel over worker processes (`--workers`).

| Script | Method | Exact? | Needs coordinates? |
| --- | --- | --- | --- |
| `methods/hop_filter.py` | Hop windows. The target is split into cores of nearby nodes (a k-hop net). A pattern of radius r is searched within r hops of each core: cheap necessary tests first (label counts, neighbour-label counts, candidate refinement), then enumeration with LAD. Each occurrence is kept by exactly one core. | yes | no |
| `methods/exact_match.py exact` | Seeded exact. VF2 runs in the ball around each target node carrying the pattern's rarest label, with radius equal to that pattern node's eccentricity. | yes | no |
| `methods/exact_match.py filter` | Spatial windows: the same tests and enumeration as hop windows, on square tiles of the plane. | misses an occurrence that sticks out of its tile (possible only through a long edge) | yes (`x`, `y`) |
| `methods/vf2.py` | Plain whole-graph igraph VF2, deduplicated by node set. The reference. | yes | no |

**Inputs.** GML files: the target has an integer `label` per node, plus `x`, `y` for the spatial filter. Each pattern has an integer `label` per node, and patterns must be connected. A pattern list is a text file with one pattern GML path per line, relative to the list's own folder. A folder of `.gml` files also works.

**Output.** One JSON line per pattern (`occurrences`, `seconds`; with `--sets`, the node sets as sorted target node ids), plus a one-line JSON summary on stdout with the timings.

## Published matchers

Five published exact matchers, all run as monomorphism, each through a small batch driver in `benchmark/drivers/`:

| Driver (`build/bin/`) | Method | Code | Input | Reports per query |
| --- | --- | --- | --- | --- |
| `ri_batch` | RI 3.6 (Bonnici et al., BMC Bioinformatics 2013) | [InfOmics/RI](https://github.com/InfOmics/RI) | `.gfu` | mappings and distinct node sets |
| `vf3_batch_VF3` | VF3 (Carletti et al., TPAMI 2018) | [MiviaLab/vf3lib](https://github.com/MiviaLab/vf3lib) | `.vft` | mappings |
| `vf3_batch_VF3L` | VF3L, VF3 without look-ahead (Carletti et al., 2019) | [MiviaLab/vf3lib](https://github.com/MiviaLab/vf3lib) | `.vft` | mappings |
| `daf_batch` | DAF (Han et al., SIGMOD 2019) | [SNUCSE-CTA/DAF](https://github.com/SNUCSE-CTA/DAF) | `.daf` | mappings |
| `fastiso_batch` | FASTiso (Coti, Reinharz & Agbeto, bioRxiv 2025) | [cbe/fastiso](https://gitlab.info.uqam.ca/cbe/fastiso) | `.gfu` | mappings |

A mapping is one assignment of pattern nodes to target nodes. A node set with several valid assignments (from pattern symmetries) has several mappings; our methods and RI count node sets.

Their code is not stored here. `benchmark/setup_tools.sh` fetches each tool at the commit listed in `benchmark/tools.lock` and builds the drivers. Each driver:
- reads the target and all queries first;
- forks the requested number of workers, which take one query at a time from a shared counter;
- prints the wall time for all queries;
- writes one TSV row per query: query index, seconds, mappings, then (`ri_batch` only) distinct node sets, then timed out (0/1).

The drivers count but don't write out the matches.

Only matchers with public code and no cap on query size are included. GraphQL, CFL, DP-iso, CECI, VEQ, RapidMatch, ASDMatch, GuP and EGSM all have a fixed maximum query size.

## Install

```bash
pip install -r requirements.txt                                     # tested: Python 3.10, numpy 2.2, python-igraph 1.0, numba 0.67
BOOST_INCLUDE=/path/to/boost/include benchmark/setup_tools.sh       # the published matchers -> build/bin/ (optional)
```

`setup_tools.sh` needs git, a C++17 compiler (`CXX`, default `g++`) and Boost headers, which DAF uses. It was tested with GCC 14.2 and Boost 1.86.

## Example: a toy pattern search

`examples/make_example.py` builds a toy target and patterns:
- **Target:** 1,600 nodes on a jittered triangular lattice, with 5 labels from common to rare and `x`, `y` coordinates.
- **Patterns:** 40 connected patterns of 3–8 nodes grown from the target, so each occurs at least once. Every eighth has its labels shuffled, so it may not occur.

```bash
python examples/make_example.py work/
```

**Our methods:**

```bash
python methods/hop_filter.py   work/target.gml work/patterns.txt --workers 4 --output work/hop.jsonl
python methods/exact_match.py exact  work/target.gml work/patterns.txt --workers 4 --output work/seeded.jsonl
python methods/exact_match.py filter work/target.gml work/patterns.txt --workers 4 --output work/filter.jsonl --sets
python methods/vf2.py          work/target.gml work/patterns.txt --workers 4 --output work/vf2.jsonl
```

**A published matcher:** convert the same inputs to the tools' formats, then call its driver with target, query list, workers, output TSV and per-query limit in seconds:

```bash
python benchmark/convert.py work/target.gml work/patterns.txt work/input   # target.{gfu,vft,daf}, q/NNNN.*, queries_*.txt
build/bin/ri_batch       work/input/target.gfu work/input/queries_gfu.txt 4 work/ri.tsv 300
build/bin/vf3_batch_VF3  work/input/target.vft work/input/queries_vft.txt 4 work/vf3.tsv 300
build/bin/vf3_batch_VF3L work/input/target.vft work/input/queries_vft.txt 4 work/vf3l.tsv 300
build/bin/daf_batch      work/input/target.daf work/input/queries_daf.txt 4 work/daf.tsv 300
build/bin/fastiso_batch  work/input/target.gfu work/input/queries_gfu.txt 4 work/fastiso.tsv 300
```

**All of them at once, then checked against each other:**

```bash
VF2=1 benchmark/run_all.sh work/target.gml work/patterns.txt results/example 4
python benchmark/summarize.py results/example
```

```
== example
  4 workers, reference VF2: 26,580 occurrences
    exact: VF2 finished 40/40; hop 40/40; seeded 40/40; filter 40/40 (short 0); RI mappings 40/40, node sets 40/40;
           VF3 mappings 40/40; VF3L mappings 40/40; DAF mappings 40/40; FASTiso mappings 40/40
    wall s: RI 0.02, VF3 0.02, VF3L 0.01, DAF 0.01, FASTiso 0.01, hop 0.23 (+0.00), filter 0.43, seeded exact 0.15, VF2 0.24
```

**How to read it.** Every method returned exactly VF2's occurrences for all 40 patterns, and every tool's mapping count equals hop windows' count. `python benchmark/tables.py results/*` turns several runs into a markdown table of wall times. The toy example is too small for the times to mean anything; run it on your own graphs for that.

## Comparing the methods on your own graph

```bash
benchmark/run_all.sh <target.gml> <patterns.txt> <out_dir> <workers> [core]
```

This runs every method on one machine so the times compare: the five tools, `hop_filter.py`, both `exact_match.py` methods, and `vf2.py` if `VF2=1`.
- **Queries:** every query, with a per-query time limit (`LIMIT`, default 300 s).
- **Spatial filter's tile size:** `core` sets it in coordinate units; without it the filter calibrates its own. The filter is skipped when the target has no coordinates.
- **Other settings:** `BIN` sets where the drivers are, and `PYTHON` which interpreter runs our methods.

**What is timed:**
- **The tools:** they prepare the target before their clock starts, so they are timed from the first query.
- **`hop_filter.py` and seeded exact:** timed by their search alone. Hop windows' setup is shown separately in brackets.
- **The spatial filter:** its setup and search together.
- **In every case:** reading the files and starting the workers are not counted.

`python tests/test_methods.py` repeats the checks on the toy example: every method against VF2 node set by node set, and, if `build/bin` exists, every tool's counts. It also runs under `pytest`.

## Notes

- **What the tools output:** RI's driver deduplicates mappings into node sets, as our methods do. VF3, VF3L, DAF and FASTiso only count mappings, so their times include no deduplication.
- **Harness changes, all outside the tools' algorithms:**
  - DAF crashes when a query's candidate space cannot be built. That means no match exists, so its driver records 0.
  - The `.vft` and `.daf` targets carry one isolated placeholder node for each label, up to the largest, that the target lacks. No pattern node can map onto one.
  - VF3 and VF3L run each query's steps from the authors' `main.cpp` once, without its repeat loop.
- **DAF's leaf counting:** DAF counts degree-1 query nodes combinatorially instead of enumerating them.
- **Licences:** this repository is MIT. The published tools keep their own licences (`benchmark/tools.lock`): vf3lib is LGPL-3.0, DAF is Apache-2.0, RI has its own terms, and FASTiso's repository has no licence file. They are fetched from their authors, not redistributed.

## Layout

```
methods/            hop_filter.py, exact_match.py, vf2.py
benchmark/          setup_tools.sh, tools.lock, drivers/*.cpp, convert.py, run_all.sh, summarize.py, tables.py
examples/           make_example.py
tests/              test_methods.py
```
