# Benchmarking exact labelled subgraph matching

Exact labelled subgraph matching on large graphs, and a harness that runs published exact matchers on the same inputs so all methods can be checked against each other and timed.

**The problem.** You have a target graph whose nodes carry labels and a small, connected pattern graph with labelled nodes. The task is to find every occurrence of the pattern. An occurrence is a set of target nodes onto which the pattern maps one-to-one, with each node going to a node of the same label and each pattern edge to a target edge. Extra target edges are allowed, so this is a labelled monomorphism (non-induced subgraph matching). Each node set counts once.

## Methods in this repository

All are exact, need no node coordinates, and are Python (numpy, python-igraph, numba) running in parallel over worker processes (`--workers`).

| Script | Method |
| --- | --- |
| `methods/hop_filter.py` | Hop windows. The target is split into cores of nearby nodes (a k-hop net). A pattern of radius r is searched within r hops of each core: cheap necessary tests first (label counts, neighbour-label counts, candidate refinement), then enumeration with LAD. Each occurrence is kept by exactly one core. `--core-hops` (default 4) sets the core spacing. |
| `methods/seeded_exact.py` | Seeded exact. VF2 runs in the ball around each target node carrying the pattern's rarest label, with radius equal to that pattern node's eccentricity. |
| `methods/vf2.py` | Plain whole-graph igraph VF2, deduplicated by node set. The reference. |

`methods/common.py` holds the GML reading and the candidate-filtering kernels they share.

**Inputs.** GML files, the same for every method and for the published tools:
- **Graphs:** the target and the patterns are simple undirected graphs (no direction, self-loops or repeated edges). Patterns must be connected; a single node is allowed.
- **Labels:** every node has a `label`, either an integer of any size or a string. Labels match by value between the target and the patterns. Other node and edge attributes are ignored.
- **Patterns:** given as `.gml` files, a folder of them, or a list file with one path per line, relative to the list's own folder.
- **Rejected inputs:** anything else stops with an error saying what is wrong.

**Output.** One JSON line per pattern, with `occurrences` and `seconds`, and with `--sets` the node sets as sorted target node ids. Each script also prints a one-line JSON summary of its timings.
- `hop_filter.py` adds `mappings` (every assignment of pattern nodes, not only node sets) and `pairs` (windows enumerated).
- `vf2.py` always writes the node sets. A query still running after `--limit` seconds is marked `capped` and has no count.

## Published matchers

Five published exact matchers, all run as monomorphism, each through a small batch driver in `benchmark/drivers/`:

| Driver (`build/bin/`) | Method | Code | Input | Follows the authors' |
| --- | --- | --- | --- | --- |
| `ri_batch` | RI 3.6 (Bonnici et al., BMC Bioinformatics 2013) | [InfOmics/RI](https://github.com/InfOmics/RI) | `.gfu` | `ri3.cpp` (`mono gfu`) |
| `vf3_batch_VF3` | VF3 (Carletti et al., TPAMI 2018) | [MiviaLab/vf3lib](https://github.com/MiviaLab/vf3lib) | `.vft` | `main.cpp` (`-u -e -r 0`) |
| `vf3_batch_VF3L` | VF3L, VF3 without look-ahead (Carletti et al., 2019) | [MiviaLab/vf3lib](https://github.com/MiviaLab/vf3lib) | `.vft` | `main.cpp` (`-u -e -r 0`) |
| `daf_batch` | DAF (Han et al., SIGMOD 2019) | [SNUCSE-CTA/DAF](https://github.com/SNUCSE-CTA/DAF) | `.daf` | `main/main.cc` |
| `fastiso_batch` | FASTiso (Coti, Reinharz & Agbeto, bioRxiv 2025) | [cbe/fastiso](https://gitlab.info.uqam.ca/cbe/fastiso) | `.gfu` | `fast_iso.cpp` (`-A monomorphism -F gfu -G -S count-solution`) |

Their code is not stored here. `benchmark/setup_tools.sh` fetches each tool at the commit listed in `benchmark/tools.lock`, unmodified, and compiles the drivers against it.

**What a driver runs.** Per query, a driver calls the same functions, in the same order and with the same options, as the authors' own program in the last column does for one query. The tools' code is not changed, and apart from two guards in `daf_batch` for queries DAF's own program gets wrong (see Notes), what a driver adds is only around those calls:
- it reads the target and all queries first;
- it forks the requested number of workers, which take one query at a time from a shared counter, with a per-query time limit;
- it keeps its result rows in memory and prints the wall time for all queries, from before the workers start until the last query finishes, so no file-system write is timed;
- it writes one TSV row per query: query index, seconds, mappings, timed out (0/1).

A mapping is one assignment of pattern nodes to target nodes, which is what every tool counts. A node set with several valid assignments (from pattern symmetries) has several mappings. Our methods count node sets, and `hop_filter.py` also reports mappings so the two can be compared.

**Match output, for verification.** A sixth argument (after the time limit) names a file for the matches. The driver then writes one line per distinct matched node set of each query: the query index, then the set's target node indices, sorted. It gets the matches through each tool's own interface: a match listener for RI (in place of `ri3.cpp`'s printing one), vf3lib's stored solutions, and FASTiso's solution list (`-S print-solution`). Every driver except `daf_batch` does this. DAF counts the matches of degree-1 pattern nodes combinatorially instead of enumerating them, so it can only be checked by its counts. Without the argument no matches are kept, so timings are unaffected.

Only matchers with public code and no cap on query size are included. GraphQL, CFL, DP-iso, CECI, VEQ, RapidMatch, ASDMatch, GuP and EGSM all have a fixed maximum query size.

## Install

```bash
pip install -r requirements.txt                                     # tested: Python 3.10, numpy 2.2, python-igraph 1.0, numba 0.67
BOOST_INCLUDE=/path/to/boost/include benchmark/setup_tools.sh       # the published matchers -> build/bin/ (optional)
```

`setup_tools.sh` needs git, a C++17 compiler (`CXX`, default `g++`) and Boost headers, which DAF uses. It was tested with GCC 14.2 and Boost 1.86.

## Example: a toy pattern search

`examples/make_example.py` builds a toy target and patterns:
- **Target:** 1,600 nodes on a jittered triangular lattice, with 5 labels from common to rare.
- **Patterns:** 40 connected patterns of 3–8 nodes grown from the target, so each occurs at least once. Every eighth has its labels shuffled, so it may not occur.

```bash
python examples/make_example.py work/
```

**Our methods:**

```bash
python methods/hop_filter.py   work/target.gml work/patterns.txt --workers 4 --output work/hop.jsonl --sets
python methods/seeded_exact.py work/target.gml work/patterns.txt --workers 4 --output work/seeded.jsonl
python methods/vf2.py          work/target.gml work/patterns.txt --workers 4 --output work/vf2.jsonl
```

**A published matcher:** convert the same inputs to the tools' formats, then call its driver with target, query list, workers, output TSV and per-query limit in seconds (an optional sixth argument names a file for the matched node sets):

```bash
python benchmark/convert.py work/target.gml work/patterns.txt work/input   # target.{gfu,vft,daf}, q/NNNN.*, queries_*.txt
build/bin/ri_batch       work/input/target.gfu work/input/queries_gfu.txt 4 work/ri.tsv 300
build/bin/vf3_batch_VF3  work/input/target.vft work/input/queries_vft.txt 4 work/vf3.tsv 300
build/bin/vf3_batch_VF3L work/input/target.vft work/input/queries_vft.txt 4 work/vf3l.tsv 300
build/bin/daf_batch      work/input/target.daf work/input/queries_daf.txt 4 work/daf.tsv 300
build/bin/fastiso_batch  work/input/target.gfu work/input/queries_gfu.txt 4 work/fastiso.tsv 300
```

**All of them at once, then checked against each other** (`SETS=1` makes every method that can write its matches do so):

```bash
VF2=1 SETS=1 benchmark/run_all.sh work/target.gml work/patterns.txt results/example 4
python benchmark/summarize.py results/example
```

Output of that run:

```
== example
  host: x1000c7s1b0n0, AMD EPYC 7742 64-Core Processor, 4 workers, 2026-09-30T00:28:29-04:00
  4 workers, reference VF2: 26,580 occurrences
    exact: VF2 finished 40/40; hop 40/40; seeded 40/40; RI mappings 40/40; VF3 mappings 40/40; VF3L mappings 40/40; DAF mappings 40/40; FASTiso mappings 40/40
    node sets vs VF2, patterns identical: hop 40/40; seeded 40/40; RI 40/40; VF3 40/40; VF3L 40/40; FASTiso 40/40
    wall s: RI 0.03, VF3 0.04, VF3L 0.05, DAF 0.01, FASTiso 0.08, hop 0.28 (+0.06), seeded exact 0.15, VF2 0.24
```

**How to read it:**
- **`exact`:** for how many of the 40 patterns each method's count agrees. Our methods' node-set counts are compared with VF2's, and every tool's mapping count with `hop_filter.py`'s.
- **`node sets vs VF2`:** for how many patterns each method found exactly VF2's node sets, set by set. Here all of them did. If a method differs, the line also gives how many sets it found that VF2 lacks (false) and how many it missed.
- **`wall s`:** seconds for all 40 queries. The toy example is too small for these to mean anything, and runs that write matches are slower anyway; time your own graphs in a run without `SETS=1`. `python benchmark/tables.py results/*` turns several runs into a markdown table of wall times.

## Comparing the methods on your own graph

```bash
benchmark/run_all.sh <target.gml> <patterns> <out_dir> <workers>
```

`<patterns>` is a single argument in one of three forms:
- a folder, meaning every `*.gml` directly inside it, in file-name order;
- a list file (any file not ending in `.gml`, e.g. `patterns.txt`) with one pattern GML path per line, relative to the list's own folder;
- one `.gml` file.

Results are keyed by pattern file name, so the names must be unique.

This runs every method on one machine so the times compare: the five tools, `hop_filter.py`, `seeded_exact.py`, and `vf2.py` if `VF2=1`.
- **Queries:** every query, with a per-query time limit (`LIMIT`, default 300 s).
- **Other settings:** `BIN` sets where the drivers are, and `PYTHON` which interpreter runs our methods. `WHICH=tools` runs only the five published tools, and `WHICH=ours` only our methods and VF2.
- **Checking that all methods find the same matches:** run it with `SETS=1`. Our methods then write their node sets, and RI, VF3, VF3L and FASTiso write theirs, and `summarize.py` compares them node set by node set with VF2 (or with `hop_filter.py` where VF2 hit its limit). DAF is checked by its counts. These runs write large files and are slower, so use a separate run for timings.

**What is timed:**
- **The tools:** they prepare the target and open their output files before the clock starts, so they are timed from the first query to the last.
- **`hop_filter.py` and seeded exact:** timed by their search alone. Hop windows' setup is shown separately in brackets.
- **`vf2.py`:** its wall time. Queries stopped by `LIMIT` count at the limit, so that time is a lower bound (marked `*`).
- **In every case:** reading the files and starting the workers are not counted.

`python tests/test_methods.py` repeats the checks on the toy example: our methods against VF2 node set by node set and, if `build/bin` exists, every tool's counts, and node sets for the tools that write them. It also runs under `pytest`.

## Notes

- **What the tools output:** all five count mappings, as their own programs do. Our methods count node sets, so their times include deduplicating mappings into node sets and the tools' times do not.
- **Where a driver differs from the authors' program.** These are the only differences:
  - **DAF, pattern labels the target lacks:** the only change to what a tool computes. Such a query has no match. The authors' program crashes or reports matches that do not exist when a query label is larger than every label of the data graph, because its label table (`DataGraph::GetTransferredLabel`) is read past its end. So `daf_batch` answers such a query 0 without calling DAF; the comparison is inside the query's time. For a missing label below the largest, DAF's own check catches it: `QueryGraph::LoadAndProcessGraph` returns false, and the driver records 0 (`main.cc` ignores that return value).
  - **DAF, crash when no match exists:** when the candidate space cannot be built, destroying DAF's partly built objects crashes (`free(): invalid pointer`). The driver records 0 and leaves them alone.
  - These were found by running the authors' programs themselves on queries with a label the target lacks: DAF's crashed or reported 229 and 650 matches where none exist.
  - **Other tools and missing labels:** RI, VF3, VF3L and FASTiso get no such help. Run as the authors ship them, each already answers these queries correctly with 0.
  - **VF3 and VF3L:** each query's steps from `main.cpp` run once, as with `-r 0`; the program's repeat loop, which re-runs a query to average its time, is not used.
  - **Match output (verification runs only):** the matches come through each tool's own interface, as described above.
- **Labels:** labels are recoded as 0..L−1 over the target and the patterns together (`methods/common.py`), so every tool sees the same labels as our methods. When most queries use a label the target lacks, as on some real graphs, the times mainly measure how fast each method rejects them. For DAF that is the driver's guard rather than DAF itself, so on such targets DAF's time is not comparable with the others'.
- **Licences:** this repository is MIT. The published tools keep their own licences (`benchmark/tools.lock`): RI is MIT, vf3lib is LGPL-3.0, DAF is Apache-2.0, and FASTiso's repository has no licence file. They are fetched from their authors, not redistributed.

## Layout

```
methods/            hop_filter.py, seeded_exact.py, vf2.py, common.py
benchmark/          setup_tools.sh, tools.lock, drivers/*.cpp, convert.py, run_all.sh, summarize.py, tables.py
examples/           make_example.py
tests/              test_methods.py
```

## Citation

If you use this repository, please cite it as:

```bibtex
@software{zhong_subgraph_matching,
  author = {Zhong, Yu},
  title  = {Benchmarking Exact Labelled Subgraph Matching},
  year   = {2026},
  url    = {https://github.com/yzhong36/subgraph_matching}
}
```
