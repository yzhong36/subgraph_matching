#!/bin/bash
# Every method on one target, one worker count, on this machine (so the times compare): the published exact matchers
# (RI, VF3, VF3L, DAF, FASTiso; built by setup_tools.sh) and ours (hop_filter.py, seeded_exact.py),
# optionally plain whole-graph VF2 (VF2=1).  Every query; LIMIT seconds per query (default 300).
#
#   benchmark/run_all.sh <target.gml> <patterns> <out_dir> <workers>
#
# patterns: a folder of pattern .gml files, a list file (one pattern GML path per line, relative to the list's folder),
# or one .gml file.  Pattern file names must be unique (results are keyed by them).
# Env: BIN (default <repo>/build/bin), PYTHON (default python3), LIMIT, VF2.
# SETS=1: verification run, for benchmark/summarize.py's check that every method finds the same node sets.  Our
# methods write their node sets and the tools (except DAF, which only counts) dump theirs, so times are not comparable.
# Results in out_dir: input/ (converted files), w<W>_<method>.{tsv,wall} for the tools, w<W>_<method>.{jsonl,json} for
# ours, w<W>_host.txt.  benchmark/summarize.py reads them.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
M=$HERE/../methods
BIN=${BIN:-$HERE/../build/bin}
PY=${PYTHON:-python3}
LIMIT=${LIMIT:-300}
G=$1; P=$2; R=$3; W=$4
mkdir -p "$R"
I=$R/input
[ -f "$I/ids.txt" ] || "$PY" "$HERE/convert.py" "$G" "$P" "$I"
echo "$(hostname), $(lscpu | grep 'Model name' | sed 's/.*: *//'), $W workers, $(date -Is)" | tee "$R/w${W}_host.txt"
step() { echo "$(date +%T) w$W $1  load $(cut -d' ' -f1-3 /proc/loadavg)"; }

tool() {   # tool <label> <key> <driver> <format> <dumps matches 1|0>
    step "$1"
    local extra=()
    [ "${SETS:-0}" = 1 ] && [ "$5" = 1 ] && extra=("$R/w${W}_$2.matches")
    "$BIN/$3" "$I/target.$4" "$I/queries_$4.txt" "$W" "$R/w${W}_$2.tsv" "$LIMIT" ${extra[@]+"${extra[@]}"} > "$R/w${W}_$2.wall"
}
tool RI      ri      ri_batch       gfu 1
tool VF3     vf3     vf3_batch_VF3  vft 1
tool VF3L    vf3l    vf3_batch_VF3L vft 1
tool DAF     daf     daf_batch      daf 0
tool FASTiso fastiso fastiso_batch  gfu 1

SETARG=(); [ "${SETS:-0}" = 1 ] && SETARG=(--sets)
step hop;     "$PY" "$M/hop_filter.py" "$G" "$P" --workers "$W" ${SETARG[@]+"${SETARG[@]}"} --output "$R/w${W}_hop.jsonl" > "$R/w${W}_hop.json"
step seeded;  "$PY" "$M/seeded_exact.py" "$G" "$P" --workers "$W" ${SETARG[@]+"${SETARG[@]}"} --output "$R/w${W}_seeded.jsonl" > "$R/w${W}_seeded.json"
if [ "${VF2:-0}" = 1 ]; then
    step VF2; "$PY" "$M/vf2.py" "$G" "$P" --workers "$W" --limit "$LIMIT" --output "$R/w${W}_vf2.jsonl" > "$R/w${W}_vf2.json"
fi
echo "$(date +%T) done"
