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
[ -f "$I/queries_daf.txt" ] || "$PY" "$HERE/convert.py" "$G" "$P" "$I"
echo "$(hostname), $(lscpu | grep 'Model name' | sed 's/.*: *//'), $W workers, $(date -Is)" | tee "$R/w${W}_host.txt"
step() { echo "$(date +%T) w$W $1  load $(cut -d' ' -f1-3 /proc/loadavg)"; }

step RI;      "$BIN/ri_batch"       "$I/target.gfu" "$I/queries_gfu.txt" "$W" "$R/w${W}_ri.tsv" "$LIMIT"      > "$R/w${W}_ri.wall"
step VF3;     "$BIN/vf3_batch_VF3"  "$I/target.vft" "$I/queries_vft.txt" "$W" "$R/w${W}_vf3.tsv" "$LIMIT"     > "$R/w${W}_vf3.wall"
step VF3L;    "$BIN/vf3_batch_VF3L" "$I/target.vft" "$I/queries_vft.txt" "$W" "$R/w${W}_vf3l.tsv" "$LIMIT"    > "$R/w${W}_vf3l.wall"
step DAF;     "$BIN/daf_batch"      "$I/target.daf" "$I/queries_daf.txt" "$W" "$R/w${W}_daf.tsv" "$LIMIT"     > "$R/w${W}_daf.wall"
step FASTiso; "$BIN/fastiso_batch"  "$I/target.gfu" "$I/queries_gfu.txt" "$W" "$R/w${W}_fastiso.tsv" "$LIMIT" > "$R/w${W}_fastiso.wall"

step hop;     "$PY" "$M/hop_filter.py" "$G" "$P" --workers "$W" --output "$R/w${W}_hop.jsonl" > "$R/w${W}_hop.json"
step seeded;  "$PY" "$M/seeded_exact.py" "$G" "$P" --workers "$W" --output "$R/w${W}_seeded.jsonl" > "$R/w${W}_seeded.json"
if [ "${VF2:-0}" = 1 ]; then
    step VF2; "$PY" "$M/vf2.py" "$G" "$P" --workers "$W" --limit "$LIMIT" --output "$R/w${W}_vf2.jsonl" > "$R/w${W}_vf2.json"
fi
echo "$(date +%T) done"
