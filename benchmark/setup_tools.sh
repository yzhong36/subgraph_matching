#!/bin/bash
# Fetch the published exact matchers at the commits in tools.lock and build the batch drivers (drivers/*.cpp) that run
# them on many queries, with the flags the benchmark used.
#
#   BOOST_INCLUDE=/path/to/boost/include benchmark/setup_tools.sh [build_dir]
#
# build_dir defaults to <repo>/build.  Needs git, a C++17 compiler (CXX, default g++; the benchmark used GCC 14.2) and
# Boost headers for DAF (BOOST_INCLUDE, default /usr/include; the benchmark used Boost 1.86).  Writes
#   build/tools/<name>   the sources at the pinned commits
#   build/bin/           ri_batch, vf3_batch_VF3, vf3_batch_VF3L, daf_batch, fastiso_batch
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
B=$(mkdir -p "${1:-$HERE/../build}" && cd "${1:-$HERE/../build}" && pwd)
T=$B/tools
CXX=${CXX:-g++}
BOOST_INCLUDE=${BOOST_INCLUDE:-/usr/include}
mkdir -p "$T" "$B/bin"

fetch() {   # fetch <name> <url> <commit>
    local d=$T/$1
    if [ ! -d "$d/.git" ]; then git init -q "$d"; git -C "$d" remote add origin "$2"; fi
    if ! git -C "$d" cat-file -e "$3^{commit}" 2>/dev/null; then
        git -C "$d" fetch -q --depth 1 origin "$3" || git -C "$d" fetch -q origin
    fi
    git -C "$d" checkout -q --detach "$3"
    [ "$(git -C "$d" rev-parse HEAD)" = "$3" ] || { echo "$1: not at $3" >&2; exit 1; }
    echo "$1 at $3"
}
while read -r name url commit _rest; do
    case "$name" in ''|'#'*) continue ;; esac
    fetch "$name" "$url" "$commit"
done < "$HERE/tools.lock"

[ -f "$BOOST_INCLUDE/boost/dynamic_bitset.hpp" ] || { echo "Boost headers not found in $BOOST_INCLUDE (set BOOST_INCLUDE)" >&2; exit 1; }
D=$HERE/drivers
set -x
# RI 3.6 and vf3lib: header-only for the drivers
$CXX -O3 -std=c++17 -I "$T/RI/rilib" -I "$T/RI/include" "$D/ri_batch.cpp" -o "$B/bin/ri_batch"
for v in VF3 VF3L; do
    $CXX -std=c++14 -O3 -march=native -Wno-deprecated -D$v -I "$T/vf3lib/include" "$D/vf3_batch.cpp" -o "$B/bin/vf3_batch_$v"
done
# DAF: its sources compiled with the driver (-include cstdint for GCC >= 13)
$CXX -std=c++17 -O3 -march=native -DNDEBUG -include cstdint -I "$T/DAF" -I "$BOOST_INCLUDE" "$D/daf_batch.cpp" "$T"/DAF/src/*.cc \
    -o "$B/bin/daf_batch"
# FASTiso: the header its CMake project would configure (FASTiso 1.0), then its sources with the driver
mkdir -p "$T/fastiso/build"
sed -e 's/@PROJECT_NAME@/FASTiso/' -e 's/@FASTiso_VERSION_MAJOR@/1/' -e 's/@FASTiso_VERSION_MINOR@/0/' \
    "$T/fastiso/FASTiso_config.h.in" > "$T/fastiso/build/FASTiso_config.hpp"
$CXX -std=c++17 -O3 -march=native -I "$T/fastiso/header" -I "$T/fastiso/build" "$D/fastiso_batch.cpp" "$T"/fastiso/src/*.cpp \
    -o "$B/bin/fastiso_batch"
set +x
ls -la "$B/bin"
