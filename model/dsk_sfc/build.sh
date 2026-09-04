#!/bin/bash
# Builds the DSK simulator without cmake, which is what the measurements in
# docs/DSK_MODEL_CHANGES.md were taken from.
#
#   ./build.sh                      the model this project runs, optimised
#   ./build.sh --upstream OUTPUT    an unmodified binary at the flags upstream
#                                   ships, for tests/dsk_build_equivalence.c
#                                   to compare against
#
# --upstream builds in a scratch copy of this tree with upstream/'s versions of
# every changed file copied over its counterpart, so the reference
# binary is upstream's code at upstream's flags and nothing here is disturbed.
set -e
# Resolve the caller's output path before moving into the source directory.
caller_dir="$PWD"
cd "$(dirname "$0")"

source_dir="$PWD"
output="$PWD/dsk_SFC"
flags="-O2 -flto=$(nproc) -fno-math-errno -msse"

if [ "$1" = "--upstream" ]; then
    [ -n "$2" ] || { echo "usage: $0 --upstream OUTPUT" >&2; exit 1; }
    case "$2" in
        /*) output="$2" ;;
        *)  output="$caller_dir/$2" ;;
    esac
    # What CMAKE_BUILD_TYPE Debug and the project's own two flags amount to:
    # no optimisation at all.
    flags="-g -fno-math-errno -msse"
    source_dir="$(mktemp -d)"
    trap 'rm -rf "$source_dir"' EXIT
    cp -a "$PWD/." "$source_dir/"
    cp "$PWD/upstream/dsk_sfc_main.cpp" "$PWD/upstream/dsk_sfc_globalvars.h" "$source_dir/"
    cp "$PWD/upstream/modules/module_finance_sfc.cpp" "$PWD/upstream/modules/module_finance_sfc.h" "$source_dir/modules/"
fi

obj="$(mktemp -d)"
trap 'rm -rf "$obj" ${source_dir:+$([ "$source_dir" != "$PWD" ] && echo "$source_dir")}' EXIT

cd "$source_dir"
export OBJDIR="$obj" FLAGS="$flags"
ls newmat10/*.cpp auxiliary/*.cpp modules/*.cpp modules/WITCH_input/*.cpp dsk_sfc_main.cpp \
  | xargs -P "$(nproc)" -I@ bash -c 'g++ -std=c++11 $FLAGS -I. -c "@" -o "$OBJDIR/$(echo "@" | tr / _ | sed "s/\.cpp$/.o/")" 2>/dev/null'
g++ -std=c++11 $flags "$obj"/*.o -o "$output"

echo "built $output"
