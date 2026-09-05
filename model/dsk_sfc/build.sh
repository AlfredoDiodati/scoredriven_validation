#!/bin/bash
# Builds the DSK simulator without cmake, which is what the measurements in
# docs/DSK_MODEL_CHANGES.md were taken from.
#
#   ./build.sh                      the model this project runs, optimised
#   ./build.sh --upstream OUTPUT    an unmodified binary at the flags upstream
#                                   ships, for tests/dsk_build_equivalence.c
#                                   to compare against
#   ./build.sh --sanitize OUTPUT    the same code under AddressSanitizer and
#                                   UndefinedBehaviorSanitizer, for
#                                   tests/dsk_memory_safety.c
#
# --upstream builds in a scratch copy of this tree with upstream/'s versions of
# every changed file copied over its counterpart, so the reference
# binary is upstream's code at upstream's flags and nothing here is disturbed.
#
# --sanitize builds this project's own code, not upstream's. The speed work
# replaced array subscripts with pointer arithmetic in several hot loops, and an
# index that runs off the end of an array reads whatever is next to it: on the
# seeds the equivalence tests happen to cover that could read a plausible number
# and produce output that compares equal. Comparing outputs cannot find that
# class of mistake, so the sanitizers are asked instead.
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

if [ "$1" = "--sanitize" ]; then
    [ -n "$2" ] || { echo "usage: $0 --sanitize OUTPUT" >&2; exit 1; }
    case "$2" in
        /*) output="$2" ;;
        *)  output="$caller_dir/$2" ;;
    esac
    # -O1 rather than -O2: the sanitizers want frame pointers and readable
    # stacks, and this build is asked about correctness rather than speed.
    flags="-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1 -fno-math-errno -msse"
fi

obj="$(mktemp -d)"
trap 'rm -rf "$obj" ${source_dir:+$([ "$source_dir" != "$PWD" ] && echo "$source_dir")}' EXIT

cd "$source_dir"
export OBJDIR="$obj" FLAGS="$flags"
ls newmat10/*.cpp auxiliary/*.cpp modules/*.cpp modules/WITCH_input/*.cpp dsk_sfc_main.cpp \
  | xargs -P "$(nproc)" -I@ bash -c 'g++ -std=c++11 $FLAGS -I. -c "@" -o "$OBJDIR/$(echo "@" | tr / _ | sed "s/\.cpp$/.o/")" 2>/dev/null'
g++ -std=c++11 $flags "$obj"/*.o -o "$output"

echo "built $output"
