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
# -ffp-contract=off so that byte equality with upstream does not rest on the
# instruction set alone. Contraction of a multiply and an add into one FMA
# changes results, and -msse does not enable FMA, so today nothing is contracted
# either way; asking for it explicitly is what keeps that true if the flags ever
# gain -march=native. docs/DSK_MODEL_CHANGES.md records the 356 KB of output
# -march=native moved when it was tried.
# -Wall on this project's own headers. The vendored newmat10 and rapidjson are
# upstream's and are compiled as they are, so the warning flag goes on the
# translation units this project wrote or rewrote rather than on all of them.
common="-fno-math-errno -msse -ffp-contract=off"
flags="-O2 -flto=$(nproc) $common"
# The scratch tree --upstream builds in, removed by the one trap installed
# below. Empty for the other two modes.
scratch_source=""

if [ "$1" = "--upstream" ]; then
    [ -n "$2" ] || { echo "usage: $0 --upstream OUTPUT" >&2; exit 1; }
    case "$2" in
        /*) output="$2" ;;
        *)  output="$caller_dir/$2" ;;
    esac
    # What CMAKE_BUILD_TYPE Debug and the project's own two flags amount to:
    # no optimisation at all.
    flags="-g $common"
    scratch_source="$(mktemp -d)"
    source_dir="$scratch_source"
    cp -a "$PWD/." "$source_dir/"
    cp "$PWD/upstream/dsk_sfc_main.cpp" "$PWD/upstream/dsk_sfc_globalvars.h" "$source_dir/"
    cp "$PWD/upstream/modules/module_finance_sfc.cpp" "$PWD/upstream/modules/module_finance_sfc.h" \
       "$PWD/upstream/modules/module_climate_sfc.cpp" "$PWD/upstream/modules/module_macro_sfc.cpp" \
       "$source_dir/modules/"
fi

if [ "$1" = "--sanitize" ]; then
    [ -n "$2" ] || { echo "usage: $0 --sanitize OUTPUT" >&2; exit 1; }
    case "$2" in
        /*) output="$2" ;;
        *)  output="$caller_dir/$2" ;;
    esac
    # -O1 rather than -O2: the sanitizers want frame pointers and readable
    # stacks, and this build is asked about correctness rather than speed.
    flags="-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1 $common"
fi

obj="$(mktemp -d)"
trap 'rm -rf "$obj" ${scratch_source:+"$scratch_source"}' EXIT

cd "$source_dir"
export OBJDIR="$obj" FLAGS="$flags"
# dsk_sfc_main.cpp is where this project's own headers are included, so it is
# the unit -Wall is asked about; newmat10 and rapidjson are upstream's and are
# compiled as they are. It is also the longest single compile, so it starts
# first and runs alongside the rest rather than after them.
#
# Compiler output is not discarded. A build that fails has to say which file and
# why, or the link reports an undefined reference with nothing to trace it to.
g++ -std=c++11 $flags -Wall -I. -c dsk_sfc_main.cpp -o "$obj/dsk_sfc_main.o" &
main_unit=$!
ls newmat10/*.cpp auxiliary/*.cpp modules/*.cpp modules/WITCH_input/*.cpp \
  | xargs -P "$(nproc)" -I@ bash -c 'g++ -std=c++11 $FLAGS -I. -c "@" -o "$OBJDIR/$(echo "@" | tr / _ | sed "s/\.cpp$/.o/")"'
wait "$main_unit"
g++ -std=c++11 $flags "$obj"/*.o -o "$output"

echo "built $output"
