#!/bin/bash
# Builds the response cache of the model named by the first argument, detached,
# on all hardware threads, then runs the qvarma_growth_states test, which
# compares that model's cache with responses recomputed from their definition.
# Capped at 12 GB of address space, since no systemd user session is reachable
# to cap resident memory; the cache file is mapped, so its 4.2 GB count towards
# that cap.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
MODEL=${1:?give a model label}
LOG=out/montecarlo/response_cache_build/run_response_cache_build_$MODEL.log
mkdir -p out/montecarlo/response_cache_build
echo "started $(date '+%F %T')" >> "$LOG"
make bin/response_cache_build bin/qvarma_growth_states_correctness >> "$LOG" 2>&1 || exit 1
( ulimit -v 12000000
  OMP_NUM_THREADS=$(nproc) ./bin/response_cache_build "$MODEL" ) >> "$LOG" 2>&1
STATUS=$?
echo "build exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
[ $STATUS -eq 0 ] || exit $STATUS
./bin/qvarma_growth_states_correctness >> "$LOG" 2>&1
echo "test exit status $? ended $(date '+%F %T')" >> "$LOG"
