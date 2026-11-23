#!/bin/bash
# Trial of the four true configuration recovery grids on the nine response
# vectors of the three other state weights (models 5 to 13 of
# montecarlo/response_cache.h), before the full runs.
#
# 1. Seri corrected, run 0 of every configuration under the nine models: 9000
#    benchmarks. It builds the nine response caches in
#    out/montecarlo/sweep_grid_response_cache/ on its way.
# 2. Seri uncorrected, the same 9000 benchmarks.
# 3. Bootstrap corrected and uncorrected, run 0 of cop_0001 to cop_0064 under
#    the nine models: 576 benchmarks each, for the time per benchmark.
# 4. Seri corrected under lp_nl of lpirfs' weight, run 0 of cop_0001 to
#    cop_0064, to compare with the rows the finished grid holds.
#
# Each test build writes its rows to its own progress file in
# out/montecarlo/state_weight_trial/ and no result. Capped at 10 GB of address
# space, since no systemd user session is reachable to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
DIR=out/montecarlo/state_weight_trial
mkdir -p "$DIR"
LOG=$DIR/trial.log
FLAGS="-O3 -march=native -Wall -Wextra -std=c11 -DMAT_DOUBLE -fopenmp -I/usr/local/include/et_al. -I/usr/include/openblas -I."
LIBS="-lm -fopenmp -lopenblas"
STATUS=0

build_and_run() {
  local name=$1 source=$2
  shift 2
  echo "$name build $(date '+%F %T')" >> "$LOG"
  cc $FLAGS "$@" -DPROGRESS_PATH="\"$DIR/$name.csv\"" -DRESULT_PATH="\"$DIR/$name.csv.gz\"" \
     "$source" -o "bin/state_weight_trial_$name" $LIBS >> "$LOG" 2>&1 || return 1
  echo "$name started $(date '+%F %T')" >> "$LOG"
  ( ulimit -v 10000000
    systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="state_weight_trial" \
                    --why="true configuration recovery trial" --mode=block "./bin/state_weight_trial_$name" ) \
    >> "$LOG" 2>&1
  local status=$?
  echo "$name exit status $status ended $(date '+%F %T')" >> "$LOG"
  return $status
}

echo "started $(date '+%F %T')" >> "$LOG"
build_and_run seri_corrected montecarlo/seri_true_configuration_recovery.c \
  -DFIRST_MODEL=5 -DLAST_MODEL=13 -DRUN_LIMIT=1 \
&& build_and_run seri_uncorrected montecarlo/seri_uncorrected_true_configuration_recovery.c \
  -DFIRST_MODEL=5 -DLAST_MODEL=13 -DRUN_LIMIT=1 \
&& build_and_run bootstrap_corrected montecarlo/true_configuration_recovery.c \
  -DFIRST_MODEL=5 -DLAST_MODEL=13 -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64 \
&& build_and_run bootstrap_uncorrected montecarlo/uncorrected_true_configuration_recovery.c \
  -DFIRST_MODEL=5 -DLAST_MODEL=13 -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64 \
&& build_and_run seri_corrected_lpirfs_nl montecarlo/seri_true_configuration_recovery.c \
  -DFIRST_MODEL=4 -DLAST_MODEL=4 -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64
STATUS=$?
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
