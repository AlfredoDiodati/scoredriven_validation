#!/bin/bash
# Test builds of the two Seri grids over qvarma_growth_states only (model 14 of
# montecarlo/response_cache.h), run 0 of the first 128 configurations, writing
# into this folder and never touching the grids' own progress files or
# results. Measures the time per benchmark before the full run.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
DIR=out/montecarlo/qvarma_growth_states_trial
LOG=$DIR/run_qvarma_growth_states_trial.log
CFLAGS="-O3 -march=native -Wall -Wextra -std=c11 -DMAT_DOUBLE -fopenmp $(pkg-config --cflags et_al.-core) -I."
LIBS="$(pkg-config --libs et_al.-core)"
mkdir -p bin
echo "started $(date '+%F %T')" >> "$LOG"
for GRID in seri_true_configuration_recovery seri_uncorrected_true_configuration_recovery; do
  cc $CFLAGS -DFIRST_MODEL=14 -DLAST_MODEL=14 -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=128 \
     -DPROGRESS_PATH="\"$DIR/${GRID}_progress.csv\"" -DRESULT_PATH="\"$DIR/${GRID}.csv.gz\"" \
     montecarlo/$GRID.c -o bin/${GRID}_qvarma_growth_states_trial $LIBS >> "$LOG" 2>&1 || exit 1
  echo "$GRID started $(date '+%F %T')" >> "$LOG"
  ( ulimit -v 12000000
    ./bin/${GRID}_qvarma_growth_states_trial ) >> "$LOG" 2>&1
  echo "$GRID exit status $? ended $(date '+%F %T')" >> "$LOG"
done
