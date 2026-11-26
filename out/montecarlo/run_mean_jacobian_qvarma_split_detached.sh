#!/bin/bash
# Moves the rows of the t-QVARMA at the sample mean of D_t (model qvarma), which
# the grids of montecarlo/ no longer run, out of the four grids' results and into
# out/old/montecarlo_mean_jacobian_qvarma/ (bin/mean_jacobian_qvarma_rows_split).
# Then checks the grids against the parted results and rewrites every report.
#
# The check: test builds over run 0 of cop_0001 to cop_0064 of the Seri
# corrected grid (models 0 to 14) and of the bootstrap corrected grid (models 0
# to 4: the bootstrap grids have not run the growth-state t-QVARMA or both
# states under the causal HP cycle at the estimated lambda), each in its
# montecarlo/ version and its old copy. Each build must write rows of exactly
# the models its version runs, equal byte for byte, after sorting, to the rows
# the matching parted result holds for those benchmarks.
#
# Outputs in out/montecarlo/mean_jacobian_qvarma_split/. Capped at 14 GB of
# address space, since no systemd user session is reachable to cap resident
# memory; the mapped response caches count towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
export LC_ALL=C
CHECK=out/montecarlo/mean_jacobian_qvarma_split
LOG=$CHECK/run_mean_jacobian_qvarma_split.log
mkdir -p "$CHECK"
echo "started $(date '+%F %T')" >> "$LOG"
# The four results as they were before the split, kept whole.
BEFORE=out/montecarlo/before_mean_jacobian_qvarma_split
mkdir -p "$BEFORE"
for GRID in true_configuration_recovery uncorrected_true_configuration_recovery \
            seri_true_configuration_recovery seri_uncorrected_true_configuration_recovery; do
  cp -n "out/montecarlo/$GRID.csv.gz" "$BEFORE/$GRID.csv.gz" || exit 1
done
make bin/mean_jacobian_qvarma_rows_split >> "$LOG" 2>&1 || exit 1
( ulimit -v 14000000; ./bin/mean_jacobian_qvarma_rows_split ) >> "$LOG" 2>&1
STATUS=$?
echo "mean_jacobian_qvarma_rows_split exit status $STATUS $(date '+%F %T')" >> "$LOG"
[ $STATUS -ne 0 ] && { echo "exit status 1 ended $(date '+%F %T')" >> "$LOG"; exit 1; }
cat out/old/montecarlo_mean_jacobian_qvarma/mean_jacobian_qvarma_rows_split.log >> "$LOG"

FLAGS="-O3 -march=native -Wall -Wextra -std=c11 -DMAT_DOUBLE -fopenmp $(pkg-config --cflags et_al.-core) -I."
LIBS="$(pkg-config --libs et_al.-core)"
STATUS=0
check_build() {
  local name=$1 source=$2 last_model=$3 result=$4 expected_models=$5
  cc $FLAGS -DFIRST_MODEL=0 -DLAST_MODEL=$last_model -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64 \
     -DPROGRESS_PATH="\"$CHECK/$name.csv\"" -DRESULT_PATH="\"$CHECK/$name.csv.gz\"" \
     "$source" -o "bin/mean_jacobian_qvarma_split_$name" $LIBS >> "$CHECK/$name.log" 2>&1 || return 1
  ( ulimit -v 14000000; "./bin/mean_jacobian_qvarma_split_$name" ) >> "$CHECK/$name.log" 2>&1 || return 1
  tail -n +2 "$CHECK/$name.csv" | sort > "$CHECK/$name.sorted.csv"
  local models
  models=$(cut -d, -f1 "$CHECK/$name.sorted.csv" | sort -u | tr '\n' ' ')
  if [ "$models" != "$expected_models" ]; then
    echo "$name: FAILED, wrote rows of models $models, expected $expected_models" >> "$LOG"
    return 1
  fi
  zcat "$result" | awk -F, -v models="$models" 'BEGIN { n = split(models, m, " "); for (i = 1; i <= n; i++) keep[m[i]] = 1 }
    NR > 1 && ($1 in keep) && $3 == 0 && $2 >= 1 && $2 <= 64' | sort > "$CHECK/$name.expected.csv"
  if cmp -s "$CHECK/$name.sorted.csv" "$CHECK/$name.expected.csv"; then
    echo "$name: $(wc -l < "$CHECK/$name.sorted.csv") rows of models $models equal the result's byte for byte" >> "$LOG"
  else
    echo "$name: FAILED, rows of models $models differ from the result's" >> "$LOG"
    return 1
  fi
}
OLD=old/montecarlo_mean_jacobian_qvarma
check_build seri_main montecarlo/seri_true_configuration_recovery.c 14 \
  out/montecarlo/seri_true_configuration_recovery.csv.gz \
  "lp_lin lp_nl_causal_hp_estimated_lambda qvarma_growth_states " || STATUS=1
check_build seri_mean_jacobian_qvarma $OLD/mean_jacobian_qvarma_seri_true_configuration_recovery.c 14 \
  out/$OLD/seri_true_configuration_recovery.csv.gz "qvarma " || STATUS=1
check_build bootstrap_main montecarlo/true_configuration_recovery.c 4 \
  out/montecarlo/true_configuration_recovery.csv.gz "lp_lin " || STATUS=1
check_build bootstrap_mean_jacobian_qvarma $OLD/mean_jacobian_qvarma_true_configuration_recovery.c 4 \
  out/$OLD/true_configuration_recovery.csv.gz "qvarma " || STATUS=1

for PROGRAM in seri_true_configuration_recovery_report seri_against_naive_grid seri_against_naive_by_configuration \
               response_variance_by_configuration seri_corrected_against_uncorrected_grid \
               corrected_against_uncorrected_grid seri_against_bootstrap_grid identification_table; do
  for BUILD in "$PROGRAM" "mean_jacobian_qvarma_$PROGRAM"; do
    make "bin/$BUILD" >> "$LOG" 2>&1 || { STATUS=1; continue; }
    ( ulimit -v 14000000; "./bin/$BUILD" ) >> "$LOG" 2>&1
    echo "$BUILD exit status $? $(date '+%F %T')" >> "$LOG"
  done
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
