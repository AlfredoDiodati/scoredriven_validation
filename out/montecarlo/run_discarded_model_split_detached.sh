#!/bin/bash
# Moves the rows of the models the grids of montecarlo/ no longer run out of the
# four grids' results: state 1 and state 2 on their own to
# out/old/montecarlo_single_states/ (bin/single_state_rows_split), both states
# under the three state weights other than the causal HP cycle at the estimated
# lambda to out/old/montecarlo_other_state_weights/
# (bin/other_state_weights_rows_split). Then checks the grids against the
# parted results and rewrites every report.
#
# The check: test builds over run 0 of cop_0001 to cop_0064 of the Seri
# corrected grid (models 0 to 13) and of the bootstrap corrected grid (models 0
# to 4, the only ones it has run), each in its montecarlo/ version and its two
# old copies. Each build must write rows of exactly the models its version
# runs, equal byte for byte, after sorting, to the rows the matching parted
# result holds for those benchmarks.
#
# Outputs in out/montecarlo/discarded_model_split/. Capped at 7 GB of address
# space, since no systemd user session is reachable to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
export LC_ALL=C
CHECK=out/montecarlo/discarded_model_split
LOG=$CHECK/run_discarded_model_split.log
mkdir -p "$CHECK"
echo "started $(date '+%F %T')" >> "$LOG"
make bin/single_state_rows_split bin/other_state_weights_rows_split >> "$LOG" 2>&1 || exit 1
for SPLIT in single_state_rows_split other_state_weights_rows_split; do
  ( ulimit -v 7000000; "./bin/$SPLIT" ) >> "$LOG" 2>&1
  STATUS=$?
  echo "$SPLIT exit status $STATUS $(date '+%F %T')" >> "$LOG"
  [ $STATUS -ne 0 ] && { echo "exit status 1 ended $(date '+%F %T')" >> "$LOG"; exit 1; }
done
cat out/old/montecarlo_single_states/single_state_rows_split.log \
    out/old/montecarlo_other_state_weights/other_state_weights_rows_split.log >> "$LOG"

FLAGS="-O3 -march=native -Wall -Wextra -std=c11 -DMAT_DOUBLE -fopenmp -I/usr/local/include/et_al. -I/usr/include/openblas -I."
LIBS="-lm -fopenmp -lopenblas"
STATUS=0
check_build() {
  local name=$1 source=$2 last_model=$3 result=$4 expected_models=$5
  cc $FLAGS -DFIRST_MODEL=0 -DLAST_MODEL=$last_model -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64 \
     -DPROGRESS_PATH="\"$CHECK/$name.csv\"" -DRESULT_PATH="\"$CHECK/$name.csv.gz\"" \
     "$source" -o "bin/discarded_model_split_$name" $LIBS >> "$CHECK/$name.log" 2>&1 || return 1
  ( ulimit -v 7000000; "./bin/discarded_model_split_$name" ) >> "$CHECK/$name.log" 2>&1 || return 1
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
SINGLE=old/montecarlo_single_states
OTHER=old/montecarlo_other_state_weights
check_build seri_main montecarlo/seri_true_configuration_recovery.c 13 \
  out/montecarlo/seri_true_configuration_recovery.csv.gz \
  "lp_lin lp_nl_causal_hp_estimated_lambda qvarma " || STATUS=1
check_build seri_single_state $SINGLE/single_state_seri_true_configuration_recovery.c 13 \
  out/$SINGLE/seri_true_configuration_recovery.csv.gz \
  "lp_s1 lp_s1_causal_hp lp_s1_causal_hp_estimated_lambda lp_s1_estimated_lambda lp_s2 lp_s2_causal_hp lp_s2_causal_hp_estimated_lambda lp_s2_estimated_lambda " || STATUS=1
check_build seri_other_state_weights $OTHER/other_state_weights_seri_true_configuration_recovery.c 13 \
  out/$OTHER/seri_true_configuration_recovery.csv.gz \
  "lp_nl lp_nl_causal_hp lp_nl_estimated_lambda " || STATUS=1
check_build bootstrap_main montecarlo/true_configuration_recovery.c 4 \
  out/montecarlo/true_configuration_recovery.csv.gz "lp_lin qvarma " || STATUS=1
check_build bootstrap_single_state $SINGLE/single_state_true_configuration_recovery.c 4 \
  out/$SINGLE/true_configuration_recovery.csv.gz "lp_s1 lp_s2 " || STATUS=1
check_build bootstrap_other_state_weights $OTHER/other_state_weights_true_configuration_recovery.c 4 \
  out/$OTHER/true_configuration_recovery.csv.gz "lp_nl " || STATUS=1

for PROGRAM in seri_true_configuration_recovery_report seri_against_naive_grid seri_against_naive_by_configuration \
               response_variance_by_configuration seri_corrected_against_uncorrected_grid \
               corrected_against_uncorrected_grid seri_against_bootstrap_grid identification_table; do
  for BUILD in "$PROGRAM" "single_state_$PROGRAM" "other_state_weights_$PROGRAM"; do
    make "bin/$BUILD" >> "$LOG" 2>&1 || { STATUS=1; continue; }
    ( ulimit -v 7000000; "./bin/$BUILD" ) >> "$LOG" 2>&1
    echo "$BUILD exit status $? $(date '+%F %T')" >> "$LOG"
  done
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
