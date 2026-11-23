#!/bin/bash
# Once the Seri grids launched by out/montecarlo/run_seri_state_weight_grids_detached.sh
# have finished, moves the rows of state 1 and state 2 on their own out of the
# four grids' results and into out/old/montecarlo_single_states/
# (bin/single_state_rows_split), then checks the grids of montecarlo/ and their
# single-state copies against the parted results and rewrites every report.
#
# The check: test builds over run 0 of cop_0001 to cop_0064 under models 0 to 4
# (the t-QVARMA and the four projections of lpirfs' weight) of the Seri
# corrected grid and of the bootstrap corrected grid, each in its montecarlo/
# version and its single-state copy. The montecarlo/ build must write exactly
# the rows of the t-QVARMA, the linear projection and both states that the
# parted main result holds for those benchmarks, the copy exactly the rows of
# state 1 and state 2 that the parted old result holds, byte for byte after
# sorting. Rows and logs go to out/old/montecarlo_single_states/removal_check/.
#
# Capped at 7 GB of address space, since no systemd user session is reachable
# to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
DIR=out/old/montecarlo_single_states
CHECK=$DIR/removal_check
LOG=$DIR/run_single_state_split.log
mkdir -p "$CHECK"
echo "started $(date '+%F %T'), waiting for the Seri grids" >> "$LOG"
until grep -q "^exit status" out/montecarlo/seri_state_weight_grids.log 2>/dev/null; do sleep 60; done
echo "Seri grids: $(tail -1 out/montecarlo/seri_state_weight_grids.log)" >> "$LOG"
grep -q "^exit status 0" out/montecarlo/seri_state_weight_grids.log || { echo "the Seri grids failed, nothing done" >> "$LOG"; exit 1; }

( ulimit -v 7000000; ./bin/single_state_rows_split ) >> "$LOG" 2>&1 || { echo "split failed, see $DIR/single_state_rows_split.log" >> "$LOG"; exit 1; }
cat "$DIR/single_state_rows_split.log" >> "$LOG"

FLAGS="-O3 -march=native -Wall -Wextra -std=c11 -DMAT_DOUBLE -fopenmp -I/usr/local/include/et_al. -I/usr/include/openblas -I."
LIBS="-lm -fopenmp -lopenblas"
STATUS=0
check_build() {
  local name=$1 source=$2 result=$3 expected_models=$4
  cc $FLAGS -DFIRST_MODEL=0 -DLAST_MODEL=4 -DRUN_LIMIT=1 -DBENCHMARK_LIMIT=64 \
     -DPROGRESS_PATH="\"$CHECK/$name.csv\"" -DRESULT_PATH="\"$CHECK/$name.csv.gz\"" \
     "$source" -o "bin/removal_check_$name" $LIBS >> "$CHECK/$name.log" 2>&1 || return 1
  ( ulimit -v 7000000; "./bin/removal_check_$name" ) >> "$CHECK/$name.log" 2>&1 || return 1
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
check_build seri_main montecarlo/seri_true_configuration_recovery.c \
  out/montecarlo/seri_true_configuration_recovery.csv.gz "lp_lin lp_nl qvarma " || STATUS=1
check_build seri_single_state old/montecarlo_single_states/single_state_seri_true_configuration_recovery.c \
  $DIR/seri_true_configuration_recovery.csv.gz "lp_s1 lp_s2 " || STATUS=1
check_build bootstrap_main montecarlo/true_configuration_recovery.c \
  out/montecarlo/true_configuration_recovery.csv.gz "lp_lin lp_nl qvarma " || STATUS=1
check_build bootstrap_single_state old/montecarlo_single_states/single_state_true_configuration_recovery.c \
  $DIR/true_configuration_recovery.csv.gz "lp_s1 lp_s2 " || STATUS=1

make bin/seri_uncorrected_true_configuration_recovery >> "$LOG" 2>&1 || STATUS=1
for PROGRAM in seri_true_configuration_recovery_report seri_against_naive_grid seri_against_naive_by_configuration \
               response_variance_by_configuration seri_corrected_against_uncorrected_grid \
               corrected_against_uncorrected_grid seri_against_bootstrap_grid identification_table; do
  for BUILD in "$PROGRAM" "single_state_$PROGRAM"; do
    make "bin/$BUILD" >> "$LOG" 2>&1 || { STATUS=1; continue; }
    ( ulimit -v 7000000; "./bin/$BUILD" ) >> "$LOG" 2>&1
    echo "$BUILD exit status $? $(date '+%F %T')" >> "$LOG"
  done
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
