#!/bin/bash
# Moves the rows of the local projections fitted with inflation and the
# interest rate in first differences out of the two Seri results of
# out/montecarlo/ into out/old/montecarlo_differenced_lp/
# (bin/differenced_lp_rows_split, which checks that the two parts are a
# partition of each result). Then checks the rows left in each main result
# against the copy taken before the differenced models were added,
# out/montecarlo/before_differenced/: after sorting they must be equal byte for
# byte. Then rewrites the reports that read the Seri results and the
# comparison of the differenced projections with the stored ones.
#
# Capped at 7 GB of address space, since no systemd user session is reachable
# to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
export LC_ALL=C
DIR=out/old/montecarlo_differenced_lp
LOG=$DIR/run_differenced_lp_split.log
mkdir -p "$DIR"
export TMPDIR="$PWD/$DIR"
echo "started $(date '+%F %T')" >> "$LOG"
make bin/differenced_lp_rows_split bin/differenced_lp_against_stored_seri_report >> "$LOG" 2>&1 || exit 1
( ulimit -v 7000000; ./bin/differenced_lp_rows_split ) >> "$LOG" 2>&1
STATUS=$?
cat "$DIR/differenced_lp_rows_split.log" >> "$LOG"
echo "differenced_lp_rows_split exit status $STATUS $(date '+%F %T')" >> "$LOG"
[ $STATUS -eq 0 ] || { echo "exit status 1 ended $(date '+%F %T')" >> "$LOG"; exit 1; }

STATUS=0
for RESULT in seri_held_out_mean_true_configuration_recovery_100_mean_runs.csv.gz \
              seri_uncorrected_true_configuration_recovery.csv.gz; do
  zcat "out/montecarlo/$RESULT" | tail -n +2 | sort > "$DIR/check_after.csv"
  zcat "out/montecarlo/before_differenced/$RESULT" | tail -n +2 | sort > "$DIR/check_before.csv"
  if cmp -s "$DIR/check_after.csv" "$DIR/check_before.csv"; then
    echo "$RESULT: $(wc -l < "$DIR/check_after.csv") rows left equal the result before the differenced models byte for byte" >> "$LOG"
  else
    echo "$RESULT: FAILED, the rows left differ from the result before the differenced models" >> "$LOG"
    STATUS=1
  fi
  echo "$RESULT moved rows per model: $(zcat "$DIR/$RESULT" | tail -n +2 | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
done
rm -f "$DIR/check_after.csv" "$DIR/check_before.csv"

for PROGRAM in seri_true_configuration_recovery_report seri_against_naive_grid seri_against_naive_by_configuration \
               response_variance_by_configuration seri_corrected_against_uncorrected_grid \
               seri_against_bootstrap_grid identification_table differenced_lp_against_stored_seri_report; do
  make "bin/$PROGRAM" >> "$LOG" 2>&1 || { STATUS=1; continue; }
  ( ulimit -v 7000000; "./bin/$PROGRAM" ) >> "$LOG" 2>&1
  echo "$PROGRAM exit status $? $(date '+%F %T')" >> "$LOG"
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
