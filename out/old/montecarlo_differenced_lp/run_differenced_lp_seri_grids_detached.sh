#!/bin/bash
# Ran from out/montecarlo/ on 2026-11-26 to 2026-11-28, when the Seri grids of
# montecarlo/ still ran the differenced models and the test was under tests/;
# their copies in old/montecarlo_differenced_lp/ run them now.
#
# The local projections with inflation and the interest rate in first
# differences, from the fits to both Seri grids, one step after the other:
#
# 1. the fits on every run: the linear projection (lpirfs' weight) and the
#    state-dependent one under the causal HP cycle at the estimated lambda;
# 2. the response caches of the differenced linear and both-states models;
# 3. tests/lp_system_differenced_correctness.c, which stops the chain on a
#    failure;
# 4. a copy of both Seri results as they are, into
#    out/montecarlo/before_differenced/;
# 5. the adjusted Seri grid (mean from the last 100 runs), then the unadjusted
#    one, each adding the two differenced models to its result.
#
# A rerun after an interruption redoes nothing finished: stored fits are kept,
# caches are read back, a grid resumes from its progress file.
#
# Capped at 12 GB of address space, since no systemd user session is reachable
# to cap resident memory; the mapped caches count towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/montecarlo/differenced_lp_seri_grids.log
step() { echo "$1 $(date '+%F %T')" >> "$LOG"; }
step "started"
make bin/abm_system_fit_lp bin/response_cache_build bin/lp_system_differenced_correctness \
     bin/seri_held_out_mean_true_configuration_recovery bin/seri_uncorrected_true_configuration_recovery \
     >> "$LOG" 2>&1 || exit 1

for FLAGS in "--difference-cointegrated" "--difference-cointegrated --causal-hp --estimate-lambda"; do
    step "fit $FLAGS started"
    ( ulimit -v 12000000; ./bin/abm_system_fit_lp $FLAGS ) >> "$LOG" 2>&1
    STATUS=$?
    step "fit $FLAGS exit status $STATUS"
    [ $STATUS -eq 0 ] || exit $STATUS
done

for MODEL in lp_lin_differenced lp_nl_differenced_causal_hp_estimated_lambda; do
    step "cache $MODEL started"
    ( ulimit -v 12000000; OMP_NUM_THREADS=$(nproc) ./bin/response_cache_build "$MODEL" ) >> "$LOG" 2>&1
    STATUS=$?
    step "cache $MODEL exit status $STATUS"
    [ $STATUS -eq 0 ] || exit $STATUS
done

./bin/lp_system_differenced_correctness >> "$LOG" 2>&1
STATUS=$?
step "test exit status $STATUS"
[ $STATUS -eq 0 ] || exit $STATUS

BEFORE=out/montecarlo/before_differenced
mkdir -p "$BEFORE"
for RESULT in seri_held_out_mean_true_configuration_recovery_100_mean_runs.csv.gz \
              seri_uncorrected_true_configuration_recovery.csv.gz; do
    [ -e "$BEFORE/$RESULT" ] || cp -p "out/montecarlo/$RESULT" "$BEFORE/$RESULT"
done
step "results copied to $BEFORE"

for GRID in seri_held_out_mean_true_configuration_recovery seri_uncorrected_true_configuration_recovery; do
    step "$GRID started"
    ( ulimit -v 12000000
      systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="$GRID" \
                      --why="Seri grid on the differenced local projections" --mode=block \
                      "./bin/$GRID" ) >> "out/montecarlo/$GRID.log" 2>&1
    STATUS=$?
    step "$GRID exit status $STATUS"
    [ $STATUS -eq 0 ] || exit $STATUS
done
for RESULT in seri_held_out_mean_true_configuration_recovery_100_mean_runs.csv.gz \
              seri_uncorrected_true_configuration_recovery.csv.gz; do
    echo "$RESULT rows per model: $(zcat "out/montecarlo/$RESULT" | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
done
step "finished"
