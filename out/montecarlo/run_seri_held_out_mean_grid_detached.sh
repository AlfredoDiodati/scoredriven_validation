#!/bin/bash
# The Seri et al. true configuration recovery grid with each configuration's
# mean response taken from its last 100 runs and the loss matrices formed on
# runs 0 to 899 (montecarlo/seri_held_out_mean_true_configuration_recovery.c),
# for the linear projection, both states under the causal HP cycle at the
# estimated lambda and the t-QVARMA with its two growth states. Resumes from its
# progress file after an interruption.
#
# Capped at 12 GB of address space, since no systemd user session is reachable
# to cap resident memory; the 4.2 GB cache is mapped and counts towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
GRID=seri_held_out_mean_true_configuration_recovery
LOG=out/montecarlo/seri_held_out_mean_grid.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/$GRID >> "$LOG" 2>&1 || exit 1
( ulimit -v 12000000
  systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="$GRID" \
                  --why="true configuration recovery with a held-out mean" --mode=block \
                  "./bin/$GRID" ) >> "out/montecarlo/$GRID.log" 2>&1
STATUS=$?
echo "$GRID exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
RESULT=out/montecarlo/${GRID}_100_mean_runs.csv.gz
[ -e "$RESULT" ] && echo "rows per model: $(zcat "$RESULT" | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
