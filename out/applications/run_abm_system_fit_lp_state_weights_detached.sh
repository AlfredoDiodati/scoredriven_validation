#!/bin/bash
# Fits the local projections of every replicate under the three state weights
# other than lpirfs' (causal HP cycle at lambda 1600, two-sided cycle at the
# estimated lambda, causal cycle at the estimated lambda), on the full sample,
# the first 188 stored periods and stored periods 188 to 375: nine runs, one
# after the other. Each binary is capped at 4 GB of address space, since no
# systemd user session was reachable to cap resident memory.
LOG=out/applications/abm_system_fit_lp_state_weights.log
echo "started $(date '+%F %T')" >> "$LOG"
STATUS=0
for BIN in ./bin/abm_system_fit_lp ./bin/abm_system_fit_lp_t188 ./bin/abm_system_fit_lp_t188_from188; do
  for FLAGS in "--causal-hp" "--estimate-lambda" "--causal-hp --estimate-lambda"; do
    echo "$BIN $FLAGS started $(date '+%F %T')" >> "$LOG"
    ( ulimit -v 4000000
      systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="abm_system_fit_lp_state_weights" \
                      --why="nine million local projection fits" --mode=block $BIN $FLAGS ) >> "$LOG" 2>&1
    STATUS=$?
    echo "$BIN $FLAGS exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
    [ $STATUS -ne 0 ] && break 2
  done
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
