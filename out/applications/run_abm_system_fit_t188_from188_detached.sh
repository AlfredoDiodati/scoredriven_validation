#!/bin/bash
# Fits stored periods 188 to 375 (0-based) of every replicate: the local
# projections first, then the t-QVARMA. Holds the machine awake for the whole
# run, since the first-188 run took about 28.5 hours.
LOG=out/applications/abm_system_fit_t188_from188.log
echo "started $(date '+%F %T')" >> "$LOG"
STATUS=0
for BIN in ./bin/abm_system_fit_lp_t188_from188 ./bin/abm_system_fit_qvarma_t188_from188; do
  echo "$BIN" >> "$LOG"
  systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="abm_system_fit_t188_from188" \
                  --why="about a million fits, more than a day" --mode=block "$BIN" >> "$LOG" 2>&1
  STATUS=$?
  [ $STATUS -ne 0 ] && break
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
