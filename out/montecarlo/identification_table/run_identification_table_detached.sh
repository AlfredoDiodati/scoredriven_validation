#!/bin/bash
# Rebuilds out/montecarlo/identification_table/ from the four grids' results,
# detached and capped at 10 GB of address space, since no systemd user session
# is reachable to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/montecarlo/identification_table/run_identification_table.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/identification_table >> "$LOG" 2>&1 || exit 1
( ulimit -v 10000000
  ./bin/identification_table ) >> "$LOG" 2>&1
STATUS=$?
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
