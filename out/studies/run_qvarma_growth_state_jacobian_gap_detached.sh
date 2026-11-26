#!/bin/bash
# Runs studies/qvarma_growth_state_jacobian_gap.c detached on all hardware
# threads, capped at 12 GB of address space, since no systemd user session is
# reachable to cap resident memory; the two mapped response caches, 6.3 GB,
# count towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/studies/qvarma_growth_state_jacobian_gap.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/qvarma_growth_state_jacobian_gap >> "$LOG" 2>&1 || exit 1
( ulimit -v 12000000
  OMP_NUM_THREADS=$(nproc) ./bin/qvarma_growth_state_jacobian_gap ) >> "$LOG" 2>&1
echo "exit status $? ended $(date '+%F %T')" >> "$LOG"
