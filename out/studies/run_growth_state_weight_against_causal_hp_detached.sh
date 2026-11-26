#!/bin/bash
# Runs studies/growth_state_weight_against_causal_hp.c detached on all hardware
# threads, capped at 8 GB of address space, since no systemd user session is
# reachable to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/studies/growth_state_weight_against_causal_hp.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/growth_state_weight_against_causal_hp >> "$LOG" 2>&1 || exit 1
( ulimit -v 8000000
  OMP_NUM_THREADS=$(nproc) ./bin/growth_state_weight_against_causal_hp ) >> "$LOG" 2>&1
STATUS=$?
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
