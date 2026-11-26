#!/bin/bash
# Runs studies/response_variance_run_count.c detached on all hardware threads,
# capped at 12 GB of address space, since no systemd user session is reachable
# to cap resident memory; the mapped response caches count towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/studies/response_variance_run_count.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/response_variance_run_count >> "$LOG" 2>&1 || exit 1
( ulimit -v 12000000
  OMP_NUM_THREADS=$(nproc) ./bin/response_variance_run_count ) >> "$LOG" 2>&1
echo "exit status $? ended $(date '+%F %T')" >> "$LOG"
