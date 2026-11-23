#!/bin/bash
# Runs tests/response_cache_state_weight_correctness detached, capped at 8 GB of
# address space, since no systemd user session is reachable to cap resident
# memory. The cap covers the 3.2 GB both-states cache it maps read-only.
LOG=out/tests/response_cache_state_weight_correctness.log
echo "started $(date '+%F %T')" >> "$LOG"
( ulimit -v 8000000
  ./bin/response_cache_state_weight_correctness ) >> "$LOG" 2>&1
STATUS=$?
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
