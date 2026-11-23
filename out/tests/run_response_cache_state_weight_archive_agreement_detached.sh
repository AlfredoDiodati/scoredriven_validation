#!/bin/bash
# Runs tests/response_cache_state_weight_archive_agreement detached, capped at 10 GB of
# address space, since no systemd user session is reachable to cap resident
# memory. The cap covers the 6.4 GB of caches of one weight it maps read-only.
LOG=out/tests/response_cache_state_weight_archive_agreement.log
echo "started $(date '+%F %T')" >> "$LOG"
( ulimit -v 10000000
  ./bin/response_cache_state_weight_archive_agreement ) >> "$LOG" 2>&1
STATUS=$?
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
