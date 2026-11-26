#!/bin/bash
# The two Seri et al. true configuration recovery grids, corrected loss then
# uncorrected loss, over both states in one vector under the causal HP cycle at
# the estimated lambda (model 13 of montecarlo/response_cache.h), the one
# state-dependent model the grids of montecarlo/ run. Each grid unpacks its
# result into its progress file, runs only the models missing from it and
# writes the result again with every model it runs; an interrupted grid resumes
# from its progress file, and a finished one reports nothing to do.
#
# After each grid, its result as it was before the new weights were added
# (out/montecarlo/before_state_weights/), without the rows of state 1 and state
# 2 on their own and of both states under lpirfs' weight, must be the first
# bytes of the new result, unpacked; the log says whether it is.
#
# Capped at 10 GB of address space, since no systemd user session is reachable
# to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/montecarlo/seri_state_weight_grids.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/seri_true_configuration_recovery bin/seri_uncorrected_true_configuration_recovery >> "$LOG" 2>&1 || exit 1
STATUS=0
for GRID in seri_true_configuration_recovery seri_uncorrected_true_configuration_recovery; do
  echo "$GRID started $(date '+%F %T')" >> "$LOG"
  ( ulimit -v 10000000
    systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="$GRID" \
                    --why="true configuration recovery under the other state weights" --mode=block \
                    "./bin/$GRID" ) >> "out/montecarlo/$GRID.log" 2>&1
  STATUS=$?
  echo "$GRID exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
  [ $STATUS -ne 0 ] && break
  OLD=out/montecarlo/before_state_weights/$GRID.csv.gz
  OLD_BYTES=$(zcat "$OLD" | awk -F, '$1 !~ /^lp_s[12](_|$)/ && $1 != "lp_nl"' | wc -c)
  if cmp -s <(zcat "$OLD" | awk -F, '$1 !~ /^lp_s[12](_|$)/ && $1 != "lp_nl"') <(zcat "out/montecarlo/$GRID.csv.gz" | head -c "$OLD_BYTES"); then
    echo "$GRID: the $OLD_BYTES bytes of the earlier result are unchanged at the start of the new one" >> "$LOG"
  else
    echo "$GRID: FAILED, the earlier result is not the start of the new one" >> "$LOG"
    STATUS=1
    break
  fi
  echo "$GRID rows per model: $(zcat "out/montecarlo/$GRID.csv.gz" | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
