#!/bin/bash
# The two bootstrap true configuration recovery grids (Hansen et al.'s MCS,
# MCS_TR at 2000 resamples), corrected loss then uncorrected loss, over both
# states in one vector under the three state weights other than lpirfs'
# (models 7, 10 and 13 of montecarlo/response_cache.h). Each grid unpacks its
# result into its progress file, runs only the models missing from it and
# writes the result again with every model it runs; an interrupted grid resumes
# from its progress file, and a finished one reports nothing to do.
#
# After each grid, its result as it was before the new weights were added
# (out/montecarlo/before_state_weights/), without the rows of state 1 and state
# 2 on their own, which are in out/old/montecarlo_single_states/, must be the
# first bytes of the new result, unpacked; the log says whether it is.
#
# Capped at 10 GB of address space, since no systemd user session is reachable
# to cap resident memory.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/montecarlo/bootstrap_state_weight_grids.log
echo "started $(date '+%F %T')" >> "$LOG"
make bin/true_configuration_recovery bin/uncorrected_true_configuration_recovery >> "$LOG" 2>&1 || exit 1
STATUS=0
for GRID in true_configuration_recovery uncorrected_true_configuration_recovery; do
  echo "$GRID started $(date '+%F %T')" >> "$LOG"
  ( ulimit -v 10000000
    systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="$GRID" \
                    --why="true configuration recovery under the other state weights" --mode=block \
                    "./bin/$GRID" ) >> "out/montecarlo/$GRID.log" 2>&1
  STATUS=$?
  echo "$GRID exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
  [ $STATUS -ne 0 ] && break
  OLD=out/montecarlo/before_state_weights/$GRID.csv.gz
  OLD_BYTES=$(zcat "$OLD" | awk -F, '$1 !~ /^lp_s[12](_|$)/' | wc -c)
  if cmp -s <(zcat "$OLD" | awk -F, '$1 !~ /^lp_s[12](_|$)/') <(zcat "out/montecarlo/$GRID.csv.gz" | head -c "$OLD_BYTES"); then
    echo "$GRID: the $OLD_BYTES bytes of the earlier result are unchanged at the start of the new one" >> "$LOG"
  else
    echo "$GRID: FAILED, the earlier result is not the start of the new one" >> "$LOG"
    STATUS=1
    break
  fi
  echo "$GRID rows per model: $(zcat "out/montecarlo/$GRID.csv.gz" | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
