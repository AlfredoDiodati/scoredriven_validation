#!/bin/bash
# The two bootstrap true configuration recovery grids (Hansen et al.'s MCS_TR),
# corrected loss then uncorrected loss, adding qvarma_growth_states (the
# t-QVARMA with its two growth states in one vector) to the three models their
# results already hold. Each grid unpacks its result into its progress file,
# runs only the missing model and writes the result again; an interrupted grid
# resumes from its progress file.
#
# Before each grid its result is copied to out/montecarlo/before_qvarma_growth_states/.
# After it, that copy, unpacked, must be the first bytes of the new result; the
# log says whether it is, and gives the rows per model.
#
# Capped at 12 GB of address space, since no systemd user session is reachable
# to cap resident memory; the 4.2 GB cache is mapped and counts towards it.
cd /home/dioda/Documents/_py/scoredriven_validation || exit 1
LOG=out/montecarlo/bootstrap_qvarma_growth_states_grids.log
BEFORE=out/montecarlo/before_qvarma_growth_states
mkdir -p "$BEFORE"
echo "started $(date '+%F %T')" >> "$LOG"
make bin/true_configuration_recovery bin/uncorrected_true_configuration_recovery >> "$LOG" 2>&1 || exit 1
STATUS=0
for GRID in true_configuration_recovery uncorrected_true_configuration_recovery; do
  [ -e "$BEFORE/$GRID.csv.gz" ] || cp -n "out/montecarlo/$GRID.csv.gz" "$BEFORE/$GRID.csv.gz"
  echo "$GRID started $(date '+%F %T')" >> "$LOG"
  ( ulimit -v 12000000
    systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown --who="$GRID" \
                    --why="true configuration recovery of the t-QVARMA growth states" --mode=block \
                    "./bin/$GRID" ) >> "out/montecarlo/$GRID.log" 2>&1
  STATUS=$?
  echo "$GRID exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
  [ $STATUS -ne 0 ] && break
  OLD_BYTES=$(zcat "$BEFORE/$GRID.csv.gz" | wc -c)
  if cmp -s <(zcat "$BEFORE/$GRID.csv.gz") <(zcat "out/montecarlo/$GRID.csv.gz" | head -c "$OLD_BYTES"); then
    echo "$GRID: the $OLD_BYTES bytes of the earlier result are unchanged at the start of the new one" >> "$LOG"
  else
    echo "$GRID: FAILED, the earlier result is not the start of the new one" >> "$LOG"
    STATUS=1
    break
  fi
  echo "$GRID rows per model: $(zcat "out/montecarlo/$GRID.csv.gz" | cut -d, -f1 | sort | uniq -c | tr '\n' ' ')" >> "$LOG"
done
echo "exit status $STATUS ended $(date '+%F %T')" >> "$LOG"
