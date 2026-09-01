#!/bin/bash
# One 500,000-fit battery, detached, memory-capped, and timed from outside the
# binary as well as inside it. Replaces run_scale_fit_detached.sh, which is left
# in place because a running bash reads its script incrementally and editing the
# file under a live run corrupts it. Swap the names when nothing is running.
#
# Records, beside the results:
#   - which et_al the run compiled against, since the t-QVARMA likelihood was
#     corrected at a large nu on 2026-08-31 and two runs are only comparable on
#     their budgets if both used the same objective
#   - the wall time of the whole process, which the binary's own timing file
#     does not cover: that one starts its clock after the 500 folders have been
#     scanned, stops it before the manifest is written, and writes nothing if
#     the process is killed
CAP="$1"
BIN=./bin/abm_system_scale_fit_qvarma_i${CAP}
STAMP=out/abm_system_scale_fit_qvarma_i${CAP}_provenance.txt
START_EPOCH=$(date +%s)
{
  echo "et_al headers this run compiled against"
  echo "  includedir: $(pkg-config --variable=includedir et_al.-core)"
  for h in sd/qvarma.h special.h ad.h dist/mv/student.h; do
    echo "  $h  md5 $(md5sum /usr/local/include/et_al./$h | cut -d' ' -f1)  mtime $(stat -c %y /usr/local/include/et_al./$h)"
  done
  grep -q special_log1p /usr/local/include/et_al./special.h \
    && echo "  carries the large-nu correction: yes" \
    || echo "  carries the large-nu correction: no"
  echo "  binary: $BIN built $(stat -c %y $BIN)"
  echo "  host: $(uname -srm), $(nproc) hardware threads, OMP_NUM_THREADS=${OMP_NUM_THREADS:-default}"
} > "$STAMP"

# The machine suspended under the first 8000-iteration attempt at 01:58 on
# 2026-09-01, 296,337 fits in, and was rebooted before it could resume. A
# battery of this length has to hold the machine awake for its own duration.
systemd-inhibit --what=sleep:idle --who="abm_system_scale_fit_qvarma_i${CAP}" \
                --why="500,000-fit battery, hours long" $BIN
STATUS=$?
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))
{
  echo
  echo "Execution time, measured from outside the process"
  echo "  process started  $(date -d @$START_EPOCH '+%F %T')"
  echo "  process ended    $(date -d @$END_EPOCH '+%F %T')"
  printf "  wall clock       %d s (%s h, %02d:%02d:%02d)\n" \
         "$ELAPSED" "$(awk -v s="$ELAPSED" 'BEGIN{printf "%.2f", s/3600}')" \
         $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
  echo "  covers the folder scan, all 500,000 fits and the manifest write;"
  echo "  the timing file's own 'elapsed' covers the fit loop alone."
  grep -h "^elapsed" out/abm_system_scale_fit_qvarma_i${CAP}_timing.txt 2>/dev/null |
    sed 's/^/  fit loop only:   /'
  echo "  exit status      $STATUS"
} >> "$STAMP"
