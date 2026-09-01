#!/bin/bash
# Waits for a running battery to exit and appends the wall time measured from
# outside it. The binary's own timing file measures its fit loop with
# omp_get_wtime, which starts after the 500 folders have been scanned and stops
# before the manifest is written, and which records nothing at all if the
# process is killed. This covers the whole process either way.
PID="$1"; CAP="$2"; START_EPOCH="$3"
STAMP=out/abm_system_scale_fit_qvarma_i${CAP}_provenance.txt
while kill -0 "$PID" 2>/dev/null; do sleep 20; done
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))
{
  echo
  echo "Execution time, measured from outside the process"
  echo "  process started  $(date -d @$START_EPOCH '+%F %T')"
  echo "  process ended    $(date -d @$END_EPOCH '+%F %T')"
  printf "  wall clock       %d s (%.2f h, %02d:%02d:%02d)\n" \
         "$ELAPSED" "$(echo "$ELAPSED/3600" | bc -l)" \
         $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
  echo "  covers the folder scan, all 500,000 fits and the manifest write;"
  echo "  the timing file's own 'elapsed' covers the fit loop alone."
  grep -h "^elapsed" out/abm_system_scale_fit_qvarma_i${CAP}_timing.txt 2>/dev/null |
    sed 's/^/  fit loop only:   /'
} >> "$STAMP"
