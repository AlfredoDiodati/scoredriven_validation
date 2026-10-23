#!/bin/bash
# Logs every 30 seconds how many t-QVARMA fit caches the timing run has written.
log=out/applications/fit_t188_check/qvarma_t188_random32_progress.txt
start=$(date +%s)
echo "seconds fits_written" > "$log"
while true; do
    n=$(find out/applications/fit_t188_check/qvarma_t188_random32 -name '*_fit.json' 2>/dev/null | wc -l)
    echo "$(( $(date +%s) - start )) $n" >> "$log"
    sleep 30
done
