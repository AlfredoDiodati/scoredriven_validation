#!/bin/bash
# One 500,000-fit battery, detached and inside a memory-capped systemd scope so
# a runaway cannot take down the terminal that started it. Records which et_al
# the run compiled against beside the results: the t-QVARMA likelihood was
# corrected at a large nu on 2026-08-31, and a run's numbers are only comparable
# with another run's if both used the same objective.
CAP="$1"
BIN=./bin/abm_system_scale_fit_qvarma_i${CAP}
STAMP=out/abm_system_scale_fit_qvarma_i${CAP}_provenance.txt
{
  echo "et_al headers this run compiled against"
  echo "  includedir: $(pkg-config --variable=includedir et_al.-core)"
  echo "  sd/qvarma.h  md5 $(md5sum /usr/local/include/et_al./sd/qvarma.h | cut -d' ' -f1)  mtime $(stat -c %y /usr/local/include/et_al./sd/qvarma.h)"
  echo "  special.h    md5 $(md5sum /usr/local/include/et_al./special.h | cut -d' ' -f1)  mtime $(stat -c %y /usr/local/include/et_al./special.h)"
  echo "  ad.h         md5 $(md5sum /usr/local/include/et_al./ad.h | cut -d' ' -f1)"
  echo "  carries the large-nu correction: $(grep -qc special_log1p /usr/local/include/et_al./special.h && echo yes || echo no)"
  echo "  binary: $BIN built $(stat -c %y $BIN)"
} > "$STAMP"
$BIN
echo "exit $?" >> "$STAMP"
