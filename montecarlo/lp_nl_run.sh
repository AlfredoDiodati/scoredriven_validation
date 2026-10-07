#!/usr/bin/env bash
#
# The state-dependent local projection with both states in one vector (lp_nl,
# _temp/Note on non-lin LP.pdf) through every Monte Carlo experiment after the
# single benchmark, which ./montecarlo/lp_run.sh runs. Every step skips what is
# already on disk, so after an interruption, a crash included, running this
# again continues where it stopped:
#
#   bin/lp_sweep        every run of cop_0191 as the benchmark, about 15 minutes
#   bin/sweep_cops      run 0 of every configuration, about 17 minutes
#   the raw stacked vector on runs 0 to 199, docs/MONTECARLO_COMPRESSED_RESPONSE.md,
#                       about 10 minutes
#   bin/sweep_grid      every run of every configuration, about a day and a half
#   the summaries of the grid's rows
#
# NUMPY_PYTHON is a Python with numpy, scipy and pandas (default python3);
# PLOT_PYTHON one with polars, plotly and kaleido, and the figures are skipped
# when it is not set. Reads out/abm_system_fit_lp/, dataset/abm_system/ and
# out/sweep_grid_response_cache/, and rebuilds none of them.

set -euo pipefail

cd "$(dirname "$0")/.."
NUMPY_PYTHON=${NUMPY_PYTHON:-python3}
make bin/lp_sweep bin/sweep_cops bin/sweep_grid bin/compressed_response_recovery

# A crash can leave the progress file ending in zero bytes, the part of the
# last write that never reached the disk. Every line before them is complete,
# so the file is cut at the first zero byte and those benchmarks are redone.
progress=montecarlo/out/sweep_grid_progress.csv
if [ -f $progress ]; then
    first_zero=$(grep -abo -m1 -P '\x00' $progress | head -1 | cut -d: -f1 || true)
    if [ -n "$first_zero" ]; then
        echo "cutting $progress at byte $first_zero, where its zero bytes start"
        truncate -s "$first_zero" $progress
    fi
fi

echo "every run of cop_0191 as the benchmark $(date)"
./bin/lp_sweep

echo "run 0 of every configuration as the benchmark $(date)"
./bin/sweep_cops

echo "the raw stacked vector on runs 0 to 199 $(date)"
raw=out/compressed_response_cache/lp_nl_raw.f32
recovery=montecarlo/out/compressed_response/recovery_lp_nl_raw_rows200.csv
[ -f $raw ] || $NUMPY_PYTHON montecarlo/compressed_response_learn.py lp_nl raw_only
[ -f $recovery ] || ./bin/compressed_response_recovery $raw 0 4 0 200 $recovery

echo "every run of every configuration as the benchmark $(date)"
./bin/sweep_grid

echo "summaries $(date)"
$NUMPY_PYTHON montecarlo/compressed_response_grid_baseline.py
$NUMPY_PYTHON montecarlo/compressed_response_report.py
python3 montecarlo/sweep_grid_presence.py
if [ -n "${PLOT_PYTHON:-}" ]; then
    $PLOT_PYTHON montecarlo/sweep_grid_identifiability_plots.py
else
    echo "PLOT_PYTHON not set, figures skipped"
fi
echo "done $(date)"
