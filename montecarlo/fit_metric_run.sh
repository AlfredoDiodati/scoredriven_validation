#!/usr/bin/env bash
#
# The fit metric of docs/MONTECARLO_FIT_METRIC.md on the 1000 by 1000 Monte
# Carlo, every step skipping what is already on disk, so running this again
# after an interruption continues where it stopped:
#
#   bin/fit_metric_oracle      the oracle version for every benchmark of the
#                              grid, about 4 minutes
#   bin/fit_metric_bootstrap   the bootstrap version for run 0 of every
#                              configuration, 200 resamples in blocks of 8,
#                              about 3 hours on 16 threads
#   the same on the first 100 configurations at 50 resamples, in blocks of 8,
#   25, 50 and 100 periods, to see how the bootstrap variance depends on the
#   block length, about 5 minutes each
#   bin/fit_metric_report      the tables
#
# Reads montecarlo/out/sweep_grid.csv.gz, out/sweep_grid_response_cache/,
# dataset/abm_system/ and out/abm_system_fit_qvarma/, and rebuilds none of them.

set -euo pipefail

cd "$(dirname "$0")/.."
make bin/fit_metric_oracle bin/fit_metric_bootstrap bin/fit_metric_report

echo "oracle $(date)"
[ -f montecarlo/out/fit_metric_oracle.csv.gz ] || ./bin/fit_metric_oracle

echo "bootstrap, run 0 of every configuration $(date)"
./bin/fit_metric_bootstrap

for block in 8 25 50 100; do
    echo "bootstrap, first 100 configurations, blocks of $block $(date)"
    FIT_METRIC_B=50 FIT_METRIC_BLOCK=$block FIT_METRIC_CONFIGURATIONS=100 ./bin/fit_metric_bootstrap
done

echo "report $(date)"
./bin/fit_metric_report
echo "done $(date)"
