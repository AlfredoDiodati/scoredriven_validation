#!/usr/bin/env bash
#
# The Monte Carlo experiment: the main pipeline's three validation protocols,
# run against one simulated run instead of against the US data.
#
# This is a pipeline of its own. It shares no source with applications/, which
# it does not touch, and it reuses that pipeline's output rather than its code:
# above all out/abm_system_fit_qvarma/, the million cached fits, because
# nothing here estimates anything. The benchmark's own parameters are one of
# those cached fits, promoted.
#
# bin/benchmark_choice picks the run and writes montecarlo/out/benchmark.env;
# every later step reads it back, so one experiment cannot end up describing
# two benchmarks. The replicate it picks is held out of every configuration's
# column. docs/MONTECARLO_VALIDATION.md is the write-up.
#
#   ./montecarlo/run.sh           the whole experiment
#   ./montecarlo/run.sh --keep    reuse the benchmark already chosen
#
# Reads out/abm_system_fit_qvarma/ and dataset/abm_system/ and rebuilds
# neither, so it does not run on a fresh clone.

set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p montecarlo/out

if [ "${1:-}" != "--keep" ] || [ ! -f montecarlo/out/benchmark.env ]; then
    ./bin/benchmark_choice
fi
sed 's/^/  /' montecarlo/out/benchmark.env

# The impulse-response loss, from the fits already on disk. Every replicate's
# IRF is compared against the benchmark replicate's own.
echo "impulse-response loss"
./bin/irf_loss

# Both score losses in one pass: q'q and the inverse-information-weighted score
# statistic, the information taken at the benchmark's estimate on its series.
echo "score losses"
./bin/score_loss

# The headline confidence set, MCS_TR over the impulse-response loss.
echo "confidence set, impulse-response loss"
./bin/mcs

# All three losses under both MCS statistics, six confidence sets.
echo "confidence sets, both statistics over all three losses"
./bin/mcs_statistic_comparison

echo "done, results in montecarlo/out"
