#!/usr/bin/env bash
#
# Discarded from the Monte Carlo pipeline. Its single benchmark, cop_0191 run 706,
# is one of the benchmarks old/naive_mae_loss/montecarlo/sweep_grid.c runs, every run of every
# configuration, under the same loss and confidence set settings except the
# number of bootstrap resamples (10000 here, 2000 in the grid). The grid
# answers the same question on all of them. docs/OLD_EXPERIMENTS.md.
#
# The Monte Carlo experiment: the main pipeline's impulse-response validation
# protocol, run against one simulated run instead of against the US data.
#
# This is a pipeline of its own. It shares no source with applications/, which
# it does not touch, and it reuses that pipeline's output rather than its code:
# above all out/abm_system_fit_qvarma/, the million cached fits, because
# nothing here estimates anything. The benchmark's own parameters are one of
# those cached fits, promoted.
#
# bin/benchmark_choice picks the run and writes out/montecarlo/benchmark.env;
# every later step reads it back, so one experiment cannot end up describing
# two benchmarks. The replicate it picks is held out of every configuration's
# column. docs/MONTECARLO_VALIDATION.md is the write-up.
#
#   ./old/montecarlo_one_benchmark/run.sh the whole experiment
#   ./old/montecarlo_one_benchmark/run.sh --keep reuse the benchmark already chosen
#
# Reads out/abm_system_fit_qvarma/ and out/dataset/abm_system/ and rebuilds
# neither, so it does not run on a fresh clone.

set -euo pipefail

cd "$(dirname "$0")/../.."
mkdir -p out/montecarlo

if [ "${1:-}" != "--keep" ] || [ ! -f out/montecarlo/benchmark.env ]; then
    ./bin/benchmark_choice
fi
sed 's/^/  /' out/montecarlo/benchmark.env

# The impulse-response loss, from the fits already on disk. Every replicate's
# IRF is compared against the benchmark replicate's own.
echo "impulse-response loss"
./bin/irf_loss

# The headline confidence set, MCS_TR over the impulse-response loss.
echo "confidence set, impulse-response loss"
./bin/mcs

# The same loss under both MCS statistics, two confidence sets.
echo "confidence sets, both statistics over the impulse-response loss"
./bin/mcs_statistic_comparison

echo "done, results in out/montecarlo"
