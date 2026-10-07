#!/usr/bin/env bash
#
# The Monte Carlo experiment under the local projections of the collaborator's
# R pipeline: the same benchmark run as montecarlo/run.sh, the same hold-out
# and the same confidence set settings, with the local projections in place of
# the t-QVARMA as the auxiliary model.
#
# The benchmark is the one montecarlo/out/benchmark.env already names, so both
# auxiliary models are asked to recover the same known answer from the same
# series. bin/benchmark_choice is not rerun here: it chooses on grounds that
# belong to the t-QVARMA fits, and rerunning it could move the benchmark under
# results already written.
#
#   ./montecarlo/lp_run.sh        the single-benchmark run, a few minutes
#   ./bin/lp_sweep                every replicate of the benchmark configuration
#                                 as the benchmark, four models: linear, each
#                                 state, and both states as one vector
#
# Reads out/abm_system_fit_lp/ (make app-abm_system_fit_lp writes it) and
# dataset/abm_system/, and rebuilds neither. docs/MONTECARLO_LP_VALIDATION.md is
# the write-up.

set -euo pipefail

cd "$(dirname "$0")/.."
test -f montecarlo/out/benchmark.env || {
    echo "montecarlo/out/benchmark.env is missing; ./montecarlo/run.sh chooses the benchmark"
    exit 1
}
sed 's/^/  /' montecarlo/out/benchmark.env

echo "local-projection losses, linear, each state and both states"
./bin/lp_irf_loss

echo "confidence sets, both statistics over the four losses"
./bin/lp_mcs

echo "done, results in montecarlo/out/lp_*"
