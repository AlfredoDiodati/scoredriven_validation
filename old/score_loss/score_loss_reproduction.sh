#!/usr/bin/env bash
#
# Whether the two score loss programs, after their move to old/ and the
# mat_chol status argument, still write the matrices already on disk.
#
# Each program writes into out/old/score_loss_reproduction/ through
# ABM_SYSTEM_SCORE_LOSS_PATH, so the existing matrices are only read. The two
# loss tables are compared byte for byte; the manifests name their own output
# paths, so for them the differing lines are printed instead.
#
# Reads dataset/abm_system/, out/us_qvarma_spec_choice_p1q1r2_fit.json,
# out/us_system.csv, out/abm_system_fit_qvarma/ and montecarlo/out/benchmark.env.
# Writes out/old/score_loss_reproduction/.

set -u

cd "$(dirname "$0")/../.."
check_dir=out/old/score_loss_reproduction
mkdir -p "$check_dir"

compare() {
    local reference=$1 rerun=$2
    if cmp -s "$reference" "$rerun"; then
        echo "identical $reference"
    else
        echo "differs $reference"
    fi
}

compare_manifest() {
    local reference=$1 rerun=$2
    echo "manifest lines that differ, $reference against $rerun:"
    diff "$reference" "$rerun" | grep '^[<>]'
}

echo "real-data score losses"
ABM_SYSTEM_SCORE_LOSS_PATH="$check_dir/abm_system_score_loss.csv" ./bin/abm_system_score_loss \
    || echo "FAILED abm_system_score_loss"
compare out/abm_system_score_loss.csv "$check_dir/abm_system_score_loss.csv"
compare out/abm_system_score_loss_weighted.csv "$check_dir/abm_system_score_loss_weighted.csv"
compare_manifest out/abm_system_score_loss_manifest.txt "$check_dir/abm_system_score_loss_manifest.txt"

echo "Monte Carlo score losses"
ABM_SYSTEM_SCORE_LOSS_PATH="$check_dir/score_loss.csv" ./bin/score_loss \
    || echo "FAILED score_loss"
compare montecarlo/out/score_loss.csv "$check_dir/score_loss.csv"
compare montecarlo/out/score_loss_weighted.csv "$check_dir/score_loss_weighted.csv"
compare_manifest montecarlo/out/score_loss_manifest.txt "$check_dir/score_loss_manifest.txt"

echo "done"
