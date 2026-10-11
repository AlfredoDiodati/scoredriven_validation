#!/usr/bin/env bash
# Reruns the four MCS statistic comparisons after moving the score losses to
# old/ and compares every output against the committed version in git.
set -u
cd "$(dirname "$0")/../.."
for program in abm_system_mcs_statistic_comparison mcs_statistic_comparison \
               abm_system_mcs_statistic_comparison_score mcs_statistic_comparison_score; do
    echo "running $program"
    ./bin/$program || echo "FAILED $program"
done
for f in $(git ls-files out/abm_system_mcs_statistic_comparison* montecarlo/out/mcs_statistic_comparison*); do
    if git diff --quiet -- "$f"; then
        echo "identical $f"
    else
        echo "differs $f"
        git diff -U0 -- "$f" | grep '^[-+][^-+]'
    fi
done
echo "done"
