#!/usr/bin/env bash
#
# Whether moving the R-levels variant out of the live local-projection code
# changed any result, on either side of the move.
#
# The live programs, now without the --r-levels flag, and the old R-levels
# programs are each run once and every output is compared with what was on disk
# before: the tracked outputs against git, the untracked loss tables against
# copies taken here first. The loss programs check every stored fit against a
# fingerprint of the series and the state series it was fitted on, so identical
# loss tables also mean the state series is unchanged on all 1,000,000 runs.
# The confidence-set reports print their own timings, so for those the
# differing lines are printed instead of a verdict.
#
# Reads out/us_system.csv, out/abm_system_fit_lp/, out/abm_system_fit_lp_r_levels/
# and out/dataset/abm_system/. Writes out/old/lp_r_levels_reproduction/.

set -u

cd "$(dirname "$0")/../.."
check_dir=out/old/lp_r_levels_reproduction
mkdir -p "$check_dir/before"

untracked="out/old/naive_mae_loss/abm_system_lp_irf_loss_lin.csv out/old/naive_mae_loss/abm_system_lp_irf_loss_s1.csv out/old/naive_mae_loss/abm_system_lp_irf_loss_s2.csv
out/abm_system_lp_r_levels_irf_loss_lin.csv out/abm_system_lp_r_levels_irf_loss_s1.csv
out/abm_system_lp_r_levels_irf_loss_s2.csv"
for f in $untracked; do cp "$f" "$check_dir/before/"; done

for program in lp_system_transform lp_system_r_levels_transform us_lp_fit us_lp_fit_r_levels \
               abm_system_lp_irf_loss abm_system_lp_r_levels_irf_loss \
               abm_system_lp_mcs abm_system_lp_r_levels_mcs; do
    echo "running $program $(date +%T)"
    ./bin/$program || echo "FAILED $program"
done

echo "fit reports, loaded rather than refitted:"
grep -H "loaded\|fitted in this run" out/us_lp_fit.txt out/us_lp_fit_r_levels.txt

for f in $untracked; do
    if cmp -s "$f" "$check_dir/before/$(basename "$f")"; then echo "identical $f"; else echo "differs $f"; fi
done

for f in $(git ls-files out/us_lp_fit* out/abm_system_lp_*); do
    if git diff --quiet -- "$f"; then
        echo "identical $f"
    else
        echo "differs $f"
        git diff -U0 -- "$f" | grep '^[-+][^-+]' | head -20
    fi
done
echo "done"
