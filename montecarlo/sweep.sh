#!/usr/bin/env bash
#
# The Monte Carlo experiment repeated with every replicate of the benchmark
# configuration standing in as the benchmark, so that a protocol which
# identifies a configuration can be told from one that drew well once.
#
# montecarlo/run.sh answers the question for a single benchmark. This answers
# it a thousand times and reports how often each of the three losses returns
# the configuration the benchmark came from.
#
# Two programs, run in that order because the second is the long one:
#
#   bin/sweep_irf     the impulse-response protocol, about 3 hours here
#   bin/sweep_score   both score protocols, about 12 hours here
#
# Both are resumable at benchmark granularity: each finished benchmark is
# appended to its summary as it completes and skipped on a rerun, so an
# interrupted sweep continues where it stopped and partial summaries are
# already readable. Running this again after it finishes does nothing.
#
# Measured costs behind those figures, on 16 cores of this machine: one pass of
# 10^6 filter evaluations plus the 15 GB dataset read takes 46 s, one pass of
# 10^6 cached fits and impulse responses takes 574 s, and one MCS_TR over 1000
# models and 999 observations at 10000 resamples takes 10.9 s.
#
# Reads out/abm_system_fit_qvarma/ and dataset/abm_system/ and rebuilds
# neither. Writes three summaries under montecarlo/out/. Progress goes to
# stderr; the summaries are the result.

set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p montecarlo/out

echo "impulse-response protocol, every benchmark"
./bin/sweep_irf

echo "score protocols, every benchmark"
./bin/sweep_score

echo "done, summaries in montecarlo/out/sweep_*.csv"
