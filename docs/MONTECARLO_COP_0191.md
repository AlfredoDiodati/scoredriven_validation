# cop_0191 against the whole Monte Carlo grid

`cop_0191` is the configuration the impulse-response protocol keeps on the US
data and the one every single-benchmark experiment uses as the truth. This page
asks two questions of the 1000 by 1000 Monte Carlo: how often is `cop_0191`
recovered when it generated the data, and how often is it returned when another
configuration did. The third table puts both against the recovery rate of
every configuration.

## Setup

Every number comes from `montecarlo/out/sweep_grid.csv.gz`, the result of
`montecarlo/sweep_grid.c`; `docs/MONTECARLO_VALIDATION.md`, "Every run of every
configuration as the benchmark", has the full description. In short:

- Each of the 1,000,000 simulated runs (1000 configurations by 1000 runs) is
  the benchmark once, standing in for the US data, so each configuration is the
  truth for 1000 benchmarks.
- The loss of a run is the mean absolute difference between its response vector
  and the benchmark's. The benchmark's run index is held out of every
  configuration, so each confidence set compares 1000 configurations on 999
  runs. No run was dropped for a missing response.
- The confidence set is MCS_TR at level 0.05, 2000 bootstrap resamples, block
  length 1, bootstrap variance, seed 123 stream 0.
- Five models give the response vector: the t-QVARMA (qvarma, p1q1r2, 525
  entries), and the local projections of `docs/MONTECARLO_LP_VALIDATION.md`:
  linear, state 1 and state 2 (400 entries each), and both states, state 1's
  and state 2's responses stacked into one vector of 800 entries as
  `_temp/Note on non-lin LP.pdf` prescribes.

Three outcomes are counted for each benchmark, all about one named
configuration: whether the confidence set holds it, whether it is the only
configuration in the set, and whether it has the smallest mean loss of all
1000.

The tables are written by `montecarlo/sweep_grid_cop_0191.py`, in percent as
below, to `montecarlo/out/sweep_grid_cop_0191.txt` and as LaTeX tables to
`montecarlo/out/sweep_grid_cop_0191.tex`.

    python montecarlo/sweep_grid_cop_0191.py

## cop_0191 is the true configuration

Its own 1000 runs as the benchmark; share of them where `cop_0191` is:

| model | in the set | alone in the set | smallest mean loss |
| --- | --- | --- | --- |
| qvarma | 87.60% | 78.10% | 82.90% |
| linear | 91.30% | 86.40% | 88.60% |
| state 1 | 85.70% | 64.60% | 75.80% |
| state 2 | 78.00% | 39.90% | 59.20% |
| both states | 91.60% | 66.80% | 80.50% |

## cop_0191 is not the true configuration

The other 999 configurations' runs as the benchmark, 999,000 benchmarks; share
of them where `cop_0191`, which did not generate the data, is:

| model | in the set | alone in the set | smallest mean loss |
| --- | --- | --- | --- |
| qvarma | 2.42% | 1.50% | 1.89% |
| linear | 2.83% | 1.56% | 2.08% |
| state 1 | 16.62% | 5.78% | 9.63% |
| state 2 | 12.04% | 2.11% | 4.84% |
| both states | 12.71% | 3.22% | 6.07% |

The "alone in the set" column counts the benchmarks where the procedure returns
`cop_0191` as its only answer and that answer is wrong.

## Every configuration as the true configuration

All 1,000,000 benchmarks, `cop_0191`'s included; share of them where the
benchmark's own configuration is:

| model | in the set | alone in the set | smallest mean loss |
| --- | --- | --- | --- |
| qvarma | 9.52% | 3.16% | 5.39% |
| linear | 10.70% | 3.90% | 6.34% |
| state 1 | 1.53% | 0.47% | 0.79% |
| state 2 | 1.61% | 0.36% | 0.73% |
| both states | 1.17% | 0.43% | 0.67% |

Under every model `cop_0191` is recovered far more often than a configuration
in general. Under the three state-dependent models it is also kept as a wrong
answer more often, 12% to 17% of the time, than the average configuration is
kept as the right one.
