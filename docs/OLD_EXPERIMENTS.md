# Discarded experiments

Experiments that are no longer part of the main pipeline or the Monte Carlo
pipeline. Their code is under `old/`, one folder per experiment; their
write-ups stay in `docs/`. Their outputs are under `out/old/`, in the folder
that mirrors the program's own, so `old/naive_mae_loss/montecarlo/` writes to
`out/old/naive_mae_loss/montecarlo/`. Every program still builds and runs with
`make old-<stem>`, the stem being its file name, and the Makefile's
`OLD_SOURCES` lists them all.

Each entry says what the experiment was, why it is out, and where its results
are recorded.

## The score losses, `old/score_loss/`

**What it was.** A loss matrix built without fitting the auxiliary model to
the simulated data: the score of the t-QVARMA log-likelihood at the benchmark's
estimate, evaluated on each simulated series, either unweighted (`q'q`) or
weighted by the inverse information matrix, `q' [T_sim I_1]^-1 q`. It ran
against the US data and against the Monte Carlo benchmark.

**Why it is out.** It gave poor results where the right answer is known: on
the single Monte Carlo benchmark the generating configuration came 5th under
`q'q` and 76th under the weighted loss. On the 389 runs of `cop_0191` that
were scored, `q'q` kept `cop_0191` in the set for 182 but ranked it first for
only 14, through large sets, and the weighted loss could not be built for 252
and kept `cop_0191` in 5 of the other 137. And there is no theoretical
justification for using either quantity as a distance between a simulated
configuration and the data.

**Where it is recorded.** `docs/OLD_SCORE_LOSS.md`, and
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` for the Monte Carlo results.

| folder | programs |
| --- | --- |
| `old/score_loss/applications/` | `abm_system_score_loss.c`, `abm_system_mcs_statistic_comparison_score.c` |
| `old/score_loss/montecarlo/` | `score_loss.c`, `mcs_statistic_comparison_score.c`, `sweep_score.c` |
| `old/score_loss/` | `score_loss_reproduction.sh`, the check that the moved programs still write the matrices on disk byte for byte |

## Monte Carlo experiments on subsets of the benchmarks

Three folders, one per experiment, discarded for the same reason.

**What they were.** The Monte Carlo question, whether the confidence set holds
the configuration a simulated benchmark was drawn from, asked of fewer
benchmarks than the grid does:

| folder | benchmarks |
| --- | --- |
| `old/montecarlo_one_benchmark/` | `cop_0191` run 706, chosen by `benchmark_choice.c`, under the t-QVARMA (`run.sh`) and the local projections (`lp_run.sh`) |
| `old/montecarlo_every_run_of_cop_0191/` | each of the 1000 runs of `cop_0191` (`sweep_irf.c`, `lp_sweep.c`, `sweep.sh`) |
| `old/montecarlo_run_0_of_every_configuration/` | run 0 of each of the 1000 configurations (`sweep_cops.c`) |

**Why they are out.** Each set of benchmarks is part of the grid,
`old/naive_mae_loss/montecarlo/sweep_grid.c`, which makes every run of every configuration the
benchmark in turn, with the same loss, the same hold-out and the same
confidence set settings. They are special cases of it. The one difference is
the number of bootstrap resamples, 10000 in these and 2000 in the grid, and on
the benchmarks they share the two agree: whether the benchmark's own
configuration is in the set differs in 1 of 5000 benchmarks for run 0 of every
configuration, by at most 2 of 1000 for the runs of `cop_0191`, and not at all
for the single benchmark. `docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md`, "Status:
discarded", has the counts.

`old/montecarlo_one_benchmark/studies/` holds the two studies that diagnosed
these experiments' results, `qvarma_recovery_diagnostics.c` and
`lp_recovery_diagnostics.c`, with the header they share. They read the single
benchmark's choice and loss tables, so they moved with it and are no longer
part of `make study`.

**Where they are recorded.** `docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` for the
t-QVARMA and the run 0 experiment, `docs/OLD_MONTECARLO_LP_VALIDATION.md` for
the local projections on the single benchmark and on every run of `cop_0191`,
including the finding that decides how the local projections are fitted.

The drivers run with `make old-montecarlo_one_benchmark`, `make
old-montecarlo_one_benchmark_lp` and `make
old-montecarlo_every_run_of_cop_0191`.

## Local projections on cumulated levels, `old/lp_r_levels/`

**What it was.** The local projections fitted, on the US data and on every
simulated run, on levels rebuilt by cumulating the five stored series from
fixed starting values, with GDP and energy in logs, in the column order GDP,
employment, inflation, interest rate, energy. This is the data preparation of
the R code the local-projection settings were first taken from (its function
`transform_data()`); `old/lp_r_levels/lp_system_r_levels.h` is the C port and
states the arithmetic. Every experiment that remains fits the five series of
the main transformation, the ones the t-QVARMA is fitted on: GDP growth, energy
growth, employment change, inflation and the interest rate
(`docs/DATA_DOCUMENTATION.md`).

**Why it is out.** The main transformation loses less information. In the Monte
Carlo experiment, where the right answer is known, the linear local projection
on the cumulated levels put it in the set for 0 of 100 runs of `cop_0191`,
against 90 of 100 on the main transformation.

The state-dependent local projections that remain still build the series that
decides each period's state the same way, the 4-period moving average of the
rebuilt log GDP; the series they are fitted on are the main transformation's.

**Where it is recorded.** `docs/OLD_LP_R_LEVELS.md` for the US data,
`docs/OLD_MONTECARLO_LP_VALIDATION.md` for the Monte Carlo.

| folder | programs |
| --- | --- |
| `old/lp_r_levels/` | `lp_system_r_levels.h`, `us_lp_fit_r_levels.c`, `abm_system_fit_lp_r_levels.c`, `abm_system_lp_r_levels_irf_loss.c`, `abm_system_lp_r_levels_mcs.c`, `lp_system_r_levels_transform.c`, `lp_r_levels_reproduction.sh` |

## The naive loss, `old/naive_mae_loss/`

**What it was.** The validation procedure with the naive loss: the mean
absolute error between a replicate's impulse-response vector and the
benchmark's, with no correction for the estimation noise of the simulated
responses. It ran on the US data, under the t-QVARMA and the local projections,
and in the Monte Carlo grid over every run of every configuration, and the
experiments that read either run used it too: the confidence set's winner and
its impulse responses, the robustness studies of that winner, the compressed
responses, the fit metric and `cop_0191` as the truth.

**Why it is out.** The expected difference between two configurations' losses
contains the difference of their response variances, which says nothing about
how close either is to the data. The confidence set can therefore favour a
configuration whose mean response is far from the data's but whose responses
are less noisy. `docs/MCS_VALIDATION_PROCEDURE.md`, "Why the variance
correction", gives the argument and the corrected loss.

**Where it is recorded.** `docs/OLD_NAIVE_MAE_LOSS.md` for the procedure,
`docs/OLD_ABM_SYSTEM_MCS_VALIDATION.md` and
`docs/OLD_ABM_SYSTEM_LP_VALIDATION.md` for the US data,
`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md` for the Monte Carlo grid, and
`docs/OLD_MONTECARLO_COP_0191.md`,
`docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md`,
`docs/OLD_MONTECARLO_FIT_METRIC.md` and
`docs/OLD_MONTECARLO_FIT_METRIC_COMPARISON.md` for the experiments that read
the grid. `docs/QVARMA_FIT_STUDIES.md` still describes the studies of the
winner.

| folder | programs |
| --- | --- |
| `old/naive_mae_loss/applications/` | `abm_system_irf_loss.c`, `abm_system_mcs.c`, `abm_system_mcs_statistic_comparison.c`, `abm_system_lp_irf_loss.c`, `abm_system_lp_mcs.c`, `abm_system_winner_irf.c`, `abm_system_winner_irf_plots.py` |
| `old/naive_mae_loss/studies/` | `us_qvarma_nu_sensitivity.c`, `abm_system_winner_diagnostics.c`, `abm_system_winner_diagnostics_plots.py`, `abm_system_winner_nu_profile.c`, `abm_system_winner_nu_likelihood_scan.c`, `abm_system_winner_tail_comparison.c`, `abm_system_winner_normality.c`, `abm_system_seed_correlation.c` |
| `old/naive_mae_loss/montecarlo/` | `sweep_grid.c`, `lp_nl_run.sh`, `sweep_grid_presence.py`, `sweep_grid_cop_0191.py`, `sweep_grid_identifiability_plots.py`, `compressed_response_learn.py`, `compressed_response_recovery.c`, `compressed_response_grid_baseline.py`, `compressed_response_report.py`, `fit_metric.h`, `fit_metric_table.h`, `fit_metric_oracle.c`, `fit_metric_bootstrap.c`, `fit_metric_pairwise.c`, `fit_metric_report.c`, `fit_metric_run.sh`, `fit_metric_v_histograms.py` |
| `old/naive_mae_loss/tests/` | `fit_metric_correctness.c` |

The C programs build and run with `make old-<stem>`; the grid with `make
old-sweep_grid`, still at -O3. The loss and confidence-set targets keep their
prerequisites, so `make old-abm_system_winner_irf` reruns the US loss table and
confidence set first. The Python scripts and the two shell drivers run from the
repository root, for example `./old/naive_mae_loss/montecarlo/fit_metric_run.sh`.

The response caches the grid built, `out/montecarlo/sweep_grid_response_cache/`,
stay where they are: they hold the responses, not the loss, and the Monte Carlo
under the corrected loss reads them.

Checked after the move: the moved programs build without warnings, and
`bin/abm_system_mcs`, `bin/fit_metric_report` and
`old/naive_mae_loss/montecarlo/sweep_grid_presence.py`, rerun from their new
location, rewrote `out/old/naive_mae_loss/applications/abm_system_mcs.csv` and
`.txt`, `out/old/naive_mae_loss/montecarlo/fit_metric_report.txt` and `.tex` and
`out/old/naive_mae_loss/montecarlo/sweep_grid_presence.txt` byte for byte. The
other programs were not rerun.
