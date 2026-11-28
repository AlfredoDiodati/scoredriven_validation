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

## The states on their own, `old/montecarlo_single_states/`

**What it was.** The Monte Carlo grids of `montecarlo/` (bootstrap or Seri
confidence set, adjusted or naive loss) run on state 1 and state 2 of the
state-dependent local projection on their own, 400 entries each, beside both
states in one vector, 800 entries, under each of the four state weights.

**Why it is out.** Including the states on their own is strictly worse than
both states jointly. Under the adjusted loss, in every grid and under every
state weight, both states in one vector put the true configuration in the set
more often, and give it the smallest mean loss more often, than state 1 or
state 2 alone. Share of the 1,000,000 benchmarks whose set holds the true
configuration, both states against state 1 and state 2, from
`out/montecarlo/identification_table/` on 2026-11-23:

- Bootstrap set, lpirfs' weight: 58.4% against 41.4% and 36.9%.
- Seri set: 24.3% against 11.7% and 10.2% under lpirfs' weight, 27.3% against
  11.5% and 11.3% under the causal HP cycle, 34.3% against 21.2% and 20.9%
  under the estimated lambda, 38.6% against 16.4% and 16.0% under both.

The naive loss is left out of the comparison: it is biased, since without the
variance correction it favours the configurations whose responses are least
noisy over the true one, so its shares (all below 1%) measure that bias rather
than the models.

**Where it is recorded.** `out/old/montecarlo_single_states/`, where
`single_state_rows_split.c` moved the single-state rows of the four results of
`out/montecarlo/` once the Seri grids finished, with a byte-for-byte check that
the two parts are a partition of each result. The grids' results before the new
state weights were added stay whole in `out/montecarlo/before_state_weights/`.

| folder | programs |
| --- | --- |
| `old/montecarlo_single_states/` | the four grids `single_state_true_configuration_recovery.c`, `single_state_uncorrected_true_configuration_recovery.c`, `single_state_seri_true_configuration_recovery.c`, `single_state_seri_uncorrected_true_configuration_recovery.c`; the reports `single_state_seri_true_configuration_recovery_report.c`, `single_state_seri_against_naive_grid.c`, `single_state_seri_against_naive_by_configuration.c`, `single_state_response_variance_by_configuration.c`, `single_state_seri_corrected_against_uncorrected_grid.c`, `single_state_identification_table.c`; `single_state_rows_split.c` |

Each is the program of `montecarlo/` with the same name less the prefix,
restricted to the single-state models (`model_is_single_state` of
`montecarlo/response_cache.h`) and reading and writing
`out/old/montecarlo_single_states/`. The grids build at -O3 with `make
old-single_state_<grid>`; the bootstrap comparisons with `make
old-single_state_corrected_against_uncorrected_grid` and `make
old-single_state_seri_against_bootstrap_grid`. The response caches of the
states stay in `out/montecarlo/sweep_grid_response_cache/`, since both states
in one vector is stacked from them.

## The t-QVARMA at the mean of D_t, `old/montecarlo_mean_jacobian_qvarma/`

**What it was.** The Monte Carlo grids of `montecarlo/` run on the t-QVARMA
with its 525 total responses computed at D, the sample mean of the score
Jacobian D_t over all periods (model `qvarma` of
`montecarlo/response_cache.h`).

**Why it is out.** It is replaced by the t-QVARMA with its two growth states
in one vector (`qvarma_growth_states`), the same fit with its responses at the
mean of D_t over recession periods and over expansion periods, stacked. Under
the Seri set, over the 1,000,000 benchmarks, the two give the same results to
within 0.01 percentage points: true configuration in the set 10.47% against
10.48% under the adjusted loss and 4.41% against 4.41% under the naive loss,
smallest mean loss 6.24% and 2.17% for both, median final set size 2 for both.
Benchmark by benchmark they disagree on whether the true configuration is in
the set in 585 benchmarks under the adjusted loss and 178 under the naive one.
`out/studies/qvarma_growth_state_jacobian_gap_report.txt` says why: the fitted
nu is large (median 695), so every D_t is close to I / (nu - 2) and any mean of
them is nearly the same matrix.

**Where it is recorded.** `out/old/montecarlo_mean_jacobian_qvarma/`, where
`mean_jacobian_qvarma_rows_split.c` moved the `qvarma` rows of the four
results of `out/montecarlo/`, with a byte-for-byte check that the two parts are
a partition of each result. The four results before the split stay whole in
`out/montecarlo/before_mean_jacobian_qvarma_split/`.
`out/montecarlo/run_mean_jacobian_qvarma_split_detached.sh` ran the split,
checked test builds of the grids against the parted results and rewrote the
reports; log in `out/montecarlo/mean_jacobian_qvarma_split/`.

| folder | programs |
| --- | --- |
| `old/montecarlo_mean_jacobian_qvarma/` | the four grids `mean_jacobian_qvarma_true_configuration_recovery.c`, `mean_jacobian_qvarma_uncorrected_true_configuration_recovery.c`, `mean_jacobian_qvarma_seri_true_configuration_recovery.c`, `mean_jacobian_qvarma_seri_uncorrected_true_configuration_recovery.c`; the reports `mean_jacobian_qvarma_seri_true_configuration_recovery_report.c`, `mean_jacobian_qvarma_seri_against_naive_grid.c`, `mean_jacobian_qvarma_seri_against_naive_by_configuration.c`, `mean_jacobian_qvarma_response_variance_by_configuration.c`, `mean_jacobian_qvarma_seri_corrected_against_uncorrected_grid.c`, `mean_jacobian_qvarma_identification_table.c`; `mean_jacobian_qvarma_rows_split.c` |

Each is the program of `montecarlo/` with the same name less the prefix,
restricted to that model (`model_is_mean_jacobian_qvarma`) and reading and
writing `out/old/montecarlo_mean_jacobian_qvarma/`. The grids build at -O3
with `make old-mean_jacobian_qvarma_<grid>`. The model's response cache and the
code computing it stay in `out/montecarlo/sweep_grid_response_cache/qvarma.f32`
and `montecarlo/response_cache.h`, which the programs of `old/naive_mae_loss/`
also read.

## The Seri grid with the mean from all runs, `old/montecarlo_full_sample_mean_seri/`

**What it was.** The adjusted Seri grid: every run of every configuration as
the benchmark, 1,000,000 per model, the loss of replicate s of configuration i

    L*_is = (1/K) sum_h (IR_ish - IR_dh)^2 - n/(n-1) (1/K) sum_h (IR_ish - IR_bar_ih)^2

with IR_bar_ih the mean of entry h over all n kept runs of configuration i,
then seri_mcs at alpha 0.05.

**Why it is out.** IR_bar_ih uses every row of column i, so the correction of
each row depends on all the other rows of the same column and the rows of the
loss matrix are not independent, which seri_mcs assumes. It is replaced by
`montecarlo/seri_held_out_mean_true_configuration_recovery.c`, which takes
IR_bar_ih from the last S = 100 runs of each configuration, forms the losses
and the benchmarks on runs 0 to 899, and scales the correction by S/(S+1)
(docs/MONTECARLO_VALIDATION.md). Over the 900,000 benchmarks both grids hold
(runs 0 to 899), the share of sets holding the true configuration was 14.7%,
38.6% and 10.5% for the linear projection, both states under the causal HP
cycle at the estimated lambda and the t-QVARMA with its two growth states,
against 15.3%, 40.2% and 11.1% with the held-out mean; the mean set size 2.80,
27.78 and 2.60 against 3.02, 30.49 and 3.26.

**Where it is recorded.** `out/old/montecarlo_full_sample_mean_seri/`: the
result `seri_true_configuration_recovery.csv.gz`, its log, its trial, its
report, the recovery table it gave in `identification_table/`, and
`seri_full_sample_mean_against_held_out_mean_grid_report.txt`, the benchmark by
benchmark comparison with the held-out grid
(`make old-seri_full_sample_mean_against_held_out_mean_grid`).

| folder | programs |
| --- | --- |
| `old/montecarlo_full_sample_mean_seri/` | the grid `full_sample_mean_seri_true_configuration_recovery.c` (`make old-full_sample_mean_seri_true_configuration_recovery`, at -O3), its report `full_sample_mean_seri_true_configuration_recovery_report.c`, the table `full_sample_mean_identification_table.c` |

Each is the program of `montecarlo/` with the same name less the prefix,
reading and writing `out/old/montecarlo_full_sample_mean_seri/`.

## The local projections with the co-integrated pair in first differences, `old/montecarlo_differenced_lp/`

**What it was.** The linear local projection and both states under the causal
HP cycle at the estimated lambda fitted with inflation and the interest rate,
the last two of the t-QVARMA's five series (K_star = 3), which it treats as
I(1) and co-integrated, in first differences, the other three series as
stored. The difference at the first fitted period (stored period 3) uses
stored period 2, which the 4-period moving average of the state series already
drops, so the fit sees the same 397 periods as the undifferenced one; the state
series, and so the estimated lambda and the state weight, are unchanged, only
the regressions change (applications/lp_system.h,
lp_system_from_stored_differenced). Both ran through the two Seri grids,
seri_mcs at alpha 0.05 over the 1000 configurations: unadjusted (squared error,
every run of every configuration a benchmark, 1,000,000 per model) and
adjusted (own deviation from the mean of runs 900 to 999 subtracted, scaled by
S/(S+1) with S = 100, runs 0 to 899 the benchmarks, 900,000 per model).

**Why it is out.** It recovers the true configuration worse than the same
projection on the stored series. Share of benchmarks whose set holds the true
configuration, share where it has the smallest mean loss, and mean set size,
stored against differenced:

- Unadjusted, linear: 7.30% against 0.15%; 4.15% against 0.10%; 2.87 against
  15.33.
- Unadjusted, both states: 2.73% against 0.10%; 1.97% against 0.10%; 3.58
  against 8.66.
- Adjusted, linear: 15.26% against 27.59%, but the median set size goes from 2
  to 120 of the 1000 configurations; smallest mean loss 8.50% against 0.25%.
- Adjusted, both states: 40.23% against 30.87%; 6.78% against 0.46%; median
  set size 20 against 58.

**Where it is recorded.** `out/old/montecarlo_differenced_lp/`: the rows of
the two differenced models of both Seri results, moved there by
`differenced_lp_rows_split.c` with a check that the two parts are a partition
of each result, and `differenced_lp_against_stored_seri_report.txt`, which has
the numbers above and the paired counts. The two results before the
differenced models were added stay whole in
`out/montecarlo/before_differenced/`. `run_differenced_lp_seri_grids_detached.sh`
and its log are the fits, caches, test and grid runs;
`run_differenced_lp_split_detached.sh` and its log the split.

| folder | programs |
| --- | --- |
| `old/montecarlo_differenced_lp/` | the grids `differenced_lp_seri_held_out_mean_true_configuration_recovery.c` and `differenced_lp_seri_uncorrected_true_configuration_recovery.c` (`make old-differenced_lp_<grid>`, at -O3); the comparison `differenced_lp_against_stored_seri_report.c`; the test `lp_system_differenced_correctness.c`; `differenced_lp_rows_split.c` |

Each grid is the program of `montecarlo/` with the same name less the prefix,
restricted to the differenced models (`model_in_differenced_lp_grids` of
`montecarlo/response_cache.h`) and reading and writing
`out/old/montecarlo_differenced_lp/`. The fits stay in
`out/applications/abm_system_fit_lp_differenced/` and
`out/applications/abm_system_fit_lp_differenced_causal_hp_estimated_lambda/`,
from `bin/abm_system_fit_lp --difference-cointegrated` (`make
old-abm_system_fit_lp_differenced` and `make
old-abm_system_fit_lp_differenced_causal_hp_estimated_lambda`); the four
models, 15 to 18, and their caches stay in `montecarlo/response_cache.h` and
`out/montecarlo/sweep_grid_response_cache/`.
