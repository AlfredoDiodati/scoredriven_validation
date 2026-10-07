# Monte Carlo validation with local projections

## The question

`docs/MONTECARLO_VALIDATION.md` asks whether the validation procedure finds a
configuration it already knows the answer for: one simulated run of the ABM
takes the place of the US data, and the Model Confidence Set should come back
holding the configuration that run was simulated from. There the auxiliary
model is the t-QVARMA. Here it is the local projections of the collaborator's R
pipeline (`_temp/main_code.R`, `_temp/functions.R`), everything else unchanged,
so the two auxiliary models are asked the same question on the same data.

The answer: fitted on the same stationary series as the t-QVARMA, the local
projections find the right configuration about as often as the t-QVARMA does.
Fitted on the levels R's pipeline builds, they almost never do. The first
version of this experiment used R's levels; that run and the investigation
that traced its failure to the levels are recorded below, after the result.

## What is fitted

Every one of the 1,000,000 simulated runs (1000 configurations, 1000 runs each,
`dataset/abm_system/`) goes through two steps.

The series are the five the t-QVARMA is fitted on, as stored: GDP growth,
energy growth, employment change, inflation and the interest rate, in that
order (`lp_system_from_stored` in `applications/lp_system.h`). The series that
decides the state is the one R's pipeline builds, the 4-period moving average
of `100 log(100 + cumsum(GDP_growth))`, and the first 3 of the 400 periods are
dropped, where that moving average does not exist: 397 periods.
`tests/lp_system_transform.c` checks that the fitted series are the stored ones
exactly and the state series is R's exactly.

Then et_al's `lp_lin_and_nl` fits both models at R's settings: 4 lags, 15
horizons, a unit shock with recursive ordering in the order above; the
state-dependent model with the Hodrick-Prescott cycle of the moving average at
lambda 1600, a logistic weight of slope 2, lagged one period. Each model's
impulse responses are a 5 by 16 by 5 array, 400 numbers: the linear model's,
and the state-dependent model's in state 1 (periods weighted towards GDP above
trend) and state 2 (below).

The state-dependent model is also compared as one vector, following
`_temp/Note on non-lin LP.pdf`: state 1's 400 responses followed by state 2's,
800 numbers, called "both states" below (`lp_nl` in the code and in the file
names, `nl` in the `lp_*` files). Its loss, the mean absolute difference over
the 800 entries, is exactly the average of the two states' losses, since the
two halves have the same length; measured on the single-benchmark tables
below, the largest relative difference between the two is 4.8e-15.

`applications/abm_system_fit_lp.c` writes the fits to one compressed archive
per configuration, `out/abm_system_fit_lp/cop_NNNN.npz`, one row per run, with
the responses, the shock matrix, the fit diagnostics, the settings and a
checksum of the data each row was fitted on; 8.6 GB in all. No run failed:
every fit has responses, every regression at every horizon and the VAR behind
the shocks were solved by ordinary least squares at full rank (none used the
minimum-norm solution et_al falls back on), and no state series was constant
(`out/abm_system_fit_lp_manifest.txt`). Fitting all 1,000,000 runs from scratch
takes 402 seconds on this machine (AMD Ryzen 7 4800H, 16 threads); the rerun
that replaced the R-level fits took 1680 seconds, since it first read and
refused every stored archive.

## One simulated run as the real data

The run standing in for the US data is the same as in the t-QVARMA experiment,
`cop_0191` run 706 (`montecarlo/out/benchmark.env`). Run 706 of every
configuration is left out, because all configurations share the same random
seeds. The loss of a run is the mean absolute difference between its 400
response values and the stand-in's, separately for the linear model, state 1
and state 2, and over 800 values for both states. Each loss table is 999 runs
by 1000 configurations; no cell was missing.

The confidence set is computed exactly as in the t-QVARMA experiment: level
0.05, 10000 bootstrap resamples of single runs, the bootstrap variance of the
resampled mean, et_al's seed 123, statistic MCS_TR (every pair of
configurations compared).

Result (`montecarlo/out/lp_mcs_<model>.txt` and `.csv`):

- linear: the set holds `cop_0191` alone, which has the smallest mean loss
  (0.689).
- state 1: the set holds `cop_0191` alone, smallest mean loss (1.336).
- state 2: the set holds `cop_0437`, `cop_0832` and `cop_0191`; `cop_0191` is
  third by mean loss (1.431 against 1.408 for `cop_0437`), and the set stopped
  on an accepted test at p = 0.0773.
- both states: the set holds `cop_0191` alone, smallest mean loss (1.383,
  next `cop_0437` at 1.425), under MCS_TR and MCS_TMAX.

The t-QVARMA experiment, same stand-in run, same settings, keeps `cop_0191`
alone.

## Every run of cop_0191 as the real data

One stand-in run cannot separate a procedure that finds the right configuration
from one that was lucky once. `montecarlo/lp_sweep.c` repeats the whole
validation 1000 times, each time with a different run of `cop_0191` (runs 0 to
999) as the stand-in and that run's index left out of every configuration,
under the same settings as above, MCS_TR only. It counts how often the set
holds `cop_0191`.

- linear: `cop_0191` in the set in 913 of 1000 repetitions (91.3%), alone in
  864; smallest mean loss in 886.
- state 1: in the set in 858 of 1000 (85.8%), alone in 645; smallest mean
  loss in 758.
- state 2: in the set in 781 of 1000 (78.1%), alone in 398; smallest mean
  loss in 592. The median set holds 2 configurations, the largest 16.
- both states: in the set in 914 of 1000 (91.4%), alone in 668; smallest mean
  loss in 805. The median set holds 1 configuration, the largest 9.

The median rank of `cop_0191` by mean loss is 1 under all four models. No run
was dropped for a missing value in any repetition.

For comparison, the t-QVARMA version of the same 1000 repetitions
(`montecarlo/out/sweep_irf.csv`) has `cop_0191` in the set in 876 of 1000
(87.6%).

Per-repetition rows: `montecarlo/out/lp_sweep_lin.csv`, `lp_sweep_s1.csv`,
`lp_sweep_s2.csv`, `lp_sweep_nl.csv`, one row per stand-in run: whether
`cop_0191` is in the set, its rank by mean loss, its p-value, the size of the
set, whether the set was decided by an accepted test, and how many runs were
dropped. Which configurations shared the set is not recorded.

The whole chain, refit, single stand-in and 1000 repetitions of the three
original models, took 12.2 hours, one confidence set at a time. The 1000
repetitions under both states were added later by `montecarlo/lp_nl_run.sh`,
16 confidence sets side by side; their time was not kept.

The same protocol with replicate 0 of each of the 1000 configurations as the
stand-in, instead of every replicate of `cop_0191`, is
`montecarlo/sweep_cops.c`; `docs/MONTECARLO_VALIDATION.md`, "The first run of
every configuration as the benchmark", describes it for both auxiliary models.

## The first run, on R's levels

The experiment was first run on the series R's `transform_data()` builds
(`lp_system_r_transform` in `applications/lp_system.h`): log GDP as
`100 log(100 + cumsum(GDP_growth))`, employment as
`0.90 + cumsum(Employment_change)`, inflation and the interest rate unchanged,
log energy like GDP, in R's column order (GDP, employment, inflation, interest
rate, energy). That transformation was checked against R's own
`transform_data()`, run unmodified (with `zoo` and `dplyr` installed) on four
stored runs (`cop_0191` run 706, `cop_0001` run 0, `cop_0532` run 123,
`cop_1000` run 999): the largest difference, relative to the larger of 1 and
the value, was 3.7e-15. Every other setting was as above. Its outputs are kept
in `montecarlo/out/lp_r_levels/`.

- Stand-in run 706: the linear and state-2 sets held `cop_0905` alone, with
  `cop_0191` 26th and 23rd by mean loss; the state-1 set held 7
  configurations, `cop_0191` first among them.
- 1000 repetitions: `cop_0191` in the set in 5 (linear), 72 (state 1) and 42
  (state 2) of 1000.

Those 1000 repetitions first took about 3 minutes each: one confidence set on
these tables took 32 to 131 seconds, against 11 to 19 for the t-QVARMA's.
et_al's `mcs` was then changed to compute a tighter bound when deciding which
configurations it can leave out of a bootstrap draw's comparisons; the change
leaves every result byte-identical and brought one confidence set to 5 to 7
seconds.

## Why the levels fail

`studies/lp_recovery_diagnostics.c` refits all 1,000,000 runs under several
variants and checks each candidate cause. Settings: stand-in `cop_0191` run
706, MCS_TR at 0.05 with 10000 resamples, and for the linear model runs 0 to 99
of `cop_0191` each taking the stand-in's place. Numbers in
`out/lp_recovery_diagnostics_report.txt` (R's levels and the variants of them)
and `out/lp_recovery_diagnostics_growth_unit_report.txt` (the stored series).

What was ruled out, all on R's levels:

- A computational error: refitting every run reproduces the stored loss tables
  to 5.4e-13 relative.
- The shock size: with shocks of one standard deviation, `cop_0191` is in the
  linear model's set 0 times in 100.
- R's form of the log level: with log GDP and energy as the cumulated growth
  rates themselves, 0 of 100, with either shock size.
- The loss: with R's mean squared difference, the stand-in run 706 set is
  again `cop_0905` alone.

What changes the outcome is the series: on the stored stationary series the
same linear model finds `cop_0191` in 90 of 100, the result the full run above
confirms at 913 of 1000.

Why, measured on the linear model at stand-in run 706 with
`studies/recovery_noise.h`, the same code for every fit. For each response
value, the variance of the configurations' average responses is divided by
the variance of one configuration's runs around its own average:

- the median is nearly the same for every fit: 0.059 on R's levels, 0.048 on
  the stored series, 0.0575 for the t-QVARMA
  (`out/qvarma_recovery_diagnostics_report.txt`);
- the 90th percentile separates them: 0.69 on R's levels, 4.2 on the stored
  series, 3.97 for the t-QVARMA.

On the stationary series about a tenth of the response values differ clearly
between configurations; on R's levels none do, and the confidence set then
keeps the configuration whose estimates scatter least. A configuration's mean
loss has Spearman correlation 0.999 with the scatter of its own runs on R's
levels, against 0.970 on the stored series and 0.531 for the t-QVARMA; with its
average's distance from the stand-in, 0.936, 0.985 and 0.887.

The t-QVARMA's estimates are therefore not more precise than the local
projections', value by value; the hypothesis that its success came from
precision alone is not supported. The protocol does reward precise estimates, but it
recovers the right configuration as long as some responses separate
configurations.

Under R's levels with unit shocks, the interest-rate shock carried 97.6% of
`cop_0191`'s loss; on the stored series, 55%. The interest rate is the least
volatile series, its volatility varies 37.6 times across configurations, and
`cop_0905`, the levels' winner, has the most volatile interest rate of the
1000, hence the smallest responses to a unit shock to it.

Not explained: an all-zero response is closer to the stand-in than 895 of the
1000 configurations' mean losses for the t-QVARMA, and than none of them for
the local projections on either kind of series.

## What this does not show

The fits no longer match the collaborator's R pipeline, which fits the levels.
A real-data comparison against the collaborator's own real-data responses would
need those rebuilt on the stored series too.

The comparison with R is on the transformation only. The local projections
themselves were compared with `lpirfs` in et_al (`docs/LP_DOCUMENTATION.md`),
on synthetic data, not on this dataset.

The state-dependent model's state still comes from R's log-GDP level; whether a
state built from the growth rates would do better was not tested.
