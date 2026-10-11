# The Model Confidence Set validation procedure

## Overview

How this project decides which ABM configurations cannot be told apart from
the one closest to a benchmark series. A configuration, a "CoP" (Configuration
of Parameters), is one row of the parameter design in
`dataset/abm_system_design.csv`, simulated for 1000 replicates under the seeds
1 to 1000 and named `cop_0001` to `cop_1000`; `docs/ABM_SYSTEM_SIMULATION.md`
describes the design. The same procedure runs in two places, and only the
benchmark differs between them:

- the main pipeline, where the benchmark is the US data
  (`applications/abm_system_irf_loss.c`, `applications/abm_system_mcs.c`), with
  the results in `docs/ABM_SYSTEM_MCS_VALIDATION.md`;
- the Monte Carlo experiment, where every simulated run of every configuration
  is the benchmark in turn, so its configuration is known
  (`montecarlo/sweep_grid.c`), with the results in
  `docs/MONTECARLO_VALIDATION.md`.

In short: fit the t-QVARMA(1,1,2), the linear projection, and the nonlinear projection with combined states to the benchmark and to every
simulated replicate of every configuration, turn each fit into its impulse
response function, measure each replicate's distance from the benchmark's
response by mean absolute error, and run the Model Confidence Set (Hansen,
Lunde and Nason 2011) over the configurations, one model per configuration and
one observation per replicate. The sections below give the reasons and the
settings. Where a section names an `applications/` file,
`montecarlo/sweep_grid.c` does the same thing against each simulated benchmark,
at 2000 bootstrap resamples rather than 10000.

The local-projection variants, `docs/ABM_SYSTEM_LP_VALIDATION.md` on the US
data and the local-projection half of the Monte Carlo grid, use the same loss
and the same confidence set settings with the local projections as the
auxiliary model in place of the t-QVARMA.

## The procedure this replicates

M. A. F. Fabiano, "Evaluating Nonlinear Simulation Models with Model
Confidence Sets" (Pisa / Sant'Anna, thesis), Sec. 3.3, applies this same
idea to the DSK agent-based model, using state-dependent local projections
(LP) as the auxiliary model instead of QVARMA. Its protocol:

1. Fit the auxiliary model to the real data and to every simulated
   Monte Carlo run, identical specification.
2. Compute each fit's own impulse response function.
3. Vectorize and stack the IRF matrices across every horizon into one flat
   vector per fit.
4. Loss = MSE between the real-data IRF vector and each simulated fit's
   own IRF vector.
5. Average the loss across Monte Carlo runs within each CoP.
6. Run the Model Confidence Set (Hansen, Lunde and Nason 2011) over the
   CoPs, using the per-run losses (not the pre-averaged numbers - MCS
   bootstraps a variance from repeated observations, which an
   already-averaged scalar has none of).

Steps 1-6 above are what `applications/abm_system_fit_qvarma.c`,
`applications/abm_system_irf_loss.c` and `applications/abm_system_mcs.c`
implement, with QVARMA standing in for LP. Two points where this project's own
procedure deliberately differs from the thesis:

- **Loss is MAE, not MSE.** See "Absolute error, not squared error" below.
- **The equivalence test's variance is a nonparametric bootstrap of the
  resampled mean** (Hansen, Lunde and Nason's own estimator,
  `MCS_VARIANCE_BOOTSTRAP` in et_al's `mcs.h`), not the thesis's own
  Gaussian quasi-likelihood-based estimate of `L_bar_i` and `sigma_i^2`.
  The thesis's method assumes approximate normality of the per-CoP mean
  loss (a parametric, CLT-style estimate); et_al's estimates the sampling
  distribution directly from the bootstrap draws, no normality assumption.
  Nothing about the comparison object or the loss changes because of this -
  only how the test converts a set of per-replicate losses into a p-value.

Also different from the thesis, structurally rather than as a design
choice: the thesis's own nonlinearity is state dependence (a
smooth-transition LP, expansion vs. recession). QVARMA's own nonlinearity
is score-driven dynamics (a different mechanism entirely - the model is
linear in its state equations, with a t-distributed score driving a
random-walk co-integrating component), so there is no "state" dimension to
run the MCS separately over the way the thesis runs it over
expansion/recession. One MCS over all configurations is the whole comparison
here.

## Why the comparison object is the IRF, not the fitted parameters

An earlier version of `abm_system_irf_loss.c` compared each fit's
`flatten_estimated` constrained-parameter vector directly against the
real-data fit's own, via `stats_mse`/`stats_mae`. This does not answer the
question a validation measure needs answered: two QVARMA fits with
different-looking coefficients can imply nearly identical dynamics, and two
with similar-looking coefficients can imply very different ones. The IRF is
the object whose distance actually measures "do these two models behave
alike" - which is also exactly what the thesis's own protocol computes a
distance between (see "The procedure this replicates" above).

`qvarma.h` already has a working point-estimate IRF: `impulse_responses(m,
D, options)`, with `D = mean_score_jacobian(m, y)`. Both are closed-form -
`mean_score_jacobian` runs the model's own filter once over the data it is
given (real or simulated), no randomness, and `impulse_responses`
differentiates the recursion directly rather than simulating shocked vs.
unshocked paths.

The comparison object is `.total[0..20]` (`contemporaneous + stationary +
cointegrated`, `ImpulseOptions.horizon = 20`), one combined K x K response
matrix per horizon, matching the thesis's own single `IR_h` per horizon rather
than its three separate components. Stacked horizon 0 first into one
length-`K*K*(H+1)` = `5*5*21` = 525 vector (`abm_system_irf_loss.c`'s own
`flatten_total_irf`).

A fit whose own `nu <= 2` cannot have its IRF computed at all -
`impulse_responses` asserts `m->nu > 2`, since `nu` enters the impulse formula
as `1/(nu-2)` - and an IRF can also come back non-finite. Either way the cell is
counted as missing, and a replicate with a missing cell is dropped from every
configuration's column, because et_al's `mcs` requires a loss matrix with no
holes. On the design experiment no cell was missing: all 1,000,000 fits gave a
finite loss.

## Absolute error, not squared error

`stats_mse` was the original choice and produced a pipeline that could not
reject anything. The mechanism, verified directly rather than assumed: a
single badly-fit replicate can produce an IRF distance many orders of
magnitude larger than the rest (observed concretely on the
constrained-parameter version of this pipeline: one replicate's squared
loss was ~3x10^17 against a typical few-hundred-to-few-hundred-thousand
range for the other 107). Squaring that difference lets it dominate not
just the model's own mean loss but the bootstrap variance the MCS
estimates a t-statistic's denominator from - since both numerator and
denominator inflate by roughly the same outlier, the standardized statistic
can stay small even though the raw means look absurdly far apart. Every
model in a 400-model run (constrained-parameter loss) survived under MSE for
exactly this reason. `stats_mae` (mean absolute error) counts the same outlier
linearly rather than squared, which is what let the MCS reject models once
switched - see `docs/ABM_SYSTEM_MCS_VALIDATION.md`, "Last recorded result, on
the older dataset".

## Why the drift-carrying t-QVARMA is not used

An earlier stage of the project also fitted the drift-carrying variant,
t-QVARMAd. Its loss was never switched from parameter distance to IRF distance,
and a joint MCS over both families under squared parameter loss kept every
model, for the reason in the section above. Its code is not part of this
repository, and the MCS here reads only the driftless t-QVARMA losses in
`out/abm_system_irf_loss.csv`.

## Grouping: one model per configuration

A CoP is one ABM parameter configuration. Each contributes one column to
`out/abm_system_irf_loss.csv`, named `cop_NNNN_qvarma_p1q1r2`, and each
row is one replicate. `abm_system_mcs.c` runs one MCS over all columns at once:
1000 models against 1000 observations.

The pipeline fitted two auxiliary specs, p1q1r2 and p1q1r4, while the dataset
was the 10,800 replicates the .Rdata route produced, and the MCS ran jointly
over both so that a p1q1r2 configuration and a p1q1r4 one competed in the same
confidence set.

p1q1r4 lost that comparison. Of the 200 models entered, the three the data
could not separate were all p1q1r2, and no p1q1r4 configuration survived at
all. The joint run is exactly the test of which auxiliary spec belongs in the
experiment, and it answered. That is why the pipeline now fits p1q1r2 alone.
Halving the number of fits is a consequence of the decision, not a reason
for it.

Restoring it is one entry in `spec_list` in
`applications/abm_system_fit_qvarma.c` and one in
`applications/abm_system_irf_loss.c`, plus the join that used to pair the two
loss tables on `replicate`. `abm_system_mcs.c` needs no change: it reads
however many columns the file holds.

## MCS settings

- `stat = MCS_TR` ("range"): every pairwise loss differential, rejects when
  any two models look different from each other. et_al's default, `MCS_TMAX`,
  compares each model with the average of the models still in the set; it is
  run alongside as a check, see `docs/ABM_SYSTEM_MCS_VALIDATION.md`, "Does the
  statistic matter".
- `block_length = 1`: the 1000 replicates are independent Monte Carlo draws,
  not a time series, so the bootstrap is iid (`mcs.h` documents
  `block_length = 1` as the literal iid-bootstrap case, not an approximation to
  one).
- `variance = MCS_VARIANCE_BOOTSTRAP`: Hansen, Lunde and Nason's own estimator
  and et_al's default. See "The procedure this replicates" above for what this
  is instead of a HAC estimate.
- `bootstrap = 10000`, five times et_al's default of 2000, so a round's p-value
  is resolved to 1/10000. The resamples are drawn once and reused in every
  elimination round, which is what the paper and the common implementations do.
- `alpha = 0.05`, seed 123, stream 0.

At 1000 models `MCS_TR` compares 499,500 pairs, and a version of et_al's `mcs`
that stored a quantity per pair per resample could not run on this machine: its
allocations, worked out from their sizes, came to about 11 GB against 7.1 GB of
memory. The current et_al keeps one resampled mean per model and
forms each pair's value from two of them, and the run on the design experiment
takes seconds.

## What the losses are built from

The losses use every one of the 1,000,000 fits, whether or not the solver
reported convergence. The fitting stage was run eight times with the iteration
cap raised from 2,000 to 1,000,000; `applications/abm_system_fit_qvarma.c`
records what each raise bought. Where it ended, per
`out/abm_system_fit_qvarma_manifest.txt` of 2026-09-10:

| why the solver last stopped | fits |
|---|---|
| converged (function decrease below tolerance, or before the reason was stored) | 352,470 |
| the line search could not lower the objective | 621,485 |
| still at the iteration cap after the last run | 5,551 |
| not converged, reason never recorded | 20,494 |

Fits stopped by the line search were not shown to be worse. Measured on the
manifest of 2026-09-09, when 33.8% of fits had converged: within each
configuration, the mean log-likelihood of its line-search-stopped fits minus the
mean log-likelihood of its converged fits had a median of -10.2 over the 1000
configurations (10th percentile -53.2, 90th percentile +27.2), and was positive
in 357 of them. For scale, the standard deviation of the log-likelihood across
the replicates of one configuration averages 156.3.
