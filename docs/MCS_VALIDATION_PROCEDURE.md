# The Model Confidence Set validation procedure

How this project decides which ABM configurations cannot be told apart from
the one closest to a benchmark series. A configuration, a "CoP" (Configuration
of Parameters), is one row of the parameter design in
`dataset/abm_system_design.csv`, simulated for 1000 replicates under the seeds
1 to 1000 and named `cop_0001` to `cop_1000`; `docs/ABM_SYSTEM_SIMULATION.md`
describes the design. The same procedure runs in two places, and only the
benchmark differs between them:

- the empirical validation, where the benchmark is the US data
  (`applications/abm_system_irf_loss.c`, `applications/abm_system_mcs.c` for
  the t-QVARMA, `applications/abm_system_lp_irf_loss.c`,
  `applications/abm_system_lp_mcs.c` for the local projections), with the
  results in `docs/ABM_SYSTEM_MCS_VALIDATION.md` and
  `docs/ABM_SYSTEM_LP_VALIDATION.md`;
- the Monte Carlo experiment, where every simulated run of every configuration
  is the benchmark in turn, so its configuration is known
  (`montecarlo/sweep_grid.c`), with the results in
  `docs/MONTECARLO_VALIDATION.md`.

In short: fit an auxiliary model to the benchmark and to every simulated
replicate of every configuration, turn each fit into its impulse response
function, measure each replicate's distance from the benchmark's response by
mean absolute error, and run the Model Confidence Set of Hansen, Lunde and
Nason (2011) over the configurations, one model per configuration and one
observation per replicate. The auxiliary models are the t-QVARMA(1,1,2) of
Blazsek, Escribano and Licht (2023), the linear local projection, and the
state-dependent local projection with its two states' responses stacked into
one vector. The sections below give the reasons and the settings.

## Status: the empirical validation is behind the Monte Carlo

The procedure is being refined in the Monte Carlo study, and the empirical
validation on the US data will be brought up to date once that is finished.
Until then the two differ:

- The Monte Carlo runs all three auxiliary models above, plus each state of
  the state-dependent local projection on its own, five response vectors in
  all.
- The empirical validation runs the t-QVARMA as the main pipeline, and the
  linear local projection and each state separately as a side run
  (`docs/ABM_SYSTEM_LP_VALIDATION.md`). The confidence set over both states
  stacked is not computed on the US data.
- The Monte Carlo uses 2000 bootstrap resamples per confidence set, the
  empirical validation 10000.

Everything below holds for both unless it says otherwise.

## The steps

1. Fit the auxiliary model to the benchmark and to every simulated replicate,
   with the same specification.
2. Compute each fit's impulse response function.
3. Stack every response at every horizon into one vector per fit.
4. The loss of a replicate is the mean absolute difference between its vector
   and the benchmark's.
5. Run the Model Confidence Set over the configurations on the per-replicate
   losses. The losses are not averaged first: the confidence set estimates the
   sampling variance of each configuration's mean loss from its repeated
   observations, and an averaged number has none.

`applications/abm_system_fit_qvarma.c` does step 1 for the t-QVARMA and
`applications/abm_system_fit_lp.c` for the local projections. The confidence
set's test statistic is standardised by a nonparametric bootstrap of the
resampled mean loss (Hansen, Lunde and Nason's own estimator,
`MCS_VARIANCE_BOOTSTRAP` in et_al's `mcs.h`), so no normality of the mean
losses is assumed.

## Why the comparison object is the impulse response, not the fitted parameters

Two t-QVARMA fits with different coefficients can imply nearly identical
dynamics, and two with similar coefficients can imply very different ones. The
coordinates of the parameter vector also differ in units and scale by orders
of magnitude, so a distance between parameter vectors is mostly a distance in
whichever coordinate is largest. The impulse response is what says whether two
fitted models behave alike, and it is in the units of the data.

For the t-QVARMA the response is et_al's point estimate,
`impulse_responses(m, D, options)` with `D = mean_score_jacobian(m, y)`.
`mean_score_jacobian` runs the model's filter once over the series it is
given, real or simulated, and `impulse_responses` differentiates the
recursion rather than simulating shocked and unshocked paths, so neither
involves randomness.

The compared vector is the total response, `.total` (the contemporaneous,
stationary and co-integrated parts summed), at horizons 0 to 20
(`ImpulseOptions.horizon = 20`): one 5 by 5 matrix per horizon, stacked horizon
0 first into 5 x 5 x 21 = 525 numbers (`flatten_total_irf` in
`applications/abm_system_irf_loss.c`).

For the local projections it is the response of each of the five series to a
unit shock in each, recursively identified, at horizons 0 to 15: 5 x 16 x 5 =
400 numbers per model, 800 for the two states stacked.
`docs/MONTECARLO_VALIDATION.md`, "The local projections", gives every setting.

A t-QVARMA fit whose degrees of freedom `nu` are at most 2 has no impulse
response, since `nu` enters it through `1/(nu-2)` and `impulse_responses`
asserts `nu > 2`; a response can also come back non-finite. Either way the cell
is counted as missing, and a replicate with a missing cell is dropped from
every configuration's column, because et_al's `mcs` requires a loss matrix with
no holes. On the design experiment no cell was missing under any auxiliary
model.

## Absolute error, not squared error

A squared loss lets one replicate with an extreme distance dominate both its
configuration's mean loss and the bootstrap estimate of that mean's variance.
The confidence set's statistic divides the first by the square root of the
second, so both grow together, the statistic stays small, and nothing is
eliminated. The mean absolute error counts the same replicate linearly.

The evidence is from an earlier version of this pipeline, run on the older
dataset of 100 configurations and 108 replicates each, with a loss computed on
the fitted parameters rather than on impulse responses: one replicate's squared
loss was about 3 x 10^17, against a few hundred to a few hundred thousand for
the other 107 replicates of its configuration, and every one of the 400 models
entered survived the confidence set. With the absolute loss, the confidence
set on that dataset kept 3 of the 200 t-QVARMA models entered (run of
2026-08-20, 2000 resamples).

## Auxiliary specification

The t-QVARMA is used without drift, with one lag of the stationary component
(p = 1), one lag of the score in that component (q = 1) and two lags of the
score in the co-integrated component (r = 2): t-QVARMA(1,1,2), written
`p1q1r2` in file and column names. Of the five series (GDP growth, energy
demand growth, employment change, inflation, interest rate), the first three
are treated as stationary and inflation and the interest rate as co-integrated
with rank 1, the Fisher relation; `applications/us_qvarma_spec_choice.c` sets
this partition.

(1,1,2) was chosen over (1,1,4) by entering both in one confidence set, on the
older dataset above, as 200 models: the three the data could not separate were
all (1,1,2), and no (1,1,4) configuration survived.
`applications/us_qvarma_spec_choice.c` fits both to the US data and records
their residual checks in `out/us_qvarma_spec_choice.txt`.

Putting (1,1,4) back is one entry in `spec_list` in
`applications/abm_system_fit_qvarma.c` and one in
`applications/abm_system_irf_loss.c`, plus a join of the two loss tables on
`replicate`. `applications/abm_system_mcs.c` needs no change: it reads however
many columns the file holds.

The drift-carrying variant of the model is not used. It was dropped at an
earlier stage of the project, before the loss moved to impulse responses, and
has not been evaluated under this procedure.

## Grouping: one model per configuration

Each configuration contributes one column to the loss table, for the t-QVARMA
`out/abm_system_irf_loss.csv` with columns named `cop_NNNN_qvarma_p1q1r2`, and
each row is one replicate. One confidence set runs over all columns at once:
1000 models against 1000 observations.

## Confidence set settings

- `stat = MCS_TR`, the range statistic: every pairwise difference of mean
  losses, scaled by its bootstrap standard deviation; a round rejects when the
  largest of them is too large. et_al's default, `MCS_TMAX`, compares each
  model with the average of the models still in the set; it is run alongside
  as a check, `docs/ABM_SYSTEM_MCS_VALIDATION.md`, "Does the statistic
  matter".
- `block_length = 1`: the replicates are independent simulations, not a time
  series, so the bootstrap resamples single replicates.
- `variance = MCS_VARIANCE_BOOTSTRAP`, as above.
- `bootstrap = 10000` in the empirical validation, five times et_al's default
  of 2000, so a round's p-value is resolved to 1/10000; 2000 in the Monte
  Carlo. The resamples are drawn once and reused in every elimination round.
- `alpha = 0.05`, seed 123, stream 0.

At 1000 models `MCS_TR` compares 499,500 pairs. et_al's `mcs` keeps one
resampled mean per model and forms each pair's difference from two of them,
and one confidence set over the design takes seconds.

## What the t-QVARMA losses are built from

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

The local projections are least squares and have no convergence status.
