# The Model Confidence Set validation procedure

How this project decides which ABM configurations cannot be told apart from
the one closest to a benchmark series. A configuration, a "CoP" (Configuration
of Parameters), is one row of the parameter design in
`dataset/abm_system_design.csv`, simulated for 1000 replicates under the seeds
1 to 1000 and named `cop_0001` to `cop_1000`; `docs/ABM_SYSTEM_SIMULATION.md`
describes the design. The same procedure runs in two places, and only the
benchmark differs between them:

- the empirical validation, where the benchmark is the US data;
- the Monte Carlo experiment, where every simulated run of every configuration
  is the benchmark in turn, so its configuration is known, described in
  `docs/MONTECARLO_VALIDATION.md`.

In short: fit an auxiliary model to the benchmark and to every simulated
replicate of every configuration, turn each fit into its impulse response
function, and measure each replicate's distance from the benchmark's response
by the squared error at each entry of the response minus the configuration's
estimated variance of that entry, averaged over the entries. Then run the
Model Confidence Set of Hansen, Lunde and Nason (2011) over the
configurations, one model per configuration and one observation per
replicate. The auxiliary models are the t-QVARMA(1,1,2) of Blazsek, Escribano
and Licht (2023), the linear local projection, and the state-dependent local
projection with its two states' responses stacked into one vector. The
sections below give the reasons and the settings.

## Status: the empirical validation is behind the Monte Carlo

The Monte Carlo runs this procedure, `montecarlo/true_configuration_recovery.c`.
The empirical validation on the US data does not yet, and will be brought up
to date once the Monte Carlo is finished. Both ran earlier with the naive loss, the mean
absolute error between response vectors with no variance correction; that
code and its outputs are under `old/naive_mae_loss/`, and
`docs/OLD_NAIVE_MAE_LOSS.md` documents it. Beyond the loss, those two runs
differed in this:

- The Monte Carlo ran all three auxiliary models above, plus each state of
  the state-dependent local projection on its own, five response vectors in
  all.
- The empirical validation ran the t-QVARMA as the main pipeline, and the
  linear local projection and each state separately as a side run
  (`docs/OLD_ABM_SYSTEM_LP_VALIDATION.md`). The confidence set over both states
  stacked was not computed on the US data.
- The Monte Carlo used 2000 bootstrap resamples per confidence set, the
  empirical validation 10000.

Everything below holds for both unless it says otherwise.

## The steps

Notation: $d$ is the benchmark, $i$ a configuration, $s = 1, \dots, n$ the
replicates of a configuration that enter the confidence set, and
$h = 1, \dots, K$ the entries of the stacked response vector, every response
at every horizon. On the US data $n = 1000$; in the Monte Carlo $n = 999$,
since the benchmark's run index is held out of every configuration
(`docs/MONTECARLO_VALIDATION.md`, "What is held out, and why").

1. Fit the auxiliary model to the benchmark and to every simulated replicate,
   with the same specification.
2. Compute each fit's impulse response function.
3. Stack every response at every horizon into one vector per fit,
   $\widehat{IR}_{d}$ for the benchmark and $\widehat{IR}_{is}$ for replicate
   $s$ of configuration $i$.
4. For each configuration and entry, estimate the variance of the response
   across the configuration's $n$ replicates,
   $$
   \hat\sigma^2_{ih} = \frac{1}{n-1} \sum_{s=1}^{n} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2,
   \qquad
   \overline{IR}_{ih} = \frac{1}{n} \sum_{s=1}^{n} \widehat{IR}_{ish}.
   $$
   The replicates are exactly the rows of the loss matrix, so in the Monte
   Carlo the benchmark's run is never among them.
5. The loss of a replicate is the squared error minus that variance, averaged
   over the entries:
   $$
   L^*_{is} = \frac{1}{K} \sum_{h=1}^{K} L^*_{ish},
   \qquad
   L^*_{ish} = \left( \widehat{IR}_{ish} - \widehat{IR}_{dh} \right)^2 - \hat\sigma^2_{ih}.
   $$
6. Run the Model Confidence Set over the configurations on the per-replicate
   losses $L^*_{is}$. The losses are not averaged first: the confidence set
   estimates the sampling variance of each configuration's mean loss from its
   repeated observations, and an averaged number has none.

A loss $L^*_{is}$ can be negative. The confidence set reads only differences
between configurations' losses, so this changes nothing. The subtracted term
is the same for every replicate of a configuration, so it moves the
configuration's mean loss and leaves the bootstrap variance of that mean
unchanged; the bootstrap treats $\hat\sigma^2_{ih}$ as known.

`applications/abm_system_fit_qvarma.c` does step 1 for the t-QVARMA and
`applications/abm_system_fit_lp.c` for the local projections. The confidence
set's test statistic is standardised by a nonparametric bootstrap of the
resampled mean loss (Hansen, Lunde and Nason's own estimator,
`MCS_VARIANCE_BOOTSTRAP` in et_al's `mcs.h`), so no normality of the mean
losses is assumed.

## Why the variance correction

Let
$$
L_{ih} = \left( \widehat{IR}_{ih} - \widehat{IR}_{dh} \right)^2
$$
be the squared error between the IRF fitted on the validation data $d$ and
the simulated series of CoP $i$ at horizon $h$, where $h$ runs over every
entry of the stacked response vector. Here no distinction is made between
auxiliary models, since the argument is common to all of them; moreover in
practice there would also be an index for the simulation, but to reduce the
verbosity of the notation it is omitted. The issue is also there with the
mean absolute error, but the math is more tedious.

Then the Model Confidence Set builds the following statistic for two CoPs $i$
and $j$:
$$
d_{ijh} = L_{ih} - L_{jh}.
$$

Now the IRFs are estimates, therefore each of them is a random variable
centred at the true $IR_k$ with variance $\sigma^2_k$, for $k$ a dummy index.
The expectation estimated is
$$
E\left[ d_{ijh} \mid \widehat{IR}_d \right]
= IR_{ih}^2 - IR_{jh}^2 - 2\, \widehat{IR}_{dh} \left( IR_{ih} - IR_{jh} \right) + \sigma^2_{ih} - \sigma^2_{jh},
$$
for which the Model Confidence Set finds an estimator by bootstrap. Note that
since the noise on the real data cannot be removed, it has to be included,
but this is not too much of a problem: assuming the auxiliary model is
consistently estimated, $\widehat{IR}_d$ will asymptotically converge to the
true IRF.

The actual problem is in the last terms, $\sigma^2_{ih} - \sigma^2_{jh}$,
which only relate to the estimation noise of the simulations, not their
vicinity to the calibration data. Therefore, the Model Confidence Set could
favour a model that has a mean IRF that is far from the data's, but has much
less noise, giving a biased confidence set.

To solve the issue, the variance of the IRF is estimated first:
$$
\hat\sigma^2_{ih} = \frac{1}{n-1} \sum_{s=1}^{n} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2,
$$
where the index $s$ denotes the simulation, $i$ the CoP and $h$ the horizon,
and $\overline{IR}_{ih}$ is the sample mean IRF of CoP $i$. Then let
$$
L^*_{ish} = L_{ish} - \hat\sigma^2_{ih},
$$
where $L_{ish}$ is the loss $L_{ih}$ above for simulation $s$. Since
$E[L_{ish} \mid \widehat{IR}_d] = (IR_{ih} - \widehat{IR}_{dh})^2 + \sigma^2_{ih}$
and $\hat\sigma^2_{ih}$ is unbiased for $\sigma^2_{ih}$, the expectation of
$L^*_{ish}$ is the distance term $(IR_{ih} - \widehat{IR}_{dh})^2$ alone. The
rest of the procedure (aggregating the horizons, the Model Confidence Set)
follows the same pattern. Note that in the Monte Carlo setting specifically,
the validation data must not be included in the computation of
$\hat\sigma^2_{ih}$.

The documentation of the naive procedure, with the mean absolute error and no
variance correction, was moved to `docs/OLD_NAIVE_MAE_LOSS.md`.

## Squared error, not absolute error

The variance correction is easier to apply under squared error. There the
estimation noise of the simulated responses enters the expected loss as a
separate additive term, $\sigma^2_{ih}$, and subtracting its estimate removes
it. Under absolute error the noise does not separate from the distance to the
data, so subtracting an estimate of it does not remove it.

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
`old/naive_mae_loss/applications/abm_system_irf_loss.c`, the program that
built the naive loss on the US data).

For the local projections it is the response of each of the five series to a
unit shock in each, recursively identified, at horizons 0 to 15: 5 x 16 x 5 =
400 numbers per model, 800 for the two states stacked.
`docs/MONTECARLO_VALIDATION.md`, "The local projections", gives every setting.

A t-QVARMA fit whose degrees of freedom `nu` are at most 2 has no impulse
response, since `nu` enters it through `1/(nu-2)` and `impulse_responses`
asserts `nu > 2`; a response can also come back non-finite. Either way the cell
is counted as missing, and a replicate with a missing cell is dropped from
every configuration's column, because et_al's `mcs` requires a loss matrix with
no holes; it is then left out of every $\hat\sigma^2_{ih}$ as well. On the
design experiment no cell was missing under any auxiliary model.

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

Putting (1,1,4) back starts with one entry in `spec_list` in
`applications/abm_system_fit_qvarma.c`, which fits it to every simulation.

The drift-carrying variant of the model is not used. It was dropped at an
earlier stage of the project, before the loss moved to impulse responses, and
has not been evaluated under this procedure.

## Grouping: one model per configuration

Each configuration contributes one column to the loss matrix, and each row is
one replicate. One confidence set runs over all columns at once: 1000 models
against $n$ observations, 1000 on the US data and 999 in the Monte Carlo.

## Confidence set settings

- `stat = MCS_TR`, the range statistic: every pairwise difference of mean
  losses, scaled by its bootstrap standard deviation; a round rejects when the
  largest of them is too large. et_al's default, `MCS_TMAX`, compares each
  model with the average of the models still in the set; it was run alongside
  as a check under the naive loss, `docs/OLD_ABM_SYSTEM_MCS_VALIDATION.md`,
  "Does the statistic matter".
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
