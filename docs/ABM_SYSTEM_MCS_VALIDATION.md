# ABM validation via Model Confidence Set over QVARMA impulse responses

## Overview

Selects which of the 100 simulated ABM parameter configurations
(`dataset/abm_system/EstimationSeriesSample1_1` .. `_100`, each a "CoP" -
Configuration of Parameters, the ABM's own structural parameters used to
generate that sample's 108 Monte Carlo replicates) cannot be statistically
distinguished, in the dynamics their simulated data implies, from the
configuration closest to the real US data.

The comparison is done through an auxiliary model's impulse response
function, not through the auxiliary model's own fitted parameters and not
through the raw simulated series directly. The auxiliary model here is the
driftless t-QVARMA (`qvarma.h`); the real-data benchmark is
`applications/us_qvarma_employment_change.c`'s own grid, two specs
(p1q1r2, p1q1r4).

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

Steps 1-6 above are what `applications/abm_system_mse_qvarma.c` and
`applications/abm_system_mcs.c` implement, with QVARMA standing in for LP.
Two points where this project's own procedure deliberately differs from
the thesis, not oversights:

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
expansion/recession. One MCS per spec is the whole comparison here.

## Why the comparison object is the IRF, not the fitted parameters

An earlier version of `abm_system_mse_qvarma.c` compared each fit's
`flatten_estimated` constrained-parameter vector directly against the
real-data fit's own, via `stats_mse`/`stats_mae`. This does not answer the
question a validation measure needs answered: two QVARMA fits with
different-looking coefficients can imply nearly identical dynamics, and two
with similar-looking coefficients can imply very different ones. The IRF is
the object whose distance actually measures "do these two models behave
alike" - which is also exactly what the thesis's own protocol computes a
distance between (see "The procedure this replicates" above), not a
coincidence.

`qvarma.h` already has a working point-estimate IRF: `impulse_responses(m,
D, options)`, with `D = mean_score_jacobian(m, y)`. Both are closed-form -
`mean_score_jacobian` runs the model's own filter once over the data it is
given (real or simulated), no randomness, and `impulse_responses`
differentiates the recursion directly rather than simulating shocked vs.
unshocked paths. `ImpulseOptions` has no field for a sign-restriction
matrix; that mechanism exists in the model's own source paper only for
confidence *bands* around the IRF (`impulse_responses`' own header comment:
"Bands are not implemented ... needs an empirical quantile et_al does not
have"), not for the point estimate this pipeline uses.

The comparison object is `.total[0..20]` (`contemporaneous + stationary +
cointegrated`, `ImpulseOptions.horizon = 20` default), one combined K x K
response matrix per horizon, matching the thesis's own single `IR_h` per
horizon rather than its three separate components. Stacked horizon 0 first
into one length-`K*K*(H+1)` = `5*5*21` = 525 vector
(`abm_system_mse_qvarma.c`'s own `flatten_total_irf`).

A fit whose own `nu <= 2` cannot have its IRF computed at all -
`impulse_responses` asserts `m->nu > 2`, since `nu` enters the impulse
formula as `1/(nu-2)`. A non-converged optimizer run can produce this.
This project's own convention (docs/MODEL_TEMPLATE.md) is that an
infeasible value an optimizer probes returns a sentinel rather than
aborting; a cached fit already on disk is exactly that case arrived at
after the fact, so `try_compute_irf` checks and skips it (counted as
missing, same as an unreadable cache file) rather than leaving it to the
library's own assert. In practice, on the 21,600 qvarma fits this pipeline
has actually been run against, zero had `nu <= 2` - every fit's IRF was
computable, whether or not the fit itself had converged by
`is_converged`'s own criterion.

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
model in a 400-model run (both qvarma and qvarmad, constrained-parameter
loss) survived under MSE for exactly this reason. `stats_mae` (mean
absolute error) counts the same outlier linearly rather than squared, which
is what let the MCS actually reject models once switched - see "Current
result" below.

## Why t-QVARMAd is excluded

`applications/abm_system_mse.c` (t-QVARMAd, constrained-parameter distance,
not switched to IRFs) and this pipeline both write their own loss table,
but `abm_system_mcs.c` reads only `out/abm_system_mse_qvarma_joint.csv`
(qvarma). Not an oversight: both families were run through one joint MCS
once, all 400 columns together, back when both used squared
constrained-parameter loss, and every one of them survived - a result of
the MSE/outlier mechanism above, not a real inability to tell the models
apart. qvarmad's own loss has not been switched to IRF distance, so the two
are not on a comparable scale to rejoin even if it were desired.

## Grouping: 100 CoPs, one MCS per spec

A CoP is one ABM parameter configuration (one sample), not a
(sample, spec) pair. The two auxiliary specs (p1q1r2, p1q1r4) are two
different auxiliary models, each with its own real-data benchmark IRF, so
`abm_system_mse_qvarma.c` computes one loss table per spec and joins them
into a single `out/abm_system_mse_qvarma_joint.csv` (both specs' columns
side by side, `_qvarma_p1q1r2` / `_qvarma_p1q1r4` suffixed). `abm_system_mcs.c`
runs one MCS over all 200 columns of that joint table at once rather than
two separate 100-model runs - a p1q1r2 CoP and a p1q1r4 CoP compete in the
same confidence set, which is a deliberate choice (not something the
thesis's own single-spec protocol has an equivalent of) rather than a
structural necessity; splitting into two 100-model runs is a one-line
change (run `mcs()` on each spec's own columns separately) if the joint
comparison turns out not to be what is wanted.

## MCS settings

- `MCS_TR` ("range": every pairwise loss differential, rejects when any
  two models look different from each other), not et_al's own default
  `MCS_TMAX` (each model against the field average). Chosen to check the
  two statistics did not disagree - on the earlier constrained-parameter
  data they gave the identical result (0 exclusions either way). Not yet
  re-confirmed against `MCS_TMAX` on the IRF/MAE data specifically - see
  "Open questions" below.
- `block_length = 1`: the 108 replicates are independent Monte Carlo draws,
  not a time series, so the bootstrap is iid (`mcs.h` documents
  `block_length = 1` as the literal iid-bootstrap case, not an
  approximation to one).
- `variance = MCS_VARIANCE_BOOTSTRAP`: et_al's own default as of the
  `mcs.h` update that added `MCSVariance` as an option at all. See "The
  procedure this replicates" above for what this is instead of a HAC
  estimate.
- `bootstrap = 2000`: et_al's own default. An earlier run used 20000 (10x)
  for less Monte Carlo noise in the p-value; at `MCS_TR`'s O(M^2) cost per
  round and M0=200, that made a full run take multiple hours. 2000 was
  fast enough (~7 minutes wall clock, this machine) to be practical while
  still being the library's own considered default rather than an
  arbitrary reduction.
- `alpha = 0.05`.

## Current result

Run of 2026-08-20, `out/abm_system_mcs_joint.csv` /
`out/abm_system_mcs_joint.txt`: **3 of 200 models survive**, all p1q1r2.

| Model | Mean loss | MCS p-value |
|---|---|---|
| EstimationSeriesSample1_91_qvarma_p1q1r2 | 0.11940 | 1.000 |
| EstimationSeriesSample1_60_qvarma_p1q1r2 | 0.11954 | 0.785 |
| EstimationSeriesSample1_15_qvarma_p1q1r2 | 0.11995 | 0.1075 |

Loss across all 200 models ranges 0.119-0.171 - tight, no extreme outliers,
unlike every constrained-parameter-distance attempt that preceded this one.
No p1q1r4 CoP survives at all.

## Running the pipeline

In order (the second and third each require the previous step's cache/output
to exist):

```
make app-abm_system_fit_qvarma        # 10,800 simulated fits, ~21,600 cache files
make app-abm_system_mse_qvarma        # IRF computation + loss table, ~5-10 min
make app-abm_system_mcs               # the MCS itself, ~7 min at these settings
```

`applications/us_qvarma_employment_change.c`'s own grid
(`out/us_qvarma_employment_change_p1q1r2_fit.json`, `..._p1q1r4_fit.json`)
must already exist; it is not part of this chain since it does not depend
on the simulated data at all.

## Open questions / not yet done

- **MCS_TMAX vs MCS_TR on the IRF/MAE data.** Confirmed to agree on the
  old, broken (MSE, constrained-parameter) data; not re-run since the loss
  definition changed. Worth checking directly rather than assuming they
  still agree, given the new data actually has models to eliminate this
  time.
- **t-QVARMAd's own loss has not been switched to IRF distance.** If it
  is, the natural next question is whether a joint qvarma-vs-qvarmad MCS
  over IRF distance behaves differently from the all-survive result the
  squared-error version gave.
- **The Fabiano thesis's own validation score**
  (`v = 1/(1+s)`, `s = d_i / sigma^2_rw`, weighting the winning CoP's
  distance by the real-data benchmark's own estimation uncertainty) is not
  implemented here. `abm_system_mcs.c` reports the MCS set and p-values
  only.
