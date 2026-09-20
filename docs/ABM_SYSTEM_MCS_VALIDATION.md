# ABM validation via Model Confidence Set over QVARMA impulse responses

## Overview

Selects which of the 1000 simulated ABM parameter configurations cannot be
statistically distinguished, in the dynamics their simulated data implies, from
the configuration closest to the real US data. Each configuration, a "CoP"
(Configuration of Parameters), is one row of the Latin hypercube design in
`dataset/abm_system_design.csv`, and was simulated for 1000 Monte Carlo
replicates under the seeds 1 to 1000; `docs/ABM_SYSTEM_SIMULATION.md` describes
the design and the simulations. The configurations are named `cop_0001` to
`cop_1000` after their design row.

The comparison is done through an auxiliary model's impulse response function,
not through the auxiliary model's own fitted parameters and not through the raw
simulated series directly. The auxiliary model is the driftless t-QVARMA(1,1,2)
(`qvarma.h`). The real-data benchmark is its fit to the US data in
`out/us_qvarma_spec_choice_p1q1r2_fit.json`, written by
`applications/us_qvarma_spec_choice.c`, which also fits p1q1r4 on the real data
because it is the record of why p1q1r2 was the one kept.

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
`applications/abm_system_mse_qvarma.c` and `applications/abm_system_mcs.c`
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

An earlier version of `abm_system_mse_qvarma.c` compared each fit's
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
length-`K*K*(H+1)` = `5*5*21` = 525 vector (`abm_system_mse_qvarma.c`'s own
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
switched - see "Last recorded result, on the older dataset" below.

## Why the drift-carrying t-QVARMA is not used

An earlier stage of the project also fitted the drift-carrying variant,
t-QVARMAd. Its loss was never switched from parameter distance to IRF distance,
and a joint MCS over both families under squared parameter loss kept every
model, for the reason in the section above. Its code is not part of this
repository, and the MCS here reads only the driftless t-QVARMA losses in
`out/abm_system_mse_qvarma_joint.csv`.

## Grouping: one model per configuration

A CoP is one ABM parameter configuration. Each contributes one column to
`out/abm_system_mse_qvarma_joint.csv`, named `cop_NNNN_qvarma_p1q1r2`, and each
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
`applications/abm_system_mse_qvarma.c`, plus the join that used to pair the two
loss tables on `replicate`. `abm_system_mcs.c` needs no change: it reads
however many columns the file holds.

## MCS settings

- `stat = MCS_TR` ("range"): every pairwise loss differential, rejects when
  any two models look different from each other. et_al's default, `MCS_TMAX`,
  compares each model with the average of the models still in the set; it is
  run alongside as a check, see "Does the statistic matter" below.
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
forms each pair's value from two of them, and the run described below takes
seconds.

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

## Result on the design experiment

Setup: `out/abm_system_mse_qvarma_joint.csv` of 2026-09-10, 1000 configurations
by 1000 replicates, no missing cell; loss is the MAE between a fit's stacked
impulse responses and the real-data fit's, horizon 20. MCS with the settings
above. Output: `out/abm_system_mcs_joint.txt` and `out/abm_system_mcs_joint.csv`
of 2026-09-11.

Mean loss per configuration runs from 0.11816 to 0.16329, average 0.13083.

**The confidence set contains one configuration, `cop_0191`.** Each round tests
whether every configuration still in the set has the same expected loss, and
drops the worst one when that is rejected. Every round rejected down to the last
two models; the last round, `cop_0191` against `cop_0148`, rejected with
p = 0.0011, meaning 11 of the 10,000 resamples produced a larger statistic than
the observed one. That p-value is the evidence against `cop_0148` having the same
expected loss as `cop_0191`, and `cop_0148` is the one dropped. With one
configuration left there is nothing to compare it with, and its MCS p-value is 1
by construction.

Hansen, Lunde and Nason's coverage result holds in this case as in any other:
the set contains the configurations with the smallest expected loss with
probability asymptotically at least `1 - alpha`, and when that configuration is
unique the set converges to it. A set of one is what the procedure returns when
the data separate the best configuration from all the others. et_al's `mcs`
reports it with `converged = 0`, which records only that no round was accepted.

The configurations eliminated last, with their losses over the 1000 replicates:

| configuration | mean loss | standard deviation | elimination round | MCS p-value |
|---|---|---|---|---|
| `cop_0191` | 0.11816 | 0.00218 | never eliminated | 1 |
| `cop_0148` | 0.11867 | 0.00424 | 999 | 0.0011 |
| `cop_0931` | 0.11931 | 0.00484 | 998 | 0.0000 |
| `cop_0429` | 0.12051 | 0.00759 | 997 | 0.0000 |
| `cop_0410` | 0.11919 | 0.00248 | 996 | 0.0000 |
| `cop_0717` | 0.11967 | 0.00387 | 995 | 0.0000 |

A p-value of 0.0000 means none of the 10,000 resamples exceeded the observed
statistic. Round 1 is the test with all 1000 configurations; round 999 is the
last.

The rejection of `cop_0148` is consistent with the size of the data rather than
an artefact: the two mean losses differ by 0.00051, and treating replicates as
independent the standard error of that difference is about
`sqrt(0.00218^2 + 0.00424^2) / sqrt(1000)`, roughly 0.00015, so the gap is about
3.4 standard errors. With 1000 replicates per configuration the test separates
mean losses that differ in the fourth decimal.

The two last configurations' design parameters, against the range the design
covers and the model's own baseline calibration in
`model/dsk_sfc/dsk_sfc_inputs.json`:

| parameter | `cop_0191` | `cop_0148` | design minimum | design maximum | baseline |
|---|---|---|---|---|---|
| Gamma | 0.20637 | 0.051444 | 0.050125 | 0.24981 | 0.194 |
| chi | -1.3917 | -1.2676 | -1.4998 | -1.2501 | -1.467 |
| psi1 | 0.4723 | 0.37852 | 0.10032 | 0.49966 | 0.113 |
| psi3 | 0.26683 | 0.34956 | 0.10036 | 0.49998 | 0.444 |
| alfa | 0.39281 | 0.20475 | 0.10025 | 0.49968 | 0.278 |
| taylor1 | 1.3263 | 1.4749 | 1 | 1.4999 | 1.186 |
| taylor2 | 0.45346 | 0.44809 | 0.00042814 | 0.49985 | 0.1 |
| taylor | 0.50534 | 0.52501 | 0.50041 | 0.94968 | 0.777 |
| kappa | 0.9419 | 0.76794 | 0.50026 | 0.94981 | 0.921 |

Both output files carry the round in which each configuration was eliminated:
`elimination_round` in the CSV, `round` in the text report, with 0 for the one
configuration no round eliminated.

### What the confidence set does not say

The Model Confidence Set is a relative procedure. It ranks the 1000
configurations against each other and returns those whose expected loss cannot
be separated from the smallest. It never tests whether the smallest is small,
and nothing else in the pipeline does either, so the result above is "this
configuration is the closest of the thousand" and not "this configuration
matches the US data".

The scale that answers the second question is already computed.
`tests/abm_system_winner_diagnostics.c` reports, in section 4 of
`out/abm_system_winner_diagnostics.txt`, the mean absolute US response over the
525 stacked elements: 0.12228. That is the loss a model whose impulse responses
were identically zero would score, since the loss is the mean absolute
difference against the US vector and a zero vector leaves the US vector itself.

| | mean loss |
|---|---:|
| `cop_0191`, the configuration kept | 0.11816 |
| a model with identically zero impulse responses | 0.12228 |
| median over the 1000 configurations | 0.13000 |
| largest over the 1000 configurations | 0.16329 |

Counted from `out/abm_system_mcs_joint.csv`: 35 of the 1000 configurations
score below 0.12228 and the other 965 score above it. The winner beats the zero
response by 3.4 per cent.

The same section gives the mechanism. The mean absolute response of a
configuration's own fits has median 0.03173 over the 1000 configurations
(5th percentile 0.02025, 95th 0.04620) against the US benchmark's 0.12228, so
the fitted simulated responses are roughly a quarter the size of the US ones
and most of the distance is the US response itself. `cop_0191` has mean
response magnitude 0.03578, rank 688 of 1000 from the smallest. The Spearman
rank correlation between a configuration's mean response magnitude and its mean
loss is 0.0694, so the loss is not simply rewarding small responses either; the
responses are small across the whole design.

The raw series say the same thing without any model in the way.
`out/abm_system_tail_origin_report.txt`, Part 1, compares unconditional moments
over the 400 stored periods against the US series over 187 quarters. The last
column counts configurations whose median over their 1000 replicates reaches
the US value:

| series | moment | US | simulated median | configurations reaching US |
|---|---|---:|---:|---:|
| GDP growth | variance | 0.5713 | 1.691 | 1000 of 1000 |
| employment change | variance | 0.1129 | 1.400 | 1000 of 1000 |
| inflation | variance | 0.6547 | 0.09712 | 0 of 1000 |
| interest rate | variance | 16.33 | 0.02168 | 0 of 1000 |
| GDP growth | skewness | -0.4202 | 0.02578 | 0 of 1000 |
| employment change | skewness | -1.678 | 0.01155 | 0 of 1000 |

GDP growth and employment change are several times too volatile at every design
point, inflation and the interest rate far too smooth, and the skewness of the
US series is not reproduced anywhere in the design. No configuration of the nine
parameters this experiment varies reproduces the second moments of the US
series, let alone the third.

None of this makes the confidence set wrong. It makes it an answer to a
narrower question than the phrase "validation" suggests, and the narrower
question is the one this document reports. An absolute measure would be the
Fabiano validation score in "Open questions" below, which weights the winning
configuration's distance by the benchmark's own estimation uncertainty; it is
not implemented here, and implementing it is what would turn this ranking into
a statement about fit.

### Does the number of resamples matter

The same run at 2000 resamples, everything else identical: the same set of one,
`cop_0191`; the last round's p-value 0.001 instead of 0.0011; every other MCS
p-value 0 in both. The elimination round changed for 918 of the 1000
configurations, by at most 29 rounds, and rounds 991 and 993 to 999 eliminated
the same configuration in both runs. Each round eliminates the configuration
whose loss gap is largest relative to its estimated uncertainty, and that
uncertainty is estimated from the resamples, so configurations with nearly equal
scores trade places when the resamples change.

### Does the statistic matter

`applications/abm_system_mcs_statistic_comparison.c` runs `MCS_TR` and
`MCS_TMAX` on the same loss table with every other setting identical, including
the seed, so both are scored on the same 10,000 resamples. Output:
`out/abm_system_mcs_statistic_comparison.txt` and `.csv` of 2026-09-11.

| | `MCS_TR` | `MCS_TMAX` |
|---|---|---|
| confidence set | `cop_0191` | `cop_0191` |
| last round's p-value | 0.0011 | 0.0011 |
| time inside `mcs` | 11.9 s | 7.9 s |

The two statistics give the same set and the same MCS p-value for every one of
the 1000 configurations. They disagree on the order of elimination: the
Spearman rank correlation between the two orders is 0.3758, and the round a
configuration left in differs by a median of 246 rounds (90th percentile 512,
largest 767, counting the configuration never eliminated as round 1000). The
last three rounds agree.

The largest disagreements are configurations with a high mean loss and a very
noisy one. `MCS_TR` keeps them almost to the end and `MCS_TMAX` drops them
hundreds of rounds earlier:

| configuration | mean loss | standard deviation | `MCS_TR` round | `MCS_TMAX` round |
|---|---|---|---|---|
| `cop_0978` | 0.13473 | 0.04421 | 994 | 898 |
| `cop_0985` | 0.13846 | 0.05031 | 992 | 837 |
| `cop_0532` | 0.14359 | 0.06294 | 991 | 712 |

For scale, the median standard deviation over the 1000 configurations is
0.00648 and its 90th percentile 0.02750. This fits how the two rules eliminate:
`MCS_TR` removes the configuration with the largest noise-scaled gap to some
other one, and a very noisy configuration rarely produces one; `MCS_TMAX`
measures each against the average of those remaining. The noise was measured;
the mechanism was not traced through et_al's code.

## Impulse responses of the winning configuration

Steps 9 and 10 of the main pipeline in the README turn the result into figures:
`applications/abm_system_winner_irf.c` computes the impulse responses of
`cop_0191`, and `applications/abm_system_winner_irf_plots.py` draws them. Both
outputs are of 2026-09-11.

### How the responses were computed

- **Parameters.** The 1000 fitted parameter sets of `cop_0191`, one per
  replicate, averaged coordinate by coordinate in the unconstrained space the
  optimizer works in, then mapped to the model's parameters once. 449 of the
  1000 fits had converged by the solver's test.
- **Score Jacobian.** The impulse responses need the average of the score
  Jacobian over the data. It was computed on each of the 1000 replicates, 400
  periods each, at the averaged parameters, and averaged.
- **At the averaged parameters:** `nu` = 2853.4, and the largest eigenvalue
  modulus of the transitory recursion is 0.396.
- **Horizon:** 0 to 20 periods. Units are percentage points: the growth rates
  and the change in employment are 100 times a log difference or a difference,
  and the interest rate is 100 times its level.
- **Bands.** 10,000,000 random rotations of the shocks, seed 20260822. A
  rotation is kept when the impact responses carry the signs of Table 1 of
  Blazsek, Escribano and Licht (2023): on the supply, demand and monetary policy
  shocks, GDP growth (+, +, -), inflation (-, +, -) and the interest rate
  (unrestricted, +, +). Energy demand growth, employment change and shocks 4
  and 5 carry no restriction. 2,467 rotations were kept, so every percentile is
  taken over 2,467 draws. The line is the median over the kept rotations and the
  band runs from the 10th to the 90th percentile.

`out/abm_system_winner_irf_manifest.txt` records all of the above for the run
that wrote the figures.

### The two sets of figures

The same responses are drawn under two ways of naming the shocks, in two
directories under `out/abm_system_winner_irf_plots/`:

- `sign_restricted/` names a shock by the sign pattern above. Shocks 4 and 5
  satisfy no restriction and have no name. These figures have bands.
- `recursive/` shows the responses without rotation, where shock number b is the
  part of series b's innovation not explained by the series before it in the
  order GDP growth, energy demand growth, employment change, inflation, interest
  rate. Each shock is named for a series, the answer depends on that order, and
  there is no band.

| figure | in | what it shows |
|---|---|---|
| `total_grid.pdf` | both | every response of the five series to all five shocks |
| `shock_supply.pdf`, `shock_demand.pdf`, `shock_monetary.pdf` | `sign_restricted/` | one named shock, its five responses at readable size |
| `shock_gdp_growth.pdf` to `shock_interest_rate.pdf`, five files | `recursive/` | one shock, its five responses at readable size |
| `cumulative.pdf` | both | the response of the level, summed over the horizon, for the four series that are differences |
| `decomposition.pdf` | both | the total response split into its transitory and co-integrated parts |
| `impact.pdf` | both | the responses at horizon 0, where the sign restrictions apply |
| `identification.pdf` | `sign_restricted/` | the median over kept rotations against the unrotated response |

In `sign_restricted/`, `cumulative.pdf`, `decomposition.pdf` and
`identification.pdf` show the three named shocks only; the other figures show
all five.

### What the figures show at impact and after 20 periods

The numbers below are read from `out/abm_system_winner_irf.csv`, sign-restricted
set, as the median with the 10th to 90th percentile in brackets.

At horizon 0:

| response | supply shock | demand shock | monetary policy shock |
|---|---|---|---|
| GDP growth | 0.457 (0.092 to 0.915) | 0.618 (0.144 to 1.073) | -0.093 (-0.232 to -0.017) |
| energy demand growth | 0.443 (0.073 to 0.899) | 0.624 (0.141 to 1.069) | -0.113 (-0.264 to 0.036) |
| employment change | 0.374 (0.036 to 0.812) | 0.593 (0.148 to 0.998) | -0.097 (-0.229 to 0.019) |
| inflation | -0.121 (-0.243 to -0.022) | 0.156 (0.036 to 0.274) | -0.017 (-0.043 to -0.003) |
| interest rate | -0.010 (-0.056 to 0.033) | 0.075 (0.032 to 0.110) | 0.005 (0.001 to 0.012) |

Every restricted entry's band lies on the required side of zero, as it must for
rotations kept on those signs.

The level after 20 periods, the cumulated response at horizon 20:

| level of | supply shock | demand shock | monetary policy shock |
|---|---|---|---|
| GDP | -0.067 (-0.284 to 0.148) | 0.345 (0.116 to 0.505) | 0.025 (-0.071 to 0.112) |
| energy demand | -0.098 (-0.346 to 0.174) | 0.340 (0.038 to 0.522) | 0.001 (-0.177 to 0.185) |
| employment | -0.110 (-0.311 to 0.048) | 0.291 (0.113 to 0.446) | -0.008 (-0.034 to 0.014) |
| prices | -0.423 (-1.175 to 0.163) | 1.089 (0.424 to 1.672) | -0.010 (-0.109 to 0.071) |

After 20 periods only the demand shock leaves all four levels higher with the
whole band above zero. Under the supply and monetary policy shocks every band at
horizon 20 contains zero.

## Last recorded result, on the older dataset

Run of 2026-08-20 over the 100 configurations and 108 replicates the .Rdata
route produced, both specs fitted, 2000 resamples. Those configurations are
named `EstimationSeriesSample1_N` and are not the design experiment's
`cop_NNNN`. **3 of 200 models survive**, all p1q1r2.

| Model | Mean loss | MCS p-value |
|---|---|---|
| EstimationSeriesSample1_91_qvarma_p1q1r2 | 0.11940 | 1.000 |
| EstimationSeriesSample1_60_qvarma_p1q1r2 | 0.11954 | 0.785 |
| EstimationSeriesSample1_15_qvarma_p1q1r2 | 0.11995 | 0.1075 |

Loss across all 200 models ranges 0.119-0.171. No p1q1r4 CoP survives at all,
which is what the grouping section above gives as the reason for dropping that
spec. It stands as the record of the .Rdata run and nothing more.

## Running the pipeline

The README's "The main pipeline" lists all ten steps, from the parameter design
to the figures, with the files each one writes. The steps this document covers
are the last five, in order. Times are from this machine (AMD Ryzen 7 4800H, 16
hardware threads, 7.1 GB of memory):

```
make app-abm_system_fit_qvarma                   # step 6: 1,000,000 fits, resumable, hours per pass
make app-abm_system_mse_qvarma                   # step 7: loss table, parallel over replicates
make app-abm_system_mcs                          # step 8: the MCS, about 14 s of it
make app-abm_system_winner_irf                   # step 9: the winner's impulse responses, about 21 s of it
python applications/abm_system_winner_irf_plots.py   # step 10: the figures, about 35 s
```

and, as a check that is not part of the pipeline:

```
make app-abm_system_mcs_statistic_comparison     # MCS_TR against MCS_TMAX, about 20 s
```

The loss table's time was measured on 4 configurations, 14 s on one thread and
2 s on 16, which projects to minutes for 1000. `make app-abm_system_mcs` and
`make app-abm_system_winner_irf` rebuild the loss table and the real-data fit
first, because the Makefile lists them as prerequisites; the "of it" in the
comments above is the step's own time. `./bin/abm_system_mcs` and
`./bin/abm_system_winner_irf` rerun one step alone on what is already on disk.

## Open questions

- **The fits that did not converge.** 64.8% of the fits did not converge by the
  solver's test, 62.1% of them because the line search could not lower the
  objective. The evidence above is that the line-search-stopped fits'
  likelihoods are not worse than the converged ones', not that their impulse
  responses are.
- **The Fabiano thesis's own validation score**
  (`v = 1/(1+s)`, `s = d_i / sigma^2_rw`, weighting the winning CoP's
  distance by the real-data benchmark's own estimation uncertainty) is not
  implemented here. `abm_system_mcs.c` reports the MCS set and p-values only.
  This is the largest of the open questions rather than one of them: it is the
  only absolute measure in the protocol, and without it the pipeline can say
  which configuration is closest but not whether it is close. See "What the
  confidence set does not say" above for what is already known on that.
- **The benchmark's own sampling error is not in the loss.** Every
  configuration is scored against one estimated US impulse response, treated as
  known. With 1000 replicates the test separates mean losses differing in the
  fourth decimal, which is far below the benchmark's own estimation
  uncertainty, so the singleton set reflects the precision of the simulation
  average and not the precision with which the US dynamics are known.
