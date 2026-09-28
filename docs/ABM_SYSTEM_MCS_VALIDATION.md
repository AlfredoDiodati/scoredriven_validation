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
switched - see "Last recorded result, on the older dataset" below.

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

Setup: `out/abm_system_irf_loss.csv` of 2026-09-10, 1000 configurations
by 1000 replicates, no missing cell; loss is the MAE between a fit's stacked
impulse responses and the real-data fit's, horizon 20. MCS with the settings
above. Output: `out/abm_system_mcs.txt` and `out/abm_system_mcs.csv`
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
`studies/abm_system_winner_diagnostics.c` reports, in section 4 of
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

Counted from `out/abm_system_mcs.csv`: 35 of the 1000 configurations
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

These are moments of runs this project's build of the simulator produced, so
the reading holds only if that build is the authors'. It is:
`tests/dsk_tail_replicate_reproduction.c` takes the replicates with the
heaviest tails in the whole experiment and reruns them in the authors'
unchanged build at the same configuration and seed, and
`out/dsk_tail_replicate_reproduction.txt` records sixteen such runs reproduced
exactly, results file, error log and stored series alike. The thinness above is
the model's.

### Is the thinness a central limit theorem

GDP here is not drawn from a distribution. It is what comes out after adding up
what 200 consumption-good firms and 20 capital-good firms each did in a
quarter, and adding up many roughly independent quantities is the situation a
central limit theorem describes. Its conclusion is a Gaussian, which has no
large rare values, so it would account for the thinness above.
`studies/abm_system_gaussian_convergence.c` asks whether the departure from a
Gaussian in a run's own time series shrinks as the sample grows.

Nothing is pooled across replicates. The system is not ergodic, so a sample
taken across replicates at a fixed date is a mixture over economies in
different states and says nothing about the aggregation inside any one of them.
Every sample here is the first n quarters of a single run.

Setup. All 1000 configurations, all 1,000,000 replicates, n = 100, 150, 200,
250, 300, 350, 400 out of each run's 400 stored quarters. The five stored
series plus the first difference of the interest rate, which has 399 values and
so stops at n = 350. Four quantities per sample: skewness, excess kurtosis
m4 / m2^2 - 3, Anderson-Darling and Jarque-Bera. Neither test statistic is read
off its asymptotic distribution: at each n the null of both is simulated from
100,000 Gaussian samples of that size, seed 20260923, and an observed statistic
becomes a p-value by where it falls in that null. Sweep 3.5 minutes on 16
threads.

Excess kurtosis, averaged over the million runs. The last row is what a
Gaussian sample of that n gives, which is not zero because the estimate is
biased downward at these sample sizes:

| series | n=100 | 150 | 200 | 250 | 300 | 350 | 400 |
|---|---:|---:|---:|---:|---:|---:|---:|
| GDP growth | 0.650 | 0.814 | 0.921 | 0.998 | 1.058 | 1.109 | 1.151 |
| energy growth | 0.694 | 0.864 | 0.973 | 1.050 | 1.110 | 1.158 | 1.198 |
| employment change | 0.666 | 0.829 | 0.932 | 1.004 | 1.060 | 1.106 | 1.143 |
| inflation | 0.719 | 0.885 | 0.985 | 1.054 | 1.108 | 1.150 | 1.184 |
| interest rate | -0.202 | -0.091 | -0.024 | 0.023 | 0.058 | 0.085 | 0.107 |
| interest rate change | 0.041 | 0.139 | 0.196 | 0.235 | 0.263 | 0.285 | |
| Gaussian sample of the same n | -0.057 | -0.038 | -0.031 | -0.022 | -0.020 | -0.017 | -0.014 |

Share of the million runs Anderson-Darling rejects at 5 per cent, which is 0.05
under the null:

| series | n=100 | 150 | 200 | 250 | 300 | 350 | 400 |
|---|---:|---:|---:|---:|---:|---:|---:|
| GDP growth | 0.150 | 0.187 | 0.225 | 0.257 | 0.291 | 0.319 | 0.347 |
| energy growth | 0.155 | 0.193 | 0.233 | 0.264 | 0.298 | 0.326 | 0.353 |
| employment change | 0.153 | 0.192 | 0.231 | 0.263 | 0.296 | 0.323 | 0.350 |
| inflation | 0.160 | 0.200 | 0.240 | 0.273 | 0.307 | 0.335 | 0.362 |
| interest rate | 0.418 | 0.432 | 0.447 | 0.454 | 0.464 | 0.470 | 0.477 |
| interest rate change | 0.150 | 0.164 | 0.179 | 0.190 | 0.204 | 0.213 | |

Skewness stays small: below 0.02 in absolute value for the four growth series
at every n, and -0.124 for inflation at n = 400. Jarque-Bera rejects more often
than Anderson-Darling in every series except the interest rate level.

The departure grows with n in every series. Excess kurtosis rises rather than
settling, and the rejection share rises with it. There is no n at which these
series look Gaussian, and none at which they are closer to Gaussian than at a
smaller n. Whatever accounts for the thin tails, it is not that the aggregation
delivers a Gaussian.

`out/abm_system_gaussian_convergence_by_sample_size.csv` holds every column
including Jarque-Bera and both null columns;
`out/abm_system_gaussian_convergence_report.txt` is the run's own report.

Two things the study does not control for. A run's quarters are dependent while
the simulated null draws them independently, so a rejection can come from the
dependence rather than from a non-Gaussian shape. And a larger n reaches later
into a run, so the sample size and the stretch of history it covers grow
together.

Anderson-Darling, Jarque-Bera, skewness and excess kurtosis are general
statistics and none of them is in et_al's `stats.h`. They are written inside
the study for now, as `sample_excess_kurtosis` is in
`studies/abm_system_tail_origin.c` and Mardia's test is in
`studies/abm_system_winner_normality.c`. They belong in et_al.

### The cross-section of firms behind GDP

`GDP_r(1) = Q1tot * dim_mach + Q2tot`, a sum over 20 K-firms and 200 C-firms,
so the cross-section a central limit theorem would act on is the 220 firms'
output. The archives hold no per-firm data.
`applications/abm_system_micro_simulate.c` reruns 5 configurations of the
design at seeds 1 to 5 with the `-m 1` switch `docs/DSK_MODEL_CHANGES.md`
records, which writes each firm's output. Those 25 runs are replicates 0 to 4
of the same configurations in `dataset/abm_system`, and the five macro series
rebuilt from one of them matched the stored replicate on all 2000 values, so
they are the same runs. `studies/abm_system_micro_clt.c` is the study.

The two firm types are kept apart. Periods 201 to 600.

The cross-section of output **levels** at a date, averaged over the periods of
a run then over the 25 runs. `H * N` is the concentration of shares: 1 when
every firm is the same size, N when one firm is everything. AD is
Anderson-Darling against a Gaussian with the null simulated at that firm count,
and would reject 5 per cent of the time on a Gaussian cross-section:

| | N | skewness | excess kurtosis | `H * N` | AD rejects |
|---|---:|---:|---:|---:|---:|
| K-firms | 20 | 2.620 | 7.102 | 5.757 | 0.938 |
| C-firms | 200 | -0.995 | 1.934 | 1.066 | 0.998 |

Average correlation over every pair of firm output changes within a run, and
the effective count `N / (1 + (N - 1) rho)` it implies:

| | N | rho | effective N |
|---|---:|---:|---:|
| K-firms | 20 | -0.0024 | 20.94 |
| C-firms | 200 | -0.0008 | 236.37 |

The levels are not the cross-section GDP growth is driven by. A firm's level
barely moves from one quarter to the next, so the spread of levels is mostly
firm size. The cross-section of **changes** is the sum GDP's change is:

| | N | excess kurtosis of change | concentration of \|change\| |
|---|---:|---:|---:|
| K-firms | 20 | 3.603 | 3.028 |
| C-firms | 200 | 14.459 | 3.103 |

GDP's change is exactly 220 times the cross-sectional mean change, so a central
limit theorem over the cross-section is a statement about that mean. If the
firms are independent at a date with cross-sectional standard deviation
sigma_t, the aggregate change has standard deviation sqrt(220) sigma_t at that
date. GDP growth over the 400 quarters of a run, averaged over the 25 runs:

| | skewness | excess kurtosis |
|---|---:|---:|
| as it stands | 0.009 | 0.571 |
| divided by sqrt(220) times that quarter's own cross-sectional spread | 0.078 | 0.187 |
| what the moving variance alone would give, 3 (E[v^2] / E[v]^2 - 1) | | 0.635 |

Dividing each quarter's growth by that quarter's own cross-sectional spread
removes about two thirds of the excess kurtosis. What is left over the time
series is mostly the cross-sectional spread moving from quarter to quarter
rather than the aggregation failing at any single date.

The three rows do not close: 0.187 and 0.635 do not sum to 0.571. Excess
kurtosis does not decompose additively, and the dispersion term overshooting
the raw figure means the variance path and the standardised shocks are not
independent of each other. How much of the overshoot is that dependence is not
measured here.

0.571 against the 1.151 the section above reports: that one is the mean over
1,000,000 replicates at n = 400, and `out/abm_system_tail_origin_report.txt`
gives 0.61 as the median over the 1000 configurations. These 5 configurations
sit at the median, and the mean is pulled up by a right tail across the design.

`dim_mach` is not read from the model's parameter file. It is recovered as
(GDP - sum of C-firm output) / (sum of K-firm output) at every period and
required to be the same number throughout, which is also the check that the
per-firm files and the aggregate file are the same run. It comes out at 40.

### Stable with index below 2, or finite variance

The section before last reports that GDP growth's sample excess kurtosis rises
with the sample, 0.650 at n = 100 to 1.151 at n = 400. That is what an infinite
fourth moment looks like. It is also what a fixed leptokurtic distribution
looks like, since sample excess kurtosis is bounded above by roughly n and
biased downward, so that measurement does not separate the two.

Stable laws are indexed by alpha in (0, 2]; alpha = 2 is the Gaussian. Every
alpha below 2 has P(|X| > x) ~ c x^-alpha and an infinite variance, so within
the family finite variance and Gaussian are the same condition and there is no
middle. `studies/abm_system_stable_tails.c` separates the cases three ways,
none of which reads a normality test:

- the slope of log sample variance on log n, over n = 50 to 400. Zero for a
  finite variance, 2/alpha - 1 for a stable alpha.
- Hill's estimator of the tail index, `alpha_hat = k / sum_i (ln |x|_(i) -
  ln |x|_(k+1))` on the largest 10, 5 and 2.5 per cent of |x|.
- the slope of log scale on log block length, summing k consecutive quarters
  for k = 1 to 32. One half for a finite variance with independent terms,
  1/alpha for a stable alpha. The scale is the interquartile range as well as
  the standard deviation, because the standard deviation of a sample from an
  infinite-variance law estimates nothing.

None of the three is read off a formula for what it should give. The same three
run first on 20,000 samples of 400 values from laws whose answer is known, seed
20260923, drawn by the Chambers-Mallows-Stuck method for the stable ones:

| law | var slope | Hill 10% | Hill 5% | Hill 2.5% | iqr slope | sd slope |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian, alpha = 2 | 0.009 | 4.76 | 6.14 | 7.73 | 0.473 | 0.488 |
| stable alpha = 1.8 | 0.119 | 2.90 | 2.86 | 2.72 | 0.530 | 0.489 |
| stable alpha = 1.5 | 0.341 | 1.77 | 1.73 | 1.76 | 0.649 | 0.490 |

Theory gives var slope 0, 0.111 and 0.333 and iqr slope 0.5, 0.556 and 0.667,
so the instrument reads correctly. Two things the calibration shows that could
not be assumed. Hill rises as k shrinks for the Gaussian, 4.76 to 7.73, and
stays flat for the stable laws, 2.90 to 2.72, so the direction separates them
more sharply than the level does; at n = 400 Hill returns 1.77 for a true 1.5.
And the sd slope is about 0.49 for all three, including alpha = 1.5, so it
carries no information at all and only the interquartile range works.

All 1000 configurations, all 1,000,000 replicates, 400 quarters each, one
sample per replicate, 1.7 minutes:

| series | var slope | Hill 10% | Hill 5% | Hill 2.5% | iqr slope | sd slope |
|---|---:|---:|---:|---:|---:|---:|
| GDP growth | 0.022 | 3.99 | 4.79 | 5.64 | 0.220 | 0.218 |
| energy growth | 0.042 | 3.98 | 4.74 | 5.52 | 0.308 | 0.306 |
| employment change | 0.025 | 3.99 | 4.75 | 5.53 | 0.174 | 0.169 |
| inflation | 0.031 | 6.95 | 8.24 | 9.65 | 0.340 | 0.337 |
| interest rate | 0.089 | 20.08 | 24.75 | 31.94 | 0.625 | 0.638 |
| interest rate change | 0.046 | 4.42 | 5.60 | 6.97 | 0.354 | 0.358 |

All three say the same thing. The variance slope of GDP growth is 0.022 against
0.009 for a Gaussian and 0.119 for alpha = 1.8, so the sample variance settles.
Hill rises as k shrinks, 3.99 to 5.64, which is the Gaussian pattern; a stable
law is flat or falling. And the block slope is 0.220, far below the 0.5 an
independent finite-variance sum gives, where a stable alpha below 2 would put
it above 0.5. None of the series is stable with an index below 2.

The block slopes sitting well below 0.5 are a separate fact: summing
consecutive quarters raises the scale of these series by much less than
independent terms would, which is negative dependence rather than anything
about the tails. The interest rate, a persistent level, is the one series above
0.5 at 0.625.

Hill on the cross-section of firm output changes at each date, 10,000
cross-sections from `dataset/abm_system_micro`:

| cross-section | firms | Hill 10% | Hill 5% | Hill 2.5% |
|---|---:|---:|---:|---:|
| K-firm changes | 20 | 3.09 | | |
| C-firm changes | 200 | 2.19 | 2.22 | 2.67 |

The K-firm cells are empty because 5 and 2.5 per cent of 20 firms is fewer than
the two order statistics the estimator needs. These are the noisiest numbers
here: 10 per cent of 200 firms is 20 values, against the 40 the n = 400 series
above give, and the calibration was run at n = 400 rather than at these sizes.

Three things none of this controls for. A run's quarters are dependent, and
dependence moves the block slope on its own, so a slope away from 0.5 is not by
itself a statement about tails. A series that is not stationary makes the
sample variance grow with n whatever its tails are. And Hill assumes the tail
is already power-law over the order statistics used, which at n = 400 and
k = 40 it need not be.

`out/abm_system_stable_tails_report.txt` and
`out/abm_system_stable_tails_by_series.csv` hold the run's own output;
`out/abm_system_micro_clt_report.txt`, `out/abm_system_micro_clt_by_run.csv`
and `out/abm_system_micro_clt_growth.csv` hold the cross-section study's.

### What the measurements leave, and which laws have those properties

Taken together the three sections above say: non-Gaussian, finite variance,
excess kurtosis 0.57 to 1.15 depending on the sample length, skewness near zero
for the four growth series and -0.124 for inflation, no stable behaviour, and
Hill below the Gaussian calibration at the same n.

"Thin tailed" is the wrong phrase for that. Excess kurtosis is positive and
Hill comes back at 3.99, 4.79 and 5.64 against the Gaussian calibration's 4.76,
6.14 and 7.73, and a lower Hill is a heavier tail. The series are heavier
tailed than a Gaussian and much lighter than the US ones, which
`out/abm_system_tail_origin_report.txt` puts at 2.87 excess kurtosis for US GDP
growth against a simulated median of 0.61.

Hill's level narrows nothing here. The Gaussian calibration returned 4.76
although a Gaussian has no tail index at all, so at n = 400 only the comparison
against that calibration carries information, not the number.

The named laws with finite variance, a tail heavier than Gaussian and no stable
behaviour fall into three groups.

Exponential-tailed, every moment finite:

| law | shape | excess kurtosis |
|---|---|---|
| Subbotin, density proportional to exp(-\|x/a\|^b) | b, with b = 2 Gaussian and b = 1 Laplace | Gamma(5/b)Gamma(1/b)/Gamma(3/b)^2 - 3, which is 0.76 at b = 1.5 and 3 at b = 1 |
| Laplace | none | 3 |
| normal inverse Gaussian, variance gamma, generalised hyperbolic | semi-heavy, exp(-alpha\|x\|) times a power | free |
| tempered stable, the truncated Levy flight | stable body, exponentially cut tail | free |

Power-law tailed with an index above 4, so the variance and the fourth moment
both exist:

| law | tail index | excess kurtosis |
|---|---|---|
| Student t, nu > 4 | nu | 6 / (nu - 4), so nu about 9 for 1.15 and about 14 for 0.57 |

Normal variance mixtures, X = sigma Z with Z Gaussian and sigma random, whose
excess kurtosis is 3 (E[sigma^4] / E[sigma^2]^2 - 1). That is the same
expression that returned 0.635 in the cross-section section above, so this
group fits the mechanism found there as well as the moments. The mixing law
picks the member: inverse-gamma gives Student t, gamma gives the variance
gamma, inverse Gaussian gives the normal inverse Gaussian, and lognormal gives
a mixture with every moment finite.

At the level of the process rather than a single draw, a GARCH-type series has
a leptokurtic unconditional law with a finite variance, and a finite kurtosis
when 3 a1^2 + 2 a1 b1 + b1^2 < 1. Time-varying dispersion is what the
cross-section section found, so this is the natural description of the series
rather than of one observation.

One fit already in this project does not settle the question.
`applications/us_qvarma_spec_choice.c` fits a Student t and gives nu near 7 on
the US data against nu in the thousands on the simulated ones, but that is the
conditional innovation after the QVARMA filter, not the unconditional marginal
these sections measure.

### Fitting those laws

`studies/abm_system_marginal_fit.c` fits all eight by maximum likelihood.
`studies/marginal_laws.h` holds them, and says why they are written there
rather than called from et_al: et_al has the Gaussian and the Student t in
their Mat forms and none of the other six, and no Bessel function of any kind
for the three that need one. All eight come through the same fitting routine on
et_al's `solver/lbfgs.h`, so what is compared is laws and not optimisers.
`tests/marginal_laws_correctness.c` is what says the densities are right, and
the section after this one is what it found.

Setup. 100 configurations of the design, 2 replicates of each, 200 series per
variable, 400 quarters each, one fit per series per law, 34.5 minutes on 16
threads. Nothing is pooled across replicates. The comparison is by the Schwarz
criterion, p log n - 2 log L, which is what charges the four- and
five-parameter laws for their parameters at n = 400. The sample is a sample
because of what the fits cost: the three laws with a Bessel function evaluate a
sixty-four-node quadrature per observation per likelihood, and a likelihood is
called a few thousand times per fit.

Share of the 200 series each law wins on the Schwarz criterion:

| law | p | GDP growth | energy | employment | inflation | interest rate | rate change |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gaussian | 2 | 0.570 | 0.560 | 0.530 | 0.535 | 0.680 | 0.745 |
| Laplace | 2 | 0.025 | 0.010 | 0.015 | 0.030 | 0.010 | 0.015 |
| Subbotin | 3 | 0.040 | 0.050 | 0.060 | 0.015 | 0.170 | 0.050 |
| Student t | 3 | 0.320 | 0.305 | 0.340 | 0.365 | 0.025 | 0.115 |
| lognormal mixture | 3 | 0.040 | 0.070 | 0.050 | 0.050 | 0.030 | 0.050 |
| normal inverse Gaussian | 4 | 0.000 | 0.000 | 0.000 | 0.000 | 0.010 | 0.000 |
| variance gamma | 4 | 0.005 | 0.005 | 0.005 | 0.005 | 0.075 | 0.025 |
| generalised hyperbolic | 5 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

Mean log-likelihood less the Gaussian's on the same series, GDP growth: Laplace
-7.87, Subbotin 4.32, Student t 6.29, lognormal mixture 5.27, normal inverse
Gaussian 6.22, variance gamma 10.90, generalised hyperbolic 6.84. Every law
with a free shape beats the Gaussian on raw likelihood and only the Laplace,
whose excess kurtosis is fixed at 3, does worse. What decides the table above
is the parameter charge, not the fit.

Mean fitted shapes for GDP growth: Subbotin b = 1.67, Student t nu = 19.4,
mixture spread 0.206, normal inverse Gaussian alpha = 426, variance gamma
kappa = 0.269.

The US series, 187 quarters, one fit each, by the same routine:

| series | lowest BIC | its shape | Gaussian BIC | winner BIC |
|---|---|---:|---:|---:|
| GDP growth | Laplace | | 436.45 | 405.92 |
| energy growth | Gaussian | | 756.64 | 756.64 |
| employment change | normal inverse Gaussian | alpha = 2.73 | 133.33 | 63.17 |
| inflation | normal inverse Gaussian | alpha = 0.859 | 461.94 | 409.04 |
| interest rate | Gaussian | | 1063.44 | 1063.44 |

The Student t fitted to the US series gives nu = 2.86 for GDP growth and
nu = 2.00 for employment change, against 19.4 and 1825 on the simulated ones.
At nu = 2 the variance does not exist, which is a different statement about the
US data than anything in the sections above, and it comes from a single sample
of 187.

Three things this does not settle, stated because the tables above look more
decisive than they are.

Convergence is not uniform across laws, and a law that did not converge is
excluded rather than counted as a loss, so a law with poor convergence is
judged on its easier cases. On the simulated series the Student t converged on
58 per cent of interest-rate fits and 87 per cent of GDP-growth fits, the
variance gamma on 79 per cent and the generalised hyperbolic on 85 per cent.
Four of the 40 US fits did not converge at all.

The mean fitted shape is not a usable summary for the Student t. Its mean nu
is 7190 for energy growth and 15040 for the interest rate, which is individual
fits running off to the Gaussian limit and dragging the mean with them. A
median would say something; the study reports a mean.

The simulated series have 400 quarters and the US series 187, and the Schwarz
penalty depends on n, so the two tables are not a like-for-like comparison of
which law wins where.

`out/abm_system_marginal_fit_report.txt`,
`out/abm_system_marginal_fit_by_law.csv` and
`out/abm_system_marginal_fit_us.csv` hold the run's own output.

### What the correctness test caught

The eight densities and their fits were new code, so
`tests/marginal_laws_correctness.c` was written before anything was fitted. It
checks log K_nu(x) against the closed forms at orders 1/2 and 3/2, against
published values at orders with none, against the recurrence
K_{nu+1} = K_{nu-1} + (2 nu / x) K_nu and against K being even in its order;
that every density integrates to 1; that the generalised hyperbolic at
lambda = -1/2 is the normal inverse Gaussian; that the scalar Gaussian and
Student t agree with et_al's Mat forms; and that maximum likelihood recovers
the parameters it was given.

It found four defects, none of which would have shown up in a fitted number as
anything other than one family fitting slightly worse than it should.

The variance gamma returned minus infinity exactly at its location. The density
there is |x - mu|^order times K_order(c |x - mu|), which is zero times infinity
numerically rather than a singularity; for a positive order the two cancel and
the value is finite. An optimiser moving the location onto an observation would
have been pushed out of a good region.

The Bessel function rejected a negative order. K is even in its order, and the
variance gamma reaches a negative one whenever its shape goes above 2.

The variance gamma's shape had to be bounded below 2, which is not a
convenience. Above 2 the density diverges at the location as a power, so with
the location free the likelihood is unbounded and there is no maximum to find.
It is the pathology the three-parameter lognormal has. Exactly at 2 the
divergence is logarithmic, which is still unbounded, and the logistic transform
reached exactly 2 in double precision, so its argument is clamped.

Three densities passed an optimiser's proposals straight into a Bessel function
that asserts on a non-positive argument. et_al's own rule is that an assert is
for programmer error and an infeasible parameter value is not one. Both routes
to a non-positive argument are reachable by ordinary steps: a scale driven
small enough to underflow, and an asymmetry driven close enough to the tail
parameter that alpha^2 - beta^2 underflows to zero while the |beta| < alpha
guard still passes. They return a sentinel now.

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

That disagreement is specific to this loss. The same comparison over the two
score losses of `docs/ABM_SYSTEM_SCORE_LOSS.md` puts the Spearman rank
correlation between the two elimination orders at 0.9992 and 0.9981, with a
median round gap of 6. Whatever makes the elimination order unstable here is a
property of the impulse-response distance, not of the confidence set.

### A loss that does not need the million fits

`docs/ABM_SYSTEM_SCORE_LOSS.md` documents a different way to reach a loss
matrix for this same confidence set. It holds the US estimate fixed and
evaluates the score of the auxiliary model's log-likelihood on each simulated
series, so nothing is fitted to simulated data and steps 6 and 7 of the
pipeline are not needed at all. Weighted by the inverse information matrix,
that loss is Rao's score statistic, it does not depend on how the auxiliary
model is parameterized, and it is the only one of the three losses under which
this confidence set stops because an equivalence test is accepted rather than
because elimination ran out of configurations. It keeps `cop_0409` and
`cop_0599`; this one keeps `cop_0191`.

The two results are not near neighbours. Under `MCS_TR` on the weighted score
loss, `cop_0191` is eliminated in round 657 of 999 at a p-value below 0.0001,
which is mid-field rather than a close miss. In the other direction the gap is
smaller: `cop_0409` and `cop_0599` survive this procedure to rounds 865 and
806. `docs/ABM_SYSTEM_SCORE_LOSS.md`, "When each winner leaves the procedures
it does not win", has the full table. Nothing in either document settles which
of the two measures should be believed.

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
make app-abm_system_irf_loss                   # step 7: loss table, parallel over replicates
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
