# The distribution of the simulated series

`docs/ABM_SYSTEM_MCS_VALIDATION.md`, "What the confidence set does not say",
reports that no configuration of the design reproduces the variance of the US
interest rate and inflation or the skewness of US GDP growth and employment
change, and that the simulated series have much lighter tails than the US
ones (excess kurtosis of GDP growth 2.87 in the US data against a median of
0.61 over the 1000 configurations, `out/abm_system_tail_origin_report.txt`).
That is a property of the model, not of this project's build of it
(`tests/dsk_tail_replicate_reproduction.c`).

This file asks where that shape comes from and what distribution describes
it: whether aggregation over firms makes the series Gaussian, whether they are
stable laws with infinite variance, and which named laws fit them. Nothing
here feeds the confidence set. Every study in it reads the stored dataset
(`dataset/abm_system/`) or the per-firm side run (`dataset/abm_system_micro/`,
README, "A side run with the per-firm output") and fits nothing of the
auxiliary model.

The five stored series are GDP growth, energy demand growth, employment
change, inflation and the interest rate (`docs/DATA_DOCUMENTATION.md`); the
studies add the first difference of the interest rate as a sixth.

## Is the thinness a central limit theorem

GDP here is not drawn from a distribution. It is what comes out after adding up
what 200 consumption-good firms and 20 capital-good firms each did in a
quarter, and adding up many roughly independent quantities is the situation a
central limit theorem describes. Its conclusion is a Gaussian, which has no
large rare values, so it would account for the thinness described at the top
of this file.
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

## The cross-section of firms behind GDP

In the model, real GDP is `GDP_r(1) = Q1tot * dim_mach + Q2tot`
(`modules/module_macro_sfc.cpp`): the total output of the 20 capital-good firms
(K-firms), counted in machines and converted by `dim_mach`, the output one
machine can produce, plus the total output of the 200 consumption-good firms
(C-firms). So the cross-section a central limit theorem would act on is the 220
firms' output. The archives hold no per-firm data.
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

0.571 against the 1.151 the first section reports: that one is the mean over
1,000,000 replicates at n = 400, and `out/abm_system_tail_origin_report.txt`
gives 0.61 as the median over the 1000 configurations. These 5 configurations
sit at the median, and the mean is pulled up by a right tail across the design.

`dim_mach` is not read from the model's parameter file. It is recovered as
(GDP - sum of C-firm output) / (sum of K-firm output) at every period and
required to be the same number throughout, which is also the check that the
per-firm files and the aggregate file are the same run. It comes out at 40.

## Stable with index below 2, or finite variance

The first section reports that GDP growth's sample excess kurtosis rises
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

## What the measurements leave, and which laws have those properties

Taken together the sections above say: non-Gaussian, finite variance,
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

## Fitting those laws

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

Three things this does not settle.

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

## What the correctness test caught

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
