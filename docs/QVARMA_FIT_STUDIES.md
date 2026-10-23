# Studies of the t-QVARMA fits

The programs under `studies/` answer a question about the pipeline's own result
and write the answer to `out/studies/`. They pass or fail nothing and gate
nothing. `make study` runs them all in the order of `STUDY_STEMS` in the
Makefile, and `make study-<stem>` runs one. Every one reads what the pipeline
already wrote to `out/` and none rebuilds it, so none runs on a fresh clone.

Seven of the studies below read the confidence set or the loss table of the
naive loss (`docs/OLD_EXPERIMENTS.md`, "The naive loss"):
`abm_system_winner_diagnostics`, `abm_system_seed_correlation`,
`us_qvarma_nu_sensitivity` and the four other studies of the winner's `nu`. They
moved with it to `old/naive_mae_loss/studies/`, their outputs to
`out/old/naive_mae_loss/studies/`, and they run with `make old-<stem>` instead
of `make study-<stem>`; `make study` no longer runs them. Their write-ups stay
here and in the documents the table names.

Where each study is written up:

| study | question | write-up |
|---|---|---|
| `abm_system_winner_diagnostics` | why `cop_0191` wins, and how its loss compares with an all-zero response | `docs/OLD_ABM_SYSTEM_MCS_VALIDATION.md`, "What the confidence set does not say"; figures below |
| `abm_system_tail_origin` | moments of the simulated series against the US ones | `docs/OLD_ABM_SYSTEM_MCS_VALIDATION.md`, same section |
| `abm_system_gaussian_convergence`, `abm_system_micro_clt`, `abm_system_stable_tails`, `abm_system_marginal_fit` | the distribution of the simulated series | `docs/ABM_SYSTEM_DISTRIBUTION.md` |
| `abm_system_seed_correlation` | whether the shared seeds correlate the configurations | `docs/ABM_SYSTEM_SIMULATION.md`, "Common random numbers" |
| `qvarma_stuck_fits`, `qvarma_conditioning`, `qvarma_convergence_criteria` | why so many fits stop without converging | this file |
| `us_qvarma_nu_sensitivity`, `abm_system_winner_nu_profile`, `abm_system_winner_nu_likelihood_scan`, `abm_system_winner_tail_comparison`, `abm_system_winner_normality` | the degrees of freedom `nu` of the winner's fits against the US fit's | this file |

`nu` is the Student t degrees of freedom of the t-QVARMA's innovations. A small
`nu` means heavy tails; as `nu` grows the Student t approaches the Gaussian.

## Why many fits stop without converging

These three were run between late August and 20 September 2026, before the
design dataset existed: the first two on the older dataset of 100
configurations of 108 replications each (named `EstimationSeriesSample1_N` in
their reports), both auxiliary specifications, the third on data simulated
from the t-QVARMA itself. They have not been rerun on the design dataset.

**`studies/qvarma_stuck_fits.c`**, `out/studies/qvarma_stuck_fits.txt`. Every
cached fit that had not converged, 20,295, probed at its own stored parameters:
the objective along the steepest-descent direction at 26 step lengths, one
decade apart, from 1e-14 to 1e11. A fit counts as stuck when no step lowers the
objective.

- 226 fits (1.1%) have no lower point on that ladder; 20,069 (98.9%) do.
- Of the 10,869 fits whose next L-BFGS iteration made no progress, 10,650 have
  a lower point on the ladder, at a median distance of 1e-10 in parameter
  space; all but one are closer than 2^-20 = 9.5e-7, the shortest step the
  line search's bisection reaches with its default of 20 halvings.
- With 60 halvings instead of 20, 94 of those 10,869 still make no progress
  and 10,717 raise the log-likelihood, by 0.00034 on average and 0.056 at most
  in one iteration.

So the fits are not blocked: a lower point exists, very close, and the gain
from reaching it is small.

**`studies/qvarma_conditioning.c`**, `out/studies/qvarma_conditioning.txt`. Five
replicates of the older dataset, each fitted from the same start at
independent budgets of 250 to 32,000 iterations, recording the log-likelihood,
the smallest residual standard deviation (the smallest diagonal of
`Omega_inv`), `nu` and the gradient norm, and the Hessian's eigenvalues at
32,000 iterations. The question is whether the likelihood has no interior
maximum, which would show as a residual standard deviation walking to zero
while the log-likelihood keeps rising, or is badly scaled, which would show as
both settling.

- In all five the smallest residual standard deviation is the interest
  rate's, and it settles: replicate 0, for example, at 0.0223 from 2000
  iterations on, while the log-likelihood stops at 713.33 from 4000 on.
- Three of the five Hessians have negative eigenvalues (19, 9 and 19 of 42),
  with condition numbers 7.6e11 to 5.6e13, and the stiffest coordinate is the
  intercept of the employment-change equation, `c[2]`. The other two are
  positive definite, with condition numbers 5.9e6 and 3.0e8.

By the study's own criterion this is a scaling problem rather than a
likelihood without a maximum; three of the five points it stops at are not
minima of the objective.

**`studies/qvarma_convergence_criteria.c`**,
`out/studies/qvarma_convergence_criteria.txt`. Which convergence test reads the
same at the same statistical distance from the optimum whatever the sample
length and the units. Data simulated from a fixed t-QVARMA(1,1,2) with `nu = 8`;
72 cells, sample length T in {200, 400, 1600} by the series multiplied by 0.01,
1 or 100 (an exact change of units of the same model) by 8 replicates. The
report's header line says the scales are {1, 100} and prints 0.01 as 0; the
program uses 0.01, 1 and 100. The reference point of each cell is the fit after
50,000 iterations, and a fit "has arrived" when every coordinate is within 0.01
of its standard error of that point. 57 of the 59 cells with a positive definite
reference Hessian arrived within 6400 iterations.

Ratio of the largest to the smallest reading of each criterion at the moment
of arrival, over those 57 cells (1 would mean it reads the same everywhere):

| criterion | max / min |
|---|---:|
| largest absolute gradient entry, what et_al's L-BFGS stops on | 36,890 |
| the same scaled by the coordinates and the objective | 1,514 |
| the same divided by T | 4,612 |
| Newton decrement, `g' H^-1 g` | 6,867,000 |

Only the Newton decrement separates the fits that have arrived (within 0.01
standard errors) from those that have not (more than 1 standard error away)
with one threshold: its largest reading on an arrived fit, 0.031, is below its
smallest on a fit that has not, 1.6.

## The degrees of freedom of the winner against the US fit

The US fit has `nu = 7.2`; the parameter set averaged over `cop_0191`'s 1000
fits has `nu = 2853`, and the individual fits put it anywhere from about 20 to
about 10^12. These studies ask whether that difference matters to the loss and
where it comes from. All ran on the design dataset; the reports of the last four
were committed on 17 September 2026. The report of `us_qvarma_nu_sensitivity`
names the confidence set's output by its earlier file name,
`out/abm_system_mcs_joint.csv`, so it predates the rename to
`out/old/naive_mae_loss/applications/abm_system_mcs.csv`; the numbers it reads
from it match the current file.

**`old/naive_mae_loss/studies/us_qvarma_nu_sensitivity.c`**,
`out/old/naive_mae_loss/studies/us_qvarma_nu_sensitivity.txt`. The US fit (`nu`
free, log-likelihood -882.14) against two variants: the same parameters with
`nu` set to the winner's 2853.4 ("swapped"), and every other parameter refitted
with `nu` held there ("refitted", log-likelihood -943.95, a likelihood ratio of
123.6). Distance of each variant's impulse response from the US benchmark's, in
the loss's own units (mean absolute difference over the 525 stacked responses):

| variant | all horizons | horizon 0 | horizons 1 to 20 |
|---|---:|---:|---:|
| swapped | 0.0323 | 0.0290 | 0.0325 |
| refitted | 0.0529 | 0.0378 | 0.0536 |

For scale, the mean losses of the 1000 configurations span 0.0451, and the two
smallest differ by 0.0005. Changing `nu` on the US side alone moves the
benchmark response by about the whole spread of the configurations.

**`old/naive_mae_loss/studies/abm_system_winner_nu_profile.c`**,
`out/old/naive_mae_loss/studies/abm_system_winner_nu_profile.txt`. A pilot on 8
of `cop_0191`'s replicates: each refitted with `nu` held at 7.2 (the US value),
10, 30, 100, 1000, 10^4, 10^5 and 10^6, every other parameter free, at most
20,000 iterations, starting from its own cached fit. Gain is the refit's
log-likelihood minus the cached fit's.

| held `nu` | median gain | replicates with a higher likelihood | mean loss |
|---:|---:|---:|---:|
| 7.2 | -258.8 | 0 of 8 | 0.11346 |
| 10 | -156.3 | 0 of 8 | 0.11441 |
| 100 | -40.2 | 0 of 8 | 0.11728 |
| 1000 | 0.7 | 4 of 8 | 0.11803 |
| 10^6 | 3.0 | 5 of 8 | 0.11799 |

Their own fits' mean loss is 0.11791. Holding `nu` at the US value costs a
median of 259 log-likelihood units and lowers the loss by 0.0045, about nine
times the gap between the two best configurations. Only 2 of the 8 refits at
`nu = 7.2` converged.

**`old/naive_mae_loss/studies/abm_system_winner_nu_likelihood_scan.c`**,
`out/old/naive_mae_loss/studies/abm_system_winner_nu_likelihood_scan.txt`. All
1000 fits of `cop_0191`, nothing refitted: the log-likelihood at the fitted
parameters with only `nu` moved. A fit is "smooth" when moving `nu` by 1% either
way changes the log-likelihood by less than 1, "knife-edge" when either move
changes it by more than 100.

- 432 fits are smooth, 141 between, 427 knife-edge. The median change for a 1%
  move is about -9.2.
- The three classes have the same mean loss to the fourth decimal (0.1182,
  0.1183, 0.1181).

**`old/naive_mae_loss/studies/abm_system_winner_tail_comparison.c`**,
`out/old/naive_mae_loss/studies/abm_system_winner_tail_comparison.txt`. The
per-period quadratic form of the residuals over the number of series, `q_t / K`
with `K = 5`, which has mean 1 and variance 0.4 under Gaussian innovations.

| residuals | variance of `q_t / K` | largest | share above 5 |
|---|---:|---:|---:|
| US, `nu` free (7.2) | 6.18 | 15.21 | 0.086 |
| US, `nu` held at 2853 | 1.35 | 7.81 | 0.021 |
| `cop_0191`, median over its 1000 fits | 0.68 | 6.48 | 0.005 |

The US values far out come from 1974 to 1981 and from 2008Q4 to 2009Q1, and
the interest rate carries the largest share of the US sum of `q_t`, 27.8%.

**`old/naive_mae_loss/studies/abm_system_winner_normality.c`**,
`out/old/naive_mae_loss/studies/abm_system_winner_normality.txt`. Mardia's tests
of whether `cop_0191`'s fitted parameter vectors, 42 coordinates on the
optimizer's unconstrained scale, are jointly normal across replicates. Joint
normality is rejected for all 1000 fits and for the 449 converged ones, with
p-values below what a double can represent, under both the skewness and the
kurtosis test. Calibrated on 1000 samples from a 42-dimensional standard normal
of the same sizes, seed 20260917: the asymptotic kurtosis test rejects 12.8% (n
= 1000) and 28.6% (n = 449) at a nominal 5%, the version with exact moments 4.4%
and 5.1%, and the skewness tests 2.4% to 7.2%.

**`old/naive_mae_loss/studies/abm_system_winner_diagnostics_plots.py`** draws
the tables of `old/naive_mae_loss/studies/abm_system_winner_diagnostics.c` into
`out/old/naive_mae_loss/studies/abm_system_winner_diagnostics_plots/`: one
histogram per fitted parameter over `cop_0191`'s fits, on the model's scale and
on the optimizer's unconstrained scale, each bar coloured by the share of its
fits that converged, with the US value, the median and the averaged parameter
set marked. It needs `polars`, `plotly` and `kaleido`.

    python old/naive_mae_loss/studies/abm_system_winner_diagnostics_plots.py
