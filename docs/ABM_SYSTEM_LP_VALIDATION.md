# Validation against the US data with local projections

The main pipeline (`docs/ABM_SYSTEM_MCS_VALIDATION.md`) asks which of the 1000
parameter configurations of the DSK model produce dynamics closest to the US
data, comparing impulse responses of the t-QVARMA. This runs the same
comparison with the collaborator's local projections as the auxiliary model,
the version `docs/MONTECARLO_LP_VALIDATION.md` checked on simulated data with a
known answer.

## What is run

- **Data.** The US series (`out/us_system.csv`) turned into the five series
  the t-QVARMA is fitted on: GDP growth, energy growth, employment change,
  inflation, interest rate (`us_system_growth_block` in
  `applications/us_data.h`). The first three quarters are dropped, where R's
  4-period moving average of log GDP, the series that decides the state, does
  not exist: 184 quarters, 1974Q1 to 2019Q4. The simulated runs, 1000 per
  configuration, 397 periods each, go through the same route
  (`lp_system_from_stored` in `applications/lp_system.h`).
- **Models.** The collaborator's settings (`_temp/estimate.R`), fitted by
  et_al's `lp_lin_and_nl`: 4 lags, horizons 0 to 15, unit shocks identified
  recursively in the order above; the state-dependent model with the HP cycle
  (lambda 1600) of the moving average through a logistic of slope 2, lagged one
  period. Three sets of responses, 400 values each: linear, state 1 (mostly
  periods with GDP above trend) and state 2 (mostly below).
- **US fit.** `applications/us_lp_fit.c`, stored in `out/us_lp_fit.npz`, with
  the fit notes in `out/us_lp_fit.txt` and every response in
  `out/us_lp_fit_irf.csv`. Every regression was solved at full rank.
- **Loss.** For each model, the mean absolute difference between the US
  responses and each simulated run's (`applications/abm_system_lp_irf_loss.c`,
  `out/abm_system_lp_irf_loss_<model>.csv`): 1000 runs by 1000
  configurations, no cell missing.
- **Confidence set.** `applications/abm_system_lp_mcs.c`: level 0.05, 10000
  bootstrap resamples of single runs, the bootstrap variance, et_al's seed 123,
  MCS_TR as the headline statistic and MCS_TMAX beside it, the t-QVARMA
  pipeline's settings (`out/abm_system_lp_mcs_<model>.txt` and `.csv`).

    make app-abm_system_lp_mcs      the US fit, the loss tables, the sets

It reads the simulated fits in `out/abm_system_fit_lp/`
(`make app-abm_system_fit_lp`) and does not rebuild them. The loss tables take
about 25 minutes, most of it reading the simulated fits and checking each
against its own series; the six confidence sets about 40 seconds.

## Result

- **Linear:** the set holds `cop_0191` alone, mean loss 0.861; the next
  configurations by mean loss are `cop_0905` (0.924) and `cop_0410` (0.963).
- **State 1:** `cop_0191` alone, mean loss 1.306; next `cop_0717` (1.341).
- **State 2:** four configurations, the elimination stopping on an accepted
  test at p = 0.164: `cop_0905` (mean loss 1.288), `cop_0191` (1.301, MCS
  p-value 0.193), `cop_0356` (1.309) and `cop_0437` (1.310).

MCS_TMAX keeps the same sets under all three models.

The t-QVARMA pipeline on the same data keeps `cop_0191` alone
(`out/abm_system_mcs_statistic_comparison.csv`). The two auxiliary models
agree.

`cop_0191` is also the configuration the Monte Carlo experiments used as their
known answer, chosen there because the t-QVARMA's US set had kept it and its
fits converged most often; nothing in the local-projection route was tuned to
it.

## The same on the collaborator's transformation

Every step above run again with `--r-levels`: the simulated runs and the US
series put through R's `transform_data` (`lp_system_r_transform`), log GDP as
`100 log(100 + cumsum(GDP_growth))`, employment as `0.90 +
cumsum(Employment_change)`, inflation and the interest rate as they are, log
energy like GDP, in R's order GDP, employment, inflation, interest rate,
energy; everything else unchanged. The collaborator's own preparation of the
US data (`rw_data_prep.R`) was not available; the thesis states the same
transformation code runs on the real and the simulated data, so the US growth
rates go through `transform_data` exactly as the simulated ones do. Outputs
carry `_r_levels`: `out/abm_system_fit_lp_r_levels/`,
`out/us_lp_fit_r_levels.*`, `out/abm_system_lp_r_levels_irf_loss_<model>.csv`,
`out/abm_system_lp_r_levels_mcs_<model>.*`. No fit failed and every regression
was at full rank, on the US data and on all 1,000,000 runs.

    make app-abm_system_fit_lp_r_levels     the simulated fits, about 7 minutes
    make app-abm_system_lp_r_levels_mcs     the US fit, the loss tables, the sets

- **Linear:** `cop_0905` alone (mean loss 2.517; next `cop_0593` 2.960);
  `cop_0191` is 26th (13.24).
- **State 1:** `cop_0905` alone (24.10; next `cop_0593` 27.07); `cop_0191`
  11th (33.14).
- **State 2:** `cop_0905` alone (4.495; next `cop_0593` 5.710); `cop_0191`
  23rd (28.36).

MCS_TMAX keeps the same single configuration under all three models.

`cop_0905` is the configuration the levels picked in the Monte Carlo
experiment when the right answer was `cop_0191`: the one whose simulated
responses scatter least, with the most volatile interest rate of the 1000
(`docs/MONTECARLO_LP_VALIDATION.md`). On R's levels the local projections
recovered a known answer in 5 (linear), 72 (state 1) and 42 (state 2) of 1000
repetitions, so this set says which configuration's estimates are most
precise rather than which is closest to the US data.

## What this does not do

- The collaborator's thesis also runs a confidence set on the two states
  together. How it combines them is not written down in the thesis or in the
  scripts received, so that set is not computed.
- The thesis's validation score, the mean distance relative to a bootstrap
  variance of the US responses, is not computed.
- The loss is the mean absolute difference, not the thesis's mean squared
  difference, and the series are the stationary ones, not the levels the
  thesis fits; `docs/MONTECARLO_LP_VALIDATION.md` records why.
