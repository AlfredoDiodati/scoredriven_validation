# Validation against the US data with local projections

The main pipeline (`docs/ABM_SYSTEM_MCS_VALIDATION.md`) asks which of the 1000
parameter configurations of the DSK model produce dynamics closest to the US
data, comparing impulse responses of the t-QVARMA. This runs the same
comparison with the collaborator's local projections as the auxiliary model,
the version the Monte Carlo grid checks on simulated data with a known answer
(`docs/MONTECARLO_VALIDATION.md`).

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

The same comparison on R's levels is discarded; `docs/OLD_LP_R_LEVELS.md`
records it and its result.

## What this does not do

- The collaborator's thesis also runs a confidence set on the two states
  together. How it combines them is not written down in the thesis or in the
  scripts received, so that set is not computed.
- The thesis's validation score, the mean distance relative to a bootstrap
  variance of the US responses, is not computed.
- The loss is the mean absolute difference, not the thesis's mean squared
  difference, and the series are the stationary ones, not the levels the
  thesis fits; `docs/OLD_MONTECARLO_LP_VALIDATION.md` records why.
