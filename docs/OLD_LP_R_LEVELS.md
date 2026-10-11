# Local projections on R's levels, discarded

## Status: discarded

The local projections were also fitted on the series the collaborator's R
function `transform_data()` builds: cumulative sums of the five stored series
started from fixed values, with GDP and energy in logs, in R's column order.
Every other experiment in this project fits the five series the t-QVARMA is
fitted on: GDP growth, energy growth, employment change, inflation and the
interest rate, `docs/DATA_DOCUMENTATION.md`.

The levels are discarded because the main transformation loses less
information. Where the right answer is known they also fail: in the Monte
Carlo experiment the linear local projection on R's levels put the right
configuration in the set for 0 of 100 runs of `cop_0191`, against 90 of 100 on
the main transformation (`docs/OLD_MONTECARLO_LP_VALIDATION.md`, "Why the
levels fail"). `docs/OLD_EXPERIMENTS.md` lists every discarded experiment.

The state-dependent model still uses R's construction for one series, the
4-period moving average of log GDP that decides each period's state;
`applications/lp_system.h` builds it, and it is the same series under both
transformations.

| program | what it does |
| --- | --- |
| `old/lp_r_levels/lp_system_r_levels.h` | R's `transform_data()`, on top of `applications/lp_system.h` |
| `old/lp_r_levels/us_lp_fit_r_levels.c` | the US fit on R's levels |
| `old/lp_r_levels/abm_system_fit_lp_r_levels.c` | every simulated run's fit on R's levels |
| `old/lp_r_levels/abm_system_lp_r_levels_irf_loss.c` | the loss tables against the US fit |
| `old/lp_r_levels/abm_system_lp_r_levels_mcs.c` | their confidence sets |
| `old/lp_r_levels/lp_system_r_levels_transform.c` | checks the header against R's arithmetic, and that the live state series equals the one R's levels build |
| `old/lp_r_levels/lp_r_levels_reproduction.sh` | checks that moving the variant out changed no result on either side |

The outputs stay where they were: `out/abm_system_fit_lp_r_levels/`,
`out/us_lp_fit_r_levels.*`, `out/abm_system_lp_r_levels_*` and, for the Monte
Carlo, `montecarlo/out/lp_r_levels/`.

## Against the US data

Every step of `docs/ABM_SYSTEM_LP_VALIDATION.md` run again on R's levels: the
simulated runs and the US series put through R's `transform_data`
(`lp_system_r_transform` in `old/lp_r_levels/lp_system_r_levels.h`), log GDP as
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

- `make old-abm_system_fit_lp_r_levels` makes the simulated fits, about 7
  minutes.
- `make old-abm_system_lp_r_levels_mcs` makes the US fit, the loss tables and
  the sets.

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
(`docs/OLD_MONTECARLO_LP_VALIDATION.md`). On R's levels the local projections
recovered a known answer in 5 (linear), 72 (state 1) and 42 (state 2) of 1000
repetitions, so this set says which configuration's estimates are most
precise rather than which is closest to the US data.
