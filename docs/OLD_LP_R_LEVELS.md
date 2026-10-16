# Local projections on cumulated levels, discarded

## Status: discarded

The local projections were also fitted on the cumulated levels: levels rebuilt
by cumulating the five stored series from fixed starting values, with GDP and
energy in logs, in the column order GDP, employment, inflation, interest rate,
energy. This is the data preparation of the R code the local-projection
settings were first taken from (its function `transform_data()`);
`old/lp_r_levels/lp_system_r_levels.h` is the C port and states the arithmetic.
Every other experiment in this project fits the five series the t-QVARMA is
fitted on: GDP growth, energy growth, employment change, inflation and the
interest rate, `docs/DATA_DOCUMENTATION.md`.

The levels are discarded because the main transformation loses less
information. Where the right answer is known they also fail: in the Monte Carlo
experiment the linear local projection on the cumulated levels put the right
configuration in the set for 0 of 100 runs of `cop_0191`, against 90 of 100 on
the main transformation (`docs/OLD_MONTECARLO_LP_VALIDATION.md`, "Why the
levels fail"). `docs/OLD_EXPERIMENTS.md` lists every discarded experiment.

The state-dependent model still uses the same construction for one series, the
4-period moving average of the rebuilt log GDP that decides each period's
state; `applications/lp_system.h` builds it, and it is the same series under
both transformations.

| program | what it does |
| --- | --- |
| `old/lp_r_levels/lp_system_r_levels.h` | the cumulated levels, on top of `applications/lp_system.h` |
| `old/lp_r_levels/us_lp_fit_r_levels.c` | the US fit on the cumulated levels |
| `old/lp_r_levels/abm_system_fit_lp_r_levels.c` | every simulated run's fit on the cumulated levels |
| `old/lp_r_levels/abm_system_lp_r_levels_irf_loss.c` | the loss tables against the US fit |
| `old/lp_r_levels/abm_system_lp_r_levels_mcs.c` | their confidence sets |
| `old/lp_r_levels/lp_system_r_levels_transform.c` | checks the header against a plain-loop transcription of the R arithmetic, and that the live state series equals the one the cumulated levels build |
| `old/lp_r_levels/lp_r_levels_reproduction.sh` | checks that moving the variant out changed no result on either side |

The outputs stay where they were: `out/abm_system_fit_lp_r_levels/`,
`out/us_lp_fit_r_levels.*`, `out/abm_system_lp_r_levels_*` and, for the Monte
Carlo, `out/montecarlo/lp_r_levels/`.

## Against the US data

Every step of `docs/OLD_ABM_SYSTEM_LP_VALIDATION.md` run again on the cumulated
levels: the simulated runs and the US series put through
`lp_system_r_transform` in `old/lp_r_levels/lp_system_r_levels.h`, log GDP as
`100 log(100 + cumsum(GDP_growth))`, employment as `0.90 +
cumsum(Employment_change)`, inflation and the interest rate as they are, log
energy like GDP, in the order GDP, employment, inflation, interest rate,
energy; everything else unchanged. The US growth rates go through the same
transformation as the simulated ones, with the same starting values. Outputs
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

`cop_0905` is the configuration the levels picked in the Monte Carlo experiment
when the right answer was `cop_0191`: the one whose simulated responses scatter
least, with the most volatile interest rate of the 1000
(`docs/OLD_MONTECARLO_LP_VALIDATION.md`). On the cumulated levels the local
projections recovered a known answer in 5 (linear), 72 (state 1) and 42 (state
2) of 1000 repetitions, so this set says which configuration's estimates are
most precise rather than which is closest to the US data.
