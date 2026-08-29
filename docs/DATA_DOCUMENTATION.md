# Data sources and the transformations used

## What is actually used

`us_real.csv` only. `applications/us_prepare_data.c` reads it and writes
`out/us_system.csv`, five columns plus a `Quarter` column, 1973Q1 to 2019Q4,
188 quarters. Every other application script reads that CSV through
`applications/us_data.h` and never touches `us_real.csv` directly.

`qvarma_data.txt` — Blazsek, Escribano and Licht's own three-series file — is
not read anywhere in that path. It is a diagnostic file only, read by
`applications/us_transformation_search.c` to check a candidate transformation
of `us_real.csv` against a published one, and by nothing downstream of it. The
background section at the end of this file is what that diagnostic
established; it explains why `qvarma_data.txt` is excluded, not how the
current variables are built.

## The five variables and their transformation

The variable set and which three are non-stationary follow
`papers/fabiano.txt` (Fabiano, "Evaluating Nonlinear Simulation Models with
Model Confidence Sets", section 4.2), which uses this same FRED-QD-sourced US
data: GDP, employment, CPI, the interest rate and energy demand.

| column in `out/us_system.csv` | `us_real.csv` source | transformation |
|---|---|---|
| `LogGDP` | `GDP` | `100 ln(GDP_t)` |
| `Employment` | `Unemployment` | `100 - Unemployment_t`, untransformed |
| `LogCPI` | `Cpi` | `100 ln(Cpi_t)` |
| `InterestRate` | `Fed_rate` | untransformed |
| `LogEnergyDemand` | `Des_Energy_demand` | `100 ln(x_t)` |

GDP, CPI and energy demand are Fabiano's own trend non-stationary,
co-integrated three, so each is a log-level rather than a growth rate or a
difference: differencing a co-integrated variable throws the relation away
(Sims et al. 1990; Cochrane 1997), so the levels are kept and the log is only
what stabilises their variance. Employment and the interest rate are
Fabiano's own stationary two and are left in their natural units,
untransformed — `static_model.md`'s own state vector, which this system
descends from, keeps `Emp_t` and `IR_t` unlogged for the same reason: nothing
about being stationary requires logging a series, and Fabiano logging all
five was a uniform-specification choice for his own VAR rather than a
requirement this project needs to inherit.

Column names carry no spaces, since they are read back by name in
`applications/us_data.h` and a name with a space in it is one more thing a
reader has to quote correctly.

Employment is not a column `us_real.csv` has; it is `100` minus
`Unemployment`, the same identity this file's background section records for
`_other/modified.csv`'s `Employement` column, and every unit root or
co-integration statistic used in this project is invariant to that affine map
— confirmed again once `docs/VARIABLE_STATUS.md`'s `Employment` numbers are
compared against the raw `Unemployment` figures in
`docs/RAW_SERIES_STATIONARITY.md`. Energy demand uses the deseasonalised
column: the raw `Energy_demand` column has a pronounced quarterly pattern a
level transformation does nothing to remove, and FRED-QD's own series are
seasonally adjusted to begin with.

## The sample

1973Q1 to 2019Q4, `us_real.csv` rows 0 to 187, 188 quarters.

Fabiano's own sample is reported as 1973:Q2 to 2019:Q4, T = 188. `us_real.csv`
row 0 is 1973Q1, and 188 quarters from there lands on 2019Q4 exactly, matching
the row count precisely; starting at 1973Q2 as the thesis text states would
give T = 187, one short. Either the thesis summary's quarter is off by one or
its own source starts a quarter later than this file — not established here —
but the row count is the fact this file can check, and 188 from 1973Q1 is what
it confirms. The end, 2019Q4, is the same choice Fabiano makes and this
project already made independently: excluding the COVID quarters, where
`Unemployment` goes from 3.8 to 13.07 and GDP growth to -10.03 in a single
quarter.

## Reproducing this

    make app-us_prepare_data

writes `out/us_system.csv`. It runs automatically before `app-us_stationarity`
and `app-us_cointegration`, and before `make applications`.

## Background: why `qvarma_data.txt` plays no part in this

Kept for the record. This is what `applications/us_transformation_search.c`
established when the project's variables were still being taken partly from
Blazsek, Escribano and Licht's own file; it is why that file was dropped from
the empirical pipeline rather than merely deprioritised.

### The two files, as the search saw them

| file | rows | span | contents |
|---|---|---|---|
| `us_real.csv` | 193 | 1973Q1 to 2021Q1 | raw levels: GDP, Consumption, Cpi, Investment, Unemployment, Energy_demand, Des_Energy_demand, Total_CO2_Emissions, Des_Total_CO2_Emissions, Fed_rate |
| `qvarma_data.txt` | 264 | 1954Q3 to 2020Q2 | the authors' own three series, tab separated, no header: GDP growth, inflation rate, effective federal funds rate |

`qvarma_data.txt` is the data Blazsek, Escribano and Licht estimated on,
supplied by them. That it is their exact dataset rather than a reconstruction
was checked against the paper's Table 2, and all six figures agreed:

| column | mean | sample sd | Table 2 |
|---|---|---|---|
| GDP growth | 2.9455 | 2.4139 | 2.9455, 2.4139 |
| inflation rate | 3.1082 | 2.1252 | 3.1082, 2.1252 |
| federal funds rate | 4.7419 | 3.5797 | 4.7419, 3.5797 |

The minima and maxima agreed too: -10.0275 and 8.7343, 0.2616 and 10.3687,
0.0600 and 17.7800. Note the standard deviation there is the `1/(n-1)` one;
`stats_var` in et_al returns the `1/n` variance, which gives 2.4093, 2.1212 and
3.5729 instead.

Row 74 of `qvarma_data.txt` is row 0 of `us_real.csv`, both being 1973Q1,
established two ways that agreed: by search, the smallest average absolute gap
anywhere in the alignment search, 0.003013, fell at offset 74; by arithmetic,
their row 0 is 1954Q3, and 1954Q3 plus 74 quarters is 1973Q1.

### Why each of their three series was rejected as a source

**GDP growth** reproduced almost exactly:
`100 (ln GDP_t - ln GDP_{t-4})` against their column, average absolute gap
0.003013, correlation 0.999861. This is not why it was dropped — it is dropped
now because the model uses the level, not the growth rate, per Fabiano's own
specification, not because `us_real.csv` could not supply the growth rate.

**The inflation rate could not be reproduced from `Cpi`.** No candidate came
close enough to be their construction:

| candidate | average gap | correlation |
|---|---|---|
| year-over-year log change of the four-quarter mean level | 0.713757 | 0.967051 |
| four-quarter mean of the year-over-year log change | 0.714957 | 0.967003 |
| `100 ln(Cpi_t / Cpi_{t-4})` | 0.760999 | 0.950423 |
| `400 ln(Cpi_t / Cpi_{t-1})` | 1.431469 | 0.806786 |

An average gap of 0.71 on a series whose own mean is 3.11 is a 23 per cent
relative error, against GDP's 0.003 from the same kind of transformation of
the same file. This is moot now that the model wants the level, `100 ln(Cpi_t)`,
built from `Cpi` directly with no reconstruction problem at all — the
inflation-rate reconstruction difficulty was specific to differencing a series
whose underlying index might be a different vintage; the level itself is not
in question.

**The federal funds rate was a level but not their column.** Against
`us_real.csv`'s `Fed_rate`:

| alignment | average gap | correlation |
|---|---|---|
| offset 74, the alignment GDP establishes | 0.697421 | 0.948107 |
| offset 75, one quarter later | 0.222434 | 0.992677 |

`Fed_rate` fit best one quarter away from the alignment the `GDP` column in
the same file establishes — a misalignment inside `us_real.csv` itself,
`Fed_rate` offset by one quarter relative to `GDP`. This project's earlier
answer was to take the rate from `qvarma_data.txt` instead and, if `Fed_rate`
were ever used, to correct the offset first. The current pipeline instead
takes `Fed_rate` as it stands, no shift, on the view that whatever internal
offset `us_real.csv` carries is a property of this file's own vintage of the
series and not something to correct by borrowing values from a different one;
`docs/VARIABLE_STATUS.md`'s own results are what say whether that choice
costs anything.

**Unemployment and energy demand were never in the authors' file at all,**
so both always came from `us_real.csv`, regardless of this question.
