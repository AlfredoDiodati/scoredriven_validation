# The US data and how it is transformed

## Source

`dataset/us_real.csv`: quarterly US series, 193 rows from 1973Q1 to 2021Q1,
with the columns `Quarter`, `GDP`, `Consumption`, `Cpi`, `Investment`,
`Unemployment`, `Energy_demand`, `Des_Energy_demand`, `Total_CO2_Emissions`,
`Des_Total_CO2_Emissions` and `Fed_rate`. `Quarter` is written as the last
month of the quarter, `1973 March` for 1973Q1. The file was received as US data
built from FRED-QD; the download date and vintage of each column are not
recorded in this repository. Five of the columns are used.

## Two steps

The US data reaches the auxiliary models in two steps, and the simulated data
is stored in the form the second step produces.

**Step 1, levels.** `applications/us_prepare_data.c` reads `us_real.csv` and
writes `out/us_system.csv`: five columns plus `Quarter`, rows 0 to 187 of the
source, 1973Q1 to 2019Q4, 188 quarters. Nothing else reads `us_real.csv`;
every other script loads `out/us_system.csv` through `load_us_system` in
`applications/us_data.h`.

| column in `out/us_system.csv` | source column | transformation |
|---|---|---|
| `LogGDP` | `GDP` | `100 ln(GDP_t)` |
| `Employment` | `Unemployment` | `100 - Unemployment_t`, the employment rate in per cent |
| `LogCPI` | `Cpi` | `100 ln(Cpi_t)` |
| `InterestRate` | `Fed_rate` | none, per cent |
| `LogEnergyDemand` | `Des_Energy_demand` | `100 ln(x_t)` |

GDP, the price level and energy demand trend upwards and are put in logs to
stabilise their variance; employment and the interest rate are already rates
and are kept in their own units. Energy demand uses the seasonally adjusted
column, since the raw `Energy_demand` has a strong quarterly pattern.

**Step 2, the five series the models are fitted on.** The first difference of
every level except the interest rate, one quarter shorter:

| series | built from step 1 | units |
|---|---|---|
| `GDP_growth` | `LogGDP_t - LogGDP_{t-1}` | 100 times the quarterly log change |
| `EN_growth` | `LogEnergyDemand_t - LogEnergyDemand_{t-1}` | 100 times the quarterly log change |
| `Employment_change` | `Employment_t - Employment_{t-1}` | percentage points |
| `Inflation` | `LogCPI_t - LogCPI_{t-1}` | 100 times the quarterly log change |
| `InterestRate` | `InterestRate_t` | per cent |

The sample is then 1973Q2 to 2019Q4, 187 quarters. This is `build_block` in
`applications/us_qvarma_spec_choice.c`, which the t-QVARMA is fitted on, and
the same arithmetic as `us_system_growth_block` in `applications/us_data.h`,
which the local projections use. The local projections drop three further
quarters, where the four-quarter moving average of log GDP that decides their
state does not exist: 184 quarters, 1974Q1 to 2019Q4.

Every simulated run is stored in this second form, by the same arithmetic
applied to the model's own levels (`applications/abm_system.h`;
`docs/ABM_SYSTEM_SIMULATION.md`, "The observables").
`tests/abm_system_layout.c` checks that the simulated route and the US route
agree.

In the t-QVARMA, the first three series are treated as stationary and
inflation and the interest rate as co-integrated with rank one, the Fisher
relation (`docs/MCS_VALIDATION_PROCEDURE.md`, "Auxiliary specification").

## The sample

The sample ends at 2019Q4 to exclude the COVID quarters, where `Unemployment`
goes from 3.8 to 13.07 and GDP growth to -10.03 in a single quarter. It starts
at the first row of `us_real.csv`, 1973Q1.

## A caveat on the interest rate

`Fed_rate` was compared with the effective federal funds rate series Blazsek,
Escribano and Licht (2023) estimated on, which the authors supplied and which
is not distributed here. Aligned on the quarter that GDP growth gives (the two
files' GDP growth agree to a mean absolute gap of 0.003 at that alignment),
the mean absolute gap between the two interest rates is 0.70 percentage points,
and one quarter later it is 0.22. `Fed_rate` therefore looks shifted by one
quarter relative to `GDP` within `us_real.csv`. The pipeline uses `Fed_rate`
as it stands, unshifted. What that costs has not been measured.

## Reproducing this

    make app-us_prepare_data

writes `out/us_system.csv`. The `make app-<name>` targets that read it depend
on `app-us_prepare_data`, directly or through the US fit, so they regenerate it
rather than read a stale copy. The studies under `studies/` do not, and read
whatever is on disk.
