# ABM validation with t-QVARMA and a Model Confidence Set

Which parameter configurations of an agent-based macroeconomic model produce
dynamics that cannot be told apart from US data. The comparison runs through an
auxiliary model fitted to both: the t-QVARMA of Blazsek, Escribano and Licht
(2023), a score-driven model that is nonlinear in its updating step while
keeping a co-integrated block in levels and Student-t innovations. The distance
between two economies is the distance between the impulse responses their fits
imply, and the Model Confidence Set of Hansen, Lunde and Nason (2011) decides
which configurations survive.

`docs/ABM_SYSTEM_MCS_VALIDATION.md` describes the validation procedure: what is
compared against what, why the comparison runs on impulse responses rather than
on fitted parameters, and the settings the confidence set is computed under.
`docs/DATA_DOCUMENTATION.md` records where the US series come from and how each
one is transformed.

## Built on et_al

This project depends on [et_al](https://github.com/AlfredoDiodati/et_al.), and
it is not optional: nothing general is implemented here. The auxiliary model
itself, the optimiser that fits it, the tests run on the residuals and the
Model Confidence Set are all et_al's:

- `et_al./sd/qvarma.h` — the t-QVARMA model: filter, link, fit, impulse
  responses and sign-restricted bands
- `et_al./solver/lbfgs.h` — the limited-memory BFGS the fit descends with
- `et_al./mcs.h` — the Model Confidence Set
- `et_al./stats.h` — sample statistics, Ljung-Box, quantiles
- `et_al./unit_root.h` — unit root tests and the Newey-West bandwidth
- `et_al./qlr_test.h` — the quasi-likelihood ratio test for the absence of
  score-driven dynamics, with the tabulated critical values compiled in
- `et_al./ad.h`, `et_al./linalg/`, `et_al./frame/`, `et_al./random.h` —
  automatic differentiation, linear algebra, dataframes with CSV and .Rdata
  readers, random number generation

What this repository carries is the application layer alone: the scripts that
prepare the data, run the fits and produce the results, plus two headers that
describe this project's own data and nothing else. Any primitive that turns out
to be general belongs in et_al rather than here.

Build flags come from `pkg-config --cflags --libs et_al.-core`. Reinstalling
et_al after a change is what makes the next `make` here recompile.

## Requirements

- A C11 compiler with OpenMP, and OpenBLAS
- et_al, installed so that `pkg-config et_al.-core` resolves
- Python with `polars`, `plotly` and `kaleido`, for the figures only; plotly
  writes the PDFs through kaleido

## Layout

    applications/     the scripts that produce results, plus us_data.h and
                      abm_system.h, which describe this project's own data
    tests/            what verifies the auxiliary model still computes what it
                      claims to
    dataset/          us_real.csv, the raw US series; the ABM's own simulated
                      output goes here too but is not tracked
    out/              every result, written here rather than printed
    docs/             the write-up and the reference documentation

## Running the pipeline

Each step reads what the previous one wrote, so the order matters. The Makefile
already encodes the dependencies, so asking for a later step builds and runs the
earlier ones.

    make app-us_prepare_data              the five US variables, 1973Q1 to 2019Q4
    make app-us_qvarma_employment_change  the auxiliary model on the real data
    make app-abm_system_extract           the simulated data into the same layout
    make app-abm_system_fit_qvarma        10,800 fits, one per simulated replicate
    make app-abm_system_mse_qvarma        impulse responses and the loss table
    make app-abm_system_mcs               the confidence set itself
    make app-abm_system_winner_irf        the surviving model's responses, with bands

    python applications/abm_system_winner_irf_plots.py

The fitting step is the long one. Every fit is cached to its own file the moment
it finishes, so an interrupted run resumes rather than starting over.

Results are written to `out/`, never printed.

## Tests

    make test         the auxiliary model computes what it claims to
    make test-stress  the same, including the slow simulation checks
    make study        parameter recovery from known truths, over sample sizes,
                      model shapes and parameter regimes

## Data not in this repository

The ABM's raw output is 229 MB of `.Rdata` under `dataset/simulated/`, one file
per parameter configuration, each holding 108 Monte Carlo replicates of 600
periods. It is not tracked. `dataset/abm_system/`, the five-variable CSVs
derived from it, and `out/abm_system_fit_qvarma/`, the cached fits, are not
tracked either: both are rebuilt from the `.Rdata` by the targets above.

`dataset/us_real.csv` is tracked, and everything derived from it is small enough
to keep, so the real-data half of the analysis reproduces from a fresh clone on
its own.
# scoredriven_validation
