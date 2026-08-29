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

## Layout

    applications/     the scripts that produce results, plus us_data.h and
                      abm_system.h, which describe this project's own data
    tests/            what verifies the auxiliary model still computes what it
                      claims to
    dataset/          us_real.csv, the raw US series; the ABM's own simulated
                      output goes here too but is not tracked
    out/              every result, written here rather than printed
    docs/             the write-up and the reference documentation

## Requirements

- A C11 compiler with OpenMP, and OpenBLAS
- et_al, installed so that `pkg-config et_al.-core` resolves
- Python with `polars`, `plotly` and `kaleido`, for the figures only; plotly
  writes the PDFs through kaleido

## Built on et_al

This project depends on [et_al](https://github.com/AlfredoDiodati/et_al.): nothing general is implemented here. The auxiliary model
itself, the optimiser that fits it, the tests run on the residuals and the
Model Confidence Set are all et_al's:

- `et_al./sd/qvarma.h` — the t-QVARMA model: filter, link, fit, impulse
  responses and sign-restricted bands
- `et_al./solver/lbfgs.h` — the limited-memory BFGS the fit descends with
- `et_al./mcs.h` — the Model Confidence Set

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

## Installing

Three things have to be in place before anything builds: OpenBLAS, et_al, and
this repository. The Python packages come last and are needed only for the
figures.

OpenBLAS first, because et_al links against it and against nothing else:

    sudo apt install libopenblas-dev        # or: sudo pacman -S openblas

Then et_al itself. It is header-only, so installing it means copying its
headers and writing the pkg-config files that carry the OpenBLAS flags, not
building a library:

    git clone https://github.com/AlfredoDiodati/et_al.
    cd et_al.
    sudo make install-model PREFIX=/usr/local

`install-model` depends on `install-core`, so that one command installs both
tiers, and it is the model tier this project needs: `sd/qvarma.h` lives there,
and a core-only install leaves the auxiliary model out. Headers go to
`$(PREFIX)/include/et_al./`, keeping et_al's own directory structure so that
its internal relative includes still resolve, and the `.pc` files go to
`$(PREFIX)/lib/pkgconfig/`. `et_al.-model.pc` declares `Requires: et_al.-core`,
so naming either one pulls in both. `PREFIX` defaults to `/usr/local`; keep it
there unless you have a reason not to, since the sources here spell their
includes `<et_al./sd/qvarma.h>`, which resolves through the compiler's default
search path rather than through et_al's own `-I`, and a prefix outside that
path then needs `CPATH` set as well as `PKG_CONFIG_PATH`. The install prints
the `export PKG_CONFIG_PATH=...` line itself when the prefix chosen is one
pkg-config does not already search.

Check that it resolves before going further:

    pkg-config --cflags --libs et_al.-core

Then this repository, and the tests, which need nothing beyond the above:

    git clone git@github.com:AlfredoDiodati/scoredriven_validation.git
    cd scoredriven_validation
    make test

The figures are the one part that leaves C:

    pip install polars plotly kaleido

`make uninstall-core` in the et_al clone, with the same `PREFIX`, reverses the
install and removes the model tier with it.

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

