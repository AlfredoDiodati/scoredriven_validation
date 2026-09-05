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
one is transformed. `docs/ABM_SYSTEM_SIMULATION.md` describes the design the
simulations run over and how they are stored, and
`docs/DSK_MODEL_CHANGES.md` records what this project changed in the simulator
and what those changes are measured and tested to leave alone.
`docs/Model_Simulation.Rmd` is not this project's writing: it is the brief the
simulation workflow was derived from, kept for reference. Everything else
written down lives in `docs/` too; no directory in this tree carries a note of
its own.

## The agent-based model

The economy simulated here is not this project's. It is the DSK
stock-flow-consistent model - Dystopian Schumpeter meeting Keynes - and the code
under `model/dsk_sfc` is the code its own authors released alongside their
peer-reviewed paper:

    https://github.com/CoMoS-SA/Reissl_2025
    commit 611ff9cb44348baa55be1bc315eefe2c117ccd44

Cite that repository and that paper for the model. The economics is theirs, and
this project starts from their program rather than from a reimplementation of
it, so that what is validated here is the published model and not a second
opinion about it.

### Made to run a million times

The validation needs 1000 parameter configurations simulated 1000 times each. At
the speed the published code runs, one 600-period simulation takes 34.3 seconds,
which puts the experiment at roughly 9800 core-hours: 73 days of the desktop it
was measured on. It now takes 0.528 seconds, 147 core-hours, 1.7 days - 65 times
faster.

None of that came from changing what the model computes. It came from how the
program stores its data and how often it repeats itself:

- **The build.** The published `CMakeLists.txt` selects a debug build and adds
  no optimisation flag, so a default compile carries none. Building it as
  release is worth 3.8x on its own, before anything in the source is touched.
- **The direction the arrays are read in.** The model keeps each firm's stock of
  machines in arrays indexed by period, then machine supplier, then firm.
  Several of its loops worked through one firm at a time, which means reading a
  single number out of each of four hundred separate places in memory, tens of
  millions of times a run. Running the firms on the inside of those loops reads
  the same numbers consecutively instead, which is the difference between a
  memory system that can anticipate the next read and one that cannot.
- **Answers that were being recomputed.** A bank pushing its non-customers to
  the end of its credit ranking rescanned the whole ranking once for each of
  them, when the answer only moves by one each time. A firm owns machines of
  very few vintages - 98 per cent of the firm-and-vintage pairs the model sweeps
  every period hold nothing - and every one of those was being visited and
  multiplied by zero.
- **Work whose result nobody reads.** Two of the arrays copied every period feed
  a single function that runs only under a shock setting this experiment does
  not use. The ages of machines nobody owns were being cleared and maintained,
  and every place that reads an age asks for it only where a machine is owned.
- **Branches the processor cannot guess.** The loop that averages each firm's
  machines over its vintages reads 52 million counts a run and acts on the 2%
  that are not zero, and which 2% is not a pattern anything can learn. Testing
  four firms at a time, and skipping all four when none holds anything, cut that
  loop by two thirds.
- **Arrays that existed to be destroyed.** Three of the model's nine
  machine-vintage arrays, 19 MB each, were whole copies kept so that one
  function could draw down a column of them; each is now a scratch of a handful
  of numbers. That took the memory a run needs from 163 MB to 74, which matters
  more than it sounds: eight runs at once had been competing for one cache, and
  the number of runs the machine finishes in a minute went up by half.

### Why the results are still the published model's

Every change has to leave the model's output identical to the published code's,
byte for byte. `make test-dsk_build_equivalence` compiles their unmodified
source, runs both programs over the same seeds, and compares all 3 MB of each
run's output; a change that moves a single digit does not go in.

Three tests do that, in the three directions a mistake could hide: the aggregate
output over five seeds, the twelve per-firm files the model writes under `-f 1`,
and four points of the parameter design away from the baseline calibration. All
three run with the model's shock channels off, which is how the experiment runs
them; `docs/DSK_MODEL_CHANGES.md` names exactly what that leaves uncovered.

That is a stronger guarantee than statistical agreement. Two runs that agree byte
for byte cannot be told apart by any test, so the question of whether the faster
version drifts away from the published one does not arise. Seventeen further
changes were tried, measured, and dropped for not being faster.
`docs/DSK_MODEL_CHANGES.md` holds the whole record: every change, every timing
and the setup it was taken under, every rejected attempt, and what the test does
and does not cover.

### Using the simulator

    make model                    builds it, no cmake needed
    make model-upstream           the unmodified reference, for the test

One run, written into an `output/` directory beside the path the executable was
invoked as:

    ./model/dsk_sfc/dsk_SFC model/dsk_sfc/dsk_sfc_inputs.json -r myrun -s 1 -f 0 -c 0 -v 0

The whole experiment goes through the driver rather than through that command:

    make app-abm_system_design            draws dataset/abm_system_design.csv
    make app-abm_system_simulate          runs the design, writes the archives

`abm_system_simulate` takes an optional first, last and replication count
(`./bin/abm_system_simulate 1 10 100` is configurations 1 to 10, 100
replications each), so one invocation is the whole design on one machine and one
configuration per job on a cluster. It writes
`dataset/abm_system/cop_NNNN/batch_NNN.npz`, ten replications to an archive, and
skips any archive already on disk, which is what makes an interrupted run
resumable. `DSK_EXECUTABLE` and `DSK_BASE_JSON` override which binary and which
baseline parameter file it uses.

## Layout

    applications/     the scripts that produce results, plus us_data.h and
                      abm_system.h, which describe this project's own data
    model/dsk_sfc/    the DSK simulator itself, a copy of upstream's source
                      with the speed work of docs/DSK_MODEL_CHANGES.md applied
    tests/            what verifies the auxiliary model still computes what it
                      claims to, and what verifies this copy of the simulator still
                      matches upstream byte for byte
    dataset/          us_real.csv, the raw US series; the ABM's own simulated
                      output goes here too but is not tracked
    out/              every result, written here rather than printed
    docs/             the write-up and the reference documentation

## Requirements

- A C11 compiler with OpenMP, and OpenBLAS
- A C++11 compiler, for the simulator under `model/dsk_sfc`; `make
  model` builds it without cmake
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
    make app-abm_system_fit_qvarma        one fit per simulated replication
    make app-abm_system_mse_qvarma        impulse responses and the loss table
    make app-abm_system_mcs               the confidence set itself
    make app-abm_system_winner_irf        the surviving model's responses, with bands

    python applications/abm_system_winner_irf_plots.py

`app-abm_system_extract` is the older route into the simulated dataset: it reads
the 108 `.Rdata` files under `dataset/simulated/` and converts them. The
Latin hypercube experiment described in `docs/ABM_SYSTEM_SIMULATION.md` replaces
it with `app-abm_system_design` and `app-abm_system_simulate`, which produce the
same archives from the simulator directly. That replacement has not been made in
the Makefile: `app-abm_system_fit_qvarma` still depends on
`app-abm_system_extract`, and `abm_system_fit_qvarma` fits every subdirectory of
`dataset/abm_system/` regardless of what wrote it, so the two datasets must not
sit there at once.

The fitting step is the long one. Every fit is cached to its own file the moment
it finishes, so an interrupted run resumes rather than starting over.

Results are written to `out/`, never printed.

## Tests

    make test         the auxiliary model computes what it claims to
    make test-stress  the same, including the slow simulation checks
    make study        parameter recovery from known truths, over sample sizes,
                      model shapes and parameter regimes

The simulator has two of its own, and they are the gate on any change
to it:

    make test-dsk_build_equivalence         every byte of the aggregate output
                                           matches upstream, over five seeds
    make test-dsk_full_output_equivalence   the same for the twelve per-firm
                                           files the model writes under -f 1
    make test-dsk_design_equivalence        the same away from the baseline, at
                                           four points of the design
    make test-dsk_long_path                 the filename bug stays fixed

`test-dsk_build_equivalence` builds the reference itself, so it takes about
three minutes: the unmodified binary carries no optimisation and one of its runs
is 35 seconds.

