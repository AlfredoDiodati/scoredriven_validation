# ABM validation with t-QVARMA and a Model Confidence Set

Which parameter configurations of an agent-based macroeconomic model sit
closest to US data, and which of them cannot be told apart from each other at
that distance. The comparison runs through an auxiliary model fitted to both:
the t-QVARMA of Blazsek, Escribano and Licht (2023), a score-driven model that
is nonlinear in its updating step while keeping a co-integrated block in levels
and Student-t innovations. The distance between two economies is the distance
between the impulse responses their fits imply, and the Model Confidence Set of
Hansen, Lunde and Nason (2011) decides which configurations survive.

**Status.** The distance is the squared error between impulse responses minus
each configuration's estimated response variance
(`docs/MCS_VALIDATION_PROCEDURE.md`). Nothing runs it yet: the Monte Carlo
study is being written first, and the empirical validation on the US data will
follow once it is finished. Every result produced so far used the naive loss,
the mean absolute error with no variance correction, which biases the
confidence set towards configurations whose responses are less noisy; its
code, its outputs and the studies built on its winner are under
`old/naive_mae_loss/` (`docs/OLD_EXPERIMENTS.md`, "The naive loss").

The six steps that produce the data and the fits the procedure compares, from
simulation to the fitted auxiliary models, are listed under "The main pipeline"
below.

The documentation, all of it under `docs/`:

- `docs/MCS_VALIDATION_PROCEDURE.md`: the validation procedure, the same for
  the US data and for the Monte Carlo: what is compared against what, why the
  comparison runs on impulse responses, and the confidence set's settings.
- `docs/MONTECARLO_VALIDATION.md`: the check that the procedure recovers an
  answer it already knows, with every simulated run of every configuration in
  turn standing in for the US data.
- `docs/DATA_DOCUMENTATION.md`: the US series and how each is transformed.
- `docs/ABM_SYSTEM_SIMULATION.md`: the parameter design, the simulations and
  how they are stored. `docs/ABM_SYSTEM_DISTRIBUTION.md`: the distribution of
  the simulated series.
- `docs/DSK_MODEL_CHANGES.md` and `docs/DSK_LONG_HORIZON.md`: what this project
  changed in the simulator, and the tests that show the changes leave its
  output alone.
- `docs/QVARMA_FIT_STUDIES.md`: the studies of the t-QVARMA fits.
  `docs/BENCHMARKS.md`: how long a fit takes, and the throughput study.
- `docs/OLD_EXPERIMENTS.md`: the discarded experiments, each with its reason
  and its own `docs/OLD_*.md` write-up.
- `docs/Model_Simulation.Rmd` is not this project's writing: it is the R
  Markdown brief the simulation workflow was derived from, kept for reference.

Every reference cited is listed under "References" at the end.

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

Those are one machine's numbers - an 8-core AMD Ryzen 7 4800H with 7.5 GB, built
with g++ 15.2.1 - and they do not carry to another without being measured again.
`docs/DSK_MODEL_CHANGES.md` gives the full specification and says which parts of
it each figure depends on; the cache size in particular decides how much of the
speedup survives when several runs share a machine.

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
byte for byte. `make model-upstream` compiles the authors' unmodified source,
the tests run both programs on the same seeds, and the files they write must
match exactly. A change that moves a single digit does not go in.

Twelve tests, each answering one question; "Tests" below lists them by name.
Four came with the long-horizon work of `docs/DSK_LONG_HORIZON.md`, which added
mechanisms that fire only past horizons no 600-period run reaches and so cannot
be checked by comparison against upstream alone.

- **Does the economy-wide output match?** The 83-column results file, compared
  byte for byte. The unmodified program is also run twice per seed and required
  to agree with itself, since a reference that varied between runs would make
  every other comparison meaningless.
- **Does the per-firm output match?** Totals can stay the same while individual
  firms change. Under `-f 1` the model writes the state of every firm and bank,
  14 files and 28 MB a run, and all of it is compared.
- **Does it still match away from the default parameters?** Four points of the
  thousand-row parameter design, including the corners of the box it covers.
  The parameters decide which branches the code takes.
- **Does the faster version read memory it does not own?** The same source
  compiled under AddressSanitizer and UndefinedBehaviorSanitizer, at three
  parameter settings. Comparing output cannot find an index that runs one place
  past the end of an array and happens to read a plausible number.
- **How small a difference would those comparisons actually catch?** The files
  are printed to ten decimals, so in principle a smaller difference could hide
  in the rounding. Each input parameter is nudged by smaller and smaller
  amounts until the output changes. For eight of the nine
  parameters, changing them by a factor of 1.000000000000001 is enough, and the
  change reaches the output within a few simulated periods. That is the last
  digit a double-precision number can hold, so a difference in the arithmetic
  would show in the files.
- **Does the filename bug stay fixed?** The published code writes
  output filenames into 64-byte buffers and corrupts memory past about 26
  characters of directory path, which is shorter than a cluster scratch
  directory.
- **Does it still reproduce the dataset the experiment was run on?** A sample of
  the design, run again against both builds and against the stored archives:
  every one of the 5 x 400 numbers must equal the stored one exactly, not to a
  tolerance.
- **Do the three long-horizon rescalings leave the economy alone?** One test
  each for counting money, machines and the consumption good in bigger units.
  Redenominating money has to be exact; the machine lots round, so that one is
  required to be small rather than exact and the bound is stated.
- **Is the bulk cancellation draw the distribution it claims to be?** It calls
  the header and the model's own generators directly rather than running the
  simulator.

All of them run with the model's shock channels off, which is how the experiment
runs it. `docs/DSK_MODEL_CHANGES.md` explains each test in full and names
exactly what none of them covers, and `docs/DSK_LONG_HORIZON.md` covers the last
five.

Byte equality is a stronger guarantee than statistical agreement. Two runs that
write the same file cannot be told apart by any test, so the question of whether
the faster version drifts away from the published one does not arise. Seventeen
further changes were tried, measured, and dropped for not being faster. The same
document holds the whole record: every change, every timing and the setup it was
taken under, and every rejected attempt.

### Using the simulator

    make model                    builds it, no cmake needed
    make model-upstream           the unmodified reference, for the test

One run, written into an `output/` directory beside the path the executable was
invoked as:

    ./model/dsk_sfc/dsk_SFC model/dsk_sfc/dsk_sfc_inputs.json -r myrun -s 1 -f 0 -c 0 -v 0

`-f 1` adds twelve per-firm files. `-m 1` adds two more holding each firm's own
output quantity, the cross-section GDP is the sum of; it is this project's
addition and `docs/DSK_MODEL_CHANGES.md` records what it changed and what says
it changed no number.

The whole experiment goes through the driver rather than through that command;
these are steps 1 and 3 of "The main pipeline" below:

    make app-abm_system_design            draws dataset/abm_system_design.csv
    make bin/abm_system_simulate          builds the driver, does not run it

    ./applications/abm_system_simulate_all.sh          the whole design, 8 at a time
    ./applications/abm_system_simulate_all.sh 8 1 16 20   a pilot, 16 configurations

Building the driver and running it are separate because with no arguments it
simulates the whole design, a million model runs and about two days of one
machine. That is not something a make target should start.

`abm_system_simulate` takes an optional first, last and replication count
(`./bin/abm_system_simulate 1 10 100` is configurations 1 to 10, 100
replications each), so one invocation is the whole design on one machine and one
configuration per job on a cluster. It writes
`dataset/abm_system/cop_NNNN/batch_NNN.npz`, ten replications to a compressed
archive, and skips any archive already on disk, which is what makes an
interrupted run resumable and a re-submitted job harmless. `DSK_EXECUTABLE`,
`DSK_BASE_JSON` and `ABM_SYSTEM_OUTPUT_DIR` override which binary, which
baseline parameter file and which output directory it uses.

`abm_system_simulate_all.sh` runs several of those over disjoint stretches of
the design, because one process is serial and this machine finishes most runs
per minute with eight at once. It holds the machine awake for the duration and
writes what it used and what it produced to
`out/abm_system_simulate_all_provenance.txt`. The full experiment has been run:
1,000,000 replications, none rejected, 41.6 hours and 15 GB at 400 runs a
minute. `docs/ABM_SYSTEM_SIMULATION.md` has the setup behind those numbers.

## Layout

    applications/     the scripts that produce results, plus us_data.h and
                      abm_system.h, which describe this project's own data
    model/dsk_sfc/    the DSK simulator itself, a copy of upstream's source
                      with the speed work of docs/DSK_MODEL_CHANGES.md applied
    montecarlo/       the Monte Carlo pipeline, with its own out/
    tests/            pass or fail. A failure stops make test and the build
    benchmarks/       how long something takes. Never part of make test: a
                      function that returns the wrong answer quickly is not fast
    studies/          a question about behaviour, answered into out/. No verdict,
                      so nothing here gates anything
    dataset/          us_real.csv, the raw US series, and the parameter design.
                      The simulated output goes here too but is not tracked:
                      abm_system/ is what the fitting stage reads,
                      abm_system_rdata/ is the older .Rdata-derived copy kept out
                      of the way so only one of the two is ever fitted, and
                      throughput/ is the separate dataset the throughput study
                      times fits on
    old/              experiments discarded from the main pipeline or the
                      Monte Carlo, one folder per experiment. Still buildable
                      with make old-<stem>; docs/OLD_EXPERIMENTS.md says why
                      each was discarded, and their write-ups stay in docs/
    out/              every result, written here rather than printed
    bin/              the compiled programs, rebuilt by make, not tracked
    docs/             the write-up and the reference documentation

A script's directory says what kind of thing it is, and its name says what it
does rather than how it was once done. Three conventions follow:

- A name states the question, not a method that may change. The loss table is
  `abm_system_irf_loss`, not `abm_system_mae_qvarma`.
- A suffix that distinguished two things is removed when only one is left.
- A file under `tests/` does not repeat the word test in its name; the directory
  already says that. The name states what is being verified, specifically enough
  that a sibling covering a different aspect of the same subject cannot collide
  with it, and a file a script writes is named after the script that wrote it.

## Requirements

- A C11 compiler with OpenMP, and OpenBLAS
- A C++11 compiler, for the simulator under `model/dsk_sfc`; `make
  model` builds it without cmake
- et_al, installed so that `pkg-config et_al.-core` resolves
- Python, only to rerun the discarded experiments under `old/`: `polars`,
  `plotly` and `kaleido` for their figures (plotly writes the PDFs through
  kaleido), and `numpy`, `scipy` and `pandas` for the compressed-response
  prototype of `docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md`

### The platform everything here was run on

One machine, and nothing in this repository except the benchmarks of
2026-08-29 (`docs/BENCHMARKS.md`) has been run on any other:

    $ uname -srmo
    Linux 6.19.14-arch1-1 x86_64 GNU/Linux
    $ cat /etc/os-release | head -2
    NAME="EndeavourOS"
    PRETTY_NAME="EndeavourOS"
    $ g++ --version | head -1
    g++ (GCC) 15.2.1 20260209
    $ gcc --version | head -1
    gcc (GCC) 15.2.1 20260209
    $ ldd --version | head -1
    ldd (GNU libc) 2.43
    $ make --version | head -1
    GNU Make 4.4.1
    $ bash --version | head -1
    GNU bash, version 5.3.9(1)-release (x86_64-pc-linux-gnu)

The hardware is an AMD Ryzen 7 4800H, 8 cores with 2 threads each, 7.5 GB of
memory; `docs/ABM_SYSTEM_SIMULATION.md` gives the cache layout, which is what
the throughput figures turn on.

Every timing, every equivalence test and the million-run experiment itself were
produced there. macOS and Windows are untested. What follows is what would have
to change, read off the code rather than tried, so treat it as a starting point
and not as instructions known to work.

### macOS

The simulator's own source already expects to be built there: it selects the
three-argument `mkdir` under `__APPLE__` exactly as it does under `__linux__`.
One thing it does not do on macOS is time a run out. The authors guard

    signal(SIGALRM, catchAlarm);
    alarm(T*2);

with `#ifdef __linux__` and say in their own comment that it works on Linux
only. A run that gets stuck on Linux is killed after `2T` seconds and recorded
as a failure; on macOS it would hang, and the driver waiting on it would wait
with it.

What has to change to build and run:

- `model/dsk_sfc/build.sh` passes `-msse`, which is an x86 flag. On an Apple
  Silicon machine it has to go. On an Intel Mac it can stay.
- The same script calls `nproc`, which is GNU coreutils and is not present.
  `sysctl -n hw.ncpu` is the equivalent, or install coreutils and use `gnproc`.
- `applications/abm_system_simulate_all.sh` uses `readlink -f`, `stat -c`,
  `md5sum`, `du -sh`, `date -d @epoch` and `systemd-inhibit`. Only the last is
  optional in the script, which already skips it when it is missing; the other
  five are GNU spellings that BSD userland writes differently (`stat -f`,
  `md5`, `date -r`). Installing GNU coreutils with Homebrew and putting its
  `gnubin` on the PATH is the shortest route. Without the inhibitor the machine
  can sleep mid-run, which loses at most the batch in flight but stops the
  clock, so disable sleep another way for a run of this length.
- OpenMP is not in Apple's clang. `brew install libomp` and point the compiler
  at it, or build with `brew`'s gcc.

The driver `bin/abm_system_simulate` itself uses no Linux-only interface. It is
the shell launcher around it and the build script that carry the GNU
assumptions.

### Windows

The simulator has a `#else` branch calling the one-argument `mkdir`, so the
authors intended it to compile there, but nothing else in this repository
does. The Makefile, both shell scripts and the driver's use of `symlink` and
`fork` through `system` are POSIX.

The route that needs no porting is WSL2 with a Linux distribution inside it,
where everything above applies unchanged. Running natively under MSYS2 or
Cygwin would need the Makefile and `abm_system_simulate_all.sh` rewritten and
is not something to attempt for a first run.

### If you reproduce the simulation elsewhere

Expect the same statistics, not the same bytes. The model calls `log`, `exp`
and `pow` from the C library, and those are allowed to differ in the last bit
between one library and another. The equality the tests in
`docs/DSK_MODEL_CHANGES.md` establish is between two builds on one machine; it
is not a claim about two machines. A run on another platform that gives
different digits in a late decimal place is behaving as expected, and one that
gives visibly different series is not.

## Built on et_al

This project depends on [et_al](https://github.com/AlfredoDiodati/et_al.):
nothing general is implemented here. The auxiliary model itself, the optimiser
that fits it, the tests run on the residuals and the Model Confidence Set are
all et_al's:

- `et_al./sd/qvarma.h`: the t-QVARMA model, with its filter, link, fit,
  impulse responses and sign-restricted bands
- `et_al./lp/lp.h`: the linear and state-dependent local projections
- `et_al./solver/lbfgs.h`: the limited-memory BFGS the fit descends with
- `et_al./mcs.h`: the Model Confidence Set
- `et_al./qlr_test.h`: the quasi-likelihood ratio test for the absence of
  score-driven dynamics, with the tabulated critical values compiled in
- `et_al./ad.h`, `et_al./linalg/`, `et_al./frame/`, `et_al./random.h`:
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
figures of the discarded experiments under `old/`.

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

The figures of the discarded experiments are the one part that leaves C:

    pip install polars plotly kaleido

`make uninstall-core` in the et_al clone, with the same `PREFIX`, reverses the
install and removes the model tier with it.

## The main pipeline

Six steps take the project from the parameter design to the fitted auxiliary
models the validation procedure compares. Each step reads what the steps before
it wrote. The steps after them, the loss under the variance correction and the
confidence set on the US data, are not written yet; the naive-loss versions are
under `old/naive_mae_loss/applications/`.

| step | script | command |
|---|---|---|
| 1. Draw the parameter design | `applications/abm_system_design.c` | `make app-abm_system_design` |
| 2. Build the simulator | `model/dsk_sfc/` | `make model` |
| 3. Simulate every configuration | `applications/abm_system_simulate_all.sh`, which runs `applications/abm_system_simulate.c` | `make bin/abm_system_simulate`, then `./applications/abm_system_simulate_all.sh` |
| 4. Prepare the US data | `applications/us_prepare_data.c` | `make app-us_prepare_data` |
| 5. Fit the auxiliary model to the US data | `applications/us_qvarma_spec_choice.c` | `make app-us_qvarma_spec_choice` |
| 6. Fit the auxiliary model to every simulation | `applications/abm_system_fit_qvarma.c` | `make app-abm_system_fit_qvarma` |

What each step is for:

1. Draws the 1000 parameter configurations the experiment compares, a Latin
   hypercube over nine of the model's parameters.
   `docs/ABM_SYSTEM_SIMULATION.md` describes the design.
2. Builds the DSK simulator described above.
3. Runs the simulator 1000 times for each configuration and keeps the five
   series the auxiliary model is fitted on. It took 41.6 hours here, which is
   why no make target starts it.
4. Builds the five US series from `dataset/us_real.csv`.
   `docs/DATA_DOCUMENTATION.md` describes the transformations.
5. Fits the t-QVARMA to the US data. It fits (1,1,2) and (1,1,4) and reports
   both, which is how (1,1,2) was chosen; the (1,1,2) fit is the benchmark every
   simulation is compared against.
6. Fits the t-QVARMA(1,1,2) to each of the 1,000,000 simulated replicates. It
   caches every fit as it finishes and resumes from the cache, so it can be run
   several times; it was run eight times here, hours each.

`make app-<name>` also reruns the steps listed as its prerequisites in the
Makefile: step 5 reruns 4. Steps 1, 3 and 6 are never rerun that way.

The pipeline is C. The discarded experiments under `old/` have a few Python
scripts of their own, which their documents list. Results are written to files,
never printed.

### Outputs of the main pipeline

| file | step | what it holds |
|---|---|---|
| `dataset/abm_system_design.csv` | 1 | the 1000 configurations, one row each, nine parameters |
| `out/abm_system_design.txt` | 1 | what the last run of step 1 did |
| `dataset/abm_system/cop_NNNN/batch_NNN.npz` | 3 | the five simulated series, ten replicates per compressed archive; 15 GB, not tracked by git |
| `out/abm_system_simulate_all_provenance.txt` | 3 | checksums of the model and the design, shard ranges, wall time, stored size |
| `out/abm_system_simulate_manifest.txt` | 3 | replicates completed and rejected, per configuration |
| `out/abm_system_simulate/` | 3 | one log per shard, and the reason for any rejected replicate |
| `out/us_system.csv` | 4 | the US series |
| `out/us_qvarma_spec_choice_p1q1r2_fit.json` | 5 | the US benchmark fit |
| `out/us_qvarma_spec_choice.txt` | 5 | the (1,1,2) against (1,1,4) comparison, with residual checks |
| `out/abm_system_fit_qvarma/cop_NNNN/replicate_NNN_p1q1r2_fit.json` | 6 | one fitted parameter set per replicate, which is also the cache; not tracked by git |
| `out/abm_system_fit_qvarma/cop_NNNN/lineage.txt` | 6 | which replicates took another replicate's parameters |
| `out/abm_system_fit_qvarma_manifest.txt` | 6 | per fit: log-likelihood, gradient, convergence, cumulative iterations, why the solver stopped |

### The score route, discarded

`old/score_loss/applications/abm_system_score_loss.c` reached a loss matrix
without the million fits of steps 6 and 7: it held the US estimate of step 5
fixed and evaluated the score of the auxiliary model's log-likelihood on each
simulated series, as the plain sum of squares `q'q` and as the
inverse-information-weighted score statistic. Neither recovered the known
answer in the Monte Carlo check and neither has a theoretical justification as
a validation distance, so both are out of the pipeline.
`docs/OLD_SCORE_LOSS.md` keeps the setup and the results, and lists the `make
old-<stem>` targets that still rebuild them. Their outputs stay where they
were:

| file | what it holds |
|---|---|
| `out/abm_system_score_loss.csv` | the unweighted loss table, one row per replicate and one column per configuration |
| `out/abm_system_score_loss_weighted.csv` | the score statistic, same shape |
| `out/abm_system_score_loss_manifest.txt` | the information matrix's spectrum and conditioning, how the statistic splits over eigen-directions, and every missing cell |
| `out/abm_system_mcs_statistic_comparison_score*.*` | the confidence sets under both MCS statistics |

### A side run with the per-firm output

A separate side run, not part of the route to the confidence set:

    make app-abm_system_micro_simulate                5 configurations x 5 seeds
                                                      with the per-firm output on

It reruns configurations of the design at seeds the archives already hold and
keeps each firm's output, which `studies/abm_system_micro_clt.c` and
`studies/abm_system_stable_tails.c` read. About half a minute and 261 MB.

### Does the procedure recover an answer it already knows

The procedure cannot be checked on the US data, because nobody knows which
configuration is really closest to the United States. `montecarlo/` asks the
same question with the answer known: every simulated run of every
configuration is promoted in turn to the role the US data plays, every loss is
measured against it, and the confidence set should come back holding the
configuration that run was drawn from.

It is a pipeline of its own under `montecarlo/`, sharing no source with
`applications/`, and it fits nothing: every run's parameters come from the fit
cache of step 6. When run n is the benchmark, run n of every configuration is
held out, because every configuration shares the same seeds. It runs under the
t-QVARMA and under local projections (et_al's `lp/lp.h`, which follows the R
package `lpirfs`): linear, each of two states of a state-dependent model, and
both states' responses stacked into one vector. The local-projection fits come
from their own step, `make app-abm_system_fit_lp`, all 1,000,000 runs in a few
minutes.

The grid under the variance-corrected loss is not written yet.
`docs/MONTECARLO_VALIDATION.md` describes the design it shares with the
discarded naive-loss grid, whose results and follow-ups are recorded in
`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md` and the documents it names. Three
earlier Monte Carlo experiments used subsets of the grid's benchmarks;
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` and
`docs/OLD_MONTECARLO_LP_VALIDATION.md` record them.

`applications/lp_system.h` builds the series the fits use, the stored series
the t-QVARMA is fitted on plus a 4-quarter moving average of log GDP to set the
state, and defines the archive layout the fits are stored in, one compressed
`out/abm_system_fit_lp/cop_NNNN.npz` per configuration. `applications/us_lp_fit.c`
fits the same local projections to the US data, `out/us_lp_fit.npz`.
`docs/MONTECARLO_VALIDATION.md`, "The local projections", describes the fits.
The same local projections fitted on levels rebuilt by cumulating the stored
series are discarded, `docs/OLD_LP_R_LEVELS.md`.

### Not part of the main pipeline

Scripts in `applications/` that the six steps do not use:

- `abm_system_convert_rdata.c` is the older way to fill a dataset directory: it
  converts the 108 `.Rdata` files under `dataset/simulated/` instead of running
  the simulator. Use one or the other, never both: step 6 fits every
  subdirectory of the dataset directory regardless of what wrote it, so two
  datasets sitting there at once would be fitted together with nothing in the
  results to say so. It writes `dataset/abm_system_rdata/`, which is where its
  output was moved when the design runs took `dataset/abm_system/` over, and it
  refuses to start if the directory it is about to write into already holds
  `cop_*` directories. `ABM_SYSTEM_RDATA_DIR` overrides where it writes.
- `throughput_dataset.c`, `throughput_fit.c` and
  `throughput_iteration_budget.c` are a throughput study of the fitting step on
  a separate 500 by 1000 dataset, `docs/BENCHMARKS.md`.

Files in `out/` the six steps do not write:

- `out/abm_system_score_loss*.csv`, `out/abm_system_score_loss_manifest.txt`
  and `out/abm_system_mcs_statistic_comparison_score*.*`, from the discarded
  score route above.
- `out/abm_system_convert_rdata_manifest.txt`, from the older route above.
- `out/throughput_*`, `out/fit_speedup_options.txt` and the shell scripts
  `out/hold_awake.sh`, `out/record_run_walltime.sh` and
  `out/run_throughput_fit_detached.sh`, from the throughput study
  (`docs/BENCHMARKS.md`). The scripts are the tooling that study was run with.
- `out/dsk_dual_layout_experiment.txt`, the run behind one of the rejected
  changes in `docs/DSK_MODEL_CHANGES.md`. Nothing writes it; it is kept because
  the document cites it.
- `out/dsk_*` and the other `*_report.txt` files of the tests below.
- `out/qvarma_*` and `out/fit_contention_source.txt`, from the benchmarks
  (`docs/BENCHMARKS.md`) and the studies (`docs/QVARMA_FIT_STUDIES.md`), and
  `out/abm_system_*` files not listed above, from the studies.
- `out/refactor_equality_check/` and
  `old/naive_mae_loss/out/winner_normality_rerun/`: checks that moving code to
  `old/` changed no result. The first reruns the four statistic-comparison
  programs and diffs their outputs against git: the CSVs are identical and the
  text reports differ only in the time column. The second reruns
  `old/naive_mae_loss/studies/abm_system_winner_normality.c`; its report is
  identical apart from a title line added since.

## Tests

    make test         every test script
    make test-stress  the same, including the slow simulation checks

The auxiliary model is et_al's and is tested there, not here: `make test` in an
et_al clone runs its correctness suite, including the cache and resume
behaviour `applications/abm_system_fit_qvarma.c` relies on.

What `make test` runs is four tests of this project's own code and data and
twelve of the simulator. The four:

    make test-abm_system_layout             applications/abm_system.h stores and
                                            returns what it says it does: the
                                            transformation against its closed
                                            form, agreement with the route the
                                            US data takes, the burn-in, and an
                                            archive round trip
    make test-abm_system_dataset_finiteness every value in dataset/abm_system is a
                                            finite number, every block has the
                                            shape the layout promises, and every
                                            configuration holds each replication
                                            once. Needs the experiment's archives
    make test-lp_system_transform           applications/lp_system.h builds the
                                            series the local projections are
                                            fitted on and the state series
                                            exactly
    make test-marginal_laws_correctness     the densities and fits of
                                            studies/marginal_laws.h,
                                            docs/ABM_SYSTEM_DISTRIBUTION.md

The simulator's twelve are the gate on any change to it:

    make test-dsk_build_equivalence         every byte of the economy-wide output
                                            matches upstream, over three seeds
    make test-dsk_full_output_equivalence   the same for the twelve per-firm
                                            files the model writes under -f 1
    make test-dsk_design_equivalence        the same away from the default
                                            parameters, at four points of the design
    make test-dsk_dataset_reproduction      the same where the experiment actually
                                            ran, against the stored archives
    make test-dsk_tail_replicate_reproduction
                                            the same on the replicates that carry
                                            the tails, which is where a change
                                            acting only in extreme states would
                                            show. Needs the experiment's archives
    make test-dsk_memory_safety             every run under AddressSanitizer and
                                            UndefinedBehaviorSanitizer is clean
    make test-dsk_ulp_sensitivity           how small a difference in the
                                            arithmetic those comparisons can see
    make test-dsk_long_path                 the filename bug stays fixed
    make test-dsk_redenomination_invariance counting money in a bigger unit leaves
                                            the economy where it was
    make test-dsk_machine_lot_rebase        the same for counting machines in lots
    make test-dsk_good_unit_invariance      the same for counting the consumption
                                            good in a bigger unit
    make test-dsk_bulk_cancellation_distribution
                                            the bulk cancellation draw is the
                                            distribution it claims to be

The last four come from the long-horizon work; `docs/DSK_LONG_HORIZON.md`
explains what each mechanism is for.

    make study        the robustness checks on the pipeline's own result
    make study-robustness   the same target, under the name that says what it is
    make study-<stem>       one of them

Every study reads what the pipeline wrote to `out/` and `dataset/` and none
rebuilds it, so none runs on a fresh clone. `make study` checks for the dataset
and the fit cache first and says which is missing.
`docs/QVARMA_FIT_STUDIES.md` lists every study and where it is written up.

    make bench        the benchmarks, docs/BENCHMARKS.md

`test-dsk_build_equivalence` builds the reference itself and then runs it seven
times for three seeds - once to check the scratch path works, twice per seed to
check it repeats itself - so it takes about four and a half minutes. The
unmodified binary carries no optimisation and one of its runs is 35 seconds.

## References

- S. Blazsek, A. Escribano, A. Licht (2023). Co-integration with score-driven
  models: an application to US real GDP growth, US inflation rate, and
  effective federal funds rate. *Macroeconomic Dynamics* 27, 203-223.
- P. R. Hansen, A. Lunde, J. M. Nason (2011). The Model Confidence Set.
  *Econometrica* 79(2), 453-497.
- P. Adaemmer (2019). lpirfs: an R package to estimate impulse response
  functions by local projections. *The R Journal* 11(2), 421-438.
- M. D. McKay, R. J. Beckman, W. J. Conover (1979). A comparison of three
  methods for selecting values of input variables in the analysis of output
  from a computer code. *Technometrics* 21(2), 239-245.
- M. Stein (1987). Large sample properties of simulations using Latin hypercube
  sampling. *Technometrics* 29(2), 143-151.
- K. V. Mardia (1970). Measures of multivariate skewness and kurtosis with
  applications. *Biometrika* 57(3), 519-530.
- K. V. Mardia (1974). Applications of some measures of multivariate skewness
  and kurtosis in testing normality and robustness studies. *Sankhya B* 36(2),
  115-128.
- J. M. Chambers, C. L. Mallows, B. W. Stuck (1976). A method for simulating
  stable random variables. *Journal of the American Statistical Association*
  71(354), 340-344.
- B. M. Hill (1975). A simple general approach to inference about the tail of a
  distribution. *The Annals of Statistics* 3(5), 1163-1174.

The DSK model itself: the repository and commit under "The agent-based model"
above.
