# Compressing the auxiliary model's output before comparing it

## The answer

The 1000 by 1000 Monte Carlo (`montecarlo/sweep_grid.c`) returns the
benchmark's own configuration in about one confidence set in ten, under the
t-QVARMA and under the linear local projection alike. The reason is not the
auxiliary models themselves but what is compared: the loss averages the
absolute difference over every entry of the impulse-response vector, 525 for
the t-QVARMA and 400 for the local projection, and almost all of those entries
are run-to-run noise. The information that separates configurations sits in a
handful of linear combinations of the entries.

Replacing each run's response vector with its coordinates along those few
combinations, learned from simulated runs the benchmark is not among, and
changing nothing else, raises the share of confidence sets that hold the right
configuration from 9.2% to 79.7% for the t-QVARMA and to 94.0% for the two
models stacked, at the full experiment's scale. The sets stay small, 1.03 to
1.10 configurations on average, and when the procedure is wrong the
configuration it picks is usually a near neighbour in parameter space rather
than an unrelated one.

The same holds for the nonlinear local projections, which start lower (1% to
2% at full scale) and reach 85% to 86% per state, 93% with the two states
stacked; "The nonlinear local projections" below.

The loss (mean absolute difference between two vectors), the confidence set
and every one of its settings are those of `montecarlo/sweep_grid.c`. Only the
vector each run is mapped to changes.

## Why the raw responses fail

All measurements in this section are on the response caches
`montecarlo/sweep_grid.c` keeps, `out/sweep_grid_response_cache/qvarma.f32`
and `lp_lin.f32`, and on its result rows for those two models, one row per
benchmark, 1,000,000 benchmarks per model. They were read from
`montecarlo/out/sweep_grid_qvarma.csv.gz` and `sweep_grid_lp_lin.csv.gz`, which
an earlier version of `sweep_grid.c` wrote and which are byte for byte the
t-QVARMA and linear rows of `montecarlo/out/sweep_grid.csv.gz`; the setup and
the full results of that experiment are in `docs/MONTECARLO_VALIDATION.md`,
"Every run of every configuration as the benchmark".

**Most entries are noise.** For each entry of the response vector, the
variance of the 1000 configurations' average responses divided by the variance
of one configuration's runs around its own average (pooled over
configurations, runs 500 to 999): the median entry has a ratio of 0.06 for the
t-QVARMA and 0.05 for the local projection. A single benchmark run is one draw
around its configuration's average, so on the median entry the benchmark's own
noise is about twenty times larger than the difference between configurations.

**The signal is concentrated in few directions.** The same ratio maximised
over linear combinations of the entries (generalised eigenvalues of the
between-configuration covariance against the pooled within-configuration
covariance) is above one in only 8 directions for the t-QVARMA, 11 for the
local projection and 14 for the two stacked
(`montecarlo/out/compressed_response/learn_*.txt`). The leading ratios are
large: 87, 41, 20 for the t-QVARMA, 74, 30, 24 for the local projection. A
mean absolute error over all entries gives those few directions the same
weight as the hundreds of noise directions.

**When the procedure is wrong, it is not wrong towards a neighbour.** Over
the benchmarks where another configuration has the smallest mean loss, that
configuration is a median 266th (t-QVARMA) and 202nd (local projection)
nearest to the right one out of 999, with the nine parameters scaled to the
unit interval. Within the ten nearest: 5.5% and 6.6% of the time, against 1%
for a configuration drawn at random.

**It rewards configurations whose runs scatter little.** Under the t-QVARMA,
400 configurations collect the 946,000 wrong smallest-loss slots and 50 of them
collect 68%. A configuration's count of wrong wins has Spearman correlation
-0.47 with the scatter of its own runs (sum over entries of the median absolute
deviation from its median). A configuration whose responses are tightly
concentrated near the middle of the field is close, on average, to many
benchmarks' noisy draws.

## What was changed

Each run's response vector is mapped to a short vector before the loss is
formed. Two kinds of map were tried, both learned from runs 500 to 999 of all
1000 configurations, so runs 0 to 499 are new data to them:

- **discriminant_k**: the coordinates of the response vector along the k
  directions above with the largest ratio, scaled so that within one
  configuration each coordinate has variance one. Entries that no run moves
  (10 of the t-QVARMA's 525, 15 of the local projection's 400, which the unit
  shock fixes) are dropped first, and the within-configuration covariance gets
  a ridge of 1e-6 times its average diagonal.
- **parameters**: a ridge regression of the nine design parameters,
  standardised, on the response vector, fitted on 40 runs per configuration
  (10 for the stacked models, for memory), with the penalty chosen on a held
  out tenth of the configurations. The map is the nine fitted values.

"Stacked" means the t-QVARMA and the local projection vectors of the same run
side by side, 925 entries, before the map is learned. This is the ensemble of
the auxiliary models: one output that uses both.

`montecarlo/compressed_response_learn.py` learns the maps and writes the mapped
vectors of all 1,000,000 runs to `out/compressed_response_cache/`, in the
layout of the response cache. `montecarlo/compressed_response_recovery.c` runs
the confidence set on any such file. `montecarlo/compressed_response_report.py`
writes `montecarlo/out/compressed_response/report.txt`.

## How it was measured

**Main setup.** For every one of the 1000 configurations, runs 0 to 4 each
stand in for the data in turn: 5000 benchmarks. The rows of the loss matrix
are runs 0 to 199 of every configuration minus the benchmark's own run, 199
rows by 1000 configurations, none of them used to learn any map. Loss: mean
absolute difference between the row's vector and the benchmark's. Confidence
set: MCS_TR, alpha 0.05, 2000 bootstrap resamples, block length 1, bootstrap
variance, seed 123 stream 0, as in `montecarlo/sweep_grid.c`. Files:
`montecarlo/out/compressed_response/recovery_<vector>_rows200.csv`, and
`rank_<vector>_rows200.csv` where only the rank by mean loss was computed.

**Full-scale setup.** Benchmarks runs 0 and 1 of every configuration, 2000
benchmarks; rows all 1000 runs minus the benchmark's, as in the 1000 by 1000
experiment. Here the rows include runs 500 to 999 the maps were learned on; the
benchmark itself never is. Files `recovery_<vector>_rows1000.csv`; the raw
baselines `grid_<model>_raw_rows1000.csv` are the rows of
`montecarlo/out/sweep_grid.csv.gz` whose benchmark is run 0 or 1, copied by
`montecarlo/compressed_response_grid_baseline.py`, not a rerun.

Measured: whether the benchmark's own configuration is in the confidence set
(in set), whether it has the smallest mean loss (first), the mean set size.
Nothing was excluded: no run has a missing response under either model.

## Results

Full scale, 2000 benchmarks:

| output vector | in set | first | mean set size |
| --- | --- | --- | --- |
| t-QVARMA, raw (525 entries) | 9.2% | 4.8% | 2.97 |
| local projection, raw (400 entries) | 10.4% | 6.1% | 3.10 |
| t-QVARMA, discriminant_8 | 79.7% | 76.9% | 1.10 |
| local projection, parameters | 93.4% | 91.8% | 1.04 |
| stacked, discriminant_12 | 94.0% | 92.5% | 1.03 |
| stacked, parameters | 94.4% | 93.0% | 1.05 |

Main setup, 5000 benchmarks, 199 rows; fewer rows make the sets larger, so the
in-set shares are higher across the board:

| output vector | in set | first | mean set size | configurations never in set (of 1000) |
| --- | --- | --- | --- | --- |
| t-QVARMA, raw | 24.7% | 5.1% | 18.95 | 466 |
| local projection, raw | 20.4% | 6.1% | 12.47 | 493 |
| t-QVARMA, discriminant_8 | 83.6% | 76.8% | 1.28 | 12 |
| local projection, discriminant_20 | 90.5% | 84.7% | 1.24 | 0 |
| local projection, parameters | 94.8% | 91.1% | 1.13 | 0 |
| stacked, discriminant_12 | 95.2% | 93.2% | 1.07 | 1 |

How many directions, by the share first, same setup (rank only):

| directions | 3 | 5 | 8 | 12 | 20 | 40 |
| --- | --- | --- | --- | --- | --- | --- |
| t-QVARMA | 55.1% | 71.3% | 76.8% | 75.4% | 69.3% | 54.3% |
| local projection | 41.7% | 62.7% | 75.0% | 81.3% | 84.7% | 76.2% |
| stacked | 54.9% | 77.7% | 91.4% | 93.2% | 93.8% | not run |

Recovery rises while the added directions carry more signal than noise and
falls after, which is what the ratios predict: the t-QVARMA peaks near its 8
directions with ratio above one. The rule "keep the directions with ratio above
one" needs only the training runs and gives 8, 11 and 14; the counts in the
tables above were picked after seeing these rank results, for the local
projection and the stack, so treat the differences between neighbouring counts
as within selection noise.

**The misses become near misses.** Among benchmarks where another
configuration is first, that configuration is one of the ten nearest to the
right one in parameter space 5% of the time for the raw t-QVARMA, 20% for
discriminant_8, 28% for the stacked discriminant_12 and 39% for the local
projection's parameters, and more often still when distance is measured over
the seven parameters the responses identify (below). Full table in
`montecarlo/out/compressed_response/report.txt`.

**The local projection carries more information than the t-QVARMA.**
Held-out R2 of the parameter regression, per parameter, in the order Gamma,
chi, psi1, psi3, alfa, taylor1, taylor2, taylor, kappa:

- t-QVARMA: 0.00, 0.21, 0.17, 0.93, 0.55, 0.58, 0.84, 0.89, 0.35
- local projection: 0.35, 0.32, 0.76, 1.00, 0.70, 0.43, 0.93, 0.95, 0.88
- stacked: 0.43, 0.37, 0.74, 1.00, 0.78, 0.64, 0.95, 0.98, 0.89

Gamma and chi are poorly identified by any of the three: two configurations
differing mainly in those would not be told apart by any vector tried here.
Recovery is high anyway because neighbouring configurations of the Latin
hypercube also differ in the identified parameters.

## The nonlinear local projections

The state-dependent local projection's responses in state 1 and state 2
(`docs/MONTECARLO_LP_VALIDATION.md`, 400 entries each, cached as
`out/sweep_grid_response_cache/lp_s1.f32` and `lp_s2.f32`) went through the
same learning step and the same two setups. Stacks: the two states together
(800 entries), and all four models together (1725 entries). The regression used
20 training runs per configuration for the two states and 10 for all four, for
memory.

On the raw responses they are the worst of the four. Over all 1,000,000
benchmarks of `montecarlo/out/sweep_grid.csv.gz` the true configuration is in
the set 1.5% of the time in state 1 and 1.6% in state 2, against 9.5% for the
t-QVARMA and 10.7% for the linear projection. The cause is the same and
stronger: 7 directions with ratio above one in each state, against 11 for the
linear projection, with leading ratios of 23 and 22 against 74.

Full scale, 2000 benchmarks, setup as above:

| output vector | in set | first | mean set size |
| --- | --- | --- | --- |
| state 1, raw | 1.0% | 0.4% | 5.34 |
| state 2, raw | 1.8% | 0.7% | 6.57 |
| state 1, parameters | 86.2% | 82.6% | 1.11 |
| state 2, parameters | 85.3% | 81.4% | 1.12 |
| both states stacked, parameters | 93.4% | 91.6% | 1.05 |
| all four stacked, discriminant_12 | 94.9% | 94.0% | 1.02 |

Main setup, 5000 benchmarks, 199 rows:

| output vector | in set | first | mean set size |
| --- | --- | --- | --- |
| state 1, raw | 4.9% | 0.7% | 34.16 |
| state 2, raw | 5.1% | 0.7% | 45.47 |
| state 1, parameters | 89.7% | 82.1% | 1.29 |
| state 2, parameters | 88.9% | 81.4% | 1.31 |
| both states stacked, parameters | 94.7% | 90.7% | 1.14 |
| all four stacked, discriminant_12 | 96.2% | 94.2% | 1.06 |

Each state alone is weaker than the linear projection alone, and the two
states together recover what the linear projection recovers by itself (93.4%
in both cases at full scale). Adding both states to the t-QVARMA and linear
stack moves it from 94.0% to 94.9%: the nonlinear projections carry little that
the linear one does not. The rank-only results for every number of directions
are in `montecarlo/out/compressed_response/rank_lp_s*_rows200.csv` and
`report.txt`.

## What this does not show

**It says nothing about absolute fit to the US data, and may hide its
absence.** The directions are those along which the simulated configurations
differ from each other. A feature of the US data that every configuration gets
wrong in the same way lies outside them and contributes nothing to the loss.
The README's finding that 965 of 1000 configurations are further from the US
responses than an all-zero response is a statement about the raw vectors and
should keep being measured on them; the compressed distance ranks
configurations relative to each other only.

**The configurations are the same ones the maps were learned on.** A benchmark
run is new data, which is the situation of the US data. A configuration
outside the 1000 was not tested.

**The full-scale rows include training runs.** The benchmarks never are, and
the main setup, which has no overlap at all, gives the same ordering and
similar margins.

**The learning step is a prototype in Python.** It reproduces byte for byte
the files the numbers above were computed from. A version in C would need two
general pieces that belong in et_al rather than here: the directions that
maximise between-group against within-group variance (a generalised symmetric
eigenproblem; et_al has `mat_eig_sym` and a Cholesky factorisation to build it
from, but no routine for it), and ridge regression.

## How to rerun

Reads `out/sweep_grid_response_cache/`, `montecarlo/out/sweep_grid.csv.gz` and
`dataset/abm_system_design.csv`. The learning and report scripts need numpy,
scipy and pandas; the baseline copy needs only the Python standard library.
The response cache is 2.1 GB for the t-QVARMA; the learning step streams it
one configuration at a time and peaks at about 650 MB.

Learn the maps. Each call writes every `discriminant_<k>` for the k in
DIRECTIONS (default 3, 5, 8, 12, 20, 40) and the `parameters` map to
`out/compressed_response_cache/<label>_<map>.f32`, `<label>` being the model
names joined by underscores; a single model also gets `<model>_raw.f32`, the
raw responses of runs 0 to 199. What was learned goes to
`montecarlo/out/compressed_response/learn_<label>.txt`. DIRECTIONS changes
which files are written, not their contents, since each map keeps the leading
k directions of the same decomposition. RIDGE_RUNS changes the `parameters` map.

    python montecarlo/compressed_response_learn.py qvarma
    python montecarlo/compressed_response_learn.py lp_lin
    RIDGE_RUNS=10 python montecarlo/compressed_response_learn.py qvarma+lp_lin
    python montecarlo/compressed_response_learn.py lp_s1
    python montecarlo/compressed_response_learn.py lp_s2
    RIDGE_RUNS=20 python montecarlo/compressed_response_learn.py lp_s1+lp_s2
    RIDGE_RUNS=10 python montecarlo/compressed_response_learn.py qvarma+lp_lin+lp_s1+lp_s2

The rank tables need every k. Afterwards only the vectors the confidence-set
tables use were kept on disk (36 to 420 MB each); to rebuild just those, pass
`DIRECTIONS=8` for `qvarma`, `DIRECTIONS=20` for `lp_lin` and `DIRECTIONS=12`
for the two stacks that use discriminant_12.

Score them. The arguments are the vectors, the first and last benchmark run,
the first and one past the last row run, the output, and `rank_only` to skip
the confidence set. The main setup is `0 4 0 200`, the full-scale setup
`0 1 0 1000`.

    make bin/compressed_response_recovery
    run() { ./bin/compressed_response_recovery out/compressed_response_cache/$1.f32 "${@:2}"; }
    dir=montecarlo/out/compressed_response

    # main setup, confidence sets
    for v in qvarma_raw qvarma_discriminant_8 lp_lin_raw lp_lin_discriminant_20 lp_lin_parameters \
             qvarma_lp_lin_discriminant_12 lp_s1_raw lp_s1_parameters lp_s2_raw lp_s2_parameters \
             lp_s1_lp_s2_parameters qvarma_lp_lin_lp_s1_lp_s2_discriminant_12; do
        run $v 0 4 0 200 $dir/recovery_${v}_rows200.csv
    done

    # main setup, rank only, for every other vector; a recovery file already
    # holds the rank
    for f in out/compressed_response_cache/*.f32; do
        v=$(basename $f .f32)
        [ -f $dir/recovery_${v}_rows200.csv ] || run $v 0 4 0 200 $dir/rank_${v}_rows200.csv rank_only
    done

    # full scale, confidence sets
    for v in qvarma_discriminant_8 lp_lin_parameters qvarma_lp_lin_discriminant_12 qvarma_lp_lin_parameters \
             lp_s1_parameters lp_s2_parameters lp_s1_lp_s2_parameters qvarma_lp_lin_lp_s1_lp_s2_discriminant_12; do
        run $v 0 1 0 1000 $dir/recovery_${v}_rows1000.csv
    done

    # full scale, raw baselines, copied from the 1000 by 1000 experiment
    python montecarlo/compressed_response_grid_baseline.py

    python montecarlo/compressed_response_report.py

These lists are the files in `montecarlo/out/compressed_response/`. Which
command produced each was not logged at the time, so they are reconstructed from
the file names and the program's usage. Two differences from the files on disk:
`qvarma_lp_lin_discriminant_40` has no rank file (never scored, "not run" in the
directions table), and six vectors (`lp_s1_raw`, `lp_s2_raw`, the three
`parameters` of the states and the four-model `discriminant_12`) have both a
rank file and a recovery file. A raw vector at full scale is not scored here
because `<model>_raw.f32` holds only runs 0 to 199; its baseline is the copy.

One confidence-set evaluation of 5000 benchmarks takes 9 to 11 minutes on 6
threads of this machine, with the 1000 by 1000 experiment using the rest; with
`rank_only` it takes seconds. The raw t-QVARMA baseline over 199 rows took 24
minutes. With the machine otherwise idle and 14 threads, the nonlinear
evaluations took about 5 minutes each over 199 rows and 3 at full scale.
