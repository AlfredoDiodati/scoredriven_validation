# Monte Carlo validation: does the procedure recover a known answer

## Overview

The main pipeline asks which ABM configuration sits closest to the US data and
which of them cannot be told apart from it. That question has no known answer,
so a result can be reported but not checked. This experiment asks the same
question with the answer known in advance.

A simulated run is promoted to the role the US data plays. Every loss is
measured against it, the validation procedure runs exactly as it does in the
main pipeline, and the question becomes: does the confidence set come back
holding the configuration the benchmark was drawn from? It should. The
benchmark is one of that configuration's own runs, so the configuration that
generated it is, by construction, the right answer. Every run of every
configuration takes that role in turn, 1,000,000 benchmarks per auxiliary
model, so the answer is known for every configuration and not only for one.

This is a check on the procedure, not on the ABM. A confidence set that misses
the configuration its own benchmark came from would say the procedure cannot
identify a configuration even when one is definitely there, which would make
the main pipeline's result hard to interpret. A confidence set that finds it
says the machinery works on a case where the answer is known, and says nothing
about whether the DSK model resembles the United States.

`docs/MCS_VALIDATION_PROCEDURE.md` describes the procedure. Everything it says
about settings, losses and their limitations applies here unchanged; only the
benchmark differs.

The first run of this experiment used the naive loss, the mean absolute error
between response vectors with no variance correction: the grid over every run
of every configuration and the experiments that read it, the compressed
responses, the fit metric and `cop_0191` as the truth. They are discarded, with
their code and outputs under `old/naive_mae_loss/`;
`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md` records the grid and its results.
Before it, smaller sets of benchmarks were used: one run of `cop_0191`, every
run of `cop_0191`, and run 0 of every configuration,
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md`. `docs/OLD_EXPERIMENTS.md` lists
every discarded experiment with its reason.

## The experiment: every run of every configuration as the benchmark

Every run of every configuration is the benchmark in turn, 1,000,000
benchmarks per model, five million confidence sets in all, under the five
models: the t-QVARMA (p1q1r2, total responses to horizon 20, 525 entries per
run) and the four local projections of "The local projections" below, linear,
state 1 and state 2 (400 entries each) and both states stacked (800 entries).
It answers, for each configuration, how often the procedure returns it when it
is the right answer.

- `make mc-true_configuration_recovery` builds at -O3 and runs.
- `./bin/true_configuration_recovery` does the same, once built.

Each benchmark is the procedure of `docs/MCS_VALIDATION_PROCEDURE.md` against
that run: the loss is the squared error between response vectors averaged over
the entries, minus each configuration's response variance averaged over the
entries; the benchmark's run index is held out of every column and of every
variance ("What is held out, and why" below); MCS_TR at alpha 0.05, block
length 1, bootstrap variance, 2000 resamples, seed 123 stream 0.
`montecarlo/recovery_loss.h` computes the loss, and
`tests/recovery_loss_correctness.c` checks it against a reference written from
the definitions and checks, on simulated responses, that it removes the bias
the uncorrected loss has.

Benchmarks are taken in order of run index: run 0 of every configuration, then
run 1, and so on, so a stopped run holds every configuration at the same
number of benchmarks. For each run index the configurations' variances are
computed once. Benchmarks sharing the run index then go in batches whose size
is the smallest multiple of the hardware thread count that is at least 64:
their loss matrices are formed in one pass over the response cache, then their
confidence sets run side by side, one per thread. Every hardware thread is
used, counted at run time.

Settings read from the environment:

- `TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK`: how many run indices one model does
  before the next model takes the same ones, default all 1000 (one model at a
  time, t-QVARMA first).
- `ABM_SYSTEM_LP_FIT_DIR`, `ABM_SYSTEM_INPUT_DIR`: where the local-projection
  fits and the dataset are read from, default `out/abm_system_fit_lp` and
  `dataset/abm_system`.

Compile-time overrides for a test build over part of the grid: `BOOTSTRAP`,
`FIRST_MODEL`, `LAST_MODEL`, `RUN_LIMIT`, `BENCHMARK_LIMIT`, `PROGRESS_PATH`,
`RESULT_PATH`, `RESPONSE_CACHE_DIR`. A test build leaves its rows in its
progress file and writes no result.

Resumable at benchmark granularity. Rows go to
`montecarlo/out/true_configuration_recovery_progress.csv` as each batch
finishes; a rerun skips what that file holds. When every model is done the
rows are compressed to `montecarlo/out/true_configuration_recovery.csv.gz` and
the progress file is removed. If the result exists and no progress file does,
the finished rows are unpacked back into the progress file and only the models
missing from them are run.

One row per model and benchmark: model, benchmark configuration and run,
whether the configuration is in the set, its rank by mean loss, its MCS
p-value, the set size, whether the set was decided by an accepted test, the
final p-value, runs dropped, the configuration with the smallest mean loss, the
configurations in the set. The columns are those of the naive-loss grid, so the
two can be compared row by row.

The full grid was started on 2026-10-12 as the user service
`true-configuration-recovery` (`systemd-run --user`, memory capped at 6 GB, all
16 hardware threads), log in `out/true_configuration_recovery.log`. No results
of it yet.

## The local projections

Every one of the 1,000,000 simulated runs (1000 configurations, 1000 runs each,
`dataset/abm_system/`) goes through two steps.

The series are the five the t-QVARMA is fitted on, as stored: GDP growth,
energy growth, employment change, inflation and the interest rate, in that
order (`lp_system_from_stored` in `applications/lp_system.h`). The series that
decides the state is the 4-period moving average of a log GDP level rebuilt
from the growth rates, `100 log(100 + cumsum(GDP_growth))`, and the first 3 of
the 400 periods are dropped, where that moving average does not exist: 397
periods. `tests/lp_system_transform.c` checks that the fitted series are the
stored ones exactly and that the state series equals, exactly, an independent
plain-loop implementation of the same formula.

Then et_al's `lp_lin_and_nl`, which follows the R package `lpirfs`, fits both
models: 4 lags, horizons 0 to 15, a unit shock with recursive ordering in the
order above; the state-dependent model with the Hodrick-Prescott cycle of the
moving average at lambda 1600, a logistic weight of slope 2, lagged one period.
Each model's impulse responses are a 5 by 16 by 5 array, 400 numbers: the
linear model's, and the state-dependent model's in state 1 (periods weighted
towards GDP above trend) and state 2 (below).

The state-dependent model is also compared as one vector: state 1's 400
responses followed by state 2's, 800 numbers, called "both states" below
(`lp_nl` in the code and in the file names, `nl` in the `lp_*` files). Its
loss, the mean absolute difference over the 800 entries, is exactly the average
of the two states' losses, since the two halves have the same length; measured
on the single-benchmark tables of `docs/OLD_MONTECARLO_LP_VALIDATION.md`, the
largest relative difference between the two is 4.8e-15.

`applications/abm_system_fit_lp.c` writes the fits to one compressed archive
per configuration, `out/abm_system_fit_lp/cop_NNNN.npz`, one row per run, with
the responses, the shock matrix, the fit diagnostics, the settings and a
checksum of the data each row was fitted on; 8.6 GB in all. No run failed:
every fit has responses, every regression at every horizon and the VAR behind
the shocks were solved by ordinary least squares at full rank (none used the
minimum-norm solution et_al falls back on), and no state series was constant
(`out/abm_system_fit_lp_manifest.txt`). Fitting all 1,000,000 runs from scratch
takes 402 seconds on this machine (AMD Ryzen 7 4800H, 16 threads); the rerun
that replaced the R-level fits took 1680 seconds, since it first read and
refused every stored archive.

## A separate pipeline that refits nothing

This is its own pipeline, in its own directory, writing to its own `out/`.
Nothing under `applications/` is modified or read as source: the main pipeline
stays exactly as it is, and a change here cannot affect a real-data result.

What is reused is that pipeline's **output**: `out/abm_system_fit_qvarma/`, the
million cached t-QVARMA fits, and `out/abm_system_fit_lp/`, the million
local-projection fits. No model is estimated anywhere in this experiment.

`montecarlo/response_cache.h` holds every run's response vector under one
model as float32, in `out/sweep_grid_response_cache/<model>.f32` (2.1 GB for
the t-QVARMA, 1.6 GB per local projection, 3.2 GB for both states, ignored by
git). The files were built from the fits by the naive-loss grid between
2026-09-29 and 2026-10-05, after the last local-projection refit (2026-09-27)
and the last t-QVARMA fit (2026-09-10). Both states' file holds, for every run,
its 400 state-1 entries followed by its 400 state-2 entries, bit for bit those
of the two state files. A file is mapped into memory rather than read into it,
so when memory is short the system drops parts of it and reads them back from
disk. A cache file is rebuilt only when its header (model, configurations,
runs, entries per run) disagrees with what the program expects; a change of
the fits or the dataset with the same shape is not detected, so delete the
file by hand after either changes.

## What is held out, and why

When run n of a configuration is the benchmark, run n of **every**
configuration is dropped from that benchmark's loss matrix, not just the
benchmark's own.

The benchmark's own cell is the obvious part: a run scored against itself
sits at a loss of essentially zero and would win by construction. The rest is
about the seeds. `applications/abm_system_simulate.c` runs every configuration
on the same seeds, `seed = mc + 1`, so run n of every configuration is seed
n + 1 of the model. Leaving those in would let every configuration be judged
partly on a draw the benchmark also used.

The matrices therefore have 999 rows rather than 1000, and every configuration
loses the same one, so the comparison between configurations is unchanged.
The same 999 runs are the ones each configuration's response variance is
estimated on (`docs/MCS_VALIDATION_PROCEDURE.md`, "The steps").

In practice the second part of this is a precaution rather than a correction.
The project already measured what the shared seeds buy:
`docs/ABM_SYSTEM_SIMULATION.md`, "Common random numbers", reports a mean-loss
correlation across configuration pairs of 0.0017 at the same seed against
-0.0001 at a shifted one, because the model consumes a different number of
draws per period at different parameters and the streams desynchronise. So run
n of one configuration has almost nothing in common with run n of another
beyond its label. The hold-out costs one observation in a thousand and removes
the question entirely.
