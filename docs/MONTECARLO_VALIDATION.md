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
  `out/dataset/abm_system`.

Compile-time overrides for a test build over part of the grid: `BOOTSTRAP`,
`FIRST_MODEL`, `LAST_MODEL`, `RUN_LIMIT`, `BENCHMARK_LIMIT`, `PROGRESS_PATH`,
`RESULT_PATH`, `RESPONSE_CACHE_DIR`. A test build leaves its rows in its
progress file and writes no result.

Resumable at benchmark granularity. Rows go to
`out/montecarlo/true_configuration_recovery_progress.csv` as each batch
finishes; a rerun skips what that file holds. When every model is done the
rows are compressed to `out/montecarlo/true_configuration_recovery.csv.gz` and
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
16 hardware threads), log in `out/true_configuration_recovery.log`. It stopped
when the machine rebooted at 16:32 that day, during t-QVARMA run 182, with
182,704 benchmarks in its progress file, and was not restarted. A rerun
resumes from there; at 0.125 s per benchmark it had about 168 hours left.

## The same grid under the Seri et al. procedure

`montecarlo/seri_true_configuration_recovery.c` asks the same question of the
procedure of `docs/MCS_ALTERNATIVE_CALIBRATION.md`. Everything above holds
except the loss and the confidence set:

- The loss is the squared error between response vectors averaged over the
  entries, minus each run's own squared deviation from its configuration's
  mean response averaged over the entries, scaled by $n/(n-1)$, with $n$ the
  kept runs (999). `montecarlo/own_deviation_loss.h` computes it, and
  `tests/own_deviation_loss_correctness.c` checks it against a reference
  written from the definition, checks that its column means equal those of
  `montecarlo/recovery_loss.h`, and checks on simulated responses that it
  removes the bias and shrinks the spread across runs to the predicted 2/3 of
  the common-variance loss's.
- The confidence set is et_al's `seri_mcs` at $\alpha = 0.05$, with no
  bootstrap and no seed.

- `make mc-seri_true_configuration_recovery` builds at -O3 and runs.
- `SERI_TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK` replaces
  `TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK`; the compile-time overrides are the
  same, with `ALPHA` and `MINIMUM_BATCH` in place of `BOOTSTRAP`.
- Progress goes to
  `out/montecarlo/seri_true_configuration_recovery_progress.csv`, the result
  to `out/montecarlo/seri_true_configuration_recovery.csv.gz`, with the
  columns of the grid above, so the two can be compared benchmark by
  benchmark.

The full grid was started on 2026-10-13 at 19:03 as the user service
`seri-true-configuration-recovery` (`systemd-run --user`, memory capped at
6 GB, all 16 hardware threads), log in
`out/seri_true_configuration_recovery.log`.

`make mc-seri_true_configuration_recovery_report` summarises the rows written
so far, or the finished result, into
`out/seri_true_configuration_recovery_report.txt`: per model, how often the
benchmark's configuration is in the set, alone and first by mean loss, the set
sizes, how many configurations are recovered in at least half and 90% of their
own benchmarks, and, on the t-QVARMA benchmarks the bootstrap grid above also
finished, how often each procedure's set holds the configuration. It can run
while the grid does, and reads only.

The grid finished on 2026-10-14 at 20:38, 5,000,000 benchmarks, with no error
in its log. Over all 1,000,000 benchmarks per model, alpha 0.05, from
`out/seri_true_configuration_recovery_report.txt`:

- t-QVARMA: the benchmark's configuration in the set 10.5%, alone 3.6%,
  first by mean loss 6.2%, mean set size 2.60, largest 32.
- LP linear: 14.7%, alone 5.2%, first 8.5%, mean set size 2.80, largest 39.
- LP state 1: 11.7%, alone 0.9%, first 2.5%, mean set size 12.34, largest 281.
- LP state 2: 10.1%, alone 0.8%, first 2.3%, mean set size 10.72, largest 227.
- LP both states: 24.3%, alone 1.2%, first 4.4%, mean set size 17.48,
  largest 399.
- No configuration is in its own set in 90% of its 1000 benchmarks under any
  model; the most often recovered reaches 84.5% (t-QVARMA). Under the
  t-QVARMA 122 configurations are never in their own set.
- On the 182,704 t-QVARMA benchmarks the bootstrap grid also finished: the
  configuration is in both sets 10.4%, in the Seri set only 0.03% (46), in
  the bootstrap set only 55.8%, in neither 33.8%; mean set size 2.60 against
  327.98.

`make mc-seri_against_naive_grid` compares this grid, benchmark by benchmark,
with the discarded naive-loss grid (`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md`:
mean absolute error, no correction, MCS_TR, 2000 resamples), both over the
same 1,000,000 benchmarks per model, into
`out/seri_against_naive_grid_report.txt`. Measured 2026-10-15:

- The configuration in the set, Seri against naive: t-QVARMA 10.5% and 9.5%,
  LP linear 14.7% and 10.7%, state 1 11.7% and 1.5%, state 2 10.1% and 1.6%,
  both states 24.3% and 1.2%.
- Its mean rank by mean loss, out of 1000: t-QVARMA 186 and 181, LP linear 80
  and 177, state 1 168 and 388, state 2 182 and 411, both states 112 and 407.
- Mean set size: t-QVARMA 2.60 and 2.94, LP linear 2.80 and 3.09, state 1
  12.34 and 5.29, state 2 10.72 and 6.89, both states 17.48 and 6.97.
- Median set size: t-QVARMA 2 and 2, LP linear 2 and 2, state 1 8 and 2,
  state 2 7 and 4, both states 10 and 3.
- Configurations in their own set in at least half of their benchmarks, Seri
  and naive: t-QVARMA 20 and 68, LP linear 17 and 35, state 1 7 and 3,
  state 2 10 and 6, both states 18 and 4. Never in their own set: t-QVARMA
  122 and 472, every LP model 0 under Seri against 188 to 772 under naive.
- The two losses give the same configuration the smallest mean loss in 14.7%
  of t-QVARMA benchmarks, 9.8% for LP linear and about 1% for the state
  models.

`make mc-seri_against_naive_by_configuration` computes, for both grids from
the same code, the per-configuration tables the naive grid's Python scripts
wrote (`sweep_grid_cop_0191.py`, `sweep_grid_presence.py`,
`sweep_grid_identifiability_plots.py`): `cop_0191` as the truth, as a wrong
answer and every configuration as the truth, into
`out/seri_against_naive_by_configuration_cop_0191.txt` and `.tex`; each
configuration's own and off-target rate into
`out/seri_against_naive_by_configuration_presence.txt`; and per
configuration the share of its own sets holding it and holding the most
frequent other configuration into
`out/seri_against_naive_by_configuration_identifiability.csv`, which
`montecarlo/seri_against_naive_by_configuration_plots.py` draws. Its
naive-grid presence and identifiability output equals the old scripts' files
line for line, and its naive `cop_0191` tables equal theirs to the printed
digit (checked 2026-10-15).

`make mc-response_variance_by_configuration` measures each configuration's
noise, its mean response variance over all 1000 runs averaged over the
entries (the term the variance correction subtracts), into
`out/response_variance_by_configuration.csv`, and writes
`out/response_variance_by_configuration_report.txt`: where `cop_0191` ranks by
noise, and the Spearman correlation across the 1000 configurations between
noise and the own and off-target rates of each grid. Measured 2026-10-15:

- `cop_0191` is among the least noisy: rank 47 of 1000 under the t-QVARMA,
  3, 3, 6 and 4 under LP linear, state 1, state 2 and both states, with noise
  5 to 12 times below the median configuration's.
- Noise against the off-target rate, Seri and naive: t-QVARMA -0.42 and
  -0.55, LP linear -0.27 and -0.45, state 1 -0.26 and -0.53, state 2 -0.17
  and -0.67, both states -0.80 and -0.74.
- Noise against the own rate, Seri and naive: t-QVARMA -0.56 and -0.51, LP
  linear 0.34 and -0.02, state 1 0.34 and -0.08, state 2 0.28 and -0.19,
  both states -0.05 and -0.05.

Since the column means of the two losses are equal, the rank of the
benchmark's configuration and the configuration with the smallest mean loss
are the same in both grids; only the confidence set differs.

## The same Seri grid with the loss left uncorrected

The comparison with the naive-loss grid above changes three things at once:
the loss (absolute against squared error), the correction, and the
confidence set (bootstrap MCS_TR against `seri_mcs`).
`montecarlo/seri_uncorrected_true_configuration_recovery.c` changes only the
correction. It is the grid of `montecarlo/seri_true_configuration_recovery.c`
with the loss the squared error between response vectors averaged over the
entries and nothing subtracted, $L_{is} = \frac{1}{K} \sum_h
(\widehat{IR}_{ish} - \widehat{IR}_{dh})^2$; same benchmarks, same held-out
run, same `seri_mcs` at $\alpha = 0.05$. The loss is
`montecarlo/recovery_loss.h`'s `recovery_loss_fill` with every
configuration's subtracted variance set to zero, the case
`tests/recovery_loss_correctness.c` checks as "a configuration with equal
runs has zero variance, its loss is the squared error".

- `make mc-seri_uncorrected_true_configuration_recovery` builds at -O3 and
  runs.
- `SERI_UNCORRECTED_TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK` sets the run block;
  the compile-time overrides are those of the Seri grid.
- Progress goes to
  `out/montecarlo/seri_uncorrected_true_configuration_recovery_progress.csv`,
  the result to
  `out/montecarlo/seri_uncorrected_true_configuration_recovery.csv.gz`, with
  the same columns.

Unlike the corrected grid against the common-variance one, the column means
here differ from the corrected grid's by each configuration's mean response
variance, so the rank and the smallest-loss configuration change too.

`make mc-seri_corrected_against_uncorrected_grid` compares the two grids over
the benchmarks both hold, into
`out/seri_corrected_against_uncorrected_grid_report.txt`. It reads the
uncorrected grid's progress file while the result does not exist, so it can
run while the grid does. Per model and grid: the benchmark's configuration in
the set, alone, first by mean loss, mean rank, set sizes; the two by two
table of in-set and how often the smallest-loss configuration agrees; per
configuration the own rate summaries and the correlation of the two grids'
own rates; and, against each configuration's noise from
`out/response_variance_by_configuration.csv`, the Spearman correlation of
noise with the own and off-target rates and the mean noise rank (1 the least
noisy) of the smallest-loss configuration and of the set members. The noise
ranks are the direct measure of the bias the correction is meant to remove:
with no dependence on noise they sit near the middle rank, 500.5.

### Trial over runs 0 and 1 of the t-QVARMA

Setup: a test build with `FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=2`, so runs 0
and 1 of all 1000 configurations as the benchmark under the t-QVARMA, 2000
benchmarks, every other setting as in the full grid; the comparison run
against the corrected grid's finished result on those 2000 benchmarks.
Machine: AMD Ryzen 7 4800H, 16 hardware threads, nothing else heavy running,
under a 4 GB memory cap. Rows, log and report in
`out/seri_uncorrected_true_configuration_recovery_trial/`, measured
2026-10-15.

- Time: 0.0152 s per benchmark, wall clock over the whole trial including
  reading the response cache.
- The corrected side equals the corrected grid's own trial above: in the set
  9.9%, first by mean loss 6.2%, mean set size 2.56, largest 22.
- Corrected against uncorrected: the benchmark's configuration in the set
  9.9% and 4.1%, alone 4.1% and 0.9%, first by mean loss 6.2% and 1.8%, mean
  rank 187 and 306 out of 1000, mean set size 2.56 and 2.79.
- In both sets 2.5%, corrected only 7.4%, uncorrected only 1.6%, neither
  88.5%. Same smallest-loss configuration in 11.7% of benchmarks.
- Mean noise rank of the smallest-loss configuration: 388 corrected, 164
  uncorrected; of the set members 426 and 190.
- Two runs per configuration is too few for the per-configuration rates to
  mean much; they are in the report but not summarised here.

The full grid was started on 2026-10-15 at 12:01 as the user service
`seri-uncorrected-true-configuration-recovery` (`systemd-run --user`, memory
capped at 5 GB, all 16 hardware threads), log in
`out/seri_uncorrected_true_configuration_recovery.log`, with about 26 hours
estimated at its start.

### Trial over runs 0 and 1 of the t-QVARMA

Setup: a test build with `FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=2`, so runs 0
and 1 of all 1000 configurations as the benchmark under the t-QVARMA, 2000
benchmarks, every other setting as in the full grid. Machine: AMD Ryzen 7
4800H, 16 hardware threads, nothing else heavy running, under a 6 GB memory
cap. Rows in `out/seri_true_configuration_recovery_trial/progress.csv`, log in
`out/seri_true_configuration_recovery_trial/trial.log`, measured 2026-10-13.

- Time: 0.0215 s per benchmark, wall clock over the whole trial including
  reading the response cache, against 0.125 s per benchmark for the grid
  above. The loss pass now dominates, not the confidence set.
- Rank and smallest-loss configuration: equal to the grid above's on all 2000
  benchmarks, as they must be. The benchmark's configuration has the smallest
  mean loss in 6.2% of them.
- The benchmark's configuration is in the set in 9.9% of benchmarks, against
  64.8% under the grid above on the same 2000 benchmarks.
- Set size: mean 2.56, smallest 1, largest 22, against a mean of 319 under
  the grid above. The set was decided by an accepted test in 58.9% of
  benchmarks.

## The local projections

Every one of the 1,000,000 simulated runs (1000 configurations, 1000 runs each,
`out/dataset/abm_system/`) goes through two steps.

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
