# Monte Carlo validation under the naive loss, discarded

## Status: discarded

This is the Monte Carlo grid under the naive loss, the mean absolute error
between impulse-response vectors with no variance correction, which biases
the confidence set towards configurations whose responses are less noisy
(`docs/MCS_VALIDATION_PROCEDURE.md`, "Why the variance correction"). The
design sections that still hold, the local projections and what is held out,
are kept in `docs/MONTECARLO_VALIDATION.md`. Its code and outputs are under
`old/naive_mae_loss/`, and `docs/OLD_EXPERIMENTS.md` lists every discarded
experiment. The text below is the write-up as it stood when it was discarded.

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

Earlier experiments used smaller sets of benchmarks: one run of `cop_0191`,
every run of `cop_0191`, and run 0 of every configuration. Each of those sets
is part of the grid below, so they are discarded and kept under `old/`;
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` records them and their results, and
`docs/OLD_EXPERIMENTS.md` lists every discarded experiment with its reason.

## The experiments

The grid is the experiment. In every one of its benchmarks the benchmark's run
index is held out of every configuration's column, the loss is measured
against the benchmark, and the question is whether the Model Confidence Set
(range statistic `MCS_TR`, which compares every pair of configurations; alpha
0.05, block length 1, bootstrap variance, seed 123 stream 0;
`docs/MCS_VALIDATION_PROCEDURE.md`, "Confidence set settings") holds the
benchmark's own configuration. The other two rows reuse its
responses or its confidence sets for a different question.

| experiment | benchmarks | losses and auxiliary models | resamples | program | output | written up in |
| --- | --- | --- | --- | --- | --- | --- |
| every run of every configuration | 1,000,000 per model | impulse response, five models | 2000 | `bin/sweep_grid` | `out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz` | this file, "Every run of every configuration" |
| compressed responses | 5000 or 2000 per vector | the same responses mapped to a few coordinates first | 2000 | `old/naive_mae_loss/montecarlo/compressed_response_*` | `out/old/naive_mae_loss/montecarlo/compressed_response/` | `docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md` |
| fit metric | 1,000,000 per model (oracle), 1000 per model (bootstrap) | the best configuration's mean loss against the responses' sampling standard deviation, five models | none | `old/naive_mae_loss/montecarlo/fit_metric_run.sh` | `out/old/naive_mae_loss/montecarlo/fit_metric_*` | `docs/OLD_MONTECARLO_FIT_METRIC.md` |

The five models are the t-QVARMA (p1q1r2, total responses to horizon 20, 525
entries per run) and the four local projections of "The local projections"
below: linear, state 1 and state 2, 400 entries per run, and both states, the
state-dependent model's two states stacked into one vector of 800 entries, so
that one confidence set compares the configurations on both states at once.

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

`applications/abm_system_fit_lp.c` writes the fits to one compressed archive per
configuration, `out/applications/abm_system_fit_lp/cop_NNNN.npz`, one row per
run, with the responses, the shock matrix, the fit diagnostics, the settings and
a checksum of the data each row was fitted on; 8.6 GB in all. No run failed:
every fit has responses, every regression at every horizon and the VAR behind
the shocks were solved by ordinary least squares at full rank (none used the
minimum-norm solution et_al falls back on), and no state series was constant
(`out/applications/abm_system_fit_lp_manifest.txt`). Fitting all 1,000,000 runs
from scratch takes 402 seconds on this machine (AMD Ryzen 7 4800H, 16 threads);
the rerun that replaced the R-level fits took 1680 seconds, since it first read
and refused every stored archive.

## A separate pipeline that refits nothing

This is its own pipeline, in its own directory, writing to its own
`out/old/naive_mae_loss/montecarlo/`. Nothing under `applications/` is modified
or read as source: the main pipeline stays exactly as it is, and a change here
cannot affect a real-data result.

What is reused is that pipeline's **output**:
`out/applications/abm_system_fit_qvarma/`, the million cached t-QVARMA fits, and
`out/applications/abm_system_fit_lp/`, the million local-projection fits. No
model is estimated anywhere in this experiment.

| file | what it is |
| --- | --- |
| `montecarlo/response_cache.h` | every run's response vector under one model, held as float32 and kept on disk |
| `old/naive_mae_loss/montecarlo/sweep_grid.c` | every run of every configuration as the benchmark, five models |
| `old/naive_mae_loss/montecarlo/lp_nl_run.sh` | runs the grid and the experiments that read it for both states, resuming where it stopped |
| `old/naive_mae_loss/montecarlo/sweep_grid_presence.py` | how often each configuration sits in sets that are not its own |
| `old/naive_mae_loss/montecarlo/sweep_grid_cop_0191.py` | `cop_0191` as the truth, as a wrong answer, and every configuration as the truth, `docs/OLD_MONTECARLO_COP_0191.md` |
| `old/naive_mae_loss/montecarlo/sweep_grid_identifiability_plots.py` | per-configuration recovery figures from the same rows |
| `old/naive_mae_loss/montecarlo/fit_metric.h`, `fit_metric_oracle.c`, `fit_metric_bootstrap.c`, `fit_metric_report.c`, `fit_metric_run.sh`, `fit_metric_v_histograms.py` | the fit metric, oracle and bootstrap versions, `docs/OLD_MONTECARLO_FIT_METRIC.md`; the planned comparison across methods, `docs/OLD_MONTECARLO_FIT_METRIC_COMPARISON.md` |
| `old/naive_mae_loss/montecarlo/compressed_response_learn.py`, `compressed_response_recovery.c`, `compressed_response_grid_baseline.py`, `compressed_response_report.py` | `docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md` |

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

In practice the second part of this is a precaution rather than a correction.
The project already measured what the shared seeds buy:
`docs/ABM_SYSTEM_SIMULATION.md`, "Common random numbers", reports a mean-loss
correlation across configuration pairs of 0.0017 at the same seed against
-0.0001 at a shifted one, because the model consumes a different number of
draws per period at different parameters and the streams desynchronise. So run
n of one configuration has almost nothing in common with run n of another
beyond its label. The hold-out costs one observation in a thousand and removes
the question entirely.

## Every run of every configuration as the benchmark

Every run of every configuration is the benchmark in turn, 1,000,000
benchmarks per model, five million confidence sets in all. It answers, for each
configuration, how often the procedure returns it when it is the right answer.

- `make old-sweep_grid` builds at -O3 and runs.
- `./bin/sweep_grid` does the same, once built.

`sweep_grid` is not in `MONTECARLO_STEMS`: it has its own rule in the Makefile
because it is compiled at -O3 rather than -O2. The Makefile records that -O3
changes no result, and its time measured on 8 sets side by side on 999 by 1000
tables at 2000 resamples: 0.13 s per set at -O2, 0.08 s at -O3.

Each benchmark is the procedure of `docs/MCS_VALIDATION_PROCEDURE.md` against
that run, at 2000 bootstrap resamples rather than the 10000 the real-data run
uses. The loss is the mean absolute difference between response vectors, the
benchmark's run index is held out of every column, a run with a missing
response in any column is dropped, MCS_TR at alpha 0.05, block length 1,
bootstrap variance, seed 123 stream 0.

Benchmarks are taken in order of run index: run 0 of every configuration, then
run 1, and so on, so a stopped run holds every configuration at the same number
of benchmarks. Benchmarks sharing a run index are processed 64 at a time, their
loss matrices formed in one pass over the response cache, then their confidence
sets run side by side, one per thread.

Settings read from the environment:

- `SWEEP_GRID_WORKERS`: how many confidence sets run side by side, default
  every hardware thread (16 on this machine).
- `SWEEP_GRID_RUN_BLOCK`: how many run indices one model does before the next
  model takes the same ones, default all 1000 (one model at a time, t-QVARMA
  first).
- `ABM_SYSTEM_LP_FIT_DIR`, `ABM_SYSTEM_INPUT_DIR`: where the local-projection
  fits and the dataset are read from, default
  `out/applications/abm_system_fit_lp` and `out/dataset/abm_system`.

Compile-time overrides for a test build over part of the grid: `BOOTSTRAP`,
`FIRST_MODEL`, `RUN_LIMIT`, `BENCHMARK_LIMIT`, `PROGRESS_PATH`, `RESULT_PATH`,
`RESPONSE_CACHE_DIR`.

It reads `out/applications/abm_system_fit_qvarma/`,
`out/applications/abm_system_fit_lp/` and `out/dataset/abm_system/` and rebuilds
none of them. Each model's responses are kept in
`out/montecarlo/sweep_grid_response_cache/<model>.f32` (2.1 GB for the t-QVARMA,
1.6 GB per local projection, 3.2 GB for both states, ignored by git), built on
the first run in about 10 minutes for the t-QVARMA and 20 per local projection,
then read back in seconds. Both states' file is built from the two state files,
copying each run's 400 state-1 entries followed by its 400 state-2 entries, in
seconds; for every run, the file on disk holds bit for bit that run's entries in
the two state files. The file is mapped into memory rather than read into it, so
when memory is short the system can drop parts of it and read them back from
disk; that is slower but does not stop the run. A cache file is rebuilt only
when its header (model, configurations, runs, entries per run) disagrees with
what the program expects; a change of the fits or the dataset with the same
shape is not detected, so delete the file by hand after either changes. The
caches on disk were built between 2026-09-29 and 2026-10-03, after the last
local-projection refit (2026-09-27) and the last t-QVARMA fit (2026-09-10); both
states' on 2026-10-05.

Resumable at benchmark granularity. Rows go to
`out/old/naive_mae_loss/montecarlo/sweep_grid_progress.csv` as each batch of 64 finishes; a rerun
skips what that file holds. When every model is done the rows are compressed to
`out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz` and the progress file is removed. If
`sweep_grid.csv.gz` exists and no progress file does, the finished rows are
unpacked back into the progress file and only the models missing from them are
run; if none is missing, the program does nothing. A crash can leave the
progress file ending in zero bytes, the part of the last write that never
reached the disk; `old/naive_mae_loss/montecarlo/lp_nl_run.sh` cuts the file at the first of them
before resuming.

Cost: the first four models spanned several sessions, the last finished on
2026-10-04. That last session did 1,187,552 benchmarks (the rest of state 1 and
all of state 2) at 0.0956 s of wall-clock time per benchmark, 8 workers on the
machine the README describes (AMD Ryzen 7 4800H, 16 threads). The earlier
sessions' timing was not kept. Both states ran from 2026-10-05 to 2026-10-07.
A crash of the machine stopped them after 471,768 benchmarks; the session after
it did the remaining 528,232 at 0.1196 s per benchmark, 16 workers, comparing
800 entries per run rather than 400. The timing before the crash was not kept.
Both timings were read from the runs' console logs, which are not tracked by
git. After both states were added, the t-QVARMA and linear rows of the new
`sweep_grid.csv.gz` were compared with `sweep_grid_qvarma.csv.gz` and
`sweep_grid_lp_lin.csv.gz` and are identical; the state 1 and state 2 rows have
no separate copy to compare with.

### Outputs

| file | what it holds |
| --- | --- |
| `out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz` | one row per model and benchmark: model, benchmark configuration and run, whether the configuration is in the set, its rank by mean loss, its MCS p-value, the set size, whether the set was decided by an accepted test, the final p-value, runs dropped, the configuration with the smallest mean loss, the configurations in the set |
| `out/old/naive_mae_loss/montecarlo/sweep_grid_qvarma.csv.gz`, `sweep_grid_lp_lin.csv.gz` | written by an earlier version of `sweep_grid.c` that split the rows per model; byte for byte the t-QVARMA and linear rows of `sweep_grid.csv.gz`, header included. The current program does not write them |
| `out/old/naive_mae_loss/montecarlo/sweep_grid_cop_0191.txt`, `.tex` | `old/naive_mae_loss/montecarlo/sweep_grid_cop_0191.py`, in plain text and as LaTeX tables: per model, in percent, how often `cop_0191` is in the set, alone and first by mean loss over its own 1000 benchmarks and over the other 999,000, and how often the true configuration is over all 1,000,000; `docs/OLD_MONTECARLO_COP_0191.md` |
| `out/old/naive_mae_loss/montecarlo/sweep_grid_presence.txt` | `old/naive_mae_loss/montecarlo/sweep_grid_presence.py`: per model, how often each configuration is in the sets of benchmarks that are not its own |
| `out/old/naive_mae_loss/montecarlo/sweep_grid_identifiability.csv`, `sweep_grid_identifiability_plots/` | `old/naive_mae_loss/montecarlo/sweep_grid_identifiability_plots.py`: per model and configuration, the share of its own 1000 sets that hold it, and the share holding the other configuration found in them most often, with two bar charts per model. The files on disk cover the t-QVARMA and the linear local projection only; see below |

    python old/naive_mae_loss/montecarlo/sweep_grid_presence.py
    python old/naive_mae_loss/montecarlo/sweep_grid_identifiability_plots.py

The first needs only the Python standard library and covers all five models.
The second needs `polars`, `plotly` and `kaleido`; its table and figures on disk
were drawn on 2026-10-02, when only the t-QVARMA and the linear local
projection were done, and have not been redrawn since, so they cover those two
models. Rerunning the script adds the other three.

### Results

No response was missing under any model, so every loss matrix is 999 runs by
1000 configurations and no run was dropped.

| model | own configuration in the set | alone in the set | smallest mean loss | mean set size | largest set |
| --- | --- | --- | --- | --- | --- |
| t-QVARMA | 95,151 (9.5%) | 31,557 | 53,926 | 2.94 | 124 |
| LP linear | 107,025 (10.7%) | 38,964 | 63,433 | 3.09 | 75 |
| LP state 1 | 15,318 (1.5%) | 4,699 | 7,854 | 5.29 | 295 |
| LP state 2 | 16,060 (1.6%) | 3,648 | 7,305 | 6.89 | 336 |
| LP both states | 11,653 (1.2%) | 4,300 | 6,746 | 6.97 | 548 |

Counts out of 1,000,000 benchmarks per model.

Stacking the two states does not help on the raw responses: both states
recover the right configuration less often than either state alone (1.2%
against 1.5% and 1.6%), with larger sets on average (6.97 against 5.29 and
6.89). The loss over the stacked vector is the average of the two states'
losses.

Per configuration, out of its own 1000 benchmarks:

| model | configurations recovered in at least half | in at least 90% | never |
| --- | --- | --- | --- |
| t-QVARMA | 68 | 4 | 472 |
| LP linear | 35 | 1 | 188 |
| LP state 1 | 3 | 0 | 660 |
| LP state 2 | 6 | 0 | 621 |
| LP both states | 4 | 1 | 772 |

Recovery is concentrated: a few dozen configurations are found most of the
time, and about half are never found under the t-QVARMA.

When the procedure is wrong, a few configurations collect the wrong answers
(`sweep_grid_presence.txt`). The off-target rate of a configuration is the share
of the other configurations' 999,000 benchmarks whose set holds it. Under the
t-QVARMA the highest are `cop_0799` (6.85%) and `cop_0505` (6.84%), and the ten
highest take 15.7% of all wrong set slots; under the linear projection
`cop_0410` (4.6%) and `cop_0914` (4.5%), ten highest 11.0%; under state 1
`cop_0437` (25.8%) and `cop_0410` (19.1%), ten highest 24.3%; under state 2
`cop_0717` (25.5%) and `cop_0410` (23.6%), ten highest 22.6%; under both states
`cop_0717` (28.2%), `cop_0437` (26.9%) and `cop_0410` (26.9%), ten highest
23.2%. `cop_0191` is among the ten highest under every local projection. In 245
configurations under the t-QVARMA, 162 under the linear projection, 770 under
state 1, 768 under state 2 and 870 under both states, the off-target rate
exceeds the configuration's own recovery rate.
`docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md` traces this to the loss averaging over
hundreds of response entries that are mostly noise.

**2000 against 10000 resamples.** Run 0 of every configuration is a benchmark
both here (2000 resamples) and in the discarded `sweep_cops.csv.gz` (10000,
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md`), on the same loss matrices. Over
those 1000 benchmarks per model, whether the own configuration is
in the set differs once in 5000 (linear local projection, 93 against 94); the
set size differs for 51 (t-QVARMA), 27 (linear), 68 (state 1), 99 (state 2) and
94 (both states) of 1000. The lower resample count changes the set membership of
the right answer essentially never.

### Recovery by whether the benchmark's fit converged

Only the t-QVARMA has a convergence status; the local projections are least
squares. Convergence read from
`out/applications/abm_system_fit_qvarma_manifest.txt`, whose fits are unchanged
since 2026-09-10, before every experiment here.

Over the grid's t-QVARMA benchmarks, the own configuration is in the set for
35,764 of the 352,470 converged ones (10.1%) and 59,387 of the 647,530
unconverged ones (9.2%). Using unconverged benchmarks is not what keeps the
grid's recovery near one in ten.

## What this experiment cannot show

The impulse-response protocols have now been asked about every configuration
("Every run of every configuration as the benchmark"), and the answer for
`cop_0191` does not generalise. `cop_0191` is recovered from 876 of its own 1000
runs under the t-QVARMA, which only 4 configurations exceed, and from more of
its own runs than any other configuration under each of the four local
projections (913, 857, 780 and 916 of 1000 under linear, state 1, state 2 and
both states). Over all configurations the rate is about one in ten. Results
that rest on `cop_0191` alone, the discarded experiments of
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md`, describe one of the easiest
configurations to identify.

`cop_0191` was the impulse-response protocol's winner on the US data; the grid
shows its impulse responses are unusually distinctive, which explains why it is
recovered here without the procedure being good in general. The configurations
the discarded score losses favoured on the US data (`docs/OLD_SCORE_LOSS.md`)
are hard to recover under the impulse-response loss: `cop_0409` from 45 of its
1000 runs under the t-QVARMA and at most 1 under the local projections,
`cop_0931` from none under the t-QVARMA and from 27 to 238 under the local
projections (95 under both states).

On its own it cannot separate the loss from the auxiliary model: a protocol that
fails to recover a configuration may use a poor loss, or its auxiliary model may
not distinguish these configurations' dynamics at all.
`docs/OLD_MONTECARLO_COMPRESSED_RESPONSE.md` separates the two for the
impulse-response protocols. Keeping the same responses and changing only the
vector the loss compares raises recovery over the grid from about one in ten to
80% to 95%, so most of the failure is the loss.

It inherits every limitation of the protocols it runs. In particular the fits
it reads are the main pipeline's, about a third of which converged, and
`docs/MCS_VALIDATION_PROCEDURE.md` records the measurement that the
unconverged ones' log-likelihoods are not distinguishable from the converged
ones'. Most benchmarks of the grid are unconverged fits, and "Recovery by
whether the benchmark's fit converged" above measures what that changes.
