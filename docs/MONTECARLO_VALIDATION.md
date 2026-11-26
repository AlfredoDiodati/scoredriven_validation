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

The grids run three models: the linear projection, both states stacked under
the causal HP cycle at the estimated lambda ("The four grids under the other
three state weights" below) and the t-QVARMA with its two growth states stacked
("The t-QVARMA conditional on the growth state" below). The other models are
discarded: state 1 and state 2 on their own, in `old/montecarlo_single_states/`
and `out/old/montecarlo_single_states/`, both states under the other three
state weights, in `old/montecarlo_other_state_weights/` and
`out/old/montecarlo_other_state_weights/`, and the t-QVARMA at the sample mean
of D_t, in `old/montecarlo_mean_jacobian_qvarma/` and
`out/old/montecarlo_mean_jacobian_qvarma/`. `docs/OLD_EXPERIMENTS.md` says why.
The results below that quote them were measured before.

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
  fits and the dataset are read from, default
  `out/applications/abm_system_fit_lp` and `out/dataset/abm_system`.

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
16 hardware threads), log in `out/montecarlo/true_configuration_recovery.log`.
It stopped when the machine rebooted at 16:32 that day, during t-QVARMA run 182,
with 182,704 benchmarks in its progress file. A rerun resumes from there; at
0.125 s per benchmark it had about 168 hours left.

It was restarted on 2026-10-16 at 19:48 under the same unit name, memory capped
at 5 GB, all 16 hardware threads, log appended to the same file. Before the
restart, a test build of the current source and the et_al installed that day
(`FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=1 BENCHMARK_LIMIT=64`, so run 0 of
`cop_0001` to `cop_0064` under the t-QVARMA) wrote rows byte for byte equal to
the first 64 rows of the progress file, so the resumed rows and the rows from
2026-10-12 come from the same computation
(`out/montecarlo/true_configuration_recovery_resume_check/`). Its first minutes
ran at 0.093 s per benchmark, about 124 hours for what is left.

`make mc-seri_against_bootstrap_grid` compares this grid with the Seri grid
below, benchmark by benchmark, over the benchmarks both hold: the same loss
(squared error minus the variance correction, so the same mean losses, ranks and
smallest-loss configuration), different confidence sets (MCS_TR with 2000
bootstrap resamples against `seri_mcs`, both at $\alpha = 0.05$). It is
`montecarlo/seri_corrected_against_uncorrected_grid.c` built with the two grids'
files and labels, into `out/montecarlo/seri_against_bootstrap_grid_report.txt`,
and reads this grid's progress file while it runs. Measured 2026-10-19, with the
t-QVARMA and LP linear finished (1,000,000 benchmarks each) and LP state 1 at
295,384, Seri against bootstrap:

- The configuration in the set: t-QVARMA 10.5% and 66.3%, LP linear 14.7%
  and 26.0%, LP state 1 11.6% and 41.2%. Alone in the set: 3.6% and 0.0%, 5.2%
  and 3.0%, 0.9% and 0.4%.
- Mean set size: t-QVARMA 2.60 and 328.5 (largest 32 and 573), LP linear 2.80
  and 9.64 (39 and 111), LP state 1 12.3 and 119.0 (281 and 698).
- The Seri set never holds the configuration when the bootstrap set does not
  (0.0% of benchmarks under all three models).
- Configurations in their own set in at least half of their benchmarks:
  t-QVARMA 20 and 712, LP linear 17 and 31, LP state 1 7 and 129. The
  correlation across configurations of the two own rates: 0.22, 0.98, 0.65.
- Mean noise rank of the set members (1 the least noisy, middle 500.5):
  t-QVARMA 410 and 686, LP linear 464 and 466, LP state 1 457 and 534.

## The bootstrap grid with the loss left uncorrected

Designed, built, not run yet.
`montecarlo/uncorrected_true_configuration_recovery.c` is to the grid above what
the uncorrected Seri grid below is to the corrected Seri grid: the loss is the
squared error between response vectors averaged over the entries with nothing
subtracted, $L_{is} = \frac{1}{K} \sum_h (\widehat{IR}_{ish} -
\widehat{IR}_{dh})^2$, computed by `montecarlo/recovery_loss.h`'s
`recovery_loss_fill` with every configuration's subtracted variance set to zero.
Everything else is the grid above's: the same 1,000,000 benchmarks per model
under all five models, the benchmark's run index held out of every column,
MCS_TR at $\alpha = 0.05$, block length 1, bootstrap variance, 2000 resamples,
seed 123 stream 0. The source differs from
`montecarlo/true_configuration_recovery.c` only in the zero correction, the file
names and the environment variable.

What it answers: whether subtracting each configuration's mean response
variance changes which configurations the bootstrap confidence set keeps,
with the loss and the confidence set held fixed. The naive-loss grid
(`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md`) cannot answer this, since it also
uses the absolute rather than the squared error.

Since the confidence set reseeds at seed 123 stream 0 on every call, both
grids draw the same bootstrap resamples for the same benchmark, so the paired
comparison is not blurred by resampling noise. Unlike the Seri pair, the
correction here is one constant per column, and MCS_TR is not unchanged by
it: the mean loss differences between configurations move by the difference
of their variances, while the bootstrap spread of those differences does not.

- `make mc-uncorrected_true_configuration_recovery` builds at -O3 and runs.
- `UNCORRECTED_TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK` sets the run block; the
  compile-time overrides are those of the grid above.
- Progress goes to
  `out/montecarlo/uncorrected_true_configuration_recovery_progress.csv`, the
  result to `out/montecarlo/uncorrected_true_configuration_recovery.csv.gz`,
  with the same columns.
- `make mc-corrected_against_uncorrected_grid` builds
  `montecarlo/seri_corrected_against_uncorrected_grid.c` pointed at the two
  bootstrap grids and writes
  `out/montecarlo/corrected_against_uncorrected_grid_report.txt`, the same
  tables as the Seri comparison. It reads the corrected grid's finished result
  only, so it can run once the grid above is done, while this one is still
  going.

Plan, none of it run:

1. After the grid above finishes, a trial over runs 0 and 1 of the t-QVARMA
   (`FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=2`, 2000 benchmarks, about 4 minutes)
   into `out/montecarlo/uncorrected_true_configuration_recovery_trial/`, and the
   comparison report on those 2000 benchmarks.
2. Then the full grid as the user service
   `uncorrected-true-configuration-recovery`, all 16 threads, memory capped at
   5 GB. The confidence set is the cost, not the loss, so expect about what the
   grid above takes: 0.09 to 0.125 s per benchmark, 125 to 175 hours.
3. It cannot run beside the grid above: each takes every thread.

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
`out/montecarlo/seri_true_configuration_recovery.log`.

`make mc-seri_true_configuration_recovery_report` summarises the rows written so
far, or the finished result, into
`out/montecarlo/seri_true_configuration_recovery_report.txt`: per model, how
often the benchmark's configuration is in the set, alone and first by mean loss,
the set sizes, how many configurations are recovered in at least half and 90% of
their own benchmarks, and, on the t-QVARMA benchmarks the bootstrap grid above
also finished, how often each procedure's set holds the configuration. It can
run while the grid does, and reads only.

The grid finished on 2026-10-14 at 20:38, 5,000,000 benchmarks, with no error
in its log. Over all 1,000,000 benchmarks per model, alpha 0.05, from
`out/montecarlo/seri_true_configuration_recovery_report.txt`:

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
`out/montecarlo/seri_against_naive_grid_report.txt`. Measured 2026-10-15:

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

`make mc-seri_against_naive_by_configuration` computes, for both grids from the
same code, the per-configuration tables the naive grid's Python scripts wrote
(`sweep_grid_cop_0191.py`, `sweep_grid_presence.py`,
`sweep_grid_identifiability_plots.py`): `cop_0191` as the truth, as a wrong
answer and every configuration as the truth, into
`out/montecarlo/seri_against_naive_by_configuration_cop_0191.txt` and `.tex`;
each configuration's own and off-target rate into
`out/montecarlo/seri_against_naive_by_configuration_presence.txt`; and per
configuration the share of its own sets holding it and holding the most frequent
other configuration into
`out/montecarlo/seri_against_naive_by_configuration_identifiability.csv`, which
`montecarlo/seri_against_naive_by_configuration_plots.py` draws. Its naive-grid
presence and identifiability output equals the old scripts' files line for line,
and its naive `cop_0191` tables equal theirs to the printed digit (checked
2026-10-15).

`make mc-response_variance_by_configuration` measures each configuration's
noise, its mean response variance over all 1000 runs averaged over the entries
(the term the variance correction subtracts), into
`out/montecarlo/response_variance_by_configuration.csv`, and writes
`out/montecarlo/response_variance_by_configuration_report.txt`: where `cop_0191`
ranks by noise, and the Spearman correlation across the 1000 configurations
between noise and the own and off-target rates of each grid. Measured
2026-10-15:

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
`out/montecarlo/seri_corrected_against_uncorrected_grid_report.txt`. It reads
the uncorrected grid's progress file while the result does not exist, so it can
run while the grid does. Per model and grid: the benchmark's configuration in
the set, alone, first by mean loss, mean rank, set sizes; the two by two table
of in-set and how often the smallest-loss configuration agrees; per
configuration the own rate summaries and the correlation of the two grids' own
rates; and, against each configuration's noise from
`out/montecarlo/response_variance_by_configuration.csv`, the Spearman
correlation of noise with the own and off-target rates and the mean noise rank
(1 the least noisy) of the smallest-loss configuration and of the set members.
The noise ranks are the direct measure of the bias the correction is meant to
remove: with no dependence on noise they sit near the middle rank, 500.5.

### Trial over runs 0 and 1 of the t-QVARMA

Setup: a test build with `FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=2`, so runs 0
and 1 of all 1000 configurations as the benchmark under the t-QVARMA, 2000
benchmarks, every other setting as in the full grid; the comparison run
against the corrected grid's finished result on those 2000 benchmarks.
Machine: AMD Ryzen 7 4800H, 16 hardware threads, nothing else heavy running,
under a 4 GB memory cap. Rows, log and report in
`out/montecarlo/seri_uncorrected_true_configuration_recovery_trial/`, measured
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
`out/montecarlo/seri_uncorrected_true_configuration_recovery.log`, with about 26
hours estimated at its start.

### Trial over runs 0 and 1 of the t-QVARMA

Setup: a test build with `FIRST_MODEL=0 LAST_MODEL=0 RUN_LIMIT=2`, so runs 0 and
1 of all 1000 configurations as the benchmark under the t-QVARMA, 2000
benchmarks, every other setting as in the full grid. Machine: AMD Ryzen 7 4800H,
16 hardware threads, nothing else heavy running, under a 6 GB memory cap. Rows
in `out/montecarlo/seri_true_configuration_recovery_trial/progress.csv`, log in
`out/montecarlo/seri_true_configuration_recovery_trial/trial.log`, measured
2026-10-13.

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

## The Seri grid with the mean taken from runs outside the loss matrix

In the corrected Seri grid every row of a configuration's column is centred on
the mean of all the column's rows, so the rows of the loss matrix are not
independent, which `seri_mcs` assumes.
`montecarlo/seri_held_out_mean_true_configuration_recovery.c` splits each
configuration's 1000 runs: the last S = 100 (runs 900 to 999) give the mean
response and nothing else, the first 900 give the benchmarks and the rows of
the loss matrices. For loss run s of configuration i against benchmark d,

    L*_is = (1/K) sum_h (IR_ish - IR_dh)^2 - S/(S+1) (1/K) sum_h (IR_ish - IR_bar_ih)^2,

with IR_bar_ih the mean of entry h over the 100 mean runs. The factor c is the
one that makes the mean loss over the n loss rows, which is what `seri_mcs`
compares, have expectation (IR_ih - IR_dh)^2 with no sigma2_ih left. With
x_s the loss runs, xbar and s2 their mean and sample variance, m the mean of
the S mean runs and d the benchmark, per entry,

    (1/n) sum_s L*_is = (xbar - d)^2 + (1 - c) (n-1)/n s2 - c (xbar - m)^2,

whose expectation is (mu - d)^2 + sigma2 [1/n + (1 - c)(n-1)/n - c (1/n + 1/S)].
The bracket is zero at c = S/(S+1). In the corrected Seri grid m is xbar
itself, the last term is gone, and the bracket 1/n + (1 - c)(n-1)/n is zero at
c = n/(n-1): the n there is the 1/n in the expectation of (xbar - d)^2, from
averaging the rows. Keeping n/(n-1) here, with n = 899 and S = 100, would leave
the bracket at about -0.011, a loss lowered by 1.1% of each configuration's own
variance.
`own_deviation_loss_held_out_mean_deviations` in
`montecarlo/own_deviation_loss.h` computes the term, and
`tests/own_deviation_loss_held_out_mean_correctness.c` checks it against a
reference written from the definition and checks, on simulated responses
(K = 20, S = 5, 200 deviation runs, 2000 draws, seed 20261125), that its
average is the mean response variance to 0.15% while the term without
S/(S+1) is 19.8% above it (passed 2026-11-25). Everything else is the Seri
grid's: the benchmark's run index held out of every column, `seri_mcs` at
alpha 0.05, the same columns. 900 benchmarks per configuration, 900,000 per
model, for the linear projection, both states under the causal HP cycle at the
estimated lambda and the t-QVARMA with its two growth states.

### How many runs for the mean

`studies/response_variance_run_count.c` (`make study-response_variance_run_count`)
measures how many runs, taken from the last one backwards, the summed response
variance needs. For configuration i and m runs (runs 1000-m to 999), V_i(m) is
the sum over every entry of the response vector of the entry's sample variance
(divisor m - 1), and d_i(m) = |V_i(m) - V_i(1000)| / V_i(1000). A model's run
count is the smallest m from which on, for every larger m too, d_i(m) < 0.05
in at least 95% of the 1000 configurations. The runs given to the mean were
bounded at 100 per configuration. Cached responses of all 1000 runs of all
1000 configurations, no randomness. Output in
`out/studies/response_variance_run_count_report.txt` and
`out/studies/response_variance_run_count_by_size.csv`, measured 2026-11-25:

- LP linear: run count 386. At m = 100 the median d is 0.037, the 95th
  percentile 0.112, 62.7% of configurations below 0.05.
- NLP both states (causal HP cycle, estimated lambda): run count 232. At
  m = 100, median 0.026, 95th percentile 0.081, 77.8% below 0.05.
- t-QVARMA growth states: run count 939. At m = 100, median 0.126, 95th
  percentile 0.640, largest 3.45, 22.3% below 0.05; at m = 500 still only
  50.6% below.

No model meets the criterion within 100 runs, so the grid uses the bound,
S = 100, for every model.

The study also writes every configuration's V_i(m), m = 2 to 1000, to
`out/studies/response_variance_run_count_by_configuration.csv`, and
`studies/response_variance_run_count_plots.py` (polars, plotly, kaleido 1)
draws V_i(m) / V_i(1000) against m as quantile bands across configurations into
`out/studies/response_variance_run_count_plots/expanding_variance_ratio.pdf`.

### Trial and full run

Trial: a test build with `RUN_LIMIT=1 BENCHMARK_LIMIT=64`, so run 0 of
`cop_0001` to `cop_0064` under the three models, 192 benchmarks; no response
missing, 100 mean runs and 900 loss runs per model. Rows and log in
`out/montecarlo/seri_held_out_mean_trial/`.

The full grid was started on 2026-11-25 at 23:42 by
`out/montecarlo/run_seri_held_out_mean_grid_detached.sh` (detached, 12 GB
address-space cap, all 16 hardware threads), log
`out/montecarlo/seri_held_out_mean_grid.log` and
`out/montecarlo/seri_held_out_mean_true_configuration_recovery.log`. Progress
goes to `out/montecarlo/seri_held_out_mean_true_configuration_recovery_100_mean_runs_progress.csv`,
the result to `out/montecarlo/seri_held_out_mean_true_configuration_recovery_100_mean_runs.csv.gz`;
the file names carry S so a grid under another split never resumes from these
rows. Its first minutes ran at about 0.015 s per benchmark on the linear
projection.

## The four grids under the other three state weights

Besides the state-dependent projection of lpirfs' weight, the two-sided HP
cycle at lambda 1600, the four grids (bootstrap or Seri, corrected or
uncorrected loss) run state 1, state 2 and both
states under the three other weights of "The local projections" below: the
causal cycle at lambda 1600, the two-sided cycle at the estimated lambda, and
the causal cycle at the estimated lambda. That is nine more models, numbers 5
to 13 of `montecarlo/response_cache.h`, labelled `lp_s1_causal_hp`,
`lp_nl_estimated_lambda`, `lp_s2_causal_hp_estimated_lambda` and so on, each
over the same 1,000,000 benchmarks, so 14 million confidence sets per grid. All
fits are of the 400-period samples; the 188-period fits are not used, as for
lpirfs' weight. The linear projection does not depend on the weight and is not
rerun.

Each grid finds its finished result, unpacks it into its progress file, runs
only the new models and writes the result again with all of them. The results
as they were before are kept in `out/montecarlo/before_state_weights/`. The two
Seri grids ran all nine; the rows of every discarded model were then moved to
the two folders of `out/old/` above, so each grid's result holds three models,
and the bootstrap grids, when started, run only model 13, both states under
the causal HP cycle at the estimated lambda, 1 million confidence sets each.
The programs that compare grids against the naive grid
(`seri_against_naive_grid`, `seri_against_naive_by_configuration`) and the
programs under `old/` stay on the original five models.

The response caches of the nine models are in
`out/montecarlo/sweep_grid_response_cache/`, built 2026-11-19, about 11 minutes
per state model and seconds per both-states model.
`tests/response_cache_state_weight_correctness.c` checks the cache under every
weight against the archives on three configurations, and that lpirfs' caches
are unchanged byte for byte. `tests/response_cache_state_weight_archive_agreement.c`
checks every entry of the nine new caches, 1000 configurations by 1000 runs,
against the archives of its weight: equal exactly, none missing (passed
2026-11-20, after the machine went down without a clean shutdown a minute after
the last cache was written).

Unclean shutdowns leave a progress file ending in NUL bytes, sometimes after a
row cut in half. Each grid cuts its progress file at the line holding the
first NUL byte, or after its last newline, before resuming
(`montecarlo/progress_file.h`, checked by
`tests/progress_file_torn_tail_repair.c`), so those benchmarks are redone
instead of being read as finished.

### Trial

Setup: test builds of the four grids, rows in
`out/montecarlo/state_weight_trial/`, log in `trial.log` there. Machine: AMD
Ryzen 7 4800H, 16 hardware threads, nothing else heavy running, under a 10 GB
address-space cap. Run 2026-11-19 from 19:42, interrupted by the shutdown at
about 20:50 and resumed 2026-11-20 at 01:22.

- Seri corrected and Seri uncorrected, `FIRST_MODEL=5 LAST_MODEL=13
  RUN_LIMIT=1`: run 0 of all 1000 configurations under the nine models, 9000
  benchmarks each. Both hold 9000 distinct rows, 1000 per model.
- Time, Seri corrected: about 0.035 s per benchmark once the caches exist
  (over the 936 benchmarks after the first batch of the both-states models),
  which would give about 90 hours for its 9 million; the full run's opening
  rate is lower, below.
- Time, Seri uncorrected: 0.025 s per benchmark over its first 1256, 0.017 s
  over the 7744 after the restart, so 45 to 60 hours.
- Bootstrap corrected and uncorrected, also `BENCHMARK_LIMIT=64`: run 0 of
  `cop_0001` to `cop_0064` under the nine models, 576 benchmarks each, at
  0.149 and 0.132 s per benchmark wall clock including the variances of the
  run index. At the 0.09 to 0.125 s the corrected bootstrap grid ran at in
  October, the two need 450 to 750 hours together.
- Seri corrected under `lp_nl` of lpirfs' weight, run 0 of `cop_0001` to
  `cop_0064`: the 64 rows equal byte for byte the same rows of the finished
  `seri_true_configuration_recovery.csv.gz`, so the change to the model table
  leaves lpirfs' weight computing what it did.

### Full run

Two launchers, each detached, 10 GB address-space cap, all 16 hardware
threads, running its two grids one after the other:

- `out/montecarlo/run_seri_state_weight_grids_detached.sh`: Seri corrected,
  then Seri uncorrected; log `out/montecarlo/seri_state_weight_grids.log`.
- `out/montecarlo/run_bootstrap_state_weight_grids_detached.sh`: the bootstrap
  grids (Hansen et al.'s MCS), corrected, then uncorrected; log
  `out/montecarlo/bootstrap_state_weight_grids.log`.

Each grid appends to its own log (`out/montecarlo/<grid>.log`); the launcher's
log records the start and end of each grid, whether the earlier result is
unchanged at the start of the new one, and the rows per model. A launcher can
be started again after an interruption: a finished grid reports nothing to do,
an interrupted one resumes.

The Seri launcher runs since 2026-11-20; the bootstrap launcher has not been
started, so the bootstrap columns of C, E and C E in the identification table
are "-". Seri corrected ran at about 0.02 s per benchmark (0.0196 s over its
first 1.54 million) and finished on 2026-11-22.

`out/montecarlo/run_discarded_model_split_detached.sh` (run 2026-11-24, after
the Seri launcher finished) moves the rows of the discarded models out of the
four results with `bin/single_state_rows_split` and
`bin/other_state_weights_rows_split`, checks test builds of the Seri and
bootstrap corrected grids, in their `montecarlo/` version and their two old
copies, against the parted results (run 0 of `cop_0001` to `cop_0064`, models
0 to 13 for Seri and 0 to 4 for the bootstrap grid, byte for byte), and
rewrites every report of `montecarlo/` and of both old folders. Outputs and log
in `out/montecarlo/discarded_model_split/`.

## One table for the four grids

`make mc-identification_table` writes, in `out/montecarlo/identification_table/`,
`identification_table.tex`: the four grids above (Hansen et al.'s bootstrap
set or `seri_mcs`, each with the naive or the adjusted loss) as one LaTeX
table, and per grid the share of benchmarks whose set holds the true
configuration, the share where it has the smallest mean loss, and the median
final set size. Rows: LP (linear), NLP (both states in one vector under the
causal HP cycle at the estimated lambda) and QVARMA states (the t-QVARMA with
its two growth states in one vector). The discarded models have
their own tables, written by the `identification_table` copy of each old
folder into its `identification_table/` subfolder of `out/old/`. A cell is
filled only when its grid holds all 1,000,000 benchmarks of that model and is
"-" otherwise, so it can be rerun as grids finish. `identification_table_counts.txt` beside it lists the file each grid
was read from and its benchmarks per model, and `identification_table.log`
when each grid was read and how long it took.

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

Those fits build the state weight as lpirfs does, from the two-sided HP cycle
of the switching series at lambda 1600. Three other weights are fitted beside
them, each by the same script with a flag and into its own tree
(`applications/lp_system.h`, `LpSystemStateWeight`):

- `--causal-hp`: the cycle at `t` from periods `1..t` only, et_al's
  `lp_lin_and_nl_causal`, at lambda 1600; `abm_system_fit_lp_causal_hp/`.
- `--estimate-lambda`: the two-sided cycle at the lambda et_al's `hp_estimate`
  finds by maximum likelihood on each replicate's own switching series;
  `abm_system_fit_lp_estimated_lambda/`.
- both flags: the causal cycle at that estimated lambda;
  `abm_system_fit_lp_causal_hp_estimated_lambda/`.

`make app-abm_system_fit_lp<window>_<weight>` runs one, `<window>` empty for
the full sample or `_t188` or `_t188_from188`. Each row of an estimated-lambda
archive stores its own lambda and the search's flags, and the manifest gives
their quantiles. `tests/lp_system_state_weight_correctness.c` checks every
weight's rows against et_al called directly and the archives' round trip.

The US data are fitted under the same weights by `applications/us_lp_fit.c`
with the same flags, `make app-us_lp_fit_<weight>`, into
`out/applications/us_lp_fit_<weight>.npz`, `.txt` and `_irf.csv`. On the US
switching series the estimated lambda is 4.54e-05, the lower end of the
search interval (log q at its upper bound 10), under both the two-sided and
the causal cycle.

The linear projection and d do not depend on the state weight, so they are
stored once per sample, in the archives of lpirfs' weight
(`abm_system_fit_lp<window>/`, `us_lp_fit.npz`). The archives of the other
three weights hold the state-dependent responses, their fit notes and the VAR's
statuses, and no `lin_`, `lin_ols_status_` or `d_` column; a reader that needs
the linear responses of a replicate reads them from the lpirfs archive of the
same sample. Each of those archives is about two thirds of an lpirfs one.

Optimization to do: the linear projection and d are still computed once per
weight, four times per replicate and sample, since each weight's run calls
et_al's `lp_lin_and_nl` or `lp_lin_and_nl_causal`, which fit both models.
Computing them once needs a state-dependent fit in et_al that takes d and the
VAR's statuses as inputs instead of fitting the VAR itself, for the two-sided
and the causal cycle; et_al has none yet. Then one run per sample could fit
the linear projection once and the state-dependent projection under all four
weights from its d. Not measured; from the regressor counts, 21 per horizon in
the linear regressions against 41 in the state-dependent ones, the linear part
is perhaps a fifth of the 54 minutes the three other weights took over the
three samples on 2026-11-19, so about 10 minutes per full refit.

## The t-QVARMA conditional on the growth state

The t-QVARMA's response depends on the data only through D, the Jacobian of
the score with respect to the standardised shock. The model above uses the
sample mean of D_t over all periods. `qvarma_growth_states` uses two
conditional means instead, one per state of the economy, and stacks the two
responses into one vector, as the NLP stacks its two states.

The state of period t is decided by g_(t-1), the mean of GDP growth over the
four periods t-4..t-1. Period t is in recession when g_(t-1) is at or below the
20th percentile of g_(t-1) over the run's periods 5 to 400, and in expansion
otherwise: 80 recession and 316 expansion periods per run. Periods 1 to 4 have
no g_(t-1) and are in neither state. D_expansion and D_recession are the plain
means of D_t over each state's periods, and each gives the 525 total responses
of the same fit, horizons 0 to 20; the vector is the expansion responses
followed by the recession ones, 1050 entries. No weight is shared with the NLP
and nothing estimated for it is used.

The state is lagged because D_t is a function of the period's own residual: a
period of very low growth is mostly one with a large negative GDP shock, so
conditioning on growth at t would select periods with large residuals, and so
small D_t, by construction.

`./bin/response_cache_build qvarma_growth_states` builds the model's response
cache without running a confidence set: the filter of every cached t-QVARMA
fit is run again to get D_t, nothing is refitted. On 2026-11-24 it took 336 s
on 16 threads, 4.2 GB, and no run lacked a response.
`out/montecarlo/run_response_cache_build_detached.sh <model>` runs it detached
and then `tests/qvarma_growth_states_correctness.c`, which on runs 0 to 9 of
`cop_0001`, `cop_0500` and `cop_1000` recomputes the split and the two means
from their definitions and compares them, and the cached responses, exactly.

All four grids run it. The two bootstrap grids have not been run on it:
`out/montecarlo/run_bootstrap_qvarma_growth_states_grids_detached.sh` runs
them, corrected loss then uncorrected, adding the model to their results the
way the Seri launcher below does, log
`out/montecarlo/bootstrap_qvarma_growth_states_grids.log`. A test build
of each Seri grid over run 0 of the first 128 configurations
(`out/montecarlo/qvarma_growth_states_trial/`) found no missing response and
took 0.049 s per benchmark under the corrected loss and 0.024 s under the
uncorrected one. The full runs were started on 2026-11-24 at 19:13 by
`out/montecarlo/run_seri_qvarma_growth_states_grids_detached.sh`, corrected
loss first, each adding the model to the result that already holds the other
three, whose earlier copy is kept in `out/montecarlo/before_qvarma_growth_states/`
and checked to be unchanged at the start of the new result; log
`out/montecarlo/seri_qvarma_growth_states_grids.log`.

## A separate pipeline that refits nothing

This is its own pipeline, in its own directory, writing to its own
`out/montecarlo/`. Nothing under `applications/` is modified or read as source:
the main pipeline stays exactly as it is, and a change here cannot affect a
real-data result.

What is reused is that pipeline's **output**:
`out/applications/abm_system_fit_qvarma/`, the million cached t-QVARMA fits, and
`out/applications/abm_system_fit_lp/`, the million local-projection fits. No
model is estimated anywhere in this experiment.

`montecarlo/response_cache.h` holds every run's response vector under one model
as float32, in `out/montecarlo/sweep_grid_response_cache/<model>.f32` (2.1 GB
for the t-QVARMA, 1.6 GB per local projection, 3.2 GB for both states, ignored
by git). The files were built from the fits by the naive-loss grid between
2026-09-29 and 2026-10-05, after the last local-projection refit (2026-09-27)
and the last t-QVARMA fit (2026-09-10). Both states' file holds, for every run,
its 400 state-1 entries followed by its 400 state-2 entries, bit for bit those
of the two state files. A file is mapped into memory rather than read into it,
so when memory is short the system drops parts of it and reads them back from
disk. A cache file is rebuilt only when its header (model, configurations, runs,
entries per run) disagrees with what the program expects; a change of the fits
or the dataset with the same shape is not detected, so delete the file by hand
after either changes.

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
