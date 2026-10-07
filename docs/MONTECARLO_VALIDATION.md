# Monte Carlo validation: does the procedure recover a known answer

## Overview

The main pipeline asks which ABM configuration sits closest to the US data and
which of them cannot be told apart from it. That question has no known answer,
so a result can be reported but not checked. This experiment asks the same
question with the answer known in advance.

One simulated run is promoted to the role the US data plays. Every loss is
measured against it, the three validation protocols run exactly as they do in
the main pipeline, and the question becomes: does the confidence set come back
holding the configuration the benchmark was drawn from? It should. The
benchmark is one of that configuration's own replicates, so the configuration
that generated it is, by construction, the right answer.

This is a check on the procedure, not on the ABM. A confidence set that misses
the configuration its own benchmark came from would say the procedure cannot
identify a configuration even when one is definitely there, which would make
the main pipeline's result hard to interpret. A confidence set that finds it
says the machinery works on a case where the answer is known, and says nothing
about whether the DSK model resembles the United States.

`docs/ABM_SYSTEM_MCS_VALIDATION.md` describes the impulse-response protocol and
`docs/ABM_SYSTEM_SCORE_LOSS.md` the two score protocols. Everything they say
about settings, losses and their limitations applies here unchanged; only the
benchmark differs.

## The experiments

Five experiments, each a larger set of benchmarks than the one before. In every
one the benchmark's run index is held out of every configuration's column, the
loss is measured against the benchmark, and the question is whether the Model
Confidence Set (MCS_TR, alpha 0.05, block length 1, bootstrap variance, seed 123
stream 0) holds the benchmark's own configuration. The last row of the table
below asks something else of the same confidence sets: how far each set's best
configuration is from the benchmark, measured against the sampling standard
deviation of the responses.

| experiment | benchmarks | losses and auxiliary models | resamples | program | output | written up in |
| --- | --- | --- | --- | --- | --- | --- |
| one benchmark | `cop_0191` run 706 | impulse response, `q'q`, `LM` (t-QVARMA); impulse response (four local projections) | 10000 | `montecarlo/run.sh`, `montecarlo/lp_run.sh` | `irf_loss.csv` to `mcs_statistic_comparison*`, `lp_*` | this file, "Results"; `docs/MONTECARLO_LP_VALIDATION.md` |
| every run of `cop_0191` | 1000 | impulse response (t-QVARMA, four local projections); `q'q`, `LM` stopped at 389 | 10000 | `montecarlo/sweep.sh`, `bin/lp_sweep` | `sweep_irf.csv`, `sweep_score*.csv`, `lp_sweep_*.csv` | this file, "Every replicate as the benchmark"; `docs/MONTECARLO_LP_VALIDATION.md` |
| run 0 of every configuration | 1000 per model | impulse response, five models | 10000 | `bin/sweep_cops` | `sweep_cops.csv.gz` | this file, "The first run of every configuration" |
| every run of every configuration | 1,000,000 per model | impulse response, five models | 2000 | `bin/sweep_grid` | `sweep_grid.csv.gz` | this file, "Every run of every configuration" |
| compressed responses | 5000 or 2000 per vector | the same responses mapped to a few coordinates first | 2000 | `montecarlo/compressed_response_*` | `compressed_response/` | `docs/MONTECARLO_COMPRESSED_RESPONSE.md` |
| fit metric | 1,000,000 per model (oracle), 1000 per model (bootstrap) | the best configuration's mean loss against the responses' sampling standard deviation, five models | none | `montecarlo/fit_metric_run.sh` | `fit_metric_*` | `docs/MONTECARLO_FIT_METRIC.md` |

The five models are the t-QVARMA (p1q1r2, total responses to horizon 20, 525
entries per run) and the four local projections of
`docs/MONTECARLO_LP_VALIDATION.md`: linear, state 1 and state 2, 400 entries per
run, and both states, the state-dependent model's two states stacked into one
vector of 800 entries as `_temp/Note on non-lin LP.pdf` prescribes. `q'q` and
`LM` are the two score losses of `docs/ABM_SYSTEM_SCORE_LOSS.md`, the plain and
the inverse-information-weighted squared score.

## A separate pipeline that refits nothing

This is its own pipeline, in its own directory, writing to its own `out/`.
Nothing under `applications/` is modified or read as source: the main pipeline
stays exactly as it is, and a change here cannot affect a real-data result.

What is reused is that pipeline's **output**, and the important one is
`out/abm_system_fit_qvarma/`, the million cached fits. No model is estimated
anywhere in this experiment. The benchmark's own parameters are one of those
cached fits, promoted; every replicate's parameters are read from the same
cache. That is what makes the whole experiment a few minutes rather than
another full estimation run.

| file | what it is |
| --- | --- |
| `montecarlo/benchmark_choice.c` | picks which run is promoted, and records why |
| `montecarlo/benchmark.h` | reads that choice back, so every step uses the same run |
| `montecarlo/irf_loss.c` | the impulse-response loss against the benchmark |
| `montecarlo/score_loss.c` | both score losses, in one pass |
| `montecarlo/mcs.c` | the headline confidence set |
| `montecarlo/mcs_statistic_comparison.c` | all three losses under both statistics |
| `montecarlo/run.sh` | the driver |
| `montecarlo/sweep_irf.c` | every run of `cop_0191` as the benchmark, impulse-response loss |
| `montecarlo/sweep_score.c` | the same, both score losses |
| `montecarlo/sweep.sh` | runs the two sweeps above in turn |
| `montecarlo/lp_irf_loss.c`, `lp_mcs.c`, `lp_run.sh` | the one-benchmark experiment under the local projections |
| `montecarlo/lp_sweep.c` | every run of `cop_0191` as the benchmark, local projections |
| `montecarlo/response_cache.h` | every run's response vector under one model, held as float32 and kept on disk |
| `montecarlo/sweep_cops.c` | run 0 of every configuration as the benchmark, five models |
| `montecarlo/sweep_grid.c` | every run of every configuration as the benchmark, five models |
| `montecarlo/lp_nl_run.sh` | runs every experiment after the single benchmark for both states, resuming where it stopped |
| `montecarlo/sweep_grid_presence.py` | how often each configuration sits in sets that are not its own |
| `montecarlo/sweep_grid_cop_0191.py` | `cop_0191` as the truth, as a wrong answer, and every configuration as the truth, `docs/MONTECARLO_COP_0191.md` |
| `montecarlo/fit_metric.h`, `fit_metric_oracle.c`, `fit_metric_bootstrap.c`, `fit_metric_report.c`, `fit_metric_run.sh`, `fit_metric_v_histograms.py` | the fit metric, oracle and bootstrap versions, `docs/MONTECARLO_FIT_METRIC.md`; the planned comparison across methods, `docs/MONTECARLO_FIT_METRIC_COMPARISON.md` |
| `montecarlo/sweep_grid_identifiability_plots.py` | per-configuration recovery figures from the same rows |
| `montecarlo/compressed_response_learn.py`, `compressed_response_recovery.c`, `compressed_response_grid_baseline.py`, `compressed_response_report.py` | `docs/MONTECARLO_COMPRESSED_RESPONSE.md` |

The four loss and confidence-set scripts began as copies of their
`applications/` counterparts and differ from them in three things: the
benchmark is a cached simulated fit rather than the US fit, the benchmark's
replicate is held out of every column, and the paths point into
`montecarlo/out/`. They are copies, so a change to the protocol in one tree has
to be made in the other; that is the price of the two pipelines being
independent, and it is paid deliberately.

## Which run stands in for the real data

Choosing the benchmark on grounds connected to the answer would rig the
experiment, so the grounds are fixed in advance and applied by
`montecarlo/benchmark_choice.c`, which writes its reasoning to
`montecarlo/out/benchmark_choice.txt`. Three requirements, in order.

**It must come from a configuration a real-data confidence set kept.** The
union of the three procedures' surviving configurations under MCS_TR (the
`tr_in_set` column of `out/abm_system_mcs_statistic_comparison*.csv`; the
MCS_TMAX sets are not read) is the candidate pool.
Promoting a configuration none of them kept would ask whether a poor benchmark
is recovered, which is a different question.

**Among those, the configuration whose fits converged most often.** About a
third of the million fits in this project converged, which makes convergence
the property in shortest supply. A benchmark whose own fit did not converge is
a point the optimizer stopped at, not an estimate, and the whole experiment
rests on the benchmark being a real estimate.

**Within it, the converged replicate with the smallest gradient norm whose
observed information matrix is positive definite.** The gradient norm orders
the candidates; positive definiteness is a gate rather than a preference,
because the weighted score loss inverts that matrix and a benchmark it cannot
be built at would leave one of the three protocols unrunnable. Fits the
manifest marks as inherited from another replicate are excluded; the current
manifest (`out/abm_system_fit_qvarma_manifest.txt`, 2026-09-10) has none, its
p1q1r2 rows being 352,470 converged and reused, 621,485 unconverged and held,
26,045 unconverged and resumed.

These three requirements apply to this one benchmark only. The later
experiments use every run of a configuration as a benchmark with no filter, so
most of their benchmarks are unconverged fits; "Recovery by whether the
benchmark's fit converged" below measures what that changes.

### What it selected

| step | result |
| --- | --- |
| candidate pool | `cop_0191`, `cop_0931`, `cop_0409`, `cop_0599` |
| convergence rate | 44.9%, 44.1%, 22.9%, 33.2% respectively |
| configuration chosen | `cop_0191` |
| replicate chosen | 706, the first candidate tried |

`cop_0191` replicate 706: log-likelihood 124.6637, gradient norm 4.058e-03,
converged, its own fit rather than an inherited one. Its observed information
matrix has smallest eigenvalue 6.699e-01 and condition number 1.119e+06, so the
weighting is buildable and about as well conditioned as the real data's own
(1.921e+06).

For scale, the US fit that the main pipeline uses converged with gradient norm
4.713e-03, so the benchmark is an estimate of comparable quality to the one it
replaces.

`cop_0191` is also the configuration the main pipeline's own headline
confidence set keeps. That is a consequence of the convergence criterion rather
than a choice, and it does not privilege any of the three protocols here: what
each protocol is being asked is whether it recovers `cop_0191`, and the
impulse-response protocol gets no help from having picked it against different
data.

## What is held out, and why

Replicate 706 of **every** configuration is dropped from every loss matrix, not
just `cop_0191`'s.

The benchmark's own cell is the obvious part: a replicate scored against itself
sits at a loss of essentially zero and would win by construction. The rest is
about the seeds. `applications/abm_system_simulate.c` runs every configuration
on the same seeds, `seed = mc + 1`, so replicate 706 of every configuration is
seed 707 of the model. Leaving those in would let every configuration be judged
partly on a draw the benchmark also used.

The matrices therefore have 999 rows rather than 1000, and every configuration
loses the same one, so the comparison between configurations is unchanged.

In practice the second part of this is a precaution rather than a correction.
The project already measured what the shared seeds buy:
`docs/ABM_SYSTEM_SIMULATION.md`, "Common random numbers", reports a mean-loss
correlation across configuration pairs of 0.0017 at the same seed against
-0.0001 at a shifted one, because the model consumes a different number of
draws per period at different parameters and the streams desynchronise. So
replicate 706 of `cop_0413` has almost nothing in common with replicate 706 of
`cop_0191` beyond its label. The hold-out costs one observation in a thousand
and removes the question entirely, which is cheaper than arguing it away.

## Running it

    make montecarlo            the whole experiment
    ./montecarlo/run.sh        the same thing directly
    ./montecarlo/run.sh --keep reuse the benchmark already chosen

It reads `out/abm_system_fit_qvarma/` and `dataset/abm_system/` and rebuilds
neither, so it does not run on a fresh clone. The impulse-response loss is the
long step, since it reads a million cached fits; the two score losses come from
one pass over the dataset.

`make montecarlo` builds every program in `MONTECARLO_STEMS` (all of
`montecarlo/` except `sweep_grid`) and runs `run.sh` only. Each program also has
its own target, `make mc-<stem>`, which builds and runs it alone.

### Settings read from the environment

All optional; the default is what every result in this file used.

| variable | read by | default |
| --- | --- | --- |
| `ABM_SYSTEM_FIT_DIR` | `irf_loss`, `score_loss` | `out/abm_system_fit_qvarma` |
| `ABM_SYSTEM_INPUT_DIR` | `irf_loss`, `score_loss`, `lp_irf_loss`, `lp_sweep`, `sweep_cops`, `sweep_grid` | `dataset/abm_system` |
| `ABM_SYSTEM_LP_FIT_DIR` | `lp_irf_loss`, `lp_sweep`, `sweep_cops`, `sweep_grid` | `out/abm_system_fit_lp` |
| `ABM_SYSTEM_LOSS_PATH` | `irf_loss` | `montecarlo/out/irf_loss.csv` |
| `ABM_SYSTEM_SCORE_LOSS_PATH` | `score_loss` | `montecarlo/out/score_loss.csv` |
| `LP_SWEEP_WORKERS` | `lp_sweep`, confidence sets side by side | every hardware thread |
| `SWEEP_COPS_WORKERS` | `sweep_cops`, confidence sets side by side | every hardware thread |
| `SWEEP_GRID_WORKERS`, `SWEEP_GRID_RUN_BLOCK` | `sweep_grid`, see its section | every hardware thread, all runs |

The other programs take their paths as fixed constants. In particular
`sweep_irf`, `sweep_score`, `benchmark_choice` and the t-QVARMA half of
`montecarlo/response_cache.h` (used by `sweep_cops` and `sweep_grid`) always
read `out/abm_system_fit_qvarma`, whatever `ABM_SYSTEM_FIT_DIR` says.

### Outputs

Everything lands in `montecarlo/out/`.

| file | what it holds |
| --- | --- |
| `benchmark_choice.txt` | which run was promoted and why, with the numbers at each step |
| `benchmark.env` | the choice, as two shell assignments the driver reads |
| `irf_loss.csv`, `irf_loss_manifest.txt` | the impulse-response loss against the benchmark |
| `score_loss.csv`, `score_loss_weighted.csv`, `score_loss_manifest.txt` | the two score losses |
| `mcs.txt`, `mcs.csv` | the headline confidence set, MCS_TR over the impulse-response loss |
| `mcs_statistic_comparison*.txt`, `*.csv` | all three losses under both statistics, six confidence sets |

The loss matrices are 19 MB each and are ignored by git for the same reason the
main pipeline's are.

## Results

### Setup

Benchmark `cop_0191` replicate 706, 400 periods, its own cached p1q1r2 fit.
Field of 1000 configurations by 999 replicates, replicate 706 held out of every
column, no missing cells in any of the three matrices. Confidence set settings
identical to the main pipeline's and to each other: alpha = 0.05, 10000
bootstrap resamples, block length 1, bootstrap variance of the resampled mean,
seed 123, stream 0, both statistics scored on the same resamples. Every fit
read from `out/abm_system_fit_qvarma/`; nothing estimated.

### The impulse-response protocol recovers the benchmark exactly

| loss | statistic | set | is the benchmark's configuration in it |
| --- | --- | --- | --- |
| IRF | MCS_TR | `cop_0191` | yes |
| IRF | MCS_TMAX | `cop_0191` | yes |
| `q'q` | MCS_TR | `cop_0429` | no |
| `q'q` | MCS_TMAX | `cop_0429` | no |
| `LM` | MCS_TR | 7 configurations | no |
| `LM` | MCS_TMAX | 21 configurations | no |

Under the impulse-response distance `cop_0191` has the lowest mean loss of all
1000 configurations, 0.01006 against a runner-up of 0.01232 and a field median
of 0.02711, and it is the single survivor under both statistics. The procedure
was handed a question with a known answer and returned that answer.

### Neither score protocol recovers it

| loss | benchmark's rank | benchmark's mean loss | winner's | ratio | benchmark eliminated in round |
| --- | --- | --- | --- | --- | --- |
| `q'q` | 5 of 1000 | 9.922e+07 | 5.610e+07 | 1.77 | 998 of 999 |
| `LM` | 76 of 1000 | 3.785e+03 | 2.014e+03 | 1.88 | 957 of 999 |

Both put the benchmark's own configuration well behind configurations that did
not generate the data, by a factor of about 1.8 in mean loss in each case. This
is not a marginal miss.

### This reverses the reading the real-data run invited

`docs/ABM_SYSTEM_SCORE_LOSS.md` records that on the US data the weighted score
loss was the only one of the three whose confidence set stopped because an
equivalence test was accepted rather than because elimination ran out, and
treats that as the mark in its favour. It does the same here: `LM` stops at
p = 0.0783 in round 994, while both other losses eliminate to a singleton with
a final p-value of 0.0000.

On this benchmark that property is worth nothing. The loss whose set is
"better behaved" in that sense is the one that misses the right answer by 75
places, and the loss whose procedure never accepts a test is the one that lands
on the right answer alone. Stopping via an accepted test says the bootstrap
found two configurations it could not separate; it says nothing about whether
either is correct.

### Why the score losses might fail here

The following is a hypothesis consistent with the numbers, not something this
experiment establishes.

The score is evaluated at one replicate's estimate. `theta_benchmark` is the
maximiser of replicate 706's own likelihood, so it carries that replicate's
sampling noise. Another replicate of `cop_0191` is an independent draw from the
same process, and its score at `theta_benchmark` reflects the full
replicate-to-replicate variability of the estimator, which is not small: the
benchmark configuration's own mean `LM` is 3785 against a chi-square reference
of 42. The auxiliary model at one replicate's estimate does not describe
another replicate of the same configuration.

That leaves room for a configuration that did not generate the data to score
lower, if its series happen to push the estimate less from `theta_benchmark`
than `cop_0191`'s own sampling spread does. The score loss then measures
something closer to "how flat is the likelihood at this point under this
configuration" alongside "how close is this configuration", and the two are not
the same question. The impulse-response loss compares two estimates, both
noisy in the same way, and is minimised when the implied dynamics agree.

The signal is there, for what it is worth: 3785 for the benchmark's own
configuration against a field median of 17980, so being the generating
configuration lowers the score fivefold. It is simply not enough to win, and
six configurations beat it.

### The three losses rank the field differently here too

Spearman rank correlation between per-configuration mean losses:

| pair | rho, Monte Carlo benchmark | rho, US benchmark |
| --- | --- | --- |
| `LM` vs `q'q` | 0.9192 | 0.9575 |
| `LM` vs IRF | 0.7472 | 0.4919 |
| `q'q` vs IRF | 0.6285 | 0.4773 |

The two score losses agree with each other, as they did against the US data.
They agree with the impulse-response distance noticeably better here than they
did there, and still not enough to pick the same configuration.

## Every replicate as the benchmark

One benchmark cannot separate a protocol that identifies a configuration from
one that drew well once. `montecarlo/sweep.sh` repeats the whole experiment
with each of `cop_0191`'s 1000 replicates standing in as the benchmark, and
reports how often each loss returns `cop_0191`.

    ./montecarlo/sweep.sh

Every individual run is what `run.sh` does for its single benchmark: same
losses, same hold-out of the benchmark's own seed, same confidence set
settings, MCS_TR at 10000 resamples. Only the loop and the tally are added.

### What makes it affordable, and what does not

Measured on 16 cores of this machine: one pass of $10^6$ filter evaluations
plus the 15 GB dataset read takes 46 s; one pass of $10^6$ cached fits and
impulse responses takes 574 s; one MCS_TR over 1000 models and 999
observations takes 10.9 s.

The impulse-response protocol has a shortcut. A cell's response vector does not
depend on which benchmark it is compared against, only the comparison does, so
`montecarlo/sweep_irf.c` computes all $10^6$ of them once, holds them as
float32 (2.1 GB), and each benchmark's loss matrix is then one pass of mean
absolute error over the cache. Without that the sweep would be 160 hours of
recomputing the same responses; with it the whole protocol is about 3 hours,
almost all of it the thousand confidence sets.

The score protocols have no such shortcut. The score is evaluated at the
benchmark's own estimate, so every cell must be filtered again for every
benchmark: $10^9$ evaluations. `montecarlo/sweep_score.c` amortises the
archives by holding a block of benchmarks at once and scoring each cell against
all of them while its series is in hand, which leaves the filter evaluations
themselves as the cost. About 12 hours.

Both are resumable at benchmark granularity. Each finished benchmark is
appended to its summary as it completes and skipped on a rerun, so an
interrupted sweep continues where it stopped and a partial summary is already
readable.

### Outputs

| file | what it holds |
| --- | --- |
| `montecarlo/out/sweep_irf.csv` | one row per benchmark, impulse-response protocol |
| `montecarlo/out/sweep_score.csv` | the same under `q'q` |
| `montecarlo/out/sweep_score_weighted.csv` | the same under `LM` |

Each row records whether `cop_0191` ended up in the confidence set, where it
ranked by mean loss, its MCS p-value, the size of the set, whether the
procedure stopped on an accepted test, and how many replicates were dropped.

Every one of `cop_0191`'s 1000 runs is a benchmark, whether or not its own fit
converged (449 did), unlike the single benchmark above. A benchmark whose
information matrix is not positive definite has no `LM` loss; it is written to
`sweep_score_weighted.csv` with `skipped` and scored under `q'q` only.

### Results

**Impulse-response loss**, complete, 1000 benchmarks, run 2026-09-20.
`cop_0191` is in the confidence set for 876 of them and has the smallest mean
loss for 829. No run was dropped for a missing response.

**Score losses**, incomplete. `bin/sweep_score` stopped after benchmarks 0 to
388, 389 of 1000, on 2026-09-21 and has not been resumed; a rerun continues from
benchmark 389. On those 389:

- `q'q`: `cop_0191` in the set for 182, smallest mean loss for 14. Its median
  rank by mean loss is 794 of 1000; the median set holds 23 configurations, the
  mean 174, the largest all 1000.
- `LM`: skipped for 252 of 389 because the benchmark's information matrix is
  not positive definite (202 of the 211 unconverged benchmarks, 50 of the 178
  converged). Of the 137 scored, `cop_0191` is in the set for 5, its median rank
  is 145 and it is never first (best rank 8); the median set holds 12.

So the weighted score loss cannot be built at most runs of `cop_0191`, and where
it can, it almost never returns `cop_0191`. The single-benchmark result above
(`LM` misses, rank 76) is the typical outcome, not an unlucky draw. The `q'q`
result is different from the single benchmark: on that one it missed, here it
keeps `cop_0191` about half the time, but through large sets in which it ranks
near the bottom, not by ranking it first.

The impulse-response result split by the benchmark's own fit, with convergence
read from `out/abm_system_fit_qvarma_manifest.txt` (unchanged since
2026-09-10, before the sweep):

| benchmark's fit | benchmarks | `cop_0191` in set | smallest mean loss |
| --- | --- | --- | --- |
| converged | 449 | 415 (92.4%) | 400 |
| not converged | 551 | 461 (83.7%) | 429 |

## The first run of every configuration as the benchmark

The sweep above asks whether the procedure identifies `cop_0191`. This asks
whether it identifies every configuration: for each of the 1000
configurations, its replicate 0 (seed 1 of the model) stands in as the
benchmark, and the confidence set should return that configuration. One
benchmark per configuration, 1000 confidence sets per model.

    make bin/sweep_cops && ./bin/sweep_cops

It runs the impulse-response loss only, under the t-QVARMA and the four local
projections of `docs/MONTECARLO_LP_VALIDATION.md` (linear, state 1, state 2,
both states).
Every individual run is what `montecarlo/sweep_irf.c` and
`montecarlo/lp_sweep.c` do for one benchmark: mean absolute error between
response vectors, replicate 0 held out of every column, a replicate with a
missing cell in any column dropped, MCS_TR at alpha 0.05 with 10000 resamples,
block length 1, bootstrap variance, seed 123 stream 0. The row for `cop_0191`
therefore reproduces row 0 of `sweep_irf.csv` and of `lp_sweep_<model>.csv`.

The whole result is one file, `montecarlo/out/sweep_cops.csv.gz`, one row per
model and benchmark configuration: whether the configuration is in the set,
its rank by mean loss, its MCS p-value, the set size, whether the set was
decided by an accepted test, the final p-value, replicates dropped, the
configuration with the smallest mean loss, the seconds spent in the confidence
set, and the configurations in the set. While it runs, rows go to
`montecarlo/out/sweep_cops_progress.csv`, which a rerun resumes from and which
is compressed and removed when all five models are done. A model added to
the program later is run on a rerun: the finished file is unpacked back into
the progress file, the new model's rows are appended, and the old rows are kept
as they were.

### Results

The first four models were run on 2026-09-28 and 29, one confidence set at a
time on 16 cores; one confidence set took 4.8 to 5.2 seconds on average per
model, 14.4 at most. Both states were added on 2026-10-05, 16 confidence sets
side by side with one thread each, 11.2 seconds per set on average and 13.4 at
most. No response was missing under any model, so every loss matrix is 999
replicates by 1000 configurations with nothing dropped. The rows of the first
four models are byte for byte those written before both states were added.

The `cop_0191` rows reproduce row 0 of `sweep_irf.csv` and of the four
`lp_sweep_<model>.csv` exactly, including the state-1 and state-2 sets of 3 and
4 configurations decided at p = 0.0777 and 0.4187.

Out of 1000 benchmark configurations:

| model | own configuration in the set | alone in the set | smallest mean loss | mean set size |
| --- | --- | --- | --- | --- |
| t-QVARMA | 97 | 34 | 51 | 2.99 |
| LP linear | 93 | 36 | 57 | 3.09 |
| LP state 1 | 10 | 2 | 5 | 5.27 |
| LP state 2 | 24 | 4 | 10 | 6.61 |
| LP both states | 14 | 5 | 6 | 6.76 |

The 87.6% at which the t-QVARMA recovers `cop_0191` from its own replicates
does not carry to the other configurations: with one benchmark each, the
procedure returns the right configuration for fewer than one in ten of them.

A few configurations take the smallest mean loss against many benchmarks
that are not theirs: under state 1, `cop_0437` for 139 benchmarks and
`cop_0191` for 102; under state 2, `cop_0717` for 92; under the t-QVARMA,
`cop_0505` for 47 and `cop_0799` for 38; under both states, `cop_0437` for
134, `cop_0717` for 110 and `cop_0410` for 95. `cop_0191` is in the set of 16
benchmarks other than itself under the t-QVARMA, 34 under the linear local
projection, 175 under state 1, 118 under state 2 and 143 under both states.

Replicate 0's t-QVARMA fit converged for 337 of the 1000 configurations. The
own configuration is in the set for 26 of those 337 and 71 of the other 663.

## Every run of every configuration as the benchmark

The experiment above uses one benchmark per configuration, so a configuration's
result is one draw. This uses all of them: every run of every configuration is
the benchmark in turn, 1,000,000 benchmarks per model, five million confidence
sets in all. It answers, for each configuration, how often the procedure
returns it when it is the right answer.

    make mc-sweep_grid             build at -O3 and run
    ./bin/sweep_grid               the same, once built

`sweep_grid` is not in `MONTECARLO_STEMS`, so `make montecarlo` does not build
it: it has its own rule in the Makefile because it is compiled at -O3 rather
than -O2. The Makefile records that -O3 changes no result, measured on 8 sets
side by side on 999 by 1000 tables at 2000 resamples: 0.13 s per set at -O2,
0.08 s at -O3.

Each benchmark is what `montecarlo/sweep_cops.c` does for one, with one setting
changed: 2000 bootstrap resamples instead of 10000. The loss is the mean
absolute difference between response vectors, the benchmark's run index is held
out of every column, a run with a missing response in any column is dropped,
MCS_TR at alpha 0.05, block length 1, bootstrap variance, seed 123 stream 0.

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
  fits and the dataset are read from, default `out/abm_system_fit_lp` and
  `dataset/abm_system`.

Compile-time overrides for a test build over part of the grid: `BOOTSTRAP`,
`FIRST_MODEL`, `RUN_LIMIT`, `BENCHMARK_LIMIT`, `PROGRESS_PATH`, `RESULT_PATH`,
`RESPONSE_CACHE_DIR`.

It reads `out/abm_system_fit_qvarma/`, `out/abm_system_fit_lp/` and
`dataset/abm_system/` and rebuilds none of them. Each model's responses are kept
in `out/sweep_grid_response_cache/<model>.f32` (2.1 GB for the t-QVARMA, 1.6 GB
per local projection, 3.2 GB for both states, ignored by git), built on the
first run in about 10 minutes for the t-QVARMA and 20 per local projection, then
read back in seconds. Both states' file is built from the two state files,
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
`montecarlo/out/sweep_grid_progress.csv` as each batch of 64 finishes; a rerun
skips what that file holds. When every model is done the rows are compressed to
`montecarlo/out/sweep_grid.csv.gz` and the progress file is removed. If
`sweep_grid.csv.gz` exists and no progress file does, the finished rows are
unpacked back into the progress file and only the models missing from them are
run; if none is missing, the program does nothing. A crash can leave the
progress file ending in zero bytes, the part of the last write that never
reached the disk; `montecarlo/lp_nl_run.sh` cuts the file at the first of them
before resuming.

Cost: the first four models spanned several sessions, the last finished on
2026-10-04. That last session did 1,187,552 benchmarks (the rest of state 1 and
all of state 2) at 0.0956 s of wall-clock time per benchmark, 8 workers on this
machine (AMD Ryzen 7 4800H, 16 threads), from `out/sweep_grid.log`. The earlier
sessions' timing was not kept. Both states ran from 2026-10-05 to 2026-10-07.
A crash of the machine stopped them after 471,768 benchmarks; the session after
it did the remaining 528,232 at 0.1196 s per benchmark, 16 workers, from
`out/lp_nl_run.log`, comparing 800 entries per run rather than 400. The timing
before the crash was not kept. After both
states were added, the t-QVARMA and linear rows of the new
`sweep_grid.csv.gz` were compared with `sweep_grid_qvarma.csv.gz` and
`sweep_grid_lp_lin.csv.gz` and are identical; the state 1 and state 2 rows have
no separate copy to compare with.

### Outputs

| file | what it holds |
| --- | --- |
| `montecarlo/out/sweep_grid.csv.gz` | one row per model and benchmark: model, benchmark configuration and run, whether the configuration is in the set, its rank by mean loss, its MCS p-value, the set size, whether the set was decided by an accepted test, the final p-value, runs dropped, the configuration with the smallest mean loss, the configurations in the set |
| `montecarlo/out/sweep_grid_qvarma.csv.gz`, `sweep_grid_lp_lin.csv.gz` | written by an earlier version of `sweep_grid.c` that split the rows per model; byte for byte the t-QVARMA and linear rows of `sweep_grid.csv.gz`, header included. The current program does not write them |
| `montecarlo/out/sweep_grid_cop_0191.txt`, `.tex` | `montecarlo/sweep_grid_cop_0191.py`, in plain text and as LaTeX tables: per model, in percent, how often `cop_0191` is in the set, alone and first by mean loss over its own 1000 benchmarks and over the other 999,000, and how often the true configuration is over all 1,000,000; `docs/MONTECARLO_COP_0191.md` |
| `montecarlo/out/sweep_grid_presence.txt` | `montecarlo/sweep_grid_presence.py`: per model, how often each configuration is in the sets of benchmarks that are not its own |
| `montecarlo/out/sweep_grid_identifiability.csv`, `sweep_grid_identifiability_plots/` | `montecarlo/sweep_grid_identifiability_plots.py`: per model and configuration, the share of its own 1000 sets that hold it, and the share holding the other configuration found in them most often, with two bar charts per model |

    python montecarlo/sweep_grid_presence.py
    python montecarlo/sweep_grid_identifiability_plots.py

The first needs only the Python standard library. The second needs `polars`,
`plotly` and `kaleido`. The figures and the table on disk were drawn on
2026-10-02, when only the t-QVARMA and the linear local projection were done, so
they cover those two models; rerunning the script adds the three other local
projections.

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

Counts out of 1,000,000 benchmarks per model. The rates are close to those of
"The first run of every configuration" (97, 93, 10, 24 and 14 out of 1000), so
that experiment's result was not an unlucky choice of run.

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
`docs/MONTECARLO_COMPRESSED_RESPONSE.md` traces this to the loss averaging over
hundreds of response entries that are mostly noise.

**2000 against 10000 resamples.** Run 0 of every configuration is a benchmark
both here (2000 resamples) and in `sweep_cops.csv.gz` (10000), on the same loss
matrices. Over those 1000 benchmarks per model, whether the own configuration is
in the set differs once in 5000 (linear local projection, 93 against 94); the
set size differs for 51 (t-QVARMA), 27 (linear), 68 (state 1), 99 (state 2) and
94 (both states) of 1000. The lower resample count changes the set membership of
the right answer essentially never.

### Recovery by whether the benchmark's fit converged

Only the t-QVARMA has a convergence status; the local projections are least
squares. Convergence read from `out/abm_system_fit_qvarma_manifest.txt`, whose
fits are unchanged since 2026-09-10, before every experiment here.

| experiment | converged benchmarks: own configuration in set | unconverged benchmarks: own configuration in set |
| --- | --- | --- |
| every run of `cop_0191` | 415 of 449 (92.4%) | 461 of 551 (83.7%) |
| run 0 of every configuration | 26 of 337 (7.7%) | 71 of 663 (10.7%) |
| every run of every configuration | 35,764 of 352,470 (10.1%) | 59,387 of 647,530 (9.2%) |

A converged benchmark helps `cop_0191` by about nine points and makes almost no
difference over the whole grid. Using unconverged benchmarks is not what keeps
the grid's recovery near one in ten.

## What this experiment cannot show

The impulse-response protocols have now been asked about every configuration
("Every run of every configuration as the benchmark"), and the answer for
`cop_0191` does not generalise. `cop_0191` is recovered from 876 of its own 1000
runs under the t-QVARMA, which only 4 configurations exceed, and from more of
its own runs than any other configuration under each of the four local
projections (913, 857, 780 and 916 of 1000 under linear, state 1, state 2 and
both states). Over all configurations the rate is
about one in ten. Results in this project that rest on `cop_0191` alone, the
single benchmark and the sweep over its runs, describe one of the easiest
configurations to identify.

That settles the confound the single benchmark raised. `cop_0191` was the
impulse-response protocol's winner on the US data and was then chosen as the
benchmark on convergence; the grid shows its impulse responses are unusually
distinctive, which explains why it is recovered here without the protocol being
good in general. The configurations the score protocols favoured on the US data
are hard to recover under the impulse-response protocols: `cop_0409` from 45 of
its 1000 runs under the t-QVARMA and at most 1 under the local projections,
`cop_0931` from none under the t-QVARMA and from 27 to 238 under the local
projections (95 under both states).

The score protocols have been asked only about `cop_0191`, and only on 389 of
its runs. Whether they recover any other configuration has not been run; at
$10^9$ filter evaluations per thousand benchmarks it is the expensive half.

On its own it cannot separate the loss from the auxiliary model: a protocol that
fails to recover a configuration may use a poor loss, or its auxiliary model may
not distinguish these configurations' dynamics at all.
`docs/MONTECARLO_COMPRESSED_RESPONSE.md` separates the two for the
impulse-response protocols. Keeping the same responses and changing only the
vector the loss compares raises recovery over the grid from about one in ten to
80% to 95%, so most of the failure is the loss.

It inherits every limitation of the protocols it runs. In particular the fits
it reads are the main pipeline's, about a third of which converged, and
`docs/ABM_SYSTEM_MCS_VALIDATION.md` records the measurement that the
unconverged ones' log-likelihoods are not distinguishable from the converged
ones'. The single benchmark is converged; the thousand-by-thousand field it is
compared against is not, and neither are most benchmarks of the later
experiments.
