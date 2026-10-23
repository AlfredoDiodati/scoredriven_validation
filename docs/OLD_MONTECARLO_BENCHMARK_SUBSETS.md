# Monte Carlo experiments on subsets of the benchmarks, discarded

## Status: discarded

These three experiments asked the Monte Carlo question of
`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md` on fewer benchmarks than the grid does:

| experiment | benchmarks | code |
| --- | --- | --- |
| one benchmark | `cop_0191` run 706 | `old/montecarlo_one_benchmark/` |
| every run of `cop_0191` | the 1000 runs of `cop_0191` | `old/montecarlo_every_run_of_cop_0191/` |
| run 0 of every configuration | run 0 of each of the 1000 configurations | `old/montecarlo_run_0_of_every_configuration/` |

Each set of benchmarks is part of the grid, `old/naive_mae_loss/montecarlo/sweep_grid.c`, which
makes every run of every configuration the benchmark in turn with the same
loss, the same hold-out and the same confidence set settings, except for the
number of bootstrap resamples: 10000 here, 2000 in the grid. The grid answers
each of these questions as a special case, so they are out of the pipeline.

The resample count is the one difference, and it was measured: run 0 of every
configuration is a benchmark both in the grid and in the third experiment, on
the same loss matrices, and whether the benchmark's own configuration is in
the set differs in 1 of the 5000 benchmarks across the five auxiliary models
(`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md`, "2000 against 10000 resamples"). The
benchmarks of the other two experiments are rows of the grid as well, and were
compared on 2026-10-11 by counting, in `out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz`
against the old outputs, the benchmarks whose set holds `cop_0191`:

| model | grid, 2000 resamples | every run of `cop_0191`, 10000 resamples |
| --- | --- | --- |
| t-QVARMA | 876 of 1000 | 876 (`sweep_irf.csv`) |
| LP linear | 913 | 913 (`lp_sweep_lin.csv`) |
| LP state 1 | 857 | 858 (`lp_sweep_s1.csv`) |
| LP state 2 | 780 | 781 (`lp_sweep_s2.csv`) |
| LP both states | 916 | 914 (`lp_sweep_nl.csv`) |

For the single benchmark, run 706 of `cop_0191`, the grid's set holds
`cop_0191` under all five models, as the first experiment's did
(`mcs.csv`, `lp_mcs_<model>.csv`); under the t-QVARMA it is `cop_0191` alone in
both.

The single-benchmark experiment also ran the two score losses, which are
discarded for their own reason, `docs/OLD_SCORE_LOSS.md`: the unweighted one,
written `q'q`, the sum of squares of the score of the t-QVARMA log-likelihood
at the benchmark's estimate, and the weighted one, written `LM`, the same score
weighted by the inverse information matrix (the Lagrange multiplier
statistic). Both labels are used in the tables below. The impulse-response
loss is written IRF. `MCS_TR` and `MCS_TMAX` are the two statistics of the
confidence set, every pair of configurations compared and each configuration
against the average of the rest (`docs/MCS_VALIDATION_PROCEDURE.md`,
"Confidence set settings"). The two studies in
`old/montecarlo_one_benchmark/studies/` diagnosed this experiment's results and
moved with it. `docs/OLD_EXPERIMENTS.md` lists every discarded experiment.

The outputs are under `out/old/montecarlo_one_benchmark/`, its studies' under
`out/old/montecarlo_one_benchmark/studies/`. Every program builds
and runs with `make old-<stem>`, the stem being its file name, and the three
drivers with `make old-montecarlo_one_benchmark`, `make
old-montecarlo_one_benchmark_lp` and `make
old-montecarlo_every_run_of_cop_0191`.

## Which run stands in for the real data

Choosing the benchmark on grounds connected to the answer would rig the
experiment, so the grounds are fixed in advance and applied by
`old/montecarlo_one_benchmark/benchmark_choice.c`, which writes its reasoning to
`out/old/montecarlo_one_benchmark/benchmark_choice.txt`. Three requirements, in
order.

**It must come from a configuration the real-data confidence set kept.** The
configurations surviving the impulse-response procedure under MCS_TR (the
`tr_in_set` column of
`out/old/naive_mae_loss/applications/abm_system_mcs_statistic_comparison.csv`;
the MCS_TMAX set is not read) are the candidate pool, which is `cop_0191` alone.
Promoting a configuration the procedure did not keep would ask whether a poor
benchmark is recovered, which is a different question.

**Among those, the configuration whose fits converged most often.** About a
third of the million fits in this project converged, which makes convergence
the property in shortest supply. A benchmark whose own fit did not converge is
a point the optimizer stopped at, not an estimate, and the whole experiment
rests on the benchmark being a real estimate.

**Within it, the converged replicate with the smallest gradient norm.** Fits the
manifest marks as inherited from another replicate are excluded; the current
manifest (`out/applications/abm_system_fit_qvarma_manifest.txt`, 2026-09-10) has
none, its p1q1r2 rows being 352,470 converged and reused, 621,485 unconverged
and held, 26,045 unconverged and resumed.

These three requirements apply to this one benchmark only. The later
experiments use every run of a configuration as a benchmark with no filter, so
most of their benchmarks are unconverged fits; "Recovery by whether the
benchmark's fit converged" below measures what that changes.

### What it selected

| step | result |
| --- | --- |
| candidate pool | `cop_0191` |
| convergence rate | 44.9% |
| configuration chosen | `cop_0191` |
| replicate chosen | 706 |

`cop_0191` replicate 706: log-likelihood 124.6637, gradient norm 4.058e-03,
converged, its own fit rather than an inherited one.

When the score protocols were part of this experiment, the pool was the union
of the sets kept under all three real-data losses, `cop_0191`, `cop_0931`,
`cop_0409` and `cop_0599`, and the replicate also had to have a positive
definite information matrix for the weighted score loss. Both rules gave the
same benchmark, `cop_0191` replicate 706, because `cop_0191` had the highest
convergence rate of the four (44.9% against 44.1%, 22.9% and 33.2%) and
replicate 706 passed the information-matrix check on the first try.

For scale, the US fit that the main pipeline uses converged with gradient norm
4.713e-03, so the benchmark is an estimate of comparable quality to the one it
replaces.

`cop_0191` is the configuration the main pipeline's confidence set keeps on
the US data, which is why it is the only candidate under the first requirement.
That does not help the experiment towards the answer: the question is whether
the procedure recovers `cop_0191` from one of its own simulated runs, and the
US data plays no part in that.

## Results

### Setup

Benchmark `cop_0191` replicate 706, 400 periods, its own cached p1q1r2 fit.
Field of 1000 configurations by 999 replicates, replicate 706 held out of every
column, no missing cells in any of the three matrices. Confidence set settings
identical to the main pipeline's and to each other: alpha = 0.05, 10000
bootstrap resamples, block length 1, bootstrap variance of the resampled mean,
seed 123, stream 0, both statistics scored on the same resamples. Every fit
read from `out/applications/abm_system_fit_qvarma/`; nothing estimated.

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
not generate the data, by a factor of about 1.8 in mean loss in each case.

### This reverses the reading the real-data run invited

`docs/OLD_SCORE_LOSS.md` records that on the US data the weighted score
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

The score does carry the right configuration: 3785 for the benchmark's own
configuration against a field median of 17980, so being the generating
configuration lowers it fivefold. That is not enough to win: 75 configurations
have a lower mean `LM` loss.

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
one that drew well once. `old/montecarlo_every_run_of_cop_0191/sweep.sh`
repeats the impulse-response experiment with each of `cop_0191`'s 1000
replicates standing in as the benchmark, and reports how often it returns
`cop_0191`. `old/score_loss/montecarlo/sweep_score.c`, run with `make
old-sweep_score`, does the same for the two discarded score losses.

    ./old/montecarlo_every_run_of_cop_0191/sweep.sh

Every individual run is what `old/montecarlo_one_benchmark/run.sh` does for its
single benchmark: same losses, same hold-out of the benchmark's own seed, same
confidence set settings, MCS_TR at 10000 resamples. Only the loop and the tally
are added.

### What makes it affordable, and what does not

Measured with 16 threads on the machine the README describes (8 cores, 2
threads each): one pass of $10^6$ filter evaluations
plus the 15 GB dataset read takes 46 s; one pass of $10^6$ cached fits and
impulse responses takes 574 s; one MCS_TR over 1000 models and 999
observations takes 10.9 s.

The impulse-response protocol has a shortcut. A cell's response vector does not
depend on which benchmark it is compared against, only the comparison does, so
`old/montecarlo_every_run_of_cop_0191/sweep_irf.c` computes all $10^6$ of them
once, holds them as float32 (2.1 GB), and each benchmark's loss matrix is then
one pass of mean absolute error over the cache. Without that the sweep would be
160 hours of recomputing the same responses; with it the whole protocol is
about 3 hours, almost all of it the thousand confidence sets.

The score protocols have no such shortcut. The score is evaluated at the
benchmark's own estimate, so every cell must be filtered again for every
benchmark: $10^9$ evaluations. `old/score_loss/montecarlo/sweep_score.c`
amortises the archives by holding a block of benchmarks at once and scoring
each cell against all of them while its series is in hand, which leaves the
filter evaluations themselves as the cost. About 12 hours.

Both are resumable at benchmark granularity. Each finished benchmark is
appended to its summary as it completes and skipped on a rerun, so an
interrupted sweep continues where it stopped and a partial summary is already
readable.

### Outputs

| file | what it holds |
| --- | --- |
| `out/old/montecarlo_every_run_of_cop_0191/sweep_irf.csv` | one row per benchmark, impulse-response protocol |
| `out/old/score_loss/montecarlo/sweep_score.csv` | the same under `q'q` |
| `out/old/score_loss/montecarlo/sweep_score_weighted.csv` | the same under `LM` |

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
read from `out/applications/abm_system_fit_qvarma_manifest.txt` (unchanged since
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
projections of `docs/OLD_MONTECARLO_LP_VALIDATION.md` (linear, state 1, state
2, both states). Every individual run is what
`old/montecarlo_every_run_of_cop_0191/sweep_irf.c` and
`old/montecarlo_every_run_of_cop_0191/lp_sweep.c` do for one benchmark: mean
absolute error between response vectors, replicate 0 held out of every column,
a replicate with a missing cell in any column dropped, MCS_TR at alpha 0.05
with 10000 resamples, block length 1, bootstrap variance, seed 123 stream 0.
The row for `cop_0191` therefore reproduces row 0 of `sweep_irf.csv` and of
`lp_sweep_<model>.csv`.

The whole result is one file,
`out/old/montecarlo_run_0_of_every_configuration/sweep_cops.csv.gz`, one row per
model and benchmark configuration: whether the configuration is in the set, its
rank by mean loss, its MCS p-value, the set size, whether the set was decided by
an accepted test, the final p-value, replicates dropped, the configuration with
the smallest mean loss, the seconds spent in the confidence set, and the
configurations in the set. While it runs, rows go to
`out/old/montecarlo_run_0_of_every_configuration/sweep_cops_progress.csv`, which
a rerun resumes from and which is compressed and removed when all five models
are done. A model added to the program later is run on a rerun: the finished
file is unpacked back into the progress file, the new model's rows are appended,
and the old rows are kept as they were.

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

## Recovery by whether the benchmark's fit converged

Only the t-QVARMA has a convergence status; the local projections are least
squares. Convergence read from
`out/applications/abm_system_fit_qvarma_manifest.txt`, whose fits are unchanged
since 2026-09-10, before every experiment here. The last row is the grid's, from
`docs/OLD_MONTECARLO_NAIVE_MAE_GRID.md`, for comparison.

| experiment | converged benchmarks: own configuration in set | unconverged benchmarks: own configuration in set |
| --- | --- | --- |
| every run of `cop_0191` | 415 of 449 (92.4%) | 461 of 551 (83.7%) |
| run 0 of every configuration | 26 of 337 (7.7%) | 71 of 663 (10.7%) |
| every run of every configuration | 35,764 of 352,470 (10.1%) | 59,387 of 647,530 (9.2%) |

A converged benchmark helps `cop_0191` by about nine points and makes almost no
difference over the whole grid. Using unconverged benchmarks is not what keeps
the grid's recovery near one in ten.
