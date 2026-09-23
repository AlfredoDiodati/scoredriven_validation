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
union of the three procedures' surviving configurations is the candidate pool.
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
manifest marks as inherited from another replicate are excluded.

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

### Results

To be filled in when the sweep completes.

## What this experiment cannot show

The sweep above answers the "one benchmark" objection for the benchmark
configuration, and only for it. A protocol that recovers `cop_0191` from every
one of its own replicates has been shown to identify `cop_0191`; it has not
been shown to identify any other configuration. Drawing benchmarks from
`cop_0409` or `cop_0931`, which the score protocols favoured on the US data,
is the obvious next step and has not been run.

Before the sweep, the results section above rested on a single benchmark, and
that is worth keeping in mind when reading anything in this project written
against that single run.

The benchmark's configuration was the impulse-response protocol's own winner on
the US data, and that is a confound the experiment cannot rule out. `cop_0191`
was selected on convergence, not on which procedure had favoured it, and the US
data plays no part in any calculation here. But if `cop_0191` has unusually
distinctive impulse responses, that would explain both why the
impulse-response protocol picked it against the US data and why it recovers it
here, without the protocol being better in general. Repeating the experiment
with a benchmark drawn from `cop_0409` or `cop_0931`, which the score protocols
favoured, is the obvious way to separate the two and has not been done.

It cannot separate the procedure from the auxiliary model. If a protocol fails
to recover `cop_0191`, that could mean the loss is a poor discriminator or that
the t-QVARMA does not distinguish these configurations' dynamics at all.
Nothing here tells the two apart.

It inherits every limitation of the protocols it runs. In particular the fits
it reads are the main pipeline's, about a third of which converged, and
`docs/ABM_SYSTEM_MCS_VALIDATION.md` records the measurement that the
unconverged ones' log-likelihoods are not distinguishable from the converged
ones'. The benchmark itself is converged; the thousand-by-thousand field it is
compared against is not.
