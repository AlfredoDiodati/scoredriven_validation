# Discarded experiments

Experiments that are no longer part of the main pipeline or the Monte Carlo
pipeline. Their code is under `old/`, one folder per experiment; their
write-ups stay in `docs/`; their outputs stay where they were written, under
`out/` and `montecarlo/out/`. Every program still builds and runs with
`make old-<stem>`, the stem being its file name, and the Makefile's
`OLD_SOURCES` lists them all.

Each entry says what the experiment was, why it is out, and where its results
are recorded.

## The score losses, `old/score_loss/`

**What it was.** A loss matrix built without fitting the auxiliary model to
the simulated data: the score of the t-QVARMA log-likelihood at the benchmark's
estimate, evaluated on each simulated series, either unweighted (`q'q`) or
weighted by the inverse information matrix, `q' [T_sim I_1]^-1 q`. It ran
against the US data and against the Monte Carlo benchmark.

**Why it is out.** It gave poor results where the right answer is known: on
the single Monte Carlo benchmark the generating configuration came 5th under
`q'q` and 76th under the weighted loss. On the 389 runs of `cop_0191` that
were scored, `q'q` kept `cop_0191` in the set for 182 but ranked it first for
only 14, through large sets, and the weighted loss could not be built for 252
and kept `cop_0191` in 5 of the other 137. And there is no theoretical justification for using either
quantity as a distance between a simulated configuration and the data.

**Where it is recorded.** `docs/OLD_SCORE_LOSS.md`, and
`docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` for the Monte Carlo results.

| folder | programs |
| --- | --- |
| `old/score_loss/applications/` | `abm_system_score_loss.c`, `abm_system_mcs_statistic_comparison_score.c` |
| `old/score_loss/montecarlo/` | `score_loss.c`, `mcs_statistic_comparison_score.c`, `sweep_score.c` |
| `old/score_loss/` | `score_loss_reproduction.sh`, the check that the moved programs still write the matrices on disk byte for byte |

## Monte Carlo experiments on subsets of the benchmarks

Three folders, one per experiment, discarded for the same reason.

**What they were.** The Monte Carlo question, whether the confidence set holds
the configuration a simulated benchmark was drawn from, asked of fewer
benchmarks than the grid does:

| folder | benchmarks |
| --- | --- |
| `old/montecarlo_one_benchmark/` | `cop_0191` run 706, chosen by `benchmark_choice.c`, under the t-QVARMA (`run.sh`) and the local projections (`lp_run.sh`) |
| `old/montecarlo_every_run_of_cop_0191/` | each of the 1000 runs of `cop_0191` (`sweep_irf.c`, `lp_sweep.c`, `sweep.sh`) |
| `old/montecarlo_run_0_of_every_configuration/` | run 0 of each of the 1000 configurations (`sweep_cops.c`) |

**Why they are out.** Each set of benchmarks is part of the grid,
`montecarlo/sweep_grid.c`, which makes every run of every configuration the
benchmark in turn, with the same loss, the same hold-out and the same
confidence set settings. They are special cases of it. The one difference is
the number of bootstrap resamples, 10000 in these and 2000 in the grid, and on
the benchmarks they share the two agree: whether the benchmark's own
configuration is in the set differs in 1 of 5000 benchmarks for run 0 of every
configuration, by at most 2 of 1000 for the runs of `cop_0191`, and not at all
for the single benchmark. `docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md`, "Status:
discarded", has the counts.

`old/montecarlo_one_benchmark/studies/` holds the two studies that diagnosed
these experiments' results, `qvarma_recovery_diagnostics.c` and
`lp_recovery_diagnostics.c`, with the header they share. They read the single
benchmark's choice and loss tables, so they moved with it and are no longer
part of `make study`.

**Where they are recorded.** `docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md` for the
t-QVARMA and the run 0 experiment, `docs/OLD_MONTECARLO_LP_VALIDATION.md` for
the local projections on the single benchmark and on every run of `cop_0191`,
including the finding that decides how the local projections are fitted.

The drivers run with `make old-montecarlo_one_benchmark`,
`make old-montecarlo_one_benchmark_lp` and
`make old-montecarlo_every_run_of_cop_0191`.
