# Benchmarks and the throughput study

How long the t-QVARMA fit takes, where its time goes, and how it scales. None
of this gates anything: the programs under `benchmarks/` are never part of
`make test`, and the throughput study under `applications/` is run by hand.

All of it was measured in late August and September 2026 on the older dataset
of 100 configurations of 108 replications each (named
`EstimationSeriesSample1_N` in the reports), before the design dataset existed,
and against the et_al of that time. Both have changed since: et_al corrected
the t-QVARMA log-likelihood at large `nu` on 2026-08-31, and the pipeline's
iteration cap was raised from 2000 to 1,000,000
(`docs/MCS_VALIDATION_PROCEDURE.md`, "What the t-QVARMA losses are built
from"). Read the numbers below as a record of those runs, not as the current
cost of a fit.

## The benchmarks

    make bench                   build and run all seven
    make bench-<stem>            one of them
    make bench-performance ETAL_DEV=1
                                 the same, built against the development copy
                                 of et_al instead of the installed one

Each writes `out/<stem>.txt`. `out/fit_speedup_options.txt` collects the
measurements of 2026-08-29 and the setup behind them; it states that they were
taken on an Intel i5-7400 (4 cores, no hyperthreading, 3.0 GHz base), not on
the machine the README describes. The reports themselves do not record the
machine.

| benchmark | what it measures | result in its report |
|---|---|---|
| `qvarma_fit_cost` | one objective evaluation with and without the gradient, and one full fit at 2000 iterations, single threaded, both specifications | (1,1,2): 0.18 ms with gradient, 0.087 ms without, 0.21 s for the fit, not converged at 2000 iterations |
| `qvarma_fit_io` | the file handling around one fit: reading the replicate, writing the fitted parameters as JSON, reading them back | a fit 219 ms, reading the replicate 3.0 ms, writing the JSON 0.17 ms, reading it 0.09 ms |
| `qvarma_iteration_budget` | for 4 replicates and both specifications, the iteration at which the log-likelihood first comes within 1, 0.1 and 0.01 of its value at iteration 2000 | within 1 at 1228 to 1956 iterations; within 0.01 only after 1823 or more, and not at all for two (1,1,4) fits |
| `qvarma_taped_vs_fused` | et_al's analytic log-likelihood against the automatic-differentiation path it replaced, after checking that both give the same value (identical) and gradient (worst relative gap 6.1e-13) | the analytic path 7.9 times faster on the value alone, 5.6 times with the gradient, single threaded |
| `qvarma_thread_scaling` | 8 fits at 1, 2 and 4 OpenMP threads | speedup 1.83 at 2 threads, 3.40 at 4 |
| `qvarma_process_scaling` | the same work split over 1, 2 and 4 separate processes, one thread each | speedup 1.97 at 2 processes, 3.14 at 4 |
| `fit_contention_source` | which shared resource parallel fits contend on, away from et_al: small `cblas_dgemm` calls, `malloc` and `free` of 64 KiB, and plain arithmetic, at 1 and 4 threads | 4 threads: `cblas_dgemm` 0.45 times the single-thread throughput, `malloc` 3.59 times, arithmetic 3.75 times; OpenBLAS's small products are what does not scale |

## The throughput study

How long a battery of 500,000 fits takes, the scale an indirect-inference loop
over the ABM would need, and what a larger iteration cap buys.

- `applications/throughput_dataset.c` builds the data, `make
  app-throughput_dataset`: 500 folders of 1000 series under
  `out/dataset/throughput/`, about 21 GB, not tracked. Each of the older dataset's
  100 files fills five folders: its 108 replications once, then perturbed
  copies of them, each with a level shift per series and noise per period, both
  normal with standard deviation 0.01, seed 20260830. The perturbation only
  makes the 500,000 fits start from different data.
- `applications/throughput_fit.c` fits all 500,000 at a fixed iteration cap,
  `make app-throughput_fit-i2000`, `-i4000`, `-i8000`. The cap is in the
  binary's name and in every file it writes, `out/throughput_fit_i<cap>/`
  (not tracked) and `out/throughput_fit_i<cap>_timing.txt`,
  `_provenance.txt` and `_manifest.txt`.
- `applications/throughput_iteration_budget.c` compares two caps over the same
  series, `make app-throughput_iteration_budget-i2000_4000` and
  `-i4000_8000`, writing `out/throughput_iteration_budget_i<a>_i<b>.txt` and
  `.csv`.

The specification is the pipeline's t-QVARMA(1,1,2), 16 threads. The
provenance files record the host as Linux 6.19.14-arch1-1 with 16 hardware
threads, the machine the README describes.

| cap | wall time | fits per second | converged |
|---:|---:|---:|---:|
| 2000 | 4.13 h | 33.6 | 65,764 of 500,000 |
| 4000 | 7.47 h | 18.6 | 246,743 of 500,000 |
| 8000 | 3.67 h, 203,663 fits estimated and 296,337 read back from the 4000 run's cache | not comparable | 366,236 of 500,000 |

From 2000 to 4000 iterations, over the 499,989 series with both fits: the
log-likelihood rose by 0.53 at the median and 1.45 on average, and not at all
for 17.1%; the parameters other than `nu` moved by a median of 0.5% of their
norm and by 29% at the 90th percentile. The comparison report counts 217,616
fits converged at 4000, where the timing report counts 246,743; the two count
convergence differently, and the difference has not been traced.

From 4000 to 8000, over 499,988 series: the mean gain was 0.41, 55% gained
nothing, and 119,489 more fits converged, 73.25% in all.

The provenance file of the 2000-iteration run was written afterwards and
records that it ran against an et_al that predates a correction of the
t-QVARMA log-likelihood at large `nu`, so its fits with a large `nu` are not
on the corrected likelihood.
