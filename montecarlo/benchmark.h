#ifndef MONTECARLO_BENCHMARK_H
#define MONTECARLO_BENCHMARK_H

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*
The simulated run promoted to the role the US data plays in the main pipeline,
as montecarlo/benchmark_choice.c chose it.

The choice is made once, by that script, and written to
montecarlo/out/benchmark.env. Every other step of one experiment reads it back
from there rather than being told which run to use, so the loss matrices and
the confidence sets cannot end up describing different benchmarks, and the
grounds for the choice stay in one place - montecarlo/out/benchmark_choice.txt,
which the same script writes beside it.

replicate is also the index held out of every configuration's column.
Replicate n of every configuration is seed n+1 of the model
(applications/abm_system_simulate.c's "seed = mc + 1", the same seeds for every
configuration), so holding out the row drops every simulation sharing the
benchmark's seed, not only the benchmark itself.
*/

#define BENCHMARK_ENV_PATH "montecarlo/out/benchmark.env"

typedef struct {
    char sample[64];
    int replicate;
} Benchmark;

static inline Benchmark benchmark_read(void) {
    FILE *f = fopen(BENCHMARK_ENV_PATH, "r");
    assert(f && "montecarlo: " BENCHMARK_ENV_PATH " is missing - bin/benchmark_choice writes it");

    Benchmark b;
    b.sample[0] = 0;
    b.replicate = -1;

    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "ABM_SYSTEM_BENCHMARK_SAMPLE=%63s", b.sample) == 1) continue;
        if (sscanf(line, "ABM_SYSTEM_BENCHMARK_REPLICATE=%d", &b.replicate) == 1) continue;
    }
    fclose(f);

    assert(b.sample[0] && "montecarlo: no benchmark configuration in " BENCHMARK_ENV_PATH);
    assert(b.replicate >= 0 && "montecarlo: no benchmark replicate in " BENCHMARK_ENV_PATH);
    return b;
}

#endif /* MONTECARLO_BENCHMARK_H */
