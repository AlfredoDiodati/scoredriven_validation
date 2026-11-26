/*
Builds one model's response cache, every run of every configuration, and runs
no confidence set. The grids build a missing cache themselves before their
first benchmark; this builds it ahead of them, so the cache can be checked
before a grid reads it.

    ./bin/response_cache_build <model label>

The label is one of response_cache.h's, for instance qvarma_growth_states. A
cache already on disk under the expected header is only mapped and counted.
Writes out/montecarlo/response_cache_build/response_cache_build_<label>.txt:
the configurations and runs, the runs with no usable response, and the wall
time.
*/

#include "montecarlo/response_cache.h"
#include <stdio.h>

#define OUTPUT_DIR "out/montecarlo/response_cache_build"

int main(int argc, char **argv) {
    assert(argc == 2 && "response_cache_build: give one model label");
    int model = -1;
    for (int m = 0; m < N_MODELS; m++)
        if (strcmp(argv[1], model_label[m]) == 0) model = m;
    assert(model >= 0 && "response_cache_build: unknown model label");

    openblas_set_num_threads(1);
    omp_set_max_active_levels(1);
    double started = omp_get_wtime();
    int n_samples = 0;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    assert(n_samples > 0 && "response_cache_build: no configurations in the LP fit cache");
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);

    ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, LP_FIT_DIR_DEFAULT,
                                                       LP_INPUT_DIR_DEFAULT, RESPONSE_CACHE_DIR_DEFAULT);

    mkdir("out/montecarlo", 0755);
    mkdir(OUTPUT_DIR, 0755);
    char path[512];
    snprintf(path, sizeof path, "%s/response_cache_build_%s.txt", OUTPUT_DIR, model_label[model]);
    FILE *report = fopen(path, "w");
    assert(report && "response_cache_build: cannot open the report");
    fprintf(report, "model %s, %d entries per run, cache %s/%s.f32\n", model_label[model], cache.dim,
            RESPONSE_CACHE_DIR_DEFAULT, model_label[model]);
    fprintf(report, "%d configurations x %d runs, runs with no usable response %ld\n", n_samples, n_replicates,
            cache.n_missing);
    fprintf(report, "wall time %.0f s on %d threads\n", omp_get_wtime() - started, omp_get_max_threads());
    fclose(report);

    response_cache_free(&cache);
    lp_free_samples(samples, n_samples);
    return 0;
}
