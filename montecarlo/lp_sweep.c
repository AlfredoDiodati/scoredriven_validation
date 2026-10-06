/*
The local-projection protocol run against every replicate of the benchmark
configuration in turn, the local-projection counterpart of
montecarlo/sweep_irf.c.

montecarlo/lp_irf_loss.c and montecarlo/lp_mcs.c ask, for one benchmark,
whether the confidence set returns the configuration it came from. One
benchmark cannot tell a protocol that identifies a configuration from one that
drew well once, so this repeats each run with every replicate of that
configuration standing in as the benchmark, for each of the four models
(linear, state 1, state 2, both states in one vector), and records how often
the answer comes back right.

Each individual run is what lp_irf_loss.c and lp_mcs.c do for its benchmark:
the same MAE over the 400 response entries (800 for both states), the
benchmark's own replicate held out of every column, a replicate with a hole
dropped, MCS_TR at 10000 resamples, block length 1, bootstrap variance.

As in sweep_irf.c, a cell's response vector does not depend on the benchmark,
so every cell's vector is read once and held as float32, and each benchmark's
loss matrix is one pass of absolute differences over that cache, accumulated
in double. The cache is montecarlo/response_cache.h's, the same files
montecarlo/sweep_grid.c keeps: 1.6 GB per model, 3.2 GB for both states. The
models are swept one after another, each cache released before the next is
opened.

LP_SWEEP_WORKERS benchmarks run side by side, one thread each, default one per
hardware thread.

Resumable at benchmark granularity per model: every finished benchmark is
appended to that model's summary as it finishes, and a rerun skips the
benchmarks already there and does not build a cache for a model whose summary
is complete.

The benchmark configuration is the one montecarlo/out/benchmark.env names; the
replicate that file names plays no part, since every replicate takes its turn.

Requires out/abm_system_fit_lp/, dataset/abm_system/ and
montecarlo/out/benchmark.env. Writes montecarlo/out/lp_sweep_<model>.csv, one
row per benchmark. Progress goes to stderr, so the summaries are the result.
*/

#include "montecarlo/response_cache.h"
#include "montecarlo/benchmark.h"
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <math.h>
#include <omp.h>

#define SUMMARY_STEM "montecarlo/out/lp_sweep"

static const char *FIT_DIR;
static const char *INPUT_DIR;
static int n_samples = 0;
static int n_replicates = 0;

/* Which benchmarks a model's summary already holds. */
static int *read_done(const char *path, int *count) {
    int *done = calloc((size_t)n_replicates, sizeof(int));
    assert(done);
    *count = 0;
    FILE *f = fopen(path, "r");
    if (!f) return done;
    char line[512];
    int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        int replicate;
        if (sscanf(line, "%d,", &replicate) == 1 && replicate >= 0 && replicate < n_replicates
            && !done[replicate]) {
            done[replicate] = 1;
            (*count)++;
        }
    }
    fclose(f);
    return done;
}

/* One benchmark's row of the summary. */
typedef struct {
    int in_set, rank, set_size, converged, n_dropped;
    double pvalue, final_pvalue;
} BenchmarkOutcome;

/* The whole protocol for one benchmark, on the calling thread alone: the loss
   matrix against that benchmark, then the confidence set over it. */
static BenchmarkOutcome run_benchmark(const ResponseCache *cache, const char *const *model_name,
                                      int benchmark_sample, int benchmark) {
    const float *reference = response_slot(cache, benchmark_sample, benchmark);
    const int dim = cache->dim;

    Mat values = mat_new(n_replicates - 1, n_samples);
    int *usable = malloc((size_t)(n_replicates - 1) * sizeof(int));
    assert(usable && "lp_sweep: out of memory for a benchmark");
    int row = 0;
    for (int replicate = 0; replicate < n_replicates; replicate++) {
        if (replicate == benchmark) continue;
        usable[row] = 1;
        for (int sample = 0; sample < n_samples; sample++) {
            const float *cell = response_slot(cache, sample, replicate);
            if (MISNAN(cell[0])) {
                AT(values, row, sample) = (mreal)NAN;
                usable[row] = 0;
                continue;
            }
            double sum = 0;
            for (int i = 0; i < dim; i++) sum += fabs((double)cell[i] - (double)reference[i]);
            AT(values, row, sample) = (mreal)(sum / dim);
        }
        row++;
    }

    int keep = 0;
    for (int r = 0; r < n_replicates - 1; r++) keep += usable[r];
    assert(keep > 0 && "lp_sweep: every replicate of a benchmark has a missing cell");

    Mat kept = values;
    if (keep < n_replicates - 1) {
        kept = mat_new(keep, n_samples);
        int at = 0;
        for (int r = 0; r < n_replicates - 1; r++) {
            if (!usable[r]) continue;
            for (int sample = 0; sample < n_samples; sample++) AT(kept, at, sample) = AT(values, r, sample);
            at++;
        }
        mat_free(values);
    }
    DataFrame losses = df_from_matrix(kept, model_name);
    mat_free(kept);
    free(usable);

    MCSOptions opt = mcs_options_default();
    opt.bootstrap = 10000;
    opt.block_length = 1;
    opt.variance = MCS_VARIANCE_BOOTSTRAP;
    opt.stat = MCS_TR;
    MCSResult res = mcs(&losses, opt);

    /* Rank 1 is the smallest mean loss. */
    double benchmark_mean = (double)stats_mean(df_col_numeric(&losses, model_name[benchmark_sample]));
    int rank = 1;
    for (int sample = 0; sample < n_samples; sample++) {
        if (sample == benchmark_sample) continue;
        if ((double)stats_mean(df_col_numeric(&losses, model_name[sample])) < benchmark_mean) rank++;
    }

    BenchmarkOutcome o = { mcs_in_set(&res, benchmark_sample), rank, res.n_surviving, res.converged,
                           n_replicates - 1 - keep, res.pvalue[benchmark_sample], res.final_pvalue };
    mcs_free(&res);
    df_free(&losses);
    return o;
}

static void sweep_model(const LpSample *samples, int benchmark_sample, int model, int n_workers) {
    const char *name = lp_model_name(model);
    char summary_path[256];
    snprintf(summary_path, sizeof summary_path, "%s_%s.csv", SUMMARY_STEM, name);

    int n_done;
    int *done = read_done(summary_path, &n_done);
    fprintf(stderr, "%s: %d of %d benchmarks already done\n", name, n_done, n_replicates);
    if (n_done == n_replicates) { free(done); return; }

    fprintf(stderr, "%s: building the response cache, %.1f GB\n", name,
            (double)n_samples * n_replicates * response_dim(model) * sizeof(float) / 1e9);
    ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, FIT_DIR,
                                                       INPUT_DIR, RESPONSE_CACHE_DIR_DEFAULT);
    fprintf(stderr, "%s: cache built, %ld cells missing\n", name, cache.n_missing);

    FILE *summary = fopen(summary_path, n_done ? "a" : "w");
    assert(summary && "lp_sweep: cannot open a summary for writing");
    if (!n_done) {
        fprintf(summary, "benchmark_replicate,benchmark_in_set,benchmark_rank,benchmark_pvalue,"
                         "set_size,decided_by_accepted_test,final_pvalue,n_dropped_replicates\n");
        fflush(summary);
    }

    const char **model_name = malloc((size_t)n_samples * sizeof(char *));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_lp_%s", samples[i].name, name);
        model_name[i] = name_buffer[i];
    }

    /* Benchmarks are independent, so several can run at once, each confidence
       set on one thread, with nested parallel regions switched off. With one
       worker the loop runs on the calling thread and each confidence set uses
       mcs's own threads instead. Rows reach the summary in the order they
       finish; the benchmark index in the first column identifies them. */
    if (n_workers > 1) omp_set_max_active_levels(1);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_workers) if (n_workers > 1)
    for (int benchmark = 0; benchmark < n_replicates; benchmark++) {
        if (done[benchmark]) continue;
        if (MISNAN(response_slot(&cache, benchmark_sample, benchmark)[0])) {
            fprintf(stderr, "%s benchmark %03d skipped, its own response is missing\n", name, benchmark);
            continue;
        }
        BenchmarkOutcome o = run_benchmark(&cache, model_name, benchmark_sample, benchmark);
        #pragma omp critical(lp_sweep_summary)
        {
            fprintf(summary, "%d,%d,%d,%.6f,%d,%s,%.6f,%d\n", benchmark, o.in_set, o.rank, o.pvalue,
                    o.set_size, o.converged ? "yes" : "no", o.final_pvalue, o.n_dropped);
            fflush(summary);
            fprintf(stderr, "%s benchmark %03d: rank %d, in set %d, set size %d\n", name,
                    benchmark, o.rank, o.in_set, o.set_size);
        }
    }

    fclose(summary);
    free(name_buffer);
    free(model_name);
    response_cache_free(&cache);
    free(done);
}

int main(void) {
    FIT_DIR = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!FIT_DIR) FIT_DIR = LP_FIT_DIR_DEFAULT;
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = LP_INPUT_DIR_DEFAULT;
    openblas_set_num_threads(1);

    Benchmark benchmark = benchmark_read();
    LpSample *samples = lp_list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0 && "lp_sweep: no configurations in the LP fit cache");
    n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);

    int benchmark_sample = -1;
    for (int i = 0; i < n_samples; i++)
        if (strcmp(samples[i].name, benchmark.sample) == 0) benchmark_sample = i;
    assert(benchmark_sample >= 0 && "lp_sweep: the benchmark configuration is not in the LP fit cache");

    /* LP_SWEEP_WORKERS confidence sets side by side, one thread each; by
       default one per hardware thread. Each set at 10000 resamples holds about
       100 MB. */
    const char *workers_text = getenv("LP_SWEEP_WORKERS");
    int n_workers = workers_text ? atoi(workers_text) : omp_get_num_procs();
    assert(n_workers >= 1 && "lp_sweep: LP_SWEEP_WORKERS must be at least 1");

    for (int model = MODEL_LP_LIN; model < N_MODELS; model++) sweep_model(samples, benchmark_sample, model, n_workers);

    lp_free_samples(samples, n_samples);
    return 0;
}
