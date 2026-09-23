/*
The impulse-response protocol run against every replicate of cop_0191 in turn,
to separate luck from structure.

montecarlo/ promotes one simulated run to the role the US data plays and asks
whether the confidence set returns the configuration that run came from. On the
one benchmark it was first run against, the impulse-response loss returned it
and both score losses did not. One benchmark cannot tell a protocol that works
from a protocol that got a lucky draw, so this repeats the whole thing with
each of the benchmark configuration's 1000 replicates standing in as the
benchmark, and reports how often the answer comes back right.

Everything about each individual run is what montecarlo/irf_loss.c and
montecarlo/mcs.c do: same loss, same hold-out, same confidence set settings.
The only thing this adds is the loop and the tally.

What makes the loop affordable. A cell's impulse response vector does not
depend on which benchmark it is being compared against - only the comparison
does. montecarlo/irf_loss.c recomputes all 1,000,000 of them on every run,
which is 574 seconds measured on this machine, and repeating that per benchmark
would be 160 hours. They are computed once here, held as float32, and every
benchmark's loss matrix is then one pass of mean absolute error over them. The
cache is 1000 x 1000 x 525 floats, about 2.1 GB.

float32 for the cache and double for the arithmetic: the loss is a mean of
absolute differences over 525 entries, accumulated in double, so the only
precision lost is in the stored responses themselves. That is the same
precision the responses are printed and plotted at, and it buys the cache
fitting in memory at all.

Resumable at benchmark granularity. Every finished benchmark is appended to the
summary as it finishes, and a rerun reads back which ones are already there and
skips them. The impulse-response cache is not persisted, so a restart pays the
574 seconds again but nothing else.

Requires out/abm_system_fit_qvarma/ and dataset/abm_system/, and rebuilds
neither. Writes montecarlo/out/sweep_irf.csv, one row per benchmark. Nothing
printed except progress, which goes to stderr so the summary is the result.
*/

#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define P 1
#define Q 1
#define SPEC_R 2
#define SPEC_LABEL "p1q1r2"
#define COLUMN_LABEL "qvarma_" SPEC_LABEL
#define HORIZON 20
#define IRF_DIM (K * K * (HORIZON + 1))

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define SUMMARY_PATH "montecarlo/out/sweep_irf.csv"
#define BENCHMARK_SAMPLE "cop_0191"

static int n_samples = 0;
static int n_replicates = 0;

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

typedef struct { char *name; int index; } SampleEntry;

static int trailing_index(const char *name) {
    const char *underscore = strrchr(name, '_');
    assert(underscore && "sweep_irf: a configuration directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "sweep_irf: cannot open the fit cache");

    SampleEntry *entries = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char path[560];
        snprintf(path, sizeof path, "%s/%s", dir, entry->d_name);
        struct stat st;
        if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            SampleEntry *grown = realloc(entries, (size_t)cap * sizeof(SampleEntry));
            assert(grown && "sweep_irf: out of memory listing configurations");
            entries = grown;
        }
        size_t len = strlen(entry->d_name);
        entries[n].name = malloc(len + 1);
        assert(entries[n].name);
        memcpy(entries[n].name, entry->d_name, len + 1);
        entries[n].index = trailing_index(entry->d_name);
        n++;
    }
    closedir(handle);
    qsort(entries, (size_t)n, sizeof(SampleEntry), compare_sample_entries);
    *count = n;
    return entries;
}

static int count_replicates(const char *sample) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    int count = 0;
    int *replicate = abm_system_list_replicates(dir, &count);
    assert(count > 0);
    int highest = replicate[0];
    for (int i = 1; i < count; i++) if (replicate[i] > highest) highest = replicate[i];
    free(replicate);
    return highest + 1;
}

/* Where one cell's response vector lives in the cache. */
static float *irf_slot(float *cache, int sample, int replicate) {
    return cache + ((size_t)sample * n_replicates + replicate) * IRF_DIM;
}

/* Fills out with the stacked total responses, horizon 0 first, and returns 1.
   Returns 0 when nu <= 2, which impulse_responses would abort on, or when a
   response comes back non-finite. */
static int cell_irf(const QvarmaParams *m, Mat y, float *out) {
    if (!(m->nu > 2)) return 0;
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses r = qvarma_impulse_responses(m, D, options);

    int at = 0, ok = 1;
    for (int h = 0; h <= r.horizon; h++)
        for (int i = 0; i < K * K; i++) {
            mreal value = r.total[h].d[i];
            if (MISNAN(value) || MISINF(value)) ok = 0;
            out[at++] = (float)value;
        }

    qvarma_impulse_responses_free(&r);
    mat_free(D);
    return ok;
}

/* Every cell's response vector, computed once. A cell that cannot produce one
   is marked by a NaN in its first entry and counted. */
static float *build_irf_cache(const SampleEntry *samples, long *n_missing_out) {
    size_t cells = (size_t)n_samples * n_replicates;
    float *cache = (float*)malloc(cells * IRF_DIM * sizeof(float));
    assert(cache && "sweep_irf: out of memory for the impulse-response cache");

    long n_missing = 0;
    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams working = spec_shape();

        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < n_samples; sample++) {
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[sample].name);

            int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
            for (int batch = 0; batch < n_batches; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);
                for (int b = 0; b < count; b++) {
                    if (replicate[b] < n_replicates) {
                        float *slot = irf_slot(cache, sample, replicate[b]);
                        char cache_path[640];
                        snprintf(cache_path, sizeof cache_path,
                                 "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR,
                                 samples[sample].name, replicate[b], SPEC_LABEL);
                        int ok = qvarma_load_params(&working, cache_path)
                              && cell_irf(&working, block[b], slot);
                        if (!ok) { slot[0] = (float)NAN; n_missing++; }
                    }
                    mat_free(block[b]);
                }
            }
            if (sample % 50 == 0) fprintf(stderr, "  irf cache: %d of %d configurations\n",
                                          sample, n_samples);
        }

        qvarma_params_free(&working);
    }

    *n_missing_out = n_missing;
    return cache;
}

/* Which benchmarks the summary already holds, so a rerun continues rather than
   repeats. */
static int *read_done(int *count) {
    int *done = (int*)calloc((size_t)n_replicates, sizeof(int));
    assert(done);
    *count = 0;
    FILE *f = fopen(SUMMARY_PATH, "r");
    if (!f) return done;

    char line[512];
    int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        int replicate;
        if (sscanf(line, "%d,", &replicate) == 1 && replicate >= 0 && replicate < n_replicates) {
            done[replicate] = 1;
            (*count)++;
        }
    }
    fclose(f);
    return done;
}

int main(void) {
    openblas_set_num_threads(1);

    SampleEntry *samples = list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0);
    n_replicates = count_replicates(samples[0].name);

    int benchmark_sample = -1;
    for (int i = 0; i < n_samples; i++)
        if (strcmp(samples[i].name, BENCHMARK_SAMPLE) == 0) benchmark_sample = i;
    assert(benchmark_sample >= 0 && "sweep_irf: the benchmark configuration is not in the cache");

    int n_done;
    int *done = read_done(&n_done);
    fprintf(stderr, "%d of %d benchmarks already done\n", n_done, n_replicates);
    if (n_done == n_replicates) { free(done); return 0; }

    fprintf(stderr, "building the impulse-response cache, %.1f GB\n",
            (double)n_samples * n_replicates * IRF_DIM * sizeof(float) / 1e9);
    long n_missing_cells;
    float *cache = build_irf_cache(samples, &n_missing_cells);
    fprintf(stderr, "cache built, %ld cells missing\n", n_missing_cells);

    FILE *summary = fopen(SUMMARY_PATH, n_done ? "a" : "w");
    assert(summary && "sweep_irf: cannot open the summary for writing");
    if (!n_done) {
        fprintf(summary, "benchmark_replicate,benchmark_in_set,benchmark_rank,benchmark_pvalue,"
                          "set_size,decided_by_accepted_test,final_pvalue,n_dropped_replicates\n");
        fflush(summary);
    }

    const char **model_name = (const char**)malloc((size_t)n_samples * sizeof(char*));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_%s", samples[i].name, COLUMN_LABEL);
        model_name[i] = name_buffer[i];
    }

    for (int benchmark = 0; benchmark < n_replicates; benchmark++) {
        if (done[benchmark]) continue;
        const float *reference = irf_slot(cache, benchmark_sample, benchmark);
        if (MISNAN(reference[0])) {
            fprintf(stderr, "benchmark %03d skipped, its own impulse response is missing\n",
                    benchmark);
            continue;
        }

        /* One row per replicate kept: the benchmark's own index is held out of
           every column, for the reason montecarlo/benchmark.h gives. */
        Mat values = mat_new(n_replicates - 1, n_samples);
        int *usable = (int*)malloc((size_t)(n_replicates - 1) * sizeof(int));
        int row = 0;
        for (int replicate = 0; replicate < n_replicates; replicate++) {
            if (replicate == benchmark) continue;
            usable[row] = 1;
            #pragma omp parallel for schedule(static)
            for (int sample = 0; sample < n_samples; sample++) {
                const float *cell = irf_slot(cache, sample, replicate);
                if (MISNAN(cell[0])) { AT(values, row, sample) = (mreal)NAN; continue; }
                double sum = 0;
                for (int i = 0; i < IRF_DIM; i++) sum += fabs((double)cell[i] - (double)reference[i]);
                AT(values, row, sample) = (mreal)(sum / IRF_DIM);
            }
            for (int sample = 0; sample < n_samples; sample++)
                if (MISNAN(AT(values, row, sample))) { usable[row] = 0; break; }
            row++;
        }

        int keep = 0;
        for (int r = 0; r < n_replicates - 1; r++) if (usable[r]) keep++;
        assert(keep > 0 && "sweep_irf: every replicate of a benchmark has a missing cell");

        /* One df_from_matrix rather than a thousand df_add_numeric_col calls:
           that function grows the frame by one column at a time, copying the
           whole block each time, which over a thousand columns and a thousand
           benchmarks is an hour of memcpy. */
        Mat kept = values;
        if (keep < n_replicates - 1) {
            kept = mat_new(keep, n_samples);
            int at = 0;
            for (int r = 0; r < n_replicates - 1; r++) {
                if (!usable[r]) continue;
                for (int sample = 0; sample < n_samples; sample++)
                    AT(kept, at, sample) = AT(values, r, sample);
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

        /* Where the benchmark's own configuration came, by mean loss and by the
           confidence set's verdict. Rank is 1 for the smallest mean loss. */
        double benchmark_mean = (double)stats_mean(df_col_numeric(&losses, model_name[benchmark_sample]));
        int rank = 1;
        for (int sample = 0; sample < n_samples; sample++) {
            if (sample == benchmark_sample) continue;
            if ((double)stats_mean(df_col_numeric(&losses, model_name[sample])) < benchmark_mean) rank++;
        }

        fprintf(summary, "%d,%d,%d,%.6f,%d,%s,%.6f,%d\n", benchmark,
                mcs_in_set(&res, benchmark_sample), rank, res.pvalue[benchmark_sample],
                res.n_surviving, res.converged ? "yes" : "no", res.final_pvalue,
                n_replicates - 1 - keep);
        fflush(summary);
        fprintf(stderr, "benchmark %03d: rank %d, in set %d, set size %d\n", benchmark, rank,
                mcs_in_set(&res, benchmark_sample), res.n_surviving);

        mcs_free(&res);
        df_free(&losses);
    }

    fclose(summary);
    free(name_buffer);
    free(model_name);
    free(cache);
    free(done);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
