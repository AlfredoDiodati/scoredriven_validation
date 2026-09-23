/*
The two score protocols run against every replicate of cop_0191 in turn, the
companion to montecarlo/sweep_irf.c.

On the single benchmark montecarlo/ was first run against, both score losses
failed to return the configuration the benchmark came from while the
impulse-response loss returned it alone. This repeats both of them with each of
that configuration's 1000 replicates standing in as the benchmark, so that a
protocol that is structurally unable to identify a configuration can be told
apart from one that drew badly once.

Each individual run is what montecarlo/score_loss.c and montecarlo/mcs.c do:
the same two losses, the same hold-out of the benchmark's own seed, the same
confidence set settings. The loop and the tally are all that is added.

Why this is the expensive half. The score is evaluated at the benchmark's own
estimate, so unlike an impulse response it cannot be computed once and reused:
every one of the 1,000,000 cells has to be filtered again for every benchmark,
which is 10^9 filter evaluations. What can be amortised is the 15 GB of
archives, by holding a block of benchmarks at once and scoring each cell
against all of them while its series is in hand. BENCHMARK_BLOCK sets how many;
the memory it costs is two loss matrices per benchmark in the block, about 16
MB each.

Resumable at benchmark granularity, like sweep_irf.c: finished benchmarks are
appended to the two summaries as they finish and skipped on a rerun. A
benchmark whose information matrix is not positive definite has no weighted
loss and is recorded as skipped there while still being scored under q'q.

Requires out/abm_system_fit_qvarma/ for every benchmark's estimate and
dataset/abm_system/ for every series, and rebuilds neither. Writes
montecarlo/out/sweep_score.csv and montecarlo/out/sweep_score_weighted.csv.
Progress goes to stderr so the summaries are the result.
*/

#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./linalg/decomp.h>
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

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define PLAIN_SUMMARY "montecarlo/out/sweep_score.csv"
#define WEIGHTED_SUMMARY "montecarlo/out/sweep_score_weighted.csv"
#define BENCHMARK_SAMPLE "cop_0191"
#define BENCHMARK_BLOCK 25

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
    assert(underscore && "sweep_score: a configuration directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "sweep_score: cannot open dataset/abm_system/");

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
            assert(grown);
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

/*
One benchmark's estimate and the weighting built at it: theta, and the Cholesky
factor of T_sim/T_benchmark times the observed information there. factor.d is
NULL when that matrix is not positive definite, which leaves the weighted loss
undefined for this benchmark and the plain one unaffected.
*/
typedef struct {
    int replicate;
    Vec theta;
    Mat factor;
    int weighting_usable;
} BenchmarkPoint;

static BenchmarkPoint benchmark_point(int replicate, int sim_periods) {
    BenchmarkPoint point;
    point.replicate = replicate;
    point.weighting_usable = 0;
    point.factor = (Mat){ 0, 0, 0, NULL };

    QvarmaParams m = spec_shape();
    char path[640];
    snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, BENCHMARK_SAMPLE,
             replicate, SPEC_LABEL);
    int loaded = qvarma_load_params(&m, path);
    assert(loaded && "sweep_score: a benchmark's cached fit will not load");

    int n = qvarma_n_theta(&m);
    point.theta = mat_new(n, 1);
    _qvarma_unlink(&m, point.theta);

    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, BENCHMARK_SAMPLE);
    Mat y = abm_system_read_replicate(dir, replicate);

    Mat information = _qvarma_hessian(point.theta, &m, y);
    Vec eigenvalues;
    Mat eigenvectors;
    if (mat_eig_sym_status(information, &eigenvalues, &eigenvectors) == 0) {
        mreal smallest = eigenvalues.d[0], largest = eigenvalues.d[0];
        for (int i = 1; i < n; i++) {
            if (eigenvalues.d[i] < smallest) smallest = eigenvalues.d[i];
            if (eigenvalues.d[i] > largest) largest = eigenvalues.d[i];
        }
        mreal floor_value = (mreal)(n * MEPS) * MABS(largest);
        if (smallest > floor_value) {
            mreal scale = (mreal)sim_periods / (mreal)y.c;
            Mat scaled = mat_copy(information);
            for (int i = 0; i < n * n; i++) scaled.d[i] *= scale;
            point.factor = mat_chol(scaled);
            point.weighting_usable = 1;
            mat_free(scaled);
        }
        mat_free(eigenvalues);
        mat_free(eigenvectors);
    }

    mat_free(information);
    mat_free(y);
    qvarma_params_free(&m);
    return point;
}

static void benchmark_point_free(BenchmarkPoint *point) {
    mat_free(point->theta);
    if (point->factor.d) mat_free(point->factor);
}

/* Which benchmarks a summary already holds, so a rerun continues. */
static int *read_done(const char *path, int *count) {
    int *done = (int*)calloc((size_t)n_replicates, sizeof(int));
    assert(done);
    *count = 0;
    FILE *f = fopen(path, "r");
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

/*
One loss matrix through MCS_TR, and where the benchmark's own configuration
landed. values is rows x n_samples and is consumed here.
*/
static void score_and_report(Mat values, const char *const *model_name, int benchmark_sample,
                             int benchmark, int n_dropped, FILE *summary) {
    DataFrame losses = df_from_matrix(values, model_name);

    MCSOptions opt = mcs_options_default();
    opt.bootstrap = 10000;
    opt.block_length = 1;
    opt.variance = MCS_VARIANCE_BOOTSTRAP;
    opt.stat = MCS_TR;
    MCSResult res = mcs(&losses, opt);

    double benchmark_mean = (double)stats_mean(df_col_numeric(&losses, model_name[benchmark_sample]));
    int rank = 1;
    for (int sample = 0; sample < n_samples; sample++) {
        if (sample == benchmark_sample) continue;
        if ((double)stats_mean(df_col_numeric(&losses, model_name[sample])) < benchmark_mean) rank++;
    }

    fprintf(summary, "%d,%d,%d,%.6f,%d,%s,%.6f,%d\n", benchmark,
            mcs_in_set(&res, benchmark_sample), rank, res.pvalue[benchmark_sample],
            res.n_surviving, res.converged ? "yes" : "no", res.final_pvalue, n_dropped);
    fflush(summary);

    mcs_free(&res);
    df_free(&losses);
}

int main(void) {
    openblas_set_num_threads(1);

    SampleEntry *samples = list_samples(INPUT_DIR, &n_samples);
    assert(n_samples > 0);
    n_replicates = count_replicates(samples[0].name);

    int benchmark_sample = -1;
    for (int i = 0; i < n_samples; i++)
        if (strcmp(samples[i].name, BENCHMARK_SAMPLE) == 0) benchmark_sample = i;
    assert(benchmark_sample >= 0 && "sweep_score: the benchmark configuration is not in the dataset");

    int n_done_plain, n_done_weighted;
    int *done_plain = read_done(PLAIN_SUMMARY, &n_done_plain);
    int *done_weighted = read_done(WEIGHTED_SUMMARY, &n_done_weighted);
    fprintf(stderr, "of %d benchmarks: %d done under q'q, %d under the weighted statistic\n",
            n_replicates, n_done_plain, n_done_weighted);

    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[0].name);
    Mat probe[ABM_SYSTEM_BATCH];
    int probe_replicate[ABM_SYSTEM_BATCH];
    int probe_count = abm_system_read_batch(dir, 0, probe, probe_replicate);
    assert(probe_count > 0);
    int sim_periods = probe[0].c;
    for (int b = 0; b < probe_count; b++) mat_free(probe[b]);

    FILE *plain_summary = fopen(PLAIN_SUMMARY, n_done_plain ? "a" : "w");
    FILE *weighted_summary = fopen(WEIGHTED_SUMMARY, n_done_weighted ? "a" : "w");
    assert(plain_summary && weighted_summary);
    const char *header = "benchmark_replicate,benchmark_in_set,benchmark_rank,benchmark_pvalue,"
                         "set_size,decided_by_accepted_test,final_pvalue,n_dropped_replicates\n";
    if (!n_done_plain) { fprintf(plain_summary, "%s", header); fflush(plain_summary); }
    if (!n_done_weighted) { fprintf(weighted_summary, "%s", header); fflush(weighted_summary); }

    const char **model_name = (const char**)malloc((size_t)n_samples * sizeof(char*));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_%s", samples[i].name, COLUMN_LABEL);
        model_name[i] = name_buffer[i];
    }

    for (int first = 0; first < n_replicates; first += BENCHMARK_BLOCK) {
        int last = first + BENCHMARK_BLOCK;
        if (last > n_replicates) last = n_replicates;

        int block_size = 0;
        BenchmarkPoint point[BENCHMARK_BLOCK];
        for (int benchmark = first; benchmark < last; benchmark++) {
            if (done_plain[benchmark] && done_weighted[benchmark]) continue;
            point[block_size++] = benchmark_point(benchmark, sim_periods);
        }
        if (block_size == 0) continue;
        fprintf(stderr, "block %d-%d, %d benchmarks to score\n", first, last - 1, block_size);

        /* Two loss matrices per benchmark in the block, filled in one pass over
           the archives: every cell is scored against every benchmark of the
           block while its series is in hand. */
        Mat *plain = malloc((size_t)block_size * sizeof(Mat));
        Mat *weighted = malloc((size_t)block_size * sizeof(Mat));
        assert(plain && weighted);
        for (int b = 0; b < block_size; b++) {
            plain[b] = mat_new(n_replicates, n_samples);
            weighted[b] = mat_new(n_replicates, n_samples);
            for (int i = 0; i < n_replicates * n_samples; i++) {
                plain[b].d[i] = (mreal)NAN;
                weighted[b].d[i] = (mreal)NAN;
            }
        }

        #pragma omp parallel
        {
            QvarmaParams shape = spec_shape();
            QvarmaAnalytic *filter = qvarma_analytic_new(&shape, sim_periods);
            Vec score = mat_new(point[0].theta.r, 1);
            Vec solved = mat_new(point[0].theta.r, 1);

            #pragma omp for schedule(dynamic)
            for (int sample = 0; sample < n_samples; sample++) {
                char sample_dir[560];
                snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, samples[sample].name);
                int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
                for (int batch = 0; batch < n_batches; batch++) {
                    Mat block[ABM_SYSTEM_BATCH];
                    int replicate[ABM_SYSTEM_BATCH];
                    int count = abm_system_read_batch(sample_dir, batch, block, replicate);
                    for (int c = 0; c < count; c++) {
                        int rep = replicate[c];
                        if (rep < n_replicates && block[c].c == sim_periods) {
                            for (int b = 0; b < block_size; b++) {
                                if (rep == point[b].replicate) continue;
                                mreal value = qvarma_analytic_log_likelihood(filter,
                                                  point[b].theta, block[c], score);
                                if (MISNAN(value) || MISINF(value)) continue;
                                double sum = 0;
                                for (int i = 0; i < score.r; i++)
                                    sum += (double)score.d[i] * (double)score.d[i];
                                if (isfinite(sum)) AT(plain[b], rep, sample) = (mreal)sum;

                                if (!point[b].weighting_usable) continue;
                                for (int i = 0; i < score.r; i++) solved.d[i] = score.d[i];
                                if (_trtrs('L', 'N', 'N', point[b].factor.r, 1, point[b].factor.d,
                                           point[b].factor.stride, solved.d, solved.stride) != 0)
                                    continue;
                                double w = 0;
                                for (int i = 0; i < solved.r; i++)
                                    w += (double)solved.d[i] * (double)solved.d[i];
                                if (isfinite(w)) AT(weighted[b], rep, sample) = (mreal)w;
                            }
                        }
                        mat_free(block[c]);
                    }
                }
            }

            mat_free(score);
            mat_free(solved);
            qvarma_analytic_free(filter);
            qvarma_params_free(&shape);
        }

        for (int b = 0; b < block_size; b++) {
            int benchmark = point[b].replicate;

            /* The benchmark's own row is never filled, and any row with a hole
               comes out with the rest. */
            int keep = 0;
            int *usable = (int*)malloc((size_t)n_replicates * sizeof(int));
            for (int rep = 0; rep < n_replicates; rep++) {
                usable[rep] = rep != benchmark;
                for (int sample = 0; sample < n_samples && usable[rep]; sample++)
                    if (MISNAN(AT(plain[b], rep, sample))) usable[rep] = 0;
                if (usable[rep]) keep++;
            }
            assert(keep > 0 && "sweep_score: a benchmark has no usable replicate");

            if (!done_plain[benchmark]) {
                Mat kept = mat_new(keep, n_samples);
                int at = 0;
                for (int rep = 0; rep < n_replicates; rep++) {
                    if (!usable[rep]) continue;
                    for (int sample = 0; sample < n_samples; sample++)
                        AT(kept, at, sample) = AT(plain[b], rep, sample);
                    at++;
                }
                score_and_report(kept, model_name, benchmark_sample, benchmark,
                                 n_replicates - 1 - keep, plain_summary);
                mat_free(kept);
            }

            if (!done_weighted[benchmark]) {
                if (!point[b].weighting_usable) {
                    fprintf(weighted_summary, "%d,-1,-1,nan,-1,skipped,nan,-1\n", benchmark);
                    fflush(weighted_summary);
                } else {
                    int weighted_keep = 0;
                    int *weighted_usable = (int*)malloc((size_t)n_replicates * sizeof(int));
                    for (int rep = 0; rep < n_replicates; rep++) {
                        weighted_usable[rep] = rep != benchmark;
                        for (int sample = 0; sample < n_samples && weighted_usable[rep]; sample++)
                            if (MISNAN(AT(weighted[b], rep, sample))) weighted_usable[rep] = 0;
                        if (weighted_usable[rep]) weighted_keep++;
                    }
                    assert(weighted_keep > 0);
                    Mat kept = mat_new(weighted_keep, n_samples);
                    int at = 0;
                    for (int rep = 0; rep < n_replicates; rep++) {
                        if (!weighted_usable[rep]) continue;
                        for (int sample = 0; sample < n_samples; sample++)
                            AT(kept, at, sample) = AT(weighted[b], rep, sample);
                        at++;
                    }
                    score_and_report(kept, model_name, benchmark_sample, benchmark,
                                     n_replicates - 1 - weighted_keep, weighted_summary);
                    mat_free(kept);
                    free(weighted_usable);
                }
            }

            fprintf(stderr, "  benchmark %03d scored\n", benchmark);
            free(usable);
            mat_free(plain[b]);
            mat_free(weighted[b]);
            benchmark_point_free(&point[b]);
        }

        free(plain);
        free(weighted);
    }

    fclose(plain_summary);
    fclose(weighted_summary);
    free(name_buffer);
    free(model_name);
    free(done_plain);
    free(done_weighted);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
