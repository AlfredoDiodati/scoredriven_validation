/*
The impulse-response protocol run with the first replicate of every
configuration standing in as the benchmark, for the t-QVARMA and the four
local projections (linear, state 1, state 2, both states in one vector).

montecarlo/sweep_irf.c and montecarlo/lp_sweep.c ask how often the confidence
set returns cop_0191 when each of cop_0191's 1000 replicates is the benchmark
in turn. That shows the procedure identifies one configuration. This asks the
same question once for each of the 1000 configurations: configuration c's
replicate 0 (seed 1 of the model) is the benchmark, and the confidence set
should return c.

Each individual run is what sweep_irf.c and lp_sweep.c do for one benchmark:
the mean absolute error between response vectors, the benchmark's replicate
index held out of every column, a replicate with a hole in any column dropped,
MCS_TR at 10000 resamples, block length 1, bootstrap variance, seed 123 stream
0. The held-out index is 0 for every benchmark, so the rows kept are the same
for every benchmark of one model and are decided once.

Response vectors are read once per model and held as float32, as in the two
sweeps this extends: 2.1 GB for the t-QVARMA (525 entries per cell), 1.6 GB for
each local projection (400 entries), 3.2 GB for both states in one vector. The
caches are montecarlo/response_cache.h's files, shared with
montecarlo/sweep_grid.c, and one model's is released before the next is
opened.

SWEEP_COPS_WORKERS benchmarks run side by side, one thread each, default one
per hardware thread. Resumable at benchmark granularity. Rows go to a plain
progress file as each benchmark finishes, in the order they finish, and a rerun skips what it holds. When every model is done
the rows are written to one gzip-compressed csv and the progress file is
removed, so a finished run leaves a single file. A rerun with the result
present and no progress file decompresses the result into the progress file,
so a model added since is run and the rows already there are kept.

Columns: model, benchmark configuration, whether it is in the set, its rank by
mean loss (1 is the smallest), its MCS p-value, the set size, whether the set
was decided by an accepted test, the final p-value, replicates dropped for a
hole, the configuration with the smallest mean loss, seconds spent in the
confidence set alone, and the configurations in the set, space separated.

Requires out/abm_system_fit_qvarma/, out/abm_system_fit_lp/ and
dataset/abm_system/, and rebuilds none of them. Progress goes to stderr.
*/

#include "montecarlo/response_cache.h"
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/frame.h>
#include <et_al./frame/gzip.h>
#include <cblas.h>

#define BENCHMARK_REPLICATE 0
#define PROGRESS_PATH "montecarlo/out/sweep_cops_progress.csv"
#define RESULT_PATH "montecarlo/out/sweep_cops.csv.gz"
#define HEADER "model,benchmark_cop,in_set,rank,pvalue,set_size,decided_by_accepted_test," \
               "final_pvalue,n_dropped_replicates,lowest_loss_cop,mcs_seconds,set_members\n"

static const char *LP_FIT_DIR;
static const char *INPUT_DIR;
static LpSample *samples;
static int n_samples = 0;
static int n_replicates = 0;

/* Which benchmarks the progress file already holds, per model. */
static int *read_done(int *n_done) {
    int *done = calloc((size_t)N_MODELS * n_samples, sizeof(int));
    assert(done);
    for (int model = 0; model < N_MODELS; model++) n_done[model] = 0;
    FILE *f = fopen(PROGRESS_PATH, "r");
    if (!f) return done;

    char line[16384];
    int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        char label[32];
        int cop;
        if (sscanf(line, "%31[^,],%d,", label, &cop) != 2) continue;
        for (int model = 0; model < N_MODELS; model++) {
            if (strcmp(label, model_label[model]) != 0) continue;
            for (int sample = 0; sample < n_samples; sample++)
                if (samples[sample].index == cop && !done[model * n_samples + sample]) {
                    done[model * n_samples + sample] = 1;
                    n_done[model]++;
                }
        }
    }
    fclose(f);
    return done;
}

/* Replicates kept in every loss matrix of one model: every index except the
   benchmark's, and none with a missing cell in any configuration, since et_al's
   mcs refuses a loss matrix with a hole. */
static int *kept_replicates(const ResponseCache *cache, int *n_kept) {
    int *kept = malloc((size_t)n_replicates * sizeof(int));
    assert(kept);
    int n = 0;
    for (int replicate = 0; replicate < n_replicates; replicate++) {
        if (replicate == BENCHMARK_REPLICATE) continue;
        int usable = 1;
        for (int sample = 0; sample < n_samples && usable; sample++)
            if (MISNAN(response_slot(cache, sample, replicate)[0])) usable = 0;
        if (usable) kept[n++] = replicate;
    }
    *n_kept = n;
    return kept;
}

static void sweep_model(int model, const int *done, int n_done, FILE *progress, int n_workers) {
    fprintf(stderr, "%s: %d of %d benchmarks already done\n", model_label[model], n_done, n_samples);
    if (n_done == n_samples) return;

    fprintf(stderr, "%s: building the response cache, %.1f GB\n", model_label[model],
            (double)n_samples * n_replicates * response_dim(model) * sizeof(float) / 1e9);
    ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, LP_FIT_DIR,
                                                       INPUT_DIR, RESPONSE_CACHE_DIR_DEFAULT);
    long n_missing_cells = cache.n_missing;
    int n_kept;
    int *kept = kept_replicates(&cache, &n_kept);
    assert(n_kept > 0 && "sweep_cops: every replicate has a missing cell");
    fprintf(stderr, "%s: cache built, %ld cells missing, %d replicates kept\n", model_label[model],
            n_missing_cells, n_kept);

    const char **model_name = malloc((size_t)n_samples * sizeof(char *));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_%s", samples[i].name, model_label[model]);
        model_name[i] = name_buffer[i];
    }

    /* Benchmarks are independent, so several run at once, each loss matrix
       and confidence set on one thread, with nested parallel regions switched
       off. With one worker the loop runs on the calling thread and each step
       uses all threads instead. Rows reach the progress file in the order they
       finish; the configuration in the second column identifies them. */
    const int dim = response_dim(model);
    if (n_workers > 1) omp_set_max_active_levels(1);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_workers) if (n_workers > 1)
    for (int benchmark_sample = 0; benchmark_sample < n_samples; benchmark_sample++) {
        if (done[model * n_samples + benchmark_sample]) continue;
        const float *reference = response_slot(&cache, benchmark_sample, BENCHMARK_REPLICATE);
        if (MISNAN(reference[0])) {
            fprintf(stderr, "%s %s skipped, its own response is missing\n", model_label[model],
                    samples[benchmark_sample].name);
            continue;
        }

        Mat values = mat_new(n_kept, n_samples);
        #pragma omp parallel for schedule(static)
        for (int row = 0; row < n_kept; row++)
            for (int sample = 0; sample < n_samples; sample++) {
                const float *cell = response_slot(&cache, sample, kept[row]);
                double sum = 0;
                for (int i = 0; i < dim; i++) sum += fabs((double)cell[i] - (double)reference[i]);
                AT(values, row, sample) = (mreal)(sum / dim);
            }
        DataFrame losses = df_from_matrix(values, model_name);
        mat_free(values);

        MCSOptions opt = mcs_options_default();
        opt.bootstrap = 10000;
        opt.block_length = 1;
        opt.variance = MCS_VARIANCE_BOOTSTRAP;
        opt.stat = MCS_TR;
        double started = omp_get_wtime();
        MCSResult res = mcs(&losses, opt);
        double mcs_seconds = omp_get_wtime() - started;

        double benchmark_mean = (double)stats_mean(df_col_numeric(&losses, model_name[benchmark_sample]));
        int rank = 1, lowest = 0;
        double lowest_mean = INFINITY;
        for (int sample = 0; sample < n_samples; sample++) {
            double mean = (double)stats_mean(df_col_numeric(&losses, model_name[sample]));
            if (mean < lowest_mean) { lowest_mean = mean; lowest = sample; }
            if (sample != benchmark_sample && mean < benchmark_mean) rank++;
        }

        #pragma omp critical(sweep_cops_progress)
        {
            fprintf(progress, "%s,%d,%d,%d,%.6f,%d,%d,%.6f,%d,%d,%.2f,", model_label[model],
                    samples[benchmark_sample].index, mcs_in_set(&res, benchmark_sample), rank,
                    res.pvalue[benchmark_sample], res.n_surviving, res.converged, res.final_pvalue,
                    n_replicates - 1 - n_kept, samples[lowest].index, mcs_seconds);
            int first_member = 1;
            for (int sample = 0; sample < n_samples; sample++)
                if (mcs_in_set(&res, sample)) {
                    fprintf(progress, first_member ? "%d" : " %d", samples[sample].index);
                    first_member = 0;
                }
            fprintf(progress, "\n");
            fflush(progress);
            fprintf(stderr, "%s %s: rank %d, in set %d, set size %d, %.1f s\n", model_label[model],
                    samples[benchmark_sample].name, rank, mcs_in_set(&res, benchmark_sample), res.n_surviving,
                    mcs_seconds);
        }

        mcs_free(&res);
        df_free(&losses);
    }

    free(name_buffer);
    free(model_name);
    free(kept);
    response_cache_free(&cache);
}

/* The result back into the progress file, so a rerun after a model was added
   skips the rows the result holds, keeps them, and does the new model. */
static void restore_progress(void) {
    FILE *f = fopen(RESULT_PATH, "rb");
    assert(f && "sweep_cops: cannot open the result");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc((size_t)length);
    assert(compressed && "sweep_cops: out of memory reading the result");
    size_t n_read = fread(compressed, 1, (size_t)length, f);
    assert(n_read == (size_t)length && "sweep_cops: short read of the result");
    fclose(f);

    size_t text_length;
    unsigned char *text = gzip_inflate(compressed, (size_t)length, &text_length);
    FILE *out = fopen(PROGRESS_PATH, "wb");
    assert(out && "sweep_cops: cannot open the progress file for writing");
    size_t n_written = fwrite(text, 1, text_length, out);
    assert(n_written == text_length && "sweep_cops: short write of the progress file");
    fclose(out);
    fprintf(stderr, "%s restored to %s, %zu bytes\n", RESULT_PATH, PROGRESS_PATH, text_length);
    free(text);
    free(compressed);
}

/* The progress file, compressed whole into the result and then removed. */
static void write_result(void) {
    FILE *f = fopen(PROGRESS_PATH, "rb");
    assert(f && "sweep_cops: the progress file is missing");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *text = malloc((size_t)length);
    assert(text && "sweep_cops: out of memory reading the progress file");
    size_t n_read = fread(text, 1, (size_t)length, f);
    assert(n_read == (size_t)length && "sweep_cops: short read of the progress file");
    fclose(f);

    size_t compressed_length;
    unsigned char *compressed = gzip_deflate_level(text, (size_t)length, 9, &compressed_length);
    FILE *out = fopen(RESULT_PATH, "wb");
    assert(out && "sweep_cops: cannot open the result for writing");
    size_t n_written = fwrite(compressed, 1, compressed_length, out);
    assert(n_written == compressed_length && "sweep_cops: short write of the result");
    fclose(out);
    remove(PROGRESS_PATH);
    fprintf(stderr, "wrote %s, %zu bytes from %ld\n", RESULT_PATH, compressed_length, length);
    free(compressed);
    free(text);
}

int main(void) {
    LP_FIT_DIR = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!LP_FIT_DIR) LP_FIT_DIR = LP_FIT_DIR_DEFAULT;
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = LP_INPUT_DIR_DEFAULT;
    openblas_set_num_threads(1);

    struct stat st;
    int progress_exists = stat(PROGRESS_PATH, &st) == 0;
    int restored = stat(RESULT_PATH, &st) == 0 && !progress_exists;

    samples = lp_list_samples(LP_FIT_DIR, &n_samples);
    assert(n_samples > 0 && "sweep_cops: no configurations in the LP fit cache");
    n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);
    for (int i = 0; i < n_samples; i++) {
        char path[640];
        snprintf(path, sizeof path, "%s/%s", RESPONSE_QVARMA_FIT_DIR, samples[i].name);
        FILE *probe = fopen(path, "r");
        assert(probe && "sweep_cops: a configuration has LP fits and no t-QVARMA fits");
        fclose(probe);
    }

    if (restored) restore_progress();
    int n_done[N_MODELS];
    int *done = read_done(n_done);
    int n_done_all = 0;
    for (int model = 0; model < N_MODELS; model++) n_done_all += n_done[model];
    if (restored && n_done_all == N_MODELS * n_samples) {
        remove(PROGRESS_PATH);
        fprintf(stderr, "%s holds every model, nothing to do\n", RESULT_PATH);
        free(done);
        lp_free_samples(samples, n_samples);
        return 0;
    }
    FILE *progress = fopen(PROGRESS_PATH, "a");
    assert(progress && "sweep_cops: cannot open the progress file");
    if (!progress_exists && !restored) {
        fputs(HEADER, progress);
        fflush(progress);
    }

    /* SWEEP_COPS_WORKERS confidence sets side by side, one thread each; by
       default one per hardware thread. Each set at 10000 resamples holds about
       100 MB. */
    const char *workers_text = getenv("SWEEP_COPS_WORKERS");
    int n_workers = workers_text ? atoi(workers_text) : omp_get_num_procs();
    assert(n_workers >= 1 && "sweep_cops: SWEEP_COPS_WORKERS must be at least 1");
    for (int model = 0; model < N_MODELS; model++)
        sweep_model(model, done, n_done[model], progress, n_workers);
    fclose(progress);

    /* A benchmark whose own response is missing is skipped rather than
       recorded, so the result is written once every other one is in. */
    write_result();

    free(done);
    lp_free_samples(samples, n_samples);
    return 0;
}
