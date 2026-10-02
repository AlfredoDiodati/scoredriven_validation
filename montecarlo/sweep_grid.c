/*
The impulse-response protocol run with every simulated run of every
configuration standing in as the benchmark: 1000 configurations times 1000
runs, for the t-QVARMA and the three local projections (linear, state 1,
state 2), four million confidence sets.

montecarlo/sweep_cops.c uses one run per configuration as the benchmark, run 0.
This uses all of them. Each individual benchmark is what sweep_cops.c does for
one: the mean absolute error between response vectors, the benchmark's run
index held out of every column, a run with a missing response in any column
dropped, MCS_TR with block length 1, bootstrap variance, seed 123 stream 0,
here at 2000 resamples rather than 10000.

The order is by run first: run 0 of every configuration, then run 1, and so on,
so a run stopped early holds every configuration at the same number of
benchmarks. SWEEP_GRID_RUN_BLOCK sets how many runs one model does before the
next model takes the same runs (default all of them, one model at a time,
t-QVARMA first); a smaller block keeps the four models level.

Benchmarks are taken 64 at a time, all with the same run index, so they share
the rows kept. Their 64 loss matrices are formed in one pass over the response
cache, each cell's responses read once and compared against the 64 benchmark
responses in turn. Each loss is the same sum, in the same order, as
sweep_cops.c forms. The 64 confidence sets then run side by side, one thread
each, which measured faster than one set at a time on every thread (0.121 s
against 0.182 s per set, 2000 resamples, 999 x 1000 tables).

Resumable at benchmark granularity. Rows go to a plain progress file as each
batch finishes, and a rerun skips what it holds. When every model is done the
rows are written to one gzip-compressed csv and the progress file is removed.

Columns: model, benchmark configuration, benchmark run, whether the benchmark's
configuration is in the set, its rank by mean loss (1 is the smallest), its MCS
p-value, the set size, whether the set was decided by an accepted test, the
final p-value, runs dropped for a missing response, the configuration with the
smallest mean loss, and the configurations in the set, space separated.

SWEEP_GRID_WORKERS sets how many confidence sets run side by side (default 8,
one per physical core here). Requires out/abm_system_fit_qvarma/,
out/abm_system_fit_lp/ and dataset/abm_system/, and rebuilds none of them.
Progress goes to stderr. Each model's response cache is kept in
out/sweep_grid_response_cache/ (1.6 to 2.1 GB each, ignored by git) so a
restart reads it in seconds.
*/

#include "montecarlo/response_cache.h"
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/frame.h>
#include <et_al./frame/gzip.h>
#include <cblas.h>

/* The rest can be overridden at compile time for a test build that
   covers only the first runs, models or configurations. */
#define BATCH 64
#ifndef BOOTSTRAP
#define BOOTSTRAP 2000
#endif
#ifndef PROGRESS_PATH
#define PROGRESS_PATH "montecarlo/out/sweep_grid_progress.csv"
#endif
#ifndef RESULT_PATH
#define RESULT_PATH "montecarlo/out/sweep_grid.csv.gz"
#endif
#ifndef RESPONSE_CACHE_DIR
#define RESPONSE_CACHE_DIR "out/sweep_grid_response_cache"
#endif
#ifndef FIRST_MODEL
#define FIRST_MODEL 0
#endif
#ifndef RUN_LIMIT
#define RUN_LIMIT 1000000
#endif
#ifndef BENCHMARK_LIMIT
#define BENCHMARK_LIMIT 1000000
#endif
#define HEADER "model,benchmark_cop,benchmark_run,in_set,rank,pvalue,set_size,decided_by_accepted_test," \
               "final_pvalue,n_dropped_replicates,lowest_loss_cop,set_members\n"

static const char *LP_FIT_DIR;
static const char *INPUT_DIR;
static LpSample *samples;
static int n_samples = 0;
static int n_replicates = 0;

static size_t done_index(int model, int sample, int replicate) {
    return ((size_t)model * n_samples + sample) * n_replicates + replicate;
}

/* Which benchmarks the progress file already holds. */
static unsigned char *read_done(long *n_done) {
    unsigned char *done = calloc((size_t)N_MODELS * n_samples * n_replicates, 1);
    assert(done);
    int *sample_of_cop = NULL;
    int highest = 0;
    for (int i = 0; i < n_samples; i++) if (samples[i].index > highest) highest = samples[i].index;
    sample_of_cop = malloc((size_t)(highest + 1) * sizeof(int));
    for (int c = 0; c <= highest; c++) sample_of_cop[c] = -1;
    for (int i = 0; i < n_samples; i++) sample_of_cop[samples[i].index] = i;

    *n_done = 0;
    FILE *f = fopen(PROGRESS_PATH, "r");
    if (f) {
        char line[16384];
        int first = 1;
        while (fgets(line, sizeof line, f)) {
            if (first) { first = 0; continue; }
            char label[32];
            int cop, replicate;
            if (sscanf(line, "%31[^,],%d,%d,", label, &cop, &replicate) != 3) continue;
            if (cop < 0 || cop > highest || sample_of_cop[cop] < 0) continue;
            if (replicate < 0 || replicate >= n_replicates) continue;
            for (int model = 0; model < N_MODELS; model++) {
                if (strcmp(label, model_label[model]) != 0) continue;
                size_t at = done_index(model, sample_of_cop[cop], replicate);
                if (!done[at]) { done[at] = 1; (*n_done)++; }
            }
        }
        fclose(f);
    }
    free(sample_of_cop);
    return done;
}

/* One benchmark's row of the result. */
typedef struct {
    int sample, replicate, in_set, rank, set_size, converged, lowest;
    double pvalue, final_pvalue;
    int *members;
} Outcome;

/* The confidence set over one benchmark's loss matrix, on the calling thread
   alone, and what the result records about it. Rank and the smallest mean
   loss use stats_mean over each column, as sweep_cops.c does. */
static Outcome run_benchmark(const DataFrame *losses, int benchmark_sample) {
    MCSOptions opt = mcs_options_default();
    opt.bootstrap = BOOTSTRAP;
    opt.block_length = 1;
    opt.variance = MCS_VARIANCE_BOOTSTRAP;
    opt.stat = MCS_TR;
    MCSResult res = mcs(losses, opt);

    Mat column = losses->numeric;
    column.c = 1;
    column.d = losses->numeric.d + benchmark_sample;
    double benchmark_mean = (double)stats_mean(column);
    int rank = 1, lowest = 0;
    double lowest_mean = INFINITY;
    for (int sample = 0; sample < n_samples; sample++) {
        column.d = losses->numeric.d + sample;
        double mean = (double)stats_mean(column);
        if (mean < lowest_mean) { lowest_mean = mean; lowest = sample; }
        if (sample != benchmark_sample && mean < benchmark_mean) rank++;
    }

    Outcome o;
    o.sample = benchmark_sample;
    o.in_set = mcs_in_set(&res, benchmark_sample);
    o.rank = rank;
    o.set_size = res.n_surviving;
    o.converged = res.converged;
    o.lowest = lowest;
    o.pvalue = res.pvalue[benchmark_sample];
    o.final_pvalue = res.final_pvalue;
    o.members = malloc((size_t)(res.n_surviving + 1) * sizeof(int));
    int at = 0;
    for (int sample = 0; sample < n_samples; sample++)
        if (mcs_in_set(&res, sample)) o.members[at++] = sample;
    o.members[at] = -1;
    mcs_free(&res);
    return o;
}

/* Progress over the whole session, for the time-left estimate. */
static double session_started;
static long session_finished = 0, session_remaining = 0;

static long count_remaining(int model, int run_from, int run_to, const unsigned char *done) {
    long remaining = 0;
    for (int sample = 0; sample < n_samples; sample++)
        for (int replicate = run_from; replicate < run_to; replicate++)
            remaining += !done[done_index(model, sample, replicate)];
    return remaining;
}

/* Every benchmark of one model whose run lies in [run_from, run_to). */
static void sweep_model(int model, int run_from, int run_to, const unsigned char *done, int n_workers,
                        FILE *progress) {
    long remaining = count_remaining(model, run_from, run_to, done);
    fprintf(stderr, "%s runs %d to %d: %ld benchmarks to do\n", model_label[model], run_from, run_to - 1,
            remaining);
    if (!remaining) return;

    fprintf(stderr, "%s: building the response cache, %.1f GB\n", model_label[model],
            (double)n_samples * n_replicates * response_dim(model) * sizeof(float) / 1e9);
    ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, LP_FIT_DIR,
                                                       INPUT_DIR, RESPONSE_CACHE_DIR);
    const int dim = cache.dim;

    /* A run with a missing response in any configuration is dropped from
       every loss matrix, since et_al's mcs refuses a matrix with a hole. */
    unsigned char *usable = malloc((size_t)n_replicates);
    int n_usable = 0;
    for (int replicate = 0; replicate < n_replicates; replicate++) {
        usable[replicate] = 1;
        for (int sample = 0; sample < n_samples && usable[replicate]; sample++)
            if (MISNAN(response_slot(&cache, sample, replicate)[0])) usable[replicate] = 0;
        n_usable += usable[replicate];
    }
    fprintf(stderr, "%s: cache built, %ld cells missing, %d of %d runs usable\n", model_label[model],
            cache.n_missing, n_usable, n_replicates);

    const char **model_name = malloc((size_t)n_samples * sizeof(char *));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_%s", samples[i].name, model_label[model]);
        model_name[i] = name_buffer[i];
    }

    /* One loss matrix per slot, reused across batches; the benchmark
       responses of a batch transposed, entry i of every slot together, and
       converted to double once. */
    DataFrame slot[BATCH];
    int slot_rows = -1;
    double *reference = calloc((size_t)dim * BATCH, sizeof(double));
    int *kept = malloc((size_t)n_replicates * sizeof(int));
    int batch_sample[BATCH];
    Outcome outcome[BATCH];

    for (int replicate = run_from; replicate < run_to && replicate < RUN_LIMIT; replicate++) {
        int n_kept = 0;
        for (int r = 0; r < n_replicates; r++) if (usable[r] && r != replicate) kept[n_kept++] = r;
        int dropped = n_replicates - 1 - n_kept;
        if (n_kept != slot_rows) {
            if (slot_rows >= 0) for (int k = 0; k < BATCH; k++) df_free(&slot[k]);
            Mat empty = mat_new(n_kept, n_samples);
            for (int k = 0; k < BATCH; k++) slot[k] = df_from_matrix(empty, model_name);
            mat_free(empty);
            slot_rows = n_kept;
        }

        int last = n_samples < BENCHMARK_LIMIT ? n_samples : BENCHMARK_LIMIT;
        for (int first = 0; first < last;) {
            int n_batch = 0;
            while (first < last && n_batch < BATCH) {
                int sample = first++;
                if (done[done_index(model, sample, replicate)]) continue;
                if (MISNAN(response_slot(&cache, sample, replicate)[0])) {
                    fprintf(stderr, "%s %s run %d skipped, its own response is missing\n",
                            model_label[model], samples[sample].name, replicate);
                    continue;
                }
                batch_sample[n_batch++] = sample;
            }
            if (!n_batch) continue;

            for (int k = 0; k < n_batch; k++) {
                const float *own = response_slot(&cache, batch_sample[k], replicate);
                for (int i = 0; i < dim; i++) reference[(size_t)i * BATCH + k] = (double)own[i];
            }

            /* Each loss is sum_i |cell_i - reference_i| over i in order,
               divided by dim, exactly as in sweep_cops.c; the batch only
               changes how many references one read of a cell serves. */
            #pragma omp parallel for schedule(static)
            for (int row = 0; row < n_kept; row++) {
                double sum[BATCH];
                for (int sample = 0; sample < n_samples; sample++) {
                    const float *cell = response_slot(&cache, sample, kept[row]);
                    for (int k = 0; k < BATCH; k++) sum[k] = 0;
                    for (int i = 0; i < dim; i++) {
                        double x = (double)cell[i];
                        const double *ref = reference + (size_t)i * BATCH;
                        for (int k = 0; k < BATCH; k++) sum[k] += fabs(x - ref[k]);
                    }
                    for (int k = 0; k < n_batch; k++) AT(slot[k].numeric, row, sample) = (mreal)(sum[k] / dim);
                }
            }

            #pragma omp parallel for schedule(dynamic, 1) num_threads(n_workers)
            for (int k = 0; k < n_batch; k++) {
                outcome[k] = run_benchmark(&slot[k], batch_sample[k]);
                outcome[k].replicate = replicate;
            }

            for (int k = 0; k < n_batch; k++) {
                const Outcome *o = &outcome[k];
                fprintf(progress, "%s,%d,%d,%d,%d,%.6f,%d,%d,%.6f,%d,%d,", model_label[model],
                        samples[o->sample].index, o->replicate, o->in_set, o->rank, o->pvalue, o->set_size,
                        o->converged, o->final_pvalue, dropped, samples[o->lowest].index);
                for (int at = 0; o->members[at] >= 0; at++)
                    fprintf(progress, at ? " %d" : "%d", samples[o->members[at]].index);
                fprintf(progress, "\n");
                free(outcome[k].members);
            }
            fflush(progress);
            session_finished += n_batch;
            double elapsed = omp_get_wtime() - session_started;
            fprintf(stderr, "%s run %d up to %s: %ld done this session, %.4f s each, %.1f h left in all\n",
                    model_label[model], replicate, samples[batch_sample[n_batch - 1]].name, session_finished,
                    elapsed / session_finished,
                    elapsed / session_finished * (session_remaining - session_finished) / 3600);
        }
    }

    if (slot_rows >= 0) for (int k = 0; k < BATCH; k++) df_free(&slot[k]);
    free(reference);
    free(kept);
    free(usable);
    free(name_buffer);
    free(model_name);
    response_cache_free(&cache);
}

/* The progress file, compressed whole into the result and then removed. */
static void write_result(void) {
    FILE *f = fopen(PROGRESS_PATH, "rb");
    assert(f && "sweep_grid: the progress file is missing");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *text = malloc((size_t)length);
    assert(text && "sweep_grid: out of memory reading the progress file");
    size_t n_read = fread(text, 1, (size_t)length, f);
    assert(n_read == (size_t)length && "sweep_grid: short read of the progress file");
    fclose(f);

    size_t compressed_length;
    unsigned char *compressed = gzip_deflate_level(text, (size_t)length, 9, &compressed_length);
    FILE *out = fopen(RESULT_PATH, "wb");
    assert(out && "sweep_grid: cannot open the result for writing");
    size_t n_written = fwrite(compressed, 1, compressed_length, out);
    assert(n_written == compressed_length && "sweep_grid: short write of the result");
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
    const char *workers_text = getenv("SWEEP_GRID_WORKERS");
    int n_workers = workers_text ? atoi(workers_text) : 8;
    assert(n_workers >= 1 && "sweep_grid: SWEEP_GRID_WORKERS must be at least 1");
    openblas_set_num_threads(1);
    /* A confidence set inside a worker runs on that worker's thread alone. */
    omp_set_max_active_levels(1);

    struct stat st;
    int progress_exists = stat(PROGRESS_PATH, &st) == 0;
    if (stat(RESULT_PATH, &st) == 0 && !progress_exists) {
        fprintf(stderr, "%s already written, nothing to do\n", RESULT_PATH);
        return 0;
    }

    samples = lp_list_samples(LP_FIT_DIR, &n_samples);
    assert(n_samples > 0 && "sweep_grid: no configurations in the LP fit cache");
    n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);
    for (int i = 0; i < n_samples; i++) {
        char path[640];
        snprintf(path, sizeof path, "%s/%s", RESPONSE_QVARMA_FIT_DIR, samples[i].name);
        int found = stat(path, &st) == 0;
        assert(found && "sweep_grid: a configuration has LP fits and no t-QVARMA fits");
        (void)found;
    }

    long n_done;
    unsigned char *done = read_done(&n_done);
    fprintf(stderr, "%ld benchmarks already in %s\n", n_done, PROGRESS_PATH);
    FILE *progress = fopen(PROGRESS_PATH, "a");
    assert(progress && "sweep_grid: cannot open the progress file");
    if (!progress_exists) {
        fputs(HEADER, progress);
        fflush(progress);
    }

    /* Runs are taken a block at a time, every model through the block before
       the next block starts, so a run stopped early holds every model at
       about the same number of runs. A model's cache is read back from disk
       at each block after the first. */
    const char *block_text = getenv("SWEEP_GRID_RUN_BLOCK");
    int run_block = block_text ? atoi(block_text) : n_replicates;
    assert(run_block >= 1 && "sweep_grid: SWEEP_GRID_RUN_BLOCK must be at least 1");
    for (int model = FIRST_MODEL; model < N_MODELS; model++)
        session_remaining += count_remaining(model, 0, n_replicates, done);
    session_started = omp_get_wtime();
    for (int run_from = 0; run_from < n_replicates; run_from += run_block) {
        int run_to = run_from + run_block < n_replicates ? run_from + run_block : n_replicates;
        for (int model = FIRST_MODEL; model < N_MODELS; model++)
            sweep_model(model, run_from, run_to, done, n_workers, progress);
    }
    fclose(progress);

    write_result();

    free(done);
    lp_free_samples(samples, n_samples);
    return 0;
}
