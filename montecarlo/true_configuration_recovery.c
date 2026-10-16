/*
Does the Model Confidence Set return the configuration a benchmark was
simulated from? Every simulated run of every configuration stands in as the
benchmark in turn: 1000 configurations times 1000 runs, for the t-QVARMA and
the four local projections (linear, state 1, state 2, both states in one
vector), five million confidence sets.

Each benchmark is the procedure of docs/MCS_VALIDATION_PROCEDURE.md against
that run: the squared error between response vectors minus each
configuration's mean response variance (montecarlo/recovery_loss.h), the
benchmark's run index held out of every column and of every variance, a run
with a missing response in any column dropped from both, MCS_TR with block
length 1, bootstrap variance, seed 123 stream 0, at 2000 resamples.

The order is by run first: run 0 of every configuration, then run 1, and so on,
so a run stopped early holds every configuration at the same number of
benchmarks. TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK sets how many runs one model
does before the next model takes the same runs (default all of them, one model
at a time, t-QVARMA first); a smaller block keeps the models level.

Every hardware thread is used, counted at run time with omp_get_num_procs.
Per run index the configurations' variances are computed once, on every
thread. Benchmarks sharing the run index are then taken in batches whose size
is the smallest multiple of the thread count that is at least 64: the batch's
loss matrices are formed in one pass over the response cache on every thread,
then its confidence sets run side by side, one per thread.

Resumable at benchmark granularity. Rows go to a plain progress file as each
batch finishes, and a rerun skips what it holds. When every model is done the
rows are written to one gzip-compressed csv and the progress file is removed.
A rerun with the result present and no progress file decompresses the result
into the progress file, so a model added since is run and the rows already
there are kept.

Columns: model, benchmark configuration, benchmark run, whether the benchmark's
configuration is in the set, its rank by mean loss (1 is the smallest), its MCS
p-value, the set size, whether the set was decided by an accepted test, the
final p-value, runs dropped for a missing response, the configuration with the
smallest mean loss, and the configurations in the set, space separated.

Requires out/abm_system_fit_qvarma/, out/abm_system_fit_lp/ and
out/dataset/abm_system/, and rebuilds none of them. Progress goes to stderr. Each
model's response cache is read from out/sweep_grid_response_cache/, and built
there when missing.
*/

#include "montecarlo/recovery_loss.h"
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/frame.h>
#include <et_al./frame/gzip.h>
#include <cblas.h>

/* The rest can be overridden at compile time for a test build that
   covers only the first runs, models or configurations. */
#define MINIMUM_BATCH 64
#ifndef BOOTSTRAP
#define BOOTSTRAP 2000
#endif
#ifndef PROGRESS_PATH
#define PROGRESS_PATH "out/montecarlo/true_configuration_recovery_progress.csv"
#endif
#ifndef RESULT_PATH
#define RESULT_PATH "out/montecarlo/true_configuration_recovery.csv.gz"
#endif
#ifndef RESPONSE_CACHE_DIR
#define RESPONSE_CACHE_DIR RESPONSE_CACHE_DIR_DEFAULT
#endif
#ifndef FIRST_MODEL
#define FIRST_MODEL 0
#endif
#ifndef LAST_MODEL
#define LAST_MODEL (N_MODELS - 1)
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
static int n_threads = 0;

static size_t done_index(int model, int sample, int replicate) {
    return ((size_t)model * n_samples + sample) * n_replicates + replicate;
}

/* Which benchmarks the progress file already holds. */
static unsigned char *read_done(long *n_done) {
    unsigned char *done = calloc((size_t)N_MODELS * n_samples * n_replicates, 1);
    assert(done);
    int highest = 0;
    for (int i = 0; i < n_samples; i++) if (samples[i].index > highest) highest = samples[i].index;
    int *sample_of_cop = malloc((size_t)(highest + 1) * sizeof(int));
    assert(sample_of_cop);
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
   loss use stats_mean over each column. */
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
    assert(o.members);
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
static void recover_model(int model, int run_from, int run_to, const unsigned char *done, FILE *progress) {
    long remaining = count_remaining(model, run_from, run_to, done);
    fprintf(stderr, "%s runs %d to %d: %ld benchmarks to do\n", model_label[model], run_from, run_to - 1,
            remaining);
    if (!remaining) return;

    ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, LP_FIT_DIR,
                                                       INPUT_DIR, RESPONSE_CACHE_DIR);
    const int dim = cache.dim;

    /* A run with a missing response in any configuration is dropped from
       every loss matrix and every variance, since et_al's mcs refuses a
       matrix with a hole. */
    unsigned char *usable = malloc((size_t)n_replicates);
    assert(usable);
    int n_usable = 0;
    for (int replicate = 0; replicate < n_replicates; replicate++) {
        usable[replicate] = 1;
        for (int sample = 0; sample < n_samples && usable[replicate]; sample++)
            if (MISNAN(response_slot(&cache, sample, replicate)[0])) usable[replicate] = 0;
        n_usable += usable[replicate];
    }
    fprintf(stderr, "%s: %ld cells missing, %d of %d runs usable\n", model_label[model], cache.n_missing,
            n_usable, n_replicates);

    const char **model_name = malloc((size_t)n_samples * sizeof(char *));
    char (*name_buffer)[128] = malloc((size_t)n_samples * sizeof *name_buffer);
    assert(model_name && name_buffer);
    for (int i = 0; i < n_samples; i++) {
        snprintf(name_buffer[i], sizeof name_buffer[i], "%s_%s", samples[i].name, model_label[model]);
        model_name[i] = name_buffer[i];
    }

    /* The smallest multiple of the thread count that is at least
       MINIMUM_BATCH, so the confidence sets of a batch fill every thread. */
    const int batch = (MINIMUM_BATCH + n_threads - 1) / n_threads * n_threads;

    /* One loss matrix per slot, reused across batches; the benchmark
       responses of a batch entry by entry, converted to double once. */
    DataFrame *slot = malloc((size_t)batch * sizeof(DataFrame));
    Mat *slot_numeric = malloc((size_t)batch * sizeof(Mat));
    int slot_rows = -1;
    double *reference = calloc((size_t)dim * batch, sizeof(double));
    double *mean_variance = malloc((size_t)n_samples * sizeof(double));
    int *kept = malloc((size_t)n_replicates * sizeof(int));
    int *batch_sample = malloc((size_t)batch * sizeof(int));
    Outcome *outcome = malloc((size_t)batch * sizeof(Outcome));
    assert(slot && slot_numeric && reference && mean_variance && kept && batch_sample && outcome);

    for (int replicate = run_from; replicate < run_to && replicate < RUN_LIMIT; replicate++) {
        int n_kept = 0;
        for (int r = 0; r < n_replicates; r++) if (usable[r] && r != replicate) kept[n_kept++] = r;
        int dropped = n_replicates - 1 - n_kept;
        if (n_kept != slot_rows) {
            if (slot_rows >= 0) for (int k = 0; k < batch; k++) df_free(&slot[k]);
            Mat empty = mat_new(n_kept, n_samples);
            for (int k = 0; k < batch; k++) {
                slot[k] = df_from_matrix(empty, model_name);
                slot_numeric[k] = slot[k].numeric;
            }
            mat_free(empty);
            slot_rows = n_kept;
        }

        int any_left = 0;
        for (int sample = 0; sample < n_samples && !any_left; sample++)
            any_left = !done[done_index(model, sample, replicate)];
        if (!any_left) continue;
        recovery_loss_mean_variance(&cache, kept, n_kept, mean_variance, n_threads);

        int last = n_samples < BENCHMARK_LIMIT ? n_samples : BENCHMARK_LIMIT;
        for (int first = 0; first < last;) {
            int n_batch = 0;
            while (first < last && n_batch < batch) {
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
                for (int h = 0; h < dim; h++) reference[(size_t)h * batch + k] = (double)own[h];
            }
            recovery_loss_fill(&cache, kept, n_kept, reference, n_batch, batch, mean_variance, slot_numeric,
                               n_threads);

            #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
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

    if (slot_rows >= 0) for (int k = 0; k < batch; k++) df_free(&slot[k]);
    free(slot);
    free(slot_numeric);
    free(reference);
    free(mean_variance);
    free(kept);
    free(batch_sample);
    free(outcome);
    free(usable);
    free(name_buffer);
    free(model_name);
    response_cache_free(&cache);
}

/* The result back into the progress file, so a rerun after a model was added
   skips the rows the result holds, keeps them, and does the new model. */
static void restore_progress(void) {
    FILE *f = fopen(RESULT_PATH, "rb");
    assert(f && "true_configuration_recovery: cannot open the result");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc((size_t)length);
    assert(compressed && "true_configuration_recovery: out of memory reading the result");
    size_t n_read = fread(compressed, 1, (size_t)length, f);
    assert(n_read == (size_t)length && "true_configuration_recovery: short read of the result");
    (void)n_read;
    fclose(f);

    size_t text_length;
    unsigned char *text = gzip_inflate(compressed, (size_t)length, &text_length);
    FILE *out = fopen(PROGRESS_PATH, "wb");
    assert(out && "true_configuration_recovery: cannot open the progress file for writing");
    size_t n_written = fwrite(text, 1, text_length, out);
    assert(n_written == text_length && "true_configuration_recovery: short write of the progress file");
    (void)n_written;
    fclose(out);
    fprintf(stderr, "%s restored to %s, %zu bytes\n", RESULT_PATH, PROGRESS_PATH, text_length);
    free(text);
    free(compressed);
}

/* The progress file, compressed whole into the result and then removed. */
static void write_result(void) {
    FILE *f = fopen(PROGRESS_PATH, "rb");
    assert(f && "true_configuration_recovery: the progress file is missing");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *text = malloc((size_t)length);
    assert(text && "true_configuration_recovery: out of memory reading the progress file");
    size_t n_read = fread(text, 1, (size_t)length, f);
    assert(n_read == (size_t)length && "true_configuration_recovery: short read of the progress file");
    (void)n_read;
    fclose(f);

    size_t compressed_length;
    unsigned char *compressed = gzip_deflate_level(text, (size_t)length, 9, &compressed_length);
    FILE *out = fopen(RESULT_PATH, "wb");
    assert(out && "true_configuration_recovery: cannot open the result for writing");
    size_t n_written = fwrite(compressed, 1, compressed_length, out);
    assert(n_written == compressed_length && "true_configuration_recovery: short write of the result");
    (void)n_written;
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
    n_threads = omp_get_num_procs();
    /* Every parallel loop here takes every hardware thread, so OpenBLAS runs
       at one thread inside each job and a confidence set inside a worker runs
       on that worker's thread alone. */
    openblas_set_num_threads(1);
    omp_set_max_active_levels(1);
    fprintf(stderr, "%d hardware threads\n", n_threads);

    struct stat st;
    int progress_exists = stat(PROGRESS_PATH, &st) == 0;
    int restored = stat(RESULT_PATH, &st) == 0 && !progress_exists;

    samples = lp_list_samples(LP_FIT_DIR, &n_samples);
    assert(n_samples > 0 && "true_configuration_recovery: no configurations in the LP fit cache");
    n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);
    for (int i = 0; i < n_samples; i++) {
        char path[640];
        snprintf(path, sizeof path, "%s/%s", RESPONSE_QVARMA_FIT_DIR, samples[i].name);
        int found = stat(path, &st) == 0;
        assert(found && "true_configuration_recovery: a configuration has LP fits and no t-QVARMA fits");
        (void)found;
    }

    if (restored) restore_progress();
    long n_done;
    unsigned char *done = read_done(&n_done);
    fprintf(stderr, "%ld benchmarks already in %s\n", n_done, PROGRESS_PATH);
    if (restored && n_done == (long)N_MODELS * n_samples * n_replicates) {
        remove(PROGRESS_PATH);
        fprintf(stderr, "%s holds every model, nothing to do\n", RESULT_PATH);
        free(done);
        lp_free_samples(samples, n_samples);
        return 0;
    }
    mkdir("out/montecarlo", 0755);
    FILE *progress = fopen(PROGRESS_PATH, "a");
    assert(progress && "true_configuration_recovery: cannot open the progress file");
    if (!progress_exists && !restored) {
        fputs(HEADER, progress);
        fflush(progress);
    }

    /* Runs are taken a block at a time, every model through the block before
       the next block starts, so a run stopped early holds every model at
       about the same number of runs. A model's cache is read back from disk
       at each block after the first. */
    const char *block_text = getenv("TRUE_CONFIGURATION_RECOVERY_RUN_BLOCK");
    int run_block = block_text ? atoi(block_text) : n_replicates;
    assert(run_block >= 1 && "true_configuration_recovery: the run block must be at least 1");
    for (int model = FIRST_MODEL; model <= LAST_MODEL; model++)
        session_remaining += count_remaining(model, 0, n_replicates, done);
    session_started = omp_get_wtime();
    for (int run_from = 0; run_from < n_replicates; run_from += run_block) {
        int run_to = run_from + run_block < n_replicates ? run_from + run_block : n_replicates;
        for (int model = FIRST_MODEL; model <= LAST_MODEL; model++)
            recover_model(model, run_from, run_to, done, progress);
    }
    fclose(progress);

    /* A test build over part of the grid leaves its rows in the progress
       file rather than writing a result that looks complete. */
    int whole_grid = FIRST_MODEL == 0 && LAST_MODEL == N_MODELS - 1 && RUN_LIMIT >= n_replicates
                  && BENCHMARK_LIMIT >= n_samples;
    if (whole_grid) write_result();
    else fprintf(stderr, "a test build over part of the grid; its rows stay in %s\n", PROGRESS_PATH);

    free(done);
    lp_free_samples(samples, n_samples);
    return 0;
}
