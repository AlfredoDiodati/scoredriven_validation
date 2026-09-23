/*
The impulse-response loss matrix of the Monte Carlo experiment: how far each
ABM configuration's simulated dynamics sit from one chosen simulated run's,
measured through the auxiliary model.

This is applications/abm_system_irf_loss.c's procedure with the benchmark
changed. There the benchmark is the US data and its fit; here it is one
simulated replicate and its own cached fit, so the right answer is known - the
configuration that replicate was drawn from - and the confidence set can be
asked whether it returns it. docs/MONTECARLO_VALIDATION.md is the write-up and
says why a copy rather than a shared source: the two pipelines are independent,
and applications/ is not touched by anything here.

Nothing is estimated. Every fit this reads, the benchmark's included, is one of
the million applications/abm_system_fit_qvarma.c already wrote to
out/abm_system_fit_qvarma/.

The benchmark is whichever run montecarlo/benchmark_choice.c picked, read back
through montecarlo/benchmark.h. Its replicate index is held out of every
configuration's column, not only its own: replicate n of every configuration is
seed n+1 of the model, so the row is every simulation sharing the benchmark's
seed. A benchmark left in its own validation set would score against itself at
a loss of zero.

Loss is the mean absolute error between the benchmark's stacked impulse
response vector and each simulated fit's own, horizon 0 first, absolute rather
than squared for the reason applications/abm_system_irf_loss.c gives: one
divergent fit's IRF would otherwise dominate both the mean and the bootstrap
variance the confidence set estimates from.

nu <= 2 fits are skipped rather than fed through impulse_responses, whose own
assert would abort on them, and counted as missing.

Requires out/abm_system_fit_qvarma/ and dataset/abm_system/, and rebuilds
neither. Writes montecarlo/out/irf_loss.csv and its manifest. Nothing printed.
*/

#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include "montecarlo/benchmark.h"
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <et_al./frame/join.h>
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
/* Replications per configuration, counted off the fit cache rather than fixed
   here: the .Rdata dataset had 108 and the design experiment has 1000. */
static int n_replicates = 0;
#define HORIZON 20
#define IRF_DIM (K * K * (HORIZON + 1))

#define FIT_DIR_DEFAULT "out/abm_system_fit_qvarma"
#define INPUT_DIR_DEFAULT "dataset/abm_system"
#define OUTPUT_PATH_DEFAULT "montecarlo/out/irf_loss.csv"

/* Overridable so the whole pass can be run over a couple of configurations,
   which is how the parallel loop was checked against the serial one. */
static const char *FIT_DIR;
static const char *INPUT_DIR;
static const char *OUTPUT_PATH;

/*
The benchmark every simulated replicate is measured against. Unset, it is the
US fit named in spec_list and the real data's own block, which is what the main
pipeline runs. Set, it is one simulated replicate's cached fit and that
replicate's own series, which is what montecarlo/ runs: the same procedure with
a known answer substituted for the real data, to see whether the confidence set
recovers the configuration the benchmark came from.

ABM_SYSTEM_BENCHMARK_REPLICATE is then also the replicate index dropped from
every column, benchmark's own configuration included. Replicate n of every
configuration is seed n+1 of the model (applications/abm_system_simulate.c's
own "seed = mc + 1", the same seeds for every configuration), so dropping the
row is exactly dropping every simulation that shares the benchmark's seed. A
benchmark left in its own validation set would be scored against itself at a
loss of zero.
*/
static const char *BENCHMARK_SAMPLE;
static int BENCHMARK_REPLICATE = -1;
#define MODEL_LABEL "qvarma"

/* label names the cache files already on disk (save_fit's own naming,
   applications/abm_system_fit_qvarma.c); column_label is what a column of
   this file's own output is called, distinct from label so it cannot
   name the spec a column was fitted under. */
typedef struct { int r; const char *label; const char *column_label; const char *real_fit_path; } Spec;
static const Spec spec_list[] = {
    { 2, "p1q1r2", MODEL_LABEL "_p1q1r2", "out/us_qvarma_spec_choice_p1q1r2_fit.json" }
};
#define N_SPECS ((int)(sizeof spec_list / sizeof spec_list[0]))

static QvarmaParams spec_shape(int r) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, r, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

/* Stacks total[0..horizon], each K x K, into one flat column vector,
   horizon 0 first - the "vectorize and stack across horizons" step.
   Indexed by the K macro (ABM_SYSTEM_K), not r->K - same value always,
   since every QvarmaParams impulse_responses is ever called on here has that
   same K, but "r->K" itself would preprocess to "r->ABM_SYSTEM_K" and fail
   to compile, K being a macro in this file the same as every other spec
   constant. */
static Vec flatten_total_irf(const QvarmaImpulseResponses *r) {
    Vec out = mat_new(K * K * (r->horizon + 1), 1);
    int at = 0;
    for (int h = 0; h <= r->horizon; h++)
        for (int i = 0; i < K * K; i++)
            out.d[at++] = r->total[h].d[i];
    return out;
}

/* Fills out and returns 1 on success; returns 0 (out untouched) when m's
   own nu <= 2, which impulse_responses' own assert would otherwise abort
   on - see this file's own header comment on why that is checked here
   instead. */
static int try_compute_irf(const QvarmaParams *m, Mat y, Vec *out) {
    if (!(m->nu > 2)) return 0;
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses r = qvarma_impulse_responses(m, D, options);
    *out = flatten_total_irf(&r);
    qvarma_impulse_responses_free(&r);
    mat_free(D);
    return 1;
}


/* How many replicates the fit cache holds for one sample, so the table's shape
   comes from the data rather than from a constant that has to be remembered. */
static int count_replicates(const char *sample, const char *label) {
    char dir[512];
    snprintf(dir, sizeof dir, "%s/%s", FIT_DIR, sample);

    DIR *handle = opendir(dir);
    assert(handle && "montecarlo/irf_loss: cannot open a sample's fit directory");

    char suffix[64];
    snprintf(suffix, sizeof suffix, "_%s_fit.json", label);

    int n = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL)
        if (strstr(entry->d_name, suffix)) n++;
    closedir(handle);

    assert(n > 0 && "montecarlo/irf_loss: a sample directory holds no fit for this spec");
    return n;
}

typedef struct { char *name; int index; } SampleEntry;

/* dataset/abm_system/EstimationSeriesSample1_<N>'s own <N>, so column order
   is 1, 2, ..., 100 rather than readdir's arbitrary order or a string
   sort's "_10" before "_2". */
static int trailing_index(const char *name) {
    const char *underscore = strrchr(name, '_');
    assert(underscore && "montecarlo/irf_loss: a sample directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "montecarlo/irf_loss: cannot open out/abm_system_fit_qvarma/ - run abm_system_fit_qvarma first");

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
            assert(grown && "montecarlo/irf_loss: out of memory listing samples");
            entries = grown;
        }
        size_t len = strlen(entry->d_name);
        entries[n].name = malloc(len + 1);
        assert(entries[n].name && "montecarlo/irf_loss: out of memory copying a sample name");
        memcpy(entries[n].name, entry->d_name, len + 1);
        entries[n].index = trailing_index(entry->d_name);
        n++;
    }
    closedir(handle);
    qsort(entries, (size_t)n, sizeof(SampleEntry), compare_sample_entries);
    *count = n;
    return entries;
}

/* One replicate's five series, out of the compressed archive holding it.
   abm_system.h's own reader is what applications/abm_system_fit_qvarma.c and
   applications/abm_system_simulate.c go through, so the three cannot disagree
   about the layout. An earlier version of this file read a per-replicate CSV,
   which is not what the dataset has ever contained. */
static void read_y(const char *sample, int replicate, Mat *y_out) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    Mat y = abm_system_read_replicate(dir, replicate);
    assert(y.r == K && "montecarlo/irf_loss: a replicate has the wrong number of series");
    *y_out = y;
}

/* One spec's loss table: "replicate" plus one suffixed column per sample -
   the suffix is what keeps this spec's columns from colliding with the
   other spec's once both are joined. real_y is applications/
   us_qvarma_spec_choice.c's own real-data block, shared across both
   specs since the data does not depend on which spec is being fit. */
static DataFrame build_spec_losses(Spec spec, Mat real_y, const SampleEntry *samples, int n_samples,
                                   FILE *manifest) {
    QvarmaParams real = spec_shape(spec.r);
    char benchmark_path[640];
    snprintf(benchmark_path, sizeof benchmark_path, "%s/%s/replicate_%03d_%s_fit.json",
             FIT_DIR, BENCHMARK_SAMPLE, BENCHMARK_REPLICATE, spec.label);
    int loaded = qvarma_load_params(&real, benchmark_path);
    assert(loaded && "montecarlo/irf_loss: could not load the benchmark fit - run "
                      "us_qvarma_spec_choice.c's grid, or abm_system_fit_qvarma, first");
    Vec real_irf;
    int real_ok = try_compute_irf(&real, real_y, &real_irf);
    assert(real_ok && "montecarlo/irf_loss: the real-data fit's own nu <= 2 - "
                       "its impulse response function cannot be computed at all");

    Mat values = mat_new(n_replicates, n_samples + 1);
    int n_missing = 0;

    /* One task per replicate, so the thousand rows spread over the machine's
       threads and every column of a row is computed by the thread that owns it.
       Each thread carries its own parameter block because qvarma_load_params
       writes into it; values is written one distinct cell at a time, and the
       manifest is the only shared sink, taken under a lock because a missing
       cell is rare enough that contention for it never arises. */
    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams working = spec_shape(spec.r);

        #pragma omp for schedule(dynamic)
        for (int row = 0; row < n_replicates; row++) {
            AT(values, row, 0) = (mreal)row;
            /* The benchmark's own seed, dropped from every configuration rather
               than computed and thrown away, since nothing downstream reads it. */
            if (row == BENCHMARK_REPLICATE) {
                for (int col = 0; col < n_samples; col++) AT(values, row, col + 1) = (mreal)NAN;
                continue;
            }
            for (int col = 0; col < n_samples; col++) {
                char cache_path[560];
                snprintf(cache_path, sizeof cache_path, "%s/%s/replicate_%03d_%s_fit.json",
                         FIT_DIR, samples[col].name, row, spec.label);
                Mat sim_y;
                Vec sim_irf;
                int ok = qvarma_load_params(&working, cache_path);
                if (ok) {
                    read_y(samples[col].name, row, &sim_y);
                    ok = try_compute_irf(&working, sim_y, &sim_irf);
                    mat_free(sim_y);
                }
                if (ok) {
                    mreal loss = stats_mae(real_irf, sim_irf);
                    mat_free(sim_irf);
                    /* An impulse response can come back with a non-finite entry
                       even from parameters that loaded and passed the nu test,
                       and one such cell is enough to make the whole confidence
                       set refuse to run: et_al's mcs asserts a finite loss
                       matrix, because a hole in one model's column loses every
                       comparison it takes part in and the p-values give no sign
                       of it. */
                    if (!MISNAN(loss) && !MISINF(loss)) {
                        AT(values, row, col + 1) = loss;
                        continue;
                    }
                }
                AT(values, row, col + 1) = (mreal)NAN;
                n_missing++;
                #pragma omp critical
                fprintf(manifest, "%s missing: %s replicate %03d\n", spec.label,
                        samples[col].name, row);
            }
        }

        qvarma_params_free(&working);
    }

    /* The confidence set needs a rectangle with no holes. et_al's mcs asserts
       that, because a model whose column carries one NaN loses every
       comparison it takes part in and the p-values report nothing unusual: a
       hole eliminates a model rather than being visible as a hole. A cell can
       only be missing here when a fit would not load or its impulse response
       came back non-finite, and both are properties of one (configuration,
       replicate) pair rather than of a whole configuration. Dropping the
       replicate keeps every configuration in the comparison and costs one
       observation from all of them, which is the direction that does not
       silently change which models are being compared. */
    int keep = 0;
    int *usable = (int*)malloc((size_t)n_replicates * sizeof(int));
    for (int row = 0; row < n_replicates; row++) {
        usable[row] = 1;
        for (int col = 0; col < n_samples; col++)
            if (MISNAN(AT(values, row, col + 1))) { usable[row] = 0; break; }
        if (usable[row]) keep++;
    }
    assert(keep > 0 && "montecarlo/irf_loss: every replicate has a missing cell");

    if (keep < n_replicates) {
        Mat kept = mat_new(keep, n_samples + 1);
        int at = 0;
        for (int row = 0; row < n_replicates; row++) {
            if (!usable[row]) continue;
            for (int col = 0; col <= n_samples; col++) AT(kept, at, col) = AT(values, row, col);
            at++;
        }
        mat_free(values);
        values = kept;
        int held_out = 1;
        fprintf(manifest, "%s: %d of %d replicates dropped for holding a missing cell\n",
                spec.label, n_replicates - keep - held_out, n_replicates);
        fprintf(manifest, "%s: replicate %03d held out of every column, the benchmark's own "
                               "seed\n", spec.label, BENCHMARK_REPLICATE);
    }
    free(usable);


    char **col_names = (char**)malloc((size_t)(n_samples + 1) * sizeof(char*));
    col_names[0] = frame_strdup("replicate");
    for (int col = 0; col < n_samples; col++) {
        char buf[128];
        snprintf(buf, sizeof buf, "%s_%s", samples[col].name, spec.column_label);
        col_names[col + 1] = frame_strdup(buf);
    }
    DataFrame df = df_from_matrix(values, (const char *const *)col_names);
    for (int col = 0; col <= n_samples; col++) free(col_names[col]);
    free(col_names);
    mat_free(values);

    /* Out of the cells actually computed, which is one row short of the table
       when a benchmark replicate is held out of every column. */
    int computed_rows = n_replicates - 1;
    fprintf(manifest, "%s: %d of %d cells missing\n\n", spec.label, n_missing,
            computed_rows * n_samples);

    qvarma_params_free(&real);
    mat_free(real_irf);
    return df;
}

int main(void) {
    FIT_DIR = getenv("ABM_SYSTEM_FIT_DIR");
    if (!FIT_DIR) FIT_DIR = FIT_DIR_DEFAULT;
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = INPUT_DIR_DEFAULT;
    OUTPUT_PATH = getenv("ABM_SYSTEM_LOSS_PATH");
    if (!OUTPUT_PATH) OUTPUT_PATH = OUTPUT_PATH_DEFAULT;
    Benchmark benchmark = benchmark_read();
    BENCHMARK_SAMPLE = benchmark.sample;
    BENCHMARK_REPLICATE = benchmark.replicate;

    int n_samples;
    SampleEntry *samples = list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0 && "montecarlo/irf_loss: no sample directories under out/abm_system_fit_qvarma/");

    /* The benchmark's own series: the US block, or the simulated replicate that
       stands in for it. Both are K x T in the same row order, which is what
       lets one procedure run against either. */
    Mat real_y;
    read_y(BENCHMARK_SAMPLE, BENCHMARK_REPLICATE, &real_y);

    n_replicates = count_replicates(samples[0].name, spec_list[0].label);

    char manifest_path[640];
    /* Named after the loss table with its extension dropped, so the pair reads
       as <stem>.csv and <stem>_manifest.txt rather than as a path with two
       extensions stuck together. */
    size_t stem = strlen(OUTPUT_PATH);
    if (stem > 4 && strcmp(OUTPUT_PATH + stem - 4, ".csv") == 0) stem -= 4;
    snprintf(manifest_path, sizeof manifest_path, "%.*s_manifest.txt", (int)stem, OUTPUT_PATH);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "montecarlo/irf_loss: cannot open the manifest path for writing");
    fprintf(manifest, "%d samples, %d replicates each, %d spec, "
                       "loss = MAE between stacked impulse response vectors (horizon %d)\n",
            n_samples, n_replicates, N_SPECS, HORIZON);
    fprintf(manifest, "benchmark: %s replicate %03d, its own cached fit and its own series, "
                       "and that replicate is held out of every column\n\n",
            BENCHMARK_SAMPLE, BENCHMARK_REPLICATE);

    /* One table per spec, and with one spec there is nothing to join it to.
       The join that stood here paired the r = 2 and r = 4 columns on the
       replicate index; if a second spec goes back into spec_list, it comes
       back with it. */
    assert(N_SPECS == 1 && "montecarlo/irf_loss: written for one spec - restore the join for more");
    DataFrame losses = build_spec_losses(spec_list[0], real_y, samples, n_samples, manifest);
    mat_free(real_y);

    assert(losses.r > 0 && "montecarlo/irf_loss: the loss table has no rows");
    df_write_csv(&losses, OUTPUT_PATH, csv_write_options_default());

    df_free(&losses);

    fclose(manifest);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
