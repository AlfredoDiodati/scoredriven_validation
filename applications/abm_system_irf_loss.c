/*
The loss matrix the Model Confidence Set runs on: how far each ABM
configuration's simulated dynamics sit from the real US data's, measured
through the auxiliary model.

The comparison object is the impulse response function, not the fitted
parameters themselves. That is the point of the whole file, and it changed
from an earlier version that compared constrained parameter vectors directly.
Parameter distance matched neither the procedure this project replicates (a
thesis using the same Model Confidence Set protocol over local-projection
IRFs - "Evaluating Nonlinear Simulation Models with Model Confidence Sets",
Fabiano, Pisa/Sant'Anna) nor what a loss between two fitted models is supposed
to measure: two QVARMA fits with different-looking coefficients can imply
nearly identical dynamics, and two with similar-looking coefficients can imply
very different ones. The IRF is the object whose distance answers "do these
two models behave alike", which raw parameter distance does not.

For the auxiliary spec p1q1r2:
  1. Fit the real-data QVARMA (out/us_qvarma_spec_choice_p1q1r2_fit.json,
     applications/us_qvarma_spec_choice.c's own grid) and compute its
     own impulse response function via qvarma.h's impulse_responses.
  2. For every simulated fit
     applications/abm_system_fit_qvarma.c already wrote to
     out/abm_system_fit_qvarma/, compute that fit's own impulse response
     function the same way.
  3. Stack every horizon's K x K response matrix (qvarma.h's own "total",
     contemporaneous + stationary + cointegrated) into one flat vector per
     model, horizon 0 first - the same "vectorize and stack across
     horizons" step the thesis's own protocol uses (its Sec 3.3, step 3).
  4. Loss is the mean absolute error between the real model's IRF vector
     and each simulated model's IRF vector, absolute rather than squared
     because a squared difference lets one badly-behaved fit's IRF
     dominate both the mean and the bootstrap variance abm_system_mcs.c
     estimates from, which is exactly what left every model
     indistinguishable from every other one before this file used absolute
     error at the parameter-vector stage - nothing about switching the
     comparison object to IRFs removes that risk, since a divergent fit
     produces a divergent IRF too.

nu <= 2 fits are skipped, not fed through impulse_responses: qvarma.h's own
impulse_responses asserts m->nu > 2 (the multivariate-t degrees of freedom
enters the impulse formula as 1/(nu-2)), which a badly non-converged
optimizer run can and does produce. The convention this project follows is
that an infeasible parameter value from an optimizer probing the space
returns a sentinel rather than aborting - a cached fit already written to
disk is exactly that case, arrived at after the fact rather than during the
fit itself, so it is checked and skipped here (counted as missing, same as
an unreadable cache file) rather than left to the library's own assert.

Requires out/abm_system_fit_qvarma/ to already hold every replicate's fit (run
applications/abm_system_fit_qvarma.c to completion first),
out/us_qvarma_spec_choice_p1q1r2_fit.json to exist
(applications/us_qvarma_spec_choice.c, run separately) and out/us_system.csv to
exist (applications/us_prepare_data.c), which load_us_system reads directly to
rebuild the real data's own K x ESTIMATION_PERIODS block the same way
applications/us_qvarma_spec_choice.c itself does.

Output: out/abm_system_irf_loss.csv, the loss matrix
applications/abm_system_mcs.c reads. One row per replicate and one column per
configuration, plus a leading "replicate" column, so over the design experiment
that is 1000 rows and 1001 columns. Each column is named for its configuration
and the spec it was fitted under, and every numeric column except "replicate"
is a model to the confidence set.

The matrix has no holes, which is a requirement rather than a nicety: et_al's
mcs asserts a finite loss matrix because one NaN in a model's column loses
every comparison that model takes part in, and the p-values that come back look
ordinary. A cell can go missing when a fit will not load, when its nu is not
above two, or when its impulse response comes back non-finite. Where that
happens the replicate is dropped from every column rather than the
configuration from every row, since a missing cell belongs to one
(configuration, replicate) pair and dropping the row is the direction that does
not change which models are being compared.
The manifest is named after the loss table rather than after this script, so
the two cannot be separated: the default is
out/abm_system_irf_loss_manifest.txt, and pointing
ABM_SYSTEM_LOSS_PATH somewhere else moves both. It names every missing cell and
counts the replicates dropped.

The fits are used as they stand, converged or not. About a third of them
converged and the rest stopped because the line search could not move;
docs/ABM_SYSTEM_MCS_VALIDATION.md records the measurement that their
log-likelihoods are not distinguishable from the converged ones.

Not part of make applications - buildable on its own via
make app-abm_system_irf_loss, same as abm_system_fit_qvarma.c itself.
Nothing printed.
*/

#include "abm_system.h"
#include <et_al./sd/qvarma.h>
#include "us_data.h"
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
#define OUTPUT_PATH_DEFAULT "out/abm_system_irf_loss.csv"

/* Overridable so the whole pass can be run over a couple of configurations,
   which is how the parallel loop was checked against the serial one. */
static const char *FIT_DIR;
static const char *INPUT_DIR;
static const char *OUTPUT_PATH;
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

/* Real data, applications/us_qvarma_spec_choice.c's own build_block:
   growth/change of GDP, energy demand and employment, inflation, and the
   interest rate in levels - identical row convention to what
   dataset/abm_system/'s own simulated series already use (both feed the
   same K=5, same row order, into the same fit code, and ROW_GDP_GROWTH
   etc. below are abm_system.h's own enum, not a fresh one), which is what
   makes comparing their IRFs meaningful at all. */
#define REAL_PERIODS (ESTIMATION_PERIODS - 1)
static Mat build_real_block(Mat original) {
    Mat y = mat_new(K, REAL_PERIODS);
    for (int t = 1; t < ESTIMATION_PERIODS; t++) {
        int c = t - 1;
        AT(y, ROW_GDP_GROWTH, c) = AT(original, LOG_GDP, t) - AT(original, LOG_GDP, t - 1);
        AT(y, ROW_EN_GROWTH, c) = AT(original, LOG_ENERGY_DEMAND, t)
                                 - AT(original, LOG_ENERGY_DEMAND, t - 1);
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = AT(original, EMPLOYMENT, t) - AT(original, EMPLOYMENT, t - 1);
        AT(y, ROW_INFLATION, c) = AT(original, LOG_CPI, t) - AT(original, LOG_CPI, t - 1);
        AT(y, ROW_INTEREST_RATE, c) = AT(original, INTEREST_RATE, t);
    }
    return y;
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
    assert(handle && "abm_system_irf_loss: cannot open a sample's fit directory");

    char suffix[64];
    snprintf(suffix, sizeof suffix, "_%s_fit.json", label);

    int n = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL)
        if (strstr(entry->d_name, suffix)) n++;
    closedir(handle);

    assert(n > 0 && "abm_system_irf_loss: a sample directory holds no fit for this spec");
    return n;
}

typedef struct { char *name; int index; } SampleEntry;

/* dataset/abm_system/EstimationSeriesSample1_<N>'s own <N>, so column order
   is 1, 2, ..., 100 rather than readdir's arbitrary order or a string
   sort's "_10" before "_2". */
static int trailing_index(const char *name) {
    const char *underscore = strrchr(name, '_');
    assert(underscore && "abm_system_irf_loss: a sample directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_irf_loss: cannot open out/abm_system_fit_qvarma/ - run abm_system_fit_qvarma first");

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
            assert(grown && "abm_system_irf_loss: out of memory listing samples");
            entries = grown;
        }
        size_t len = strlen(entry->d_name);
        entries[n].name = malloc(len + 1);
        assert(entries[n].name && "abm_system_irf_loss: out of memory copying a sample name");
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
    assert(y.r == K && "abm_system_irf_loss: a replicate has the wrong number of series");
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
    int loaded = qvarma_load_params(&real, spec.real_fit_path);
    assert(loaded && "abm_system_irf_loss: could not load the real-data fit - run "
                      "us_qvarma_spec_choice.c's grid first");
    Vec real_irf;
    int real_ok = try_compute_irf(&real, real_y, &real_irf);
    assert(real_ok && "abm_system_irf_loss: the real-data fit's own nu <= 2 - "
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
    assert(keep > 0 && "abm_system_irf_loss: every replicate has a missing cell");

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
        fprintf(manifest, "%s: %d of %d replicates dropped for holding a missing cell\n",
                spec.label, n_replicates - keep, n_replicates);
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

    fprintf(manifest, "%s: %d of %d cells missing\n\n", spec.label, n_missing,
            n_replicates * n_samples);

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

    int n_samples;
    SampleEntry *samples = list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0 && "abm_system_irf_loss: no sample directories under out/abm_system_fit_qvarma/");

    Mat original = load_us_system();
    Mat real_y = build_real_block(original);
    mat_free(original);

    n_replicates = count_replicates(samples[0].name, spec_list[0].label);

    char manifest_path[640];
    /* Named after the loss table with its extension dropped, so the pair reads
       as <stem>.csv and <stem>_manifest.txt rather than as a path with two
       extensions stuck together. */
    size_t stem = strlen(OUTPUT_PATH);
    if (stem > 4 && strcmp(OUTPUT_PATH + stem - 4, ".csv") == 0) stem -= 4;
    snprintf(manifest_path, sizeof manifest_path, "%.*s_manifest.txt", (int)stem, OUTPUT_PATH);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "abm_system_irf_loss: cannot open the manifest path for writing");
    fprintf(manifest, "%d samples, %d replicates each, %d spec, "
                       "loss = MAE between stacked impulse response vectors (horizon %d)\n\n",
            n_samples, n_replicates, N_SPECS, HORIZON);

    /* One table per spec, and with one spec there is nothing to join it to.
       The join that stood here paired the r = 2 and r = 4 columns on the
       replicate index; if a second spec goes back into spec_list, it comes
       back with it. */
    assert(N_SPECS == 1 && "abm_system_irf_loss: written for one spec - restore the join for more");
    DataFrame losses = build_spec_losses(spec_list[0], real_y, samples, n_samples, manifest);
    mat_free(real_y);

    assert(losses.r > 0 && "abm_system_irf_loss: the loss table has no rows");
    df_write_csv(&losses, OUTPUT_PATH, csv_write_options_default());

    df_free(&losses);

    fclose(manifest);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
