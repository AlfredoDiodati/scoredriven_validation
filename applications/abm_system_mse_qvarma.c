/*
The driftless t-QVARMA half of the joint loss table - see abm_system_mse.c
for the drift model's own half and for why this is two files rather than
one (qvarma.h and qvarma_d.h cannot both be included in the same
translation unit, being two independent models that happen to define
identically-named types - docs/MODEL_TEMPLATE.md's own entry 16). This file
never includes qvarma_d.h and never will.

The comparison object is the impulse response function, not the fitted
parameters themselves - this is the point of the whole file, and it changed
from an earlier version that compared flatten_estimated's constrained
parameter vectors directly. That matched neither the procedure this project
is actually replicating (a thesis using the same Model Confidence Set
protocol over local-projection IRFs - "Evaluating Nonlinear Simulation
Models with Model Confidence Sets", Fabiano, Pisa/Sant'Anna) nor what a
loss between two fitted models is supposed to measure: two QVARMA fits with
different-looking coefficients can still imply nearly identical dynamics,
and two with similar-looking coefficients can imply very different ones -
the IRF is the object whose distance actually answers "do these two models
behave alike", which raw parameter distance does not.

For each of the two auxiliary specs (p1q1r2, p1q1r4):
  1. Fit the real-data QVARMA (out/us_qvarma_employment_change_p1q1rN_fit.json,
     applications/us_qvarma_employment_change.c's own grid) and compute its
     own impulse response function via qvarma.h's impulse_responses.
  2. For every one of the 10,800 simulated fits
     applications/abm_system_fit_qvarma.c already wrote to
     out/abm_system_fit_qvarma/, compute that fit's own impulse response
     function the same way.
  3. Stack every horizon's K x K response matrix (qvarma.h's own "total",
     contemporaneous + stationary + cointegrated) into one flat vector per
     model, horizon 0 first - the same "vectorize and stack across
     horizons" step the thesis's own protocol uses (its Sec 3.3, step 3).
  4. Loss is the mean absolute error between the real model's IRF vector
     and each simulated model's IRF vector - absolute rather than squared
     for the same reason abm_system_mse.c uses stats_mae rather than
     stats_mse: a squared difference lets one badly-behaved fit's IRF
     dominate both the mean and the bootstrap variance abm_system_mcs.c
     estimates from, which is exactly what left every model
     indistinguishable from every other one before this file used absolute
     error at the parameter-vector stage - nothing about switching the
     comparison object to IRFs removes that risk, since a divergent fit
     produces a divergent IRF too.

nu <= 2 fits are skipped, not fed through impulse_responses: qvarma.h's own
impulse_responses asserts m->nu > 2 (the multivariate-t degrees of freedom
enters the impulse formula as 1/(nu-2)), which a badly non-converged
optimizer run can and does produce. This project's own convention is that
an infeasible parameter value from an optimizer probing the space returns a
sentinel rather than aborting (see docs/MODEL_TEMPLATE.md's "Implementing a
new model" policy) - a cached fit already written to disk is exactly that
case, arrived at after the fact rather than during the fit itself, so it is
checked and skipped here (counted as missing, same as an unreadable cache
file) rather than left to the library's own assert.

This file additionally reads abm_system_mse.c's own output
(out/abm_system_mse_qvarmad_joint.csv) - but no longer joins it into a
combined table, since that file's own loss is still parameter-vector MAE,
not IRF MAE, and joining two different loss definitions side by side would
mislabel them as comparable. See "Requires" below for what this file still
depends on abm_system_mse.c for.

Requires out/abm_system_fit_qvarma/ to already hold every replicate's fit
(run applications/abm_system_fit_qvarma.c to completion first), both
out/us_qvarma_employment_change_p1q1r2_fit.json and ..._p1q1r4_fit.json to
exist (applications/us_qvarma_employment_change.c's own grid, run
separately), and out/us_system.csv to exist (applications/us_prepare_data.c),
which load_us_system reads directly to rebuild the real data's own
K x ESTIMATION_PERIODS block the same way
applications/us_qvarma_employment_change.c itself does.

Output: out/abm_system_mse_qvarma_joint.csv, a 108 x 201 table - a leading
"replicate" column (0-107) followed by the 100 p1q1r2 columns then the 100
p1q1r4 columns, each labeled with qvarma and its own spec -
abm_system_mcs.c's own input. A cell is NaN when that (sample, replicate,
spec) cache is missing, unreadable, or has nu <= 2, rather than silently
skipped or aborting the whole run; out/abm_system_mse_qvarma_manifest.txt
records exactly which cells that happened for, if any, plus a per-spec
count.

Not part of make applications - buildable on its own via
make app-abm_system_mse_qvarma, same as abm_system_fit_qvarma.c itself.
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

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define P 1
#define Q 1
#define N_REPLICATES 108
#define HORIZON 20
#define IRF_DIM (K * K * (HORIZON + 1))

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define OUTPUT_PATH "out/abm_system_mse_qvarma_joint.csv"
#define MODEL_LABEL "qvarma"

/* label names the cache files already on disk (save_fit's own naming,
   applications/abm_system_fit_qvarma.c); column_label is what a column of
   this file's own output is called, distinct from label so it cannot
   collide with abm_system_mse.c's identically-shaped p1q1r2/p1q1r4
   columns should the two ever be compared side by side again. */
typedef struct { int r; const char *label; const char *column_label; const char *real_fit_path; } Spec;
static const Spec spec_list[] = {
    { 2, "p1q1r2", MODEL_LABEL "_p1q1r2", "out/us_qvarma_employment_change_p1q1r2_fit.json" },
    { 4, "p1q1r4", MODEL_LABEL "_p1q1r4", "out/us_qvarma_employment_change_p1q1r4_fit.json" }
};
#define N_SPECS ((int)(sizeof spec_list / sizeof spec_list[0]))

static QvarmaParams spec_shape(int r) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, r, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

/* Real data, applications/us_qvarma_employment_change.c's own build_block:
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

typedef struct { char *name; int index; } SampleEntry;

/* dataset/abm_system/EstimationSeriesSample1_<N>'s own <N>, so column order
   is 1, 2, ..., 100 rather than readdir's arbitrary order or a string
   sort's "_10" before "_2". */
static int trailing_index(const char *name) {
    const char *underscore = strrchr(name, '_');
    assert(underscore && "abm_system_mse_qvarma: a sample directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_mse_qvarma: cannot open out/abm_system_fit_qvarma/ - run abm_system_fit_qvarma first");

    SampleEntry *entries = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char path[560];
        snprintf(path, sizeof path, "%s/%s", dir, entry->d_name);
        struct stat st;
        if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) continue;
        if (n == cap) { cap = cap ? cap * 2 : 16; entries = realloc(entries, (size_t)cap * sizeof(SampleEntry)); }
        size_t len = strlen(entry->d_name);
        entries[n].name = malloc(len + 1);
        memcpy(entries[n].name, entry->d_name, len + 1);
        entries[n].index = trailing_index(entry->d_name);
        n++;
    }
    closedir(handle);
    qsort(entries, (size_t)n, sizeof(SampleEntry), compare_sample_entries);
    *count = n;
    return entries;
}

/* One replicate's own K x T series, same convention
   applications/abm_system_fit_qvarma_cluster.c's own read_y uses - kept as
   its own copy here rather than shared, same reason qvarma.h and
   qvarma_d.h stay two files: this project's own convention is that
   independent scripts do not import functions from one another. */
static void read_y(const char *sample, int replicate, Mat *y_out) {
    char csv_path[560];
    snprintf(csv_path, sizeof csv_path, "%s/%s/replicate_%03d.csv", INPUT_DIR, sample, replicate);
    DataFrame df = df_read_csv(csv_path, csv_read_options_default());
    Mat y = mat_new(K, df.r);
    static const char *row_name[K] = {
        "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
    };
    for (int k = 0; k < K; k++) {
        Mat column = df_col_numeric(&df, row_name[k]);
        for (int t = 0; t < df.r; t++) AT(y, k, t) = AT(column, t, 0);
    }
    df_free(&df);
    *y_out = y;
}

/* One spec's loss table: "replicate" plus one suffixed column per sample -
   the suffix is what keeps this spec's columns from colliding with the
   other spec's once both are joined. real_y is applications/
   us_qvarma_employment_change.c's own real-data block, shared across both
   specs since the data does not depend on which spec is being fit. */
static DataFrame build_spec_losses(Spec spec, Mat real_y, const SampleEntry *samples, int n_samples,
                                   FILE *manifest) {
    QvarmaParams real = spec_shape(spec.r);
    int loaded = qvarma_load_params(&real, spec.real_fit_path);
    assert(loaded && "abm_system_mse_qvarma: could not load the real-data fit - run "
                      "us_qvarma_employment_change.c's grid first");
    Vec real_irf;
    int real_ok = try_compute_irf(&real, real_y, &real_irf);
    assert(real_ok && "abm_system_mse_qvarma: the real-data fit's own nu <= 2 - "
                       "its impulse response function cannot be computed at all");

    Mat values = mat_new(N_REPLICATES, n_samples + 1);
    QvarmaParams working = spec_shape(spec.r);
    int n_missing = 0;

    for (int row = 0; row < N_REPLICATES; row++) {
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
                AT(values, row, col + 1) = stats_mae(real_irf, sim_irf);
                mat_free(sim_irf);
            } else {
                AT(values, row, col + 1) = (mreal)NAN;
                n_missing++;
                fprintf(manifest, "%s missing: %s replicate %03d\n", spec.label,
                        samples[col].name, row);
            }
        }
    }
    qvarma_params_free(&working);

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
            N_REPLICATES * n_samples);

    qvarma_params_free(&real);
    mat_free(real_irf);
    return df;
}

int main(void) {
    int n_samples;
    SampleEntry *samples = list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0 && "abm_system_mse_qvarma: no sample directories under out/abm_system_fit_qvarma/");

    Mat original = load_us_system();
    Mat real_y = build_real_block(original);
    mat_free(original);

    FILE *manifest = fopen("out/abm_system_mse_qvarma_manifest.txt", "w");
    assert(manifest && "abm_system_mse_qvarma: cannot open the manifest path for writing");
    fprintf(manifest, "%d samples, %d replicates each, %d specs, joined on replicate, "
                       "loss = MAE between stacked impulse response vectors (horizon %d)\n\n",
            n_samples, N_REPLICATES, N_SPECS, HORIZON);

    assert(N_SPECS == 2 && "abm_system_mse_qvarma: the join step below is written for exactly two specs");
    DataFrame left = build_spec_losses(spec_list[0], real_y, samples, n_samples, manifest);
    DataFrame right = build_spec_losses(spec_list[1], real_y, samples, n_samples, manifest);
    mat_free(real_y);

    DataFrame own_half = df_join(&left, &right, "replicate", JOIN_INNER);
    assert(own_half.r == N_REPLICATES &&
           "abm_system_mse_qvarma: every replicate exists in both specs' tables, so an inner "
           "join on replicate should keep all of them - something upstream disagrees");
    df_write_csv(&own_half, OUTPUT_PATH, csv_write_options_default());

    df_free(&left);
    df_free(&right);
    df_free(&own_half);

    fclose(manifest);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
