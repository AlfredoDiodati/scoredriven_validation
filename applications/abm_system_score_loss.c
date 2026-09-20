/*
Two loss matrices for the Model Confidence Set, built without fitting
anything: how strongly each simulated series pulls the real data's own fitted
auxiliary model away from where the real data put it.

applications/abm_system_irf_loss.c fits a t-QVARMA to every simulated
replicate and measures the distance between that fit's impulse responses and
the real data's. Here nothing is estimated. The real data's fitted parameter
vector theta_US (out/us_qvarma_spec_choice_p1q1r2_fit.json, the same p1q1r2
spec) is held fixed, and for configuration i and replicate n the score

    q_in = d log L(theta_US ; y_in) / d theta

is evaluated on that replicate's simulated series y_in. It is zero when the
simulated series leaves the real data's estimate exactly at its own maximum
and grows with how far the simulation would push it, so a configuration whose
dynamics match the US data scores low for the same reason the IRF distance
does, one estimation step earlier.

The score is et_al.'s own analytic gradient of the t-QVARMA log-likelihood
(qvarma_analytic_log_likelihood with a non-NULL gradient), which is exact
rather than a difference quotient.

Two losses are written, from the same pass over the data.

The plain sum of squares, q_in' q_in, out/abm_system_score_loss.csv. It gives
every coordinate of theta equal weight, which means it is not a distance in
any metric the model supplies: theta is the optimizer's own unconstrained
parameterization, its coordinates carry different units and differ in scale by
orders of magnitude, and the loss is dominated by whichever of them happens to
have the largest one. Measured on this dataset, Omega_inv[2,0] alone carries
62 per cent of it. Kept because it is the object that was asked for and
because it is what the weighted version has to be compared against.

The score test statistic, out/abm_system_score_loss_weighted.csv:

    LM_in = q_in' [T_sim I_1(theta_US)]^-1 q_in

with I_1 the information matrix per observation. This is Rao's score statistic
(the Lagrange multiplier statistic) for the hypothesis that y_in was generated
by the auxiliary model at theta_US, and it is the statistically meaningful
version of the same idea for two reasons.

It is invariant to how the model is parameterized. Under a smooth
reparameterization theta -> g(theta) both q and I transform by the same
Jacobian and the quadratic form is unchanged, so the weighted loss measures
the same thing whether it is computed in the unconstrained coordinates the
optimizer steps or in the paper's own parameters. The unweighted q'q does not:
it changes if a coordinate is rescaled, which is why one entry of the scale
matrix dominates it.

It has a null distribution. Under the hypothesis that y_in came from the
auxiliary model at theta_US, LM_in is asymptotically chi-square with
qvarma_n_theta degrees of freedom, so the number can be read against a
reference rather than only against the other cells of the matrix. The manifest
prints that reference.

I_1 is estimated once, from the real data, as the observed information at
theta_US:

    I_US = -d2 log L(theta_US ; y_US) / d theta d theta'

divided by the 187 periods it was computed over, and multiplied back up by the
400 periods a simulated replicate has. et_al.'s _qvarma_hessian differences the
negative log-likelihood, so what it returns is I_US directly.

One weighting matrix for every cell, taken from the real data rather than from
each configuration's own replicates. A per-configuration weight would give
every configuration its own metric, and losses measured in different metrics
cannot be compared, which is exactly what the Model Confidence Set has to do
with them. The cost is that the chi-square calibration holds only under the
null: where a configuration is far from the US data the variance of q_in under
that configuration's own process is not I_US, so LM_in is a distance in a
fixed, model-supplied metric rather than a statistic with a valid per-cell
p-value. That is what the confidence set needs and it is not more than that.

The weighting is applied through the Cholesky factor of T_sim/T_US * I_US:
with that matrix equal to L L', the loss is the squared norm of the solution
of L x = q, one triangular solve per cell and no inverse formed anywhere. The
factor is computed once, before the parallel loop.

Both matrices are written from one sweep because reading the archives is what
the pass costs; the second quadratic form is free beside it.

The score is a sum over periods, so it scales with the length of a replicate.
Every replicate in the dataset has the same length, so this does not vary
across the cells of either matrix, but a comparison against a matrix built
over series of a different length is not meaningful.

A replicate whose scale matrix or likelihood is not usable at theta_US comes
back from the filter as minus infinity with a zeroed gradient, which would
read as a perfect loss of zero. That case is counted as a missing cell, not as
a zero.

Requires dataset/abm_system/ (applications/abm_system_simulate.c or
applications/abm_system_convert_rdata.c),
out/us_qvarma_spec_choice_p1q1r2_fit.json
(applications/us_qvarma_spec_choice.c) and out/us_system.csv
(applications/us_prepare_data.c), the last of these because the information
matrix is computed on the real data. It does not read
out/abm_system_fit_qvarma/ at all - no fit of a simulated series takes part in
either loss.

Output: two CSVs, one row per replicate and one column per configuration plus
a leading "replicate" column, with the column names
applications/abm_system_irf_loss.c uses for the same configurations so the
matrices can be read side by side. The manifest is named after the unweighted
table, out/abm_system_score_loss_manifest.txt, covers both, and records the
information matrix's eigenvalues and condition number, since a weighting that
came out of a badly conditioned matrix is worth no more than its conditioning.

Neither matrix has holes, which et_al's mcs requires: one NaN in a model's
column loses every comparison that model takes part in and the p-values come
back looking ordinary. Where a cell is missing the replicate is dropped from
every column of both matrices, the same direction
applications/abm_system_irf_loss.c drops in, since that keeps every
configuration in the comparison.

In EXPERIMENT_STEMS, buildable on its own via make app-abm_system_score_loss.
Nothing printed.
*/

#include "abm_system.h"
#include "us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./linalg/decomp.h>
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
#define REAL_FIT_PATH "out/us_qvarma_spec_choice_p1q1r2_fit.json"
#define REAL_PERIODS (ESTIMATION_PERIODS - 1)

#define INPUT_DIR_DEFAULT "dataset/abm_system"
#define OUTPUT_PATH_DEFAULT "out/abm_system_score_loss.csv"

/* Overridable so the pass can be run over a couple of configurations before it
   is turned loose on the whole design. */
static const char *INPUT_DIR;
static const char *OUTPUT_PATH;

/* Counted off the dataset rather than fixed here: the .Rdata route produced
   108 replicates per configuration and the design experiment produces 1000. */
static int n_replicates = 0;

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

/* Real data, applications/us_qvarma_spec_choice.c's own build_block: growth or
   change of GDP, energy demand and employment, inflation, and the interest
   rate in levels. Identical row convention to the simulated series, which is
   what makes one information matrix serve both. */
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

typedef struct { char *name; int index; } SampleEntry;

/* dataset/abm_system/cop_<N>'s own <N>, so column order is 1, 2, ..., 1000
   rather than readdir's arbitrary order or a string sort's "_10" before
   "_2". */
static int trailing_index(const char *name) {
    const char *underscore = strrchr(name, '_');
    assert(underscore && "abm_system_score_loss: a configuration directory name has no trailing _<N>");
    return atoi(underscore + 1);
}

static int compare_sample_entries(const void *a, const void *b) {
    return ((const SampleEntry*)a)->index - ((const SampleEntry*)b)->index;
}

static SampleEntry *list_samples(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_score_loss: cannot open dataset/abm_system/ - "
                      "run applications/abm_system_simulate_all.sh first");

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
            assert(grown && "abm_system_score_loss: out of memory listing configurations");
            entries = grown;
        }
        size_t len = strlen(entry->d_name);
        entries[n].name = malloc(len + 1);
        assert(entries[n].name && "abm_system_score_loss: out of memory copying a configuration name");
        memcpy(entries[n].name, entry->d_name, len + 1);
        entries[n].index = trailing_index(entry->d_name);
        n++;
    }
    closedir(handle);
    qsort(entries, (size_t)n, sizeof(SampleEntry), compare_sample_entries);
    *count = n;
    return entries;
}

/* One past the highest replicate index one configuration holds, so a row of
   the table is a replicate index rather than a position in whatever order the
   archives happened to be read. */
static int count_replicates(const char *sample) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    int count = 0;
    int *replicate = abm_system_list_replicates(dir, &count);
    assert(count > 0 && "abm_system_score_loss: a configuration directory holds no replicate");
    int highest = replicate[0];
    for (int i = 1; i < count; i++) if (replicate[i] > highest) highest = replicate[i];
    free(replicate);
    return highest + 1;
}

/* One replicate, to size the information matrix's rescaling and to have
   something concrete to report the weighting's own behaviour on. The periods
   are read from the data rather than assumed, since the two routes that wrote
   dataset/abm_system/ do not agree about them. Caller must mat_free. */
static Mat first_replicate(const char *sample) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    Mat block[ABM_SYSTEM_BATCH];
    int replicate[ABM_SYSTEM_BATCH];
    int count = abm_system_read_batch(dir, 0, block, replicate);
    assert(count > 0 && "abm_system_score_loss: the first archive of a configuration is empty");
    Mat y = mat_copy(block[0]);
    for (int b = 0; b < count; b++) mat_free(block[b]);
    return y;
}

/*
Where the unweighted q'q comes from, coordinate by coordinate. This is the
defect the weighting exists to remove, so it is measured on the same replicate
the weighted split below is measured on rather than asserted.
*/
static void report_coordinate_shares(Vec q, const QvarmaParams *shape, const char *sample,
                                     FILE *manifest) {
    int n = q.r;
    double total = 0;
    double *share = (double*)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        share[i] = (double)q.d[i] * (double)q.d[i];
        total += share[i];
    }

    int *order = (int*)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) order[i] = i;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            if (share[order[b]] > share[order[a]]) { int t = order[a]; order[a] = order[b]; order[b] = t; }

    fprintf(manifest, "  how q'q splits over coordinates of theta, on %s replicate 0, total %.4e\n",
            sample, total);
    double running = 0;
    char name[64];
    for (int a = 0; a < 5 && a < n; a++) {
        running += share[order[a]] / total;
        _qvarma_theta_name(shape, order[a], name, sizeof name);
        fprintf(manifest, "    %-18s share %.4f  cumulative %.4f\n",
                name, share[order[a]] / total, running);
    }
    free(order);
    free(share);
}

/*
Where the weighted statistic comes from, direction by direction. With the
information matrix equal to sum_k lambda_k v_k v_k', the statistic splits as

    q' I^-1 q = sum_k (v_k . q)^2 / lambda_k

so a direction the real data pins down poorly (small lambda_k) counts for more
than one it pins down well. That is what the weighting is for, and it is also
how the weighting can go wrong: if one direction carries almost all of the
statistic then the loss is again effectively one number, as the unweighted
q'q is, only in a different basis. The shares are reported on one replicate so
that this can be read off rather than assumed.
*/
static void report_direction_shares(Vec q, Vec eigenvalues, Mat eigenvectors, mreal scale,
                                    const char *sample, FILE *manifest) {
    int n = q.r;
    double total = 0;
    double *share = (double*)malloc((size_t)n * sizeof(double));
    for (int k = 0; k < n; k++) {
        double projection = 0;
        for (int i = 0; i < n; i++) projection += (double)AT(eigenvectors, i, k) * (double)q.d[i];
        share[k] = projection * projection / ((double)eigenvalues.d[k] * (double)scale);
        total += share[k];
    }

    int *order = (int*)malloc((size_t)n * sizeof(int));
    for (int k = 0; k < n; k++) order[k] = k;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            if (share[order[b]] > share[order[a]]) { int t = order[a]; order[a] = order[b]; order[b] = t; }

    fprintf(manifest, "  how the statistic splits over eigen-directions, on %s replicate 0, "
                       "total %.4e\n", sample, total);
    double running = 0;
    for (int a = 0; a < 5 && a < n; a++) {
        running += share[order[a]] / total;
        fprintf(manifest, "    direction %2d  eigenvalue %.4e  share %.4f  cumulative %.4f\n",
                order[a], (double)eigenvalues.d[order[a]], share[order[a]] / total, running);
    }
    free(order);
    free(share);
}

/*
The Cholesky factor of T_sim/T_US * I_US, which is the matrix the weighted loss
solves against: with it equal to L L', q' [T_sim I_1]^-1 q is the squared norm
of the solution of L x = q. Returning the factor rather than the inverse is
what keeps the hot loop to one triangular solve and forms no inverse at all.

The eigenvalues go to the manifest before the factorization, because
mat_chol's own assert on a matrix that is not positive definite says nothing
about how close to indefinite it was, and a weighting matrix is only worth its
conditioning.
*/
static Mat build_information_factor(Vec theta, const QvarmaParams *shape, Mat real_y,
                                    Mat probe_y, const char *probe_sample, FILE *manifest) {
    int n = theta.r, sim_periods = probe_y.c;
    Mat information = _qvarma_hessian(theta, shape, real_y);

    Vec eigenvalues;
    Mat eigenvectors;
    int status = mat_eig_sym_status(information, &eigenvalues, &eigenvectors);
    assert(status == 0 && "abm_system_score_loss: the information matrix at theta_US has no "
                           "eigendecomposition - the weighting cannot be built from it");

    mreal smallest = eigenvalues.d[0], largest = eigenvalues.d[0];
    for (int i = 1; i < n; i++) {
        if (eigenvalues.d[i] < smallest) smallest = eigenvalues.d[i];
        if (eigenvalues.d[i] > largest) largest = eigenvalues.d[i];
    }
    /* Curvature this far below the largest cannot be told from zero once the
       Hessian has been differenced, the same floor qvarma_standard_errors
       reads a flat direction against. */
    mreal floor_value = (mreal)(n * MEPS) * MABS(largest);
    int n_flat = 0;
    for (int i = 0; i < n; i++) if (MABS(eigenvalues.d[i]) <= floor_value) n_flat++;

    fprintf(manifest, "weighting\n");
    fprintf(manifest, "  observed information at theta_US, -d2 logL/dtheta dtheta', "
                       "differenced from the analytic gradient over %d periods of real data\n",
            real_y.c);
    fprintf(manifest, "  eigenvalues: smallest %.6e, largest %.6e, condition %.3e\n",
            (double)smallest, (double)largest,
            smallest > 0 ? (double)(largest / smallest) : INFINITY);
    fprintf(manifest, "  positive definite %s, %d directions flat against a floor of %.3e\n",
            smallest > floor_value ? "yes" : "no", n_flat, (double)floor_value);

    assert(smallest > floor_value && "abm_system_score_loss: the information matrix at theta_US is "
                                      "not positive definite - the weighted loss would not be a "
                                      "distance and is not written");

    mreal scale = (mreal)sim_periods / (mreal)real_y.c;
    Mat scaled = mat_copy(information);
    for (int i = 0; i < n * n; i++) scaled.d[i] *= scale;
    Mat factor = mat_chol(scaled);

    fprintf(manifest, "  scaled by T_sim / T_US = %d / %d before inversion, so the loss is on the "
                       "chi-square scale\n", sim_periods, real_y.c);
    fprintf(manifest, "  under the null that a replicate came from the auxiliary model at theta_US, "
                       "chi-square with %d degrees of freedom: mean %d, sd %.2f\n",
            n, n, sqrt(2.0 * n));
    fprintf(manifest, "  that calibration holds under the null only - away from it the weighted "
                       "loss is a distance in a fixed metric, not a per-cell p-value\n");

    QvarmaAnalytic *probe_filter = qvarma_analytic_new(shape, sim_periods);
    Vec probe_score = mat_new(n, 1);
    qvarma_analytic_log_likelihood(probe_filter, theta, probe_y, probe_score);
    report_coordinate_shares(probe_score, shape, probe_sample, manifest);
    report_direction_shares(probe_score, eigenvalues, eigenvectors, scale, probe_sample, manifest);
    mat_free(probe_score);
    qvarma_analytic_free(probe_filter);
    fprintf(manifest, "\n");

    mat_free(scaled);
    mat_free(eigenvalues);
    mat_free(eigenvectors);
    mat_free(information);
    return factor;
}

/* Both quadratic forms of one replicate's score, or NaN in both when the
   filter cannot be evaluated there. The zeroed gradient the filter returns
   alongside minus infinity would otherwise read as a loss of zero, which is
   the best score in the table rather than the worst.

   solved is scratch the caller owns, so the triangular solve allocates
   nothing: at a million cells the allocation would cost more than the solve. */
static void score_losses(QvarmaAnalytic *filter, Vec theta, Mat y, Vec score, Mat factor,
                         Vec solved, mreal *plain_out, mreal *weighted_out) {
    *plain_out = (mreal)NAN;
    *weighted_out = (mreal)NAN;

    mreal value = qvarma_analytic_log_likelihood(filter, theta, y, score);
    if (MISNAN(value) || MISINF(value)) return;

    double plain = 0;
    for (int i = 0; i < score.r; i++) plain += (double)score.d[i] * (double)score.d[i];
    if (!isfinite(plain)) return;

    for (int i = 0; i < score.r; i++) solved.d[i] = score.d[i];
    if (_trtrs('L', 'N', 'N', factor.r, 1, factor.d, factor.stride, solved.d, solved.stride) != 0)
        return;

    double weighted = 0;
    for (int i = 0; i < solved.r; i++) weighted += (double)solved.d[i] * (double)solved.d[i];
    if (!isfinite(weighted)) return;

    *plain_out = (mreal)plain;
    *weighted_out = (mreal)weighted;
}

/* The workspace the analytic filter evaluates through. It is sized for one
   series length, so it is built once per thread and rebuilt only if a
   replicate of a different length turns up. */
typedef struct { QvarmaAnalytic *filter; int periods; } Workspace;

static QvarmaAnalytic *workspace_for(Workspace *w, const QvarmaParams *shape, int periods) {
    if (w->filter && w->periods == periods) return w->filter;
    if (w->filter) qvarma_analytic_free(w->filter);
    w->filter = qvarma_analytic_new(shape, periods);
    w->periods = periods;
    return w->filter;
}

/* One configuration's column of both matrices, read archive by archive rather
   than replicate by replicate - see abm_system.h's own abm_system_read_batch
   on why. */
static int fill_column(const char *sample, Vec theta, const QvarmaParams *shape, Mat factor,
                       Workspace *workspace, Mat plain, Mat weighted, int col, FILE *manifest) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);

    Vec score = mat_new(theta.r, 1);
    Vec solved = mat_new(theta.r, 1);
    for (int row = 0; row < n_replicates; row++) {
        AT(plain, row, col + 1) = (mreal)NAN;
        AT(weighted, row, col + 1) = (mreal)NAN;
    }

    int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
    for (int batch = 0; batch < n_batches; batch++) {
        Mat block[ABM_SYSTEM_BATCH];
        int replicate[ABM_SYSTEM_BATCH];
        int count = abm_system_read_batch(dir, batch, block, replicate);
        for (int b = 0; b < count; b++) {
            if (replicate[b] < n_replicates) {
                QvarmaAnalytic *filter = workspace_for(workspace, shape, block[b].c);
                score_losses(filter, theta, block[b], score, factor, solved,
                             &AT(plain, replicate[b], col + 1),
                             &AT(weighted, replicate[b], col + 1));
            }
            mat_free(block[b]);
        }
    }
    mat_free(score);
    mat_free(solved);

    int n_missing = 0;
    for (int row = 0; row < n_replicates; row++)
        if (MISNAN(AT(plain, row, col + 1))) {
            n_missing++;
            #pragma omp critical
            fprintf(manifest, "missing: %s replicate %03d\n", sample, row);
        }
    return n_missing;
}

/* The table as a frame, named for the configurations it was built over. */
static DataFrame to_frame(Mat values, const SampleEntry *samples, int n_samples) {
    char **col_names = (char**)malloc((size_t)(n_samples + 1) * sizeof(char*));
    col_names[0] = frame_strdup("replicate");
    for (int col = 0; col < n_samples; col++) {
        char buf[128];
        snprintf(buf, sizeof buf, "%s_%s", samples[col].name, COLUMN_LABEL);
        col_names[col + 1] = frame_strdup(buf);
    }
    DataFrame df = df_from_matrix(values, (const char *const *)col_names);
    for (int col = 0; col <= n_samples; col++) free(col_names[col]);
    free(col_names);
    return df;
}

/* Rows with a missing cell anywhere, dropped from both matrices together so
   the two stay the same shape and the same replicates. */
static Mat drop_incomplete(Mat values, const int *usable, int keep, int n_samples) {
    if (keep == values.r) return values;
    Mat kept = mat_new(keep, n_samples + 1);
    int at = 0;
    for (int row = 0; row < values.r; row++) {
        if (!usable[row]) continue;
        for (int col = 0; col <= n_samples; col++) AT(kept, at, col) = AT(values, row, col);
        at++;
    }
    mat_free(values);
    return kept;
}

int main(void) {
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = INPUT_DIR_DEFAULT;
    OUTPUT_PATH = getenv("ABM_SYSTEM_SCORE_LOSS_PATH");
    if (!OUTPUT_PATH) OUTPUT_PATH = OUTPUT_PATH_DEFAULT;

    /* All the parallelism belongs to the outer loop over configurations.
       Without this OpenBLAS starts its own pool inside every worker thread and
       the machine's cores are oversubscribed by a factor of the thread count,
       the same reason applications/abm_system_fit_qvarma.c sets it. */
    openblas_set_num_threads(1);

    int n_samples;
    SampleEntry *samples = list_samples(INPUT_DIR, &n_samples);
    assert(n_samples > 0 && "abm_system_score_loss: no configuration directories under dataset/abm_system/");
    n_replicates = count_replicates(samples[0].name);
    Mat probe_y = first_replicate(samples[0].name);
    int sim_periods = probe_y.c;

    /* theta_US. qvarma_load_params checks every field of the shape it is given
       against the stored one, so a cache written under a different spec does
       not load at all, and _qvarma_unlink is the exact inverse of the link the
       cache's own theta went through. */
    QvarmaParams real = spec_shape();
    int loaded = qvarma_load_params(&real, REAL_FIT_PATH);
    assert(loaded && "abm_system_score_loss: could not load the real-data fit - run "
                      "us_qvarma_spec_choice.c first");
    Vec theta = mat_new(qvarma_n_theta(&real), 1);
    _qvarma_unlink(&real, theta);

    char weighted_path[640], manifest_path[640];
    size_t stem = strlen(OUTPUT_PATH);
    if (stem > 4 && strcmp(OUTPUT_PATH + stem - 4, ".csv") == 0) stem -= 4;
    snprintf(weighted_path, sizeof weighted_path, "%.*s_weighted.csv", (int)stem, OUTPUT_PATH);
    snprintf(manifest_path, sizeof manifest_path, "%.*s_manifest.txt", (int)stem, OUTPUT_PATH);

    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "abm_system_score_loss: cannot open the manifest path for writing");
    fprintf(manifest, "%d configurations, %d replicates each, %d periods each, spec %s\n",
            n_samples, n_replicates, sim_periods, SPEC_LABEL);
    fprintf(manifest, "q is the score of the t-QVARMA log-likelihood at the real-data estimate, "
                       "evaluated on the simulated series\n");
    fprintf(manifest, "theta from %s, %d coordinates, unconstrained scale\n", REAL_FIT_PATH, theta.r);
    fprintf(manifest, "nothing is estimated here: no simulated series is fitted\n");
    fprintf(manifest, "%s holds q'q\n", OUTPUT_PATH);
    fprintf(manifest, "%s holds the score statistic q' [T_sim I_1]^-1 q\n\n", weighted_path);

    Mat original = load_us_system();
    Mat real_y = build_real_block(original);
    mat_free(original);
    Mat factor = build_information_factor(theta, &real, real_y, probe_y, samples[0].name, manifest);
    mat_free(real_y);
    mat_free(probe_y);

    Mat plain = mat_new(n_replicates, n_samples + 1);
    Mat weighted = mat_new(n_replicates, n_samples + 1);
    for (int row = 0; row < n_replicates; row++) {
        AT(plain, row, 0) = (mreal)row;
        AT(weighted, row, 0) = (mreal)row;
    }

    int n_missing = 0;
    /* One task per configuration, so each thread reads a whole configuration's
       archives and writes its own column of both tables. The manifest is the
       only shared sink and is taken under a lock, a missing cell being rare
       enough that nothing ever waits for it. */
    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams shape = spec_shape();
        Workspace workspace = { NULL, 0 };

        #pragma omp for schedule(dynamic)
        for (int col = 0; col < n_samples; col++)
            n_missing += fill_column(samples[col].name, theta, &shape, factor, &workspace,
                                     plain, weighted, col, manifest);

        if (workspace.filter) qvarma_analytic_free(workspace.filter);
        qvarma_params_free(&shape);
    }

    int keep = 0;
    int *usable = (int*)malloc((size_t)n_replicates * sizeof(int));
    for (int row = 0; row < n_replicates; row++) {
        usable[row] = 1;
        for (int col = 0; col < n_samples; col++)
            if (MISNAN(AT(plain, row, col + 1))) { usable[row] = 0; break; }
        if (usable[row]) keep++;
    }
    assert(keep > 0 && "abm_system_score_loss: every replicate has a missing cell");
    if (keep < n_replicates)
        fprintf(manifest, "%d of %d replicates dropped for holding a missing cell\n",
                n_replicates - keep, n_replicates);

    plain = drop_incomplete(plain, usable, keep, n_samples);
    weighted = drop_incomplete(weighted, usable, keep, n_samples);
    free(usable);

    DataFrame plain_frame = to_frame(plain, samples, n_samples);
    DataFrame weighted_frame = to_frame(weighted, samples, n_samples);
    assert(plain_frame.r > 0 && "abm_system_score_loss: the loss table has no rows");
    df_write_csv(&plain_frame, OUTPUT_PATH, csv_write_options_default());
    df_write_csv(&weighted_frame, weighted_path, csv_write_options_default());

    fprintf(manifest, "%d of %d cells missing\n", n_missing, n_replicates * n_samples);

    df_free(&plain_frame);
    df_free(&weighted_frame);
    mat_free(plain);
    mat_free(weighted);
    mat_free(factor);
    mat_free(theta);
    qvarma_params_free(&real);
    fclose(manifest);
    for (int i = 0; i < n_samples; i++) free(samples[i].name);
    free(samples);
    return 0;
}
