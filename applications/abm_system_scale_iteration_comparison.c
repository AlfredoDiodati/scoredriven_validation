/*
What raising the solver budget from 2000 to 4000 iterations bought, over the
same 500,000 series, and what it cost.

86.85% of the fits applications/abm_system_scale_fit_qvarma.c wrote at a cap
of 2000 stopped at the cap rather than at the gradient tolerance, so their
estimates are wherever L-BFGS happened to be at iteration 2000 rather than at
a maximum. Whether that matters is not answerable from the 2000-iteration run
alone: it needs the same 500,000 fits at a larger budget and a comparison of
the two, which is what this reads and reports.

Both budgets are on disk as whole trees, out/abm_system_scale_fit_qvarma_i2000/
and _i4000/, because the fitting script names every file it writes after its
own cap. This opens both JSONs for each series and compares them. It fits
nothing and writes nothing into either tree.

The two runs share a starting point, data and solver, and L-BFGS is
deterministic and returns the best point it saw, so the 4000-iteration run
retraces the 2000-iteration path exactly for its first 2000 iterations and
then continues. The log-likelihood gain is therefore non-negative by
construction, exactly zero wherever the fit had already converged or found
nothing better. That is not an assumption this makes but a property it checks:
a negative gain would mean the two runs did not follow the same path, and the
count of them is reported rather than passed over.

Three things are compared, because "is it worth it" is three questions:

  - the objective. How much log-likelihood the extra 2000 iterations bought,
    as a distribution rather than a mean, since a handful of fits moving a
    long way and every fit moving a little are different situations with the
    same average.
  - convergence. How many fits that stopped at 2000 met the gradient
    tolerance by 4000, and what the gradient norm did for the ones that still
    did not.
  - the estimates themselves. The objective moving is not the same as the
    estimates moving, and it is the estimates the indirect-inference step
    consumes. Movement is measured on the constrained parameters, the scale
    the paper names them on, both as one Euclidean distance per fit and as
    the largest absolute change inside each named block, so a large distance
    can be attributed to nu or to the co-integrating vector rather than left
    as a number.

One thing has to be excluded before any of that means anything. The
Student-t degrees of freedom enter through nu = exp(theta) + 2, so a search
running up that coordinate reaches nu where the coded log-likelihood stops
being computable: the density's own (nu + K)/2 log(1 + q/nu) term is a
cancellation of two large quantities, and past nu of about 3.6e9 it loses
decimal places, past about 5e11 it returns values with no relation to the
likelihood at all. Measured on
dataset/abm_system_scale/EstimationSeriesSample1_47_noise4/series_743.csv by
sweeping theta_nu at that series' own 2000-iteration estimate: the value sits
at 393.788 - the Gaussian limit the t tends to - from nu of 6.6e7 through
8.0e8, wobbles in the low decimals from 3.6e9, and by nu of 1.1e13 reads 616,
then 1216, then 1664. Values of exactly 262144 and -4.611686e18 appear along
the way, which are 2^18 and -2^62 rather than likelihoods.

L-BFGS climbs that artifact where it finds it, and a larger budget gives it
more room to. Fits carrying nu above NU_DIVERGENCE_LIMIT at either budget are
therefore counted and reported on their own, and left out of every average,
quantile and total below: five such fits are enough to put a mean
log-likelihood gain in the 1e32 range and say nothing about the other 499,995.
The limit is deliberately below where the value visibly breaks, since the
precision is already gone before that.

This is a limitation of the likelihood in et_al.'s own sd/qvarma.h, not of
this script or of the fits, and belongs there rather than in a workaround
here.

Wall time comes from the two runs' own timing files, quoted rather than
recomputed.

out/abm_system_scale_iteration_comparison.txt holds the summary and
out/abm_system_scale_iteration_comparison.csv one row per series, for
whatever plot the summary does not answer. Nothing printed.

Requires both trees to exist and to hold the same series. A series present in
one and missing from the other is counted and skipped rather than aborting the
run, since a partial second budget is a normal state to want a comparison of.
*/

#include "abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./frame/csv.h>
#include <omp.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define P 1
#define Q 1
#define RLAG 2
#define SPEC_LABEL "p1q1r2"

#ifndef BASE_ITERATIONS
#define BASE_ITERATIONS 2000
#endif
#ifndef HIGH_ITERATIONS
#define HIGH_ITERATIONS 4000
#endif

/* Above this nu the coded log-likelihood has lost precision; see this file's
   own header comment for the sweep that located it. */
#define NU_DIVERGENCE_LIMIT ((mreal)1e10)

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#define INPUT_DIR "dataset/abm_system_scale"
#define BASE_DIR "out/abm_system_scale_fit_qvarma_i" STRINGIFY(BASE_ITERATIONS)
#define HIGH_DIR "out/abm_system_scale_fit_qvarma_i" STRINGIFY(HIGH_ITERATIONS)
#define BASE_TIMING "out/abm_system_scale_fit_qvarma_i" STRINGIFY(BASE_ITERATIONS) "_timing.txt"
#define HIGH_TIMING "out/abm_system_scale_fit_qvarma_i" STRINGIFY(HIGH_ITERATIONS) "_timing.txt"
#define BASE_PROVENANCE "out/abm_system_scale_fit_qvarma_i" STRINGIFY(BASE_ITERATIONS) "_provenance.txt"
#define HIGH_PROVENANCE "out/abm_system_scale_fit_qvarma_i" STRINGIFY(HIGH_ITERATIONS) "_provenance.txt"
/* Named for the pair compared, so a second pair does not overwrite the first. */
#define REPORT_PATH "out/abm_system_scale_iteration_comparison_i" STRINGIFY(BASE_ITERATIONS) \
                    "_i" STRINGIFY(HIGH_ITERATIONS) ".txt"
#define CSV_PATH "out/abm_system_scale_iteration_comparison_i" STRINGIFY(BASE_ITERATIONS) \
                 "_i" STRINGIFY(HIGH_ITERATIONS) ".csv"

/* The named blocks of QvarmaParams, in the paper's own spelling where it has
   one. Movement is reported per block so a distance can be attributed. */
enum { BLOCK_C, BLOCK_PHI_STAR, BLOCK_PSI_STAR, BLOCK_OMEGA_INV, BLOCK_NU,
       BLOCK_ALPHA, BLOCK_BETA, N_BLOCKS };
static const char *block_name[N_BLOCKS] = {
    "c", "Phi_star", "Psi_star", "Omega_inv", "nu", "alpha", "beta"
};

static int compare_names(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static char **list_subdirs(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_scale_iteration_comparison: cannot open dataset/abm_system_scale/");
    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char path[600];
        snprintf(path, sizeof path, "%s/%s", dir, entry->d_name);
        struct stat st;
        if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) continue;
        if (n == cap) { cap = cap ? cap * 2 : 16; names = realloc(names, (size_t)cap * sizeof(char*)); }
        size_t name_len = strlen(entry->d_name);
        names[n] = malloc(name_len + 1);
        memcpy(names[n], entry->d_name, name_len + 1);
        n++;
    }
    closedir(handle);
    qsort(names, (size_t)n, sizeof(char*), compare_names);
    *count = n;
    return names;
}

static int count_series(const char *folder) {
    char path[600];
    snprintf(path, sizeof path, "%s/%s", INPUT_DIR, folder);
    DIR *handle = opendir(path);
    assert(handle && "abm_system_scale_iteration_comparison: cannot open a folder under the dataset");
    int n = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL)
        if (strncmp(entry->d_name, "series_", 7) == 0) n++;
    closedir(handle);
    return n;
}

/* Whether a run's own provenance stamp says it carries the large-nu
   correction, as "yes", "no" or "unrecorded". Two runs answering differently
   are not comparable on their budgets alone, and the report says so rather
   than leaving a reader to notice. */
static const char *carries_correction(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return "unrecorded";
    char line[512];
    const char *answer = "unrecorded";
    while (fgets(line, sizeof line, f)) {
        char *at = strstr(line, "carries the large-nu correction:");
        if (!at) continue;
        answer = strstr(at, "yes") ? "yes" : strstr(at, "no") ? "no" : "unrecorded";
    }
    fclose(f);
    return answer;
}

/* Copies a run's provenance stamp into the report, indented, so the report is
   self-contained rather than pointing at a file that may move. */
static void quote_provenance(FILE *report, const char *label, const char *path) {
    FILE *f = fopen(path, "r");
    fprintf(report, "  %s (%s)\n", label, path);
    if (!f) { fprintf(report, "    no stamp on disk\n\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, f)) fprintf(report, "    %s", line);
    fprintf(report, "\n");
    fclose(f);
}

/*
What a run actually cost, in seconds.

The timing file's own "elapsed" is the fit loop of one process. A battery
interrupted and resumed has more than one of those, and the last one covers
only the fits the resume still had to estimate - the 8000-iteration run was
suspended by the machine at 296,337 fits and finished in a second process, so
its timing file reports 3.67 h for what took 9.95. The provenance stamp beside
it carries the reconstructed total on a "one clean run" line, and that is
preferred wherever it exists, so this report cannot quote a resumed segment as
though it were the whole run.
*/
static double elapsed_from_timing(const char *timing_path, const char *provenance_path) {
    FILE *p = fopen(provenance_path, "r");
    if (p) {
        char line[512];
        double seconds = -1;
        while (fgets(line, sizeof line, p)) {
            char *at = strstr(line, "one clean run");
            if (at && sscanf(at + 13, " %lf", &seconds) == 1) { fclose(p); return seconds; }
        }
        fclose(p);
    }
    FILE *f = fopen(timing_path, "r");
    if (!f) return -1;
    char line[512];
    double seconds = -1;
    while (fgets(line, sizeof line, f))
        if (strncmp(line, "elapsed", 7) == 0) sscanf(line + 7, " %lf", &seconds);
    fclose(f);
    return seconds;
}

static int load_params(QvarmaParams *out, Mat y, const char *path, QvarmaFitResult *diagnostics) {
    *out = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    out->phi_star_bound = PHI_STAR_BOUND;
    out->mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    diagnostics->params = *out;
    if (qvarma_load_fit(diagnostics, y, path)) { *out = diagnostics->params; return 1; }
    qvarma_params_free(out);
    return 0;
}

/* Squared distance and largest absolute change between two fitted parameter
   sets, on the constrained scale, accumulated into per-block maxima. */
static void compare_params(const QvarmaParams *a, const QvarmaParams *b,
                           double *sum_squares, double *base_sum_squares,
                           double *sum_squares_no_nu, double *base_sum_squares_no_nu,
                           double *max_absolute, double *block_max) {
    *sum_squares = 0; *base_sum_squares = 0; *max_absolute = 0;
    *sum_squares_no_nu = 0; *base_sum_squares_no_nu = 0;
    for (int i = 0; i < N_BLOCKS; i++) block_max[i] = 0;

    #define ACCUMULATE(block, va, vb) do { \
        double difference = (double)(vb) - (double)(va); \
        double magnitude = difference < 0 ? -difference : difference; \
        *sum_squares += difference * difference; \
        *base_sum_squares += (double)(va) * (double)(va); \
        if ((block) != BLOCK_NU) { \
            *sum_squares_no_nu += difference * difference; \
            *base_sum_squares_no_nu += (double)(va) * (double)(va); \
        } \
        if (magnitude > *max_absolute) *max_absolute = magnitude; \
        if (magnitude > block_max[block]) block_max[block] = magnitude; \
    } while (0)

    for (int k = 0; k < K; k++) ACCUMULATE(BLOCK_C, AT(a->c, k, 0), AT(b->c, k, 0));
    for (int i = 0; i < P; i++)
        ACCUMULATE(BLOCK_PHI_STAR, AT(a->Phi_star, i, 0), AT(b->Phi_star, i, 0));
    int psi_rows = qvarma_psi_star_rows(a);
    for (int l = 0; l < Q; l++)
        for (int i = 0; i < psi_rows; i++)
            for (int j = 0; j < K; j++)
                ACCUMULATE(BLOCK_PSI_STAR, AT(a->Psi_star[l], i, j), AT(b->Psi_star[l], i, j));
    for (int i = 0; i < K; i++)
        for (int j = 0; j <= i; j++)
            ACCUMULATE(BLOCK_OMEGA_INV, AT(a->Omega_inv, i, j), AT(b->Omega_inv, i, j));
    ACCUMULATE(BLOCK_NU, a->nu, b->nu);
    int K_dag = K - K_STAR;
    for (int l = 0; l < RLAG; l++)
        for (int i = 0; i < K_dag; i++)
            for (int j = 0; j < R; j++)
                ACCUMULATE(BLOCK_ALPHA, AT(a->alpha[l], i, j), AT(b->alpha[l], i, j));
    int n_beta = qvarma_n_beta_matrices(a);
    for (int l = 0; l < n_beta; l++)
        for (int i = 0; i < R; i++)
            for (int j = 0; j < K_dag; j++)
                ACCUMULATE(BLOCK_BETA, AT(a->beta[l], i, j), AT(b->beta[l], i, j));
    #undef ACCUMULATE
}

static int compare_doubles(const void *a, const void *b) {
    double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* The value at a quantile of an already-sorted array, nearest-rank. */
static double quantile(const double *sorted, long n, double q) {
    if (n == 0) return 0;
    long index = (long)(q * (double)(n - 1) + 0.5);
    if (index < 0) index = 0;
    if (index >= n) index = n - 1;
    return sorted[index];
}

typedef struct { int folder_index, series; } Task;

int main(void) {
    int n_folders;
    char **folders = list_subdirs(INPUT_DIR, &n_folders);
    assert(n_folders > 0 && "abm_system_scale_iteration_comparison: no folders in the dataset");

    int *n_series = malloc((size_t)n_folders * sizeof(int));
    long total_tasks = 0;
    for (int f = 0; f < n_folders; f++) {
        n_series[f] = count_series(folders[f]);
        total_tasks += n_series[f];
    }

    Task *tasks = malloc((size_t)total_tasks * sizeof(Task));
    long t_idx = 0;
    for (int f = 0; f < n_folders; f++)
        for (int s = 0; s < n_series[f]; s++) tasks[t_idx++] = (Task){ f, s };

    double *gain = malloc((size_t)total_tasks * sizeof(double));
    double *base_log_lik = malloc((size_t)total_tasks * sizeof(double));
    double *base_gradient = malloc((size_t)total_tasks * sizeof(double));
    double *high_gradient = malloc((size_t)total_tasks * sizeof(double));
    double *distance = malloc((size_t)total_tasks * sizeof(double));
    double *relative = malloc((size_t)total_tasks * sizeof(double));
    double *relative_no_nu = malloc((size_t)total_tasks * sizeof(double));
    double *largest = malloc((size_t)total_tasks * sizeof(double));
    int *base_niter = malloc((size_t)total_tasks * sizeof(int));
    int *high_niter = malloc((size_t)total_tasks * sizeof(int));
    double *base_nu = malloc((size_t)total_tasks * sizeof(double));
    double *high_nu = malloc((size_t)total_tasks * sizeof(double));
    char *base_converged = malloc((size_t)total_tasks);
    char *high_converged = malloc((size_t)total_tasks);
    char *usable = malloc((size_t)total_tasks);
    char *divergent = malloc((size_t)total_tasks);
    double *block_largest = malloc((size_t)total_tasks * N_BLOCKS * sizeof(double));

    static const char *row_name[K] = {
        "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
    };

    #pragma omp parallel for schedule(dynamic, 64)
    for (long i = 0; i < total_tasks; i++) {
        Task task = tasks[i];
        const char *folder = folders[task.folder_index];
        usable[i] = 0;
        divergent[i] = 0;

        char csv_path[600];
        snprintf(csv_path, sizeof csv_path, "%s/%s/series_%03d.csv", INPUT_DIR, folder, task.series);
        DataFrame df = df_read_csv(csv_path, csv_read_options_default());
        Mat y = mat_new(K, df.r);
        for (int k = 0; k < K; k++) {
            Mat column = df_col_numeric(&df, row_name[k]);
            for (int t = 0; t < df.r; t++) AT(y, k, t) = AT(column, t, 0);
        }
        df_free(&df);

        char base_path[640], high_path[640];
        snprintf(base_path, sizeof base_path, "%s/%s/series_%03d_%s_fit.json",
                 BASE_DIR, folder, task.series, SPEC_LABEL);
        snprintf(high_path, sizeof high_path, "%s/%s/series_%03d_%s_fit.json",
                 HIGH_DIR, folder, task.series, SPEC_LABEL);

        QvarmaFitResult base_fit, high_fit;
        QvarmaParams base_params, high_params;
        if (load_params(&base_params, y, base_path, &base_fit)) {
            if (load_params(&high_params, y, high_path, &high_fit)) {
                double sum_squares, base_sum_squares, max_absolute;
                double sum_squares_no_nu, base_sum_squares_no_nu;
                compare_params(&base_params, &high_params, &sum_squares, &base_sum_squares,
                               &sum_squares_no_nu, &base_sum_squares_no_nu,
                               &max_absolute, block_largest + i * N_BLOCKS);
                gain[i] = (double)(high_fit.log_likelihood - base_fit.log_likelihood);
                base_log_lik[i] = (double)base_fit.log_likelihood;
                base_gradient[i] = (double)base_fit.gradient_norm;
                high_gradient[i] = (double)high_fit.gradient_norm;
                distance[i] = sqrt(sum_squares);
                relative[i] = base_sum_squares > 0 ? sqrt(sum_squares / base_sum_squares) : 0;
                relative_no_nu[i] = base_sum_squares_no_nu > 0
                                  ? sqrt(sum_squares_no_nu / base_sum_squares_no_nu) : 0;
                largest[i] = max_absolute;
                base_niter[i] = base_fit.niter;
                high_niter[i] = high_fit.niter;
                base_converged[i] = (char)base_fit.is_converged;
                high_converged[i] = (char)high_fit.is_converged;
                base_nu[i] = (double)base_params.nu;
                high_nu[i] = (double)high_params.nu;
                divergent[i] = (char)(base_params.nu > NU_DIVERGENCE_LIMIT ||
                                      high_params.nu > NU_DIVERGENCE_LIMIT);
                usable[i] = 1;
                qvarma_fit_result_free(&high_fit);
            }
            qvarma_fit_result_free(&base_fit);
        }
        mat_free(y);
    }

    long n_usable = 0, n_missing = 0, n_divergent = 0, n_divergent_base = 0, n_divergent_high = 0;
    long base_converged_count = 0, high_converged_count = 0, newly_converged = 0;
    long negative_gain = 0, zero_gain = 0, used_extra_budget = 0;
    double gain_total = 0, worst_negative = 0;
    for (long i = 0; i < total_tasks; i++) {
        if (!usable[i]) { n_missing++; continue; }
        if (divergent[i]) {
            n_divergent++;
            if (base_nu[i] > NU_DIVERGENCE_LIMIT) n_divergent_base++;
            if (high_nu[i] > NU_DIVERGENCE_LIMIT) n_divergent_high++;
            continue;
        }
        n_usable++;
        gain_total += gain[i];
        if (base_converged[i]) base_converged_count++;
        if (high_converged[i]) high_converged_count++;
        if (!base_converged[i] && high_converged[i]) newly_converged++;
        if (gain[i] < -1e-9) { negative_gain++; if (gain[i] < worst_negative) worst_negative = gain[i]; }
        else if (gain[i] <= 1e-9) zero_gain++;
        if (high_niter[i] > BASE_ITERATIONS) used_extra_budget++;
    }

    /* Quantiles over the fits that actually moved: a distribution dominated by
       the fits that had already converged says nothing about the ones that had
       not. Both populations are reported. */
    double *sorted_gain = malloc((size_t)n_usable * sizeof(double));
    double *sorted_relative = malloc((size_t)n_usable * sizeof(double));
    double *sorted_largest = malloc((size_t)n_usable * sizeof(double));
    double *sorted_gain_moved = malloc((size_t)n_usable * sizeof(double));
    double *sorted_relative_moved = malloc((size_t)n_usable * sizeof(double));
    double *sorted_relative_no_nu = malloc((size_t)n_usable * sizeof(double));
    long n_all = 0, n_moved = 0;
    for (long i = 0; i < total_tasks; i++) {
        if (!usable[i] || divergent[i]) continue;
        sorted_gain[n_all] = gain[i];
        sorted_relative[n_all] = relative[i];
        sorted_relative_no_nu[n_all] = relative_no_nu[i];
        sorted_largest[n_all] = largest[i];
        n_all++;
        if (high_niter[i] > BASE_ITERATIONS) {
            sorted_gain_moved[n_moved] = gain[i];
            sorted_relative_moved[n_moved] = relative[i];
            n_moved++;
        }
    }
    qsort(sorted_gain, (size_t)n_all, sizeof(double), compare_doubles);
    qsort(sorted_relative, (size_t)n_all, sizeof(double), compare_doubles);
    qsort(sorted_largest, (size_t)n_all, sizeof(double), compare_doubles);
    qsort(sorted_relative_no_nu, (size_t)n_all, sizeof(double), compare_doubles);
    qsort(sorted_gain_moved, (size_t)n_moved, sizeof(double), compare_doubles);
    qsort(sorted_relative_moved, (size_t)n_moved, sizeof(double), compare_doubles);

    double block_worst[N_BLOCKS] = {0};
    double block_sum[N_BLOCKS] = {0};
    for (long i = 0; i < total_tasks; i++) {
        if (!usable[i] || divergent[i]) continue;
        for (int b = 0; b < N_BLOCKS; b++) {
            double v = block_largest[i * N_BLOCKS + b];
            block_sum[b] += v;
            if (v > block_worst[b]) block_worst[b] = v;
        }
    }

    double base_seconds = elapsed_from_timing(BASE_TIMING, BASE_PROVENANCE);
    double high_seconds = elapsed_from_timing(HIGH_TIMING, HIGH_PROVENANCE);

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "cannot open the report path for writing");
    fprintf(report, "t-QVARMA(%d,%d,%d) at a solver budget of %d iterations against %d,\n",
            P, Q, RLAG, BASE_ITERATIONS, HIGH_ITERATIONS);
    fprintf(report, "over the same %ld series under %s.\n\n", total_tasks, INPUT_DIR);
    fprintf(report, "%s\n%s\n\n", BASE_DIR, HIGH_DIR);

    const char *base_fixed = carries_correction(BASE_PROVENANCE);
    const char *high_fixed = carries_correction(HIGH_PROVENANCE);
    if (strcmp(base_fixed, high_fixed) != 0) {
        fprintf(report,
          "READ THIS FIRST: the two runs did not use the same likelihood\n\n"
          "  large-nu correction, %d-iteration run: %s\n"
          "  large-nu correction, %d-iteration run: %s\n\n"
          "  et_al's t-QVARMA density was corrected at a large nu on 2026-08-31.\n"
          "  The two objectives agree to about 1e-8 wherever nu is moderate, so\n"
          "  most of what follows is still the budget. But an unconverged search\n"
          "  amplifies a perturbation that small: measured on this project's own\n"
          "  data, two builds agreeing to 6e-13 in the gradient moved the fitted\n"
          "  log-likelihood by a mean |difference| of 7.2 across 10,800 series,\n"
          "  against a mean gain of 1.45 for an entire doubling of the budget.\n"
          "  Everything below therefore mixes the two changes, and the build is\n"
          "  the larger of them. Treat the numbers as an upper bound on what the\n"
          "  extra iterations bought, not as a measurement of it.\n\n"
          "  Fitting both budgets against one et_al is what would separate them.\n\n",
          BASE_ITERATIONS, base_fixed, HIGH_ITERATIONS, high_fixed);
    }

    fprintf(report, "Provenance of each run, as its own stamp recorded it\n\n");
    quote_provenance(report, "base", BASE_PROVENANCE);
    quote_provenance(report, "high", HIGH_PROVENANCE);
    if (n_missing)
        fprintf(report, "%ld series had a fit in only one of the two trees and are left out "
                        "of everything below.\n\n", n_missing);

    fprintf(report, "Fits left out\n\n");
    fprintf(report, "  nu enters as exp(theta) + 2, and past nu of about 3.6e9 the coded\n"
                    "  log-likelihood loses precision - the density's (nu + K)/2 log(1 + q/nu)\n"
                    "  term is a cancellation of two large quantities - while the likelihood\n"
                    "  itself is flat there, at the Gaussian limit the t tends to. L-BFGS climbs\n"
                    "  what is left, and a larger budget gives it more room to. Fits carrying\n"
                    "  nu above %g at either budget are left out of everything below.\n\n",
            (double)NU_DIVERGENCE_LIMIT);
    fprintf(report, "  left out                    %ld of %ld (%.4f%%)\n", n_divergent,
            n_divergent + n_usable, 100.0 * (double)n_divergent / (double)(n_divergent + n_usable));
    fprintf(report, "  already over the limit at %d  %ld\n", BASE_ITERATIONS, n_divergent_base);
    fprintf(report, "  over the limit at %d         %ld\n", HIGH_ITERATIONS, n_divergent_high);
    fprintf(report, "  compared below              %ld\n\n", n_usable);

    fprintf(report, "What it cost\n\n");
    if (base_seconds > 0 && high_seconds > 0) {
        fprintf(report, "  %d iterations   %9.1f s   %.2f h\n", BASE_ITERATIONS, base_seconds,
                base_seconds / 3600);
        fprintf(report, "  %d iterations   %9.1f s   %.2f h\n", HIGH_ITERATIONS, high_seconds,
                high_seconds / 3600);
        fprintf(report, "  ratio             %8.2fx   %+.2f h, %+.4f s per fit\n\n",
                high_seconds / base_seconds, (high_seconds - base_seconds) / 3600,
                (high_seconds - base_seconds) / (double)n_usable);
    } else {
        fprintf(report, "  one of the two timing files is missing or carries no elapsed line\n\n");
    }
    fprintf(report, "  %ld of %ld fits (%.2f%%) ran past iteration %d at all; the rest had already\n"
                    "  stopped, so the larger budget cost them nothing.\n\n",
            used_extra_budget, n_usable, 100.0 * (double)used_extra_budget / (double)n_usable,
            BASE_ITERATIONS);

    fprintf(report, "What it bought, on the objective\n\n");
    fprintf(report, "  The two runs share a start, data and solver, so the larger budget retraces\n"
                    "  the smaller one's path and the gain cannot be negative. Negative gains are\n"
                    "  a check on that, not a result: %ld found", negative_gain);
    if (negative_gain) fprintf(report, ", worst %.6g", worst_negative);
    fprintf(report, ".\n\n");
    fprintf(report, "  total log-likelihood gained   %.4f over %ld fits\n", gain_total, n_usable);
    fprintf(report, "  mean gain per fit             %.6f\n", gain_total / (double)n_usable);
    fprintf(report, "  fits with no gain at all      %ld (%.2f%%)\n\n", zero_gain,
            100.0 * (double)zero_gain / (double)n_usable);
    fprintf(report, "  log-likelihood gain, quantiles\n");
    fprintf(report, "  %-14s %14s %14s\n", "quantile", "all fits", "fits that ran on");
    const double quantile_points[] = { 0.0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0 };
    const char *quantile_label[] = { "min", "p25", "median", "p75", "p90", "p99", "max" };
    for (int q = 0; q < 7; q++) {
        fprintf(report, "  %-14s %14.6f", quantile_label[q],
                quantile(sorted_gain, n_all, quantile_points[q]));
        if (n_moved) fprintf(report, " %14.6f\n", quantile(sorted_gain_moved, n_moved, quantile_points[q]));
        else fprintf(report, " %14s\n", "-");
    }

    fprintf(report, "\nWhat it bought, on convergence\n\n");
    fprintf(report, "  converged at %d      %ld of %ld (%.2f%%)\n", BASE_ITERATIONS,
            base_converged_count, n_usable, 100.0 * (double)base_converged_count / (double)n_usable);
    fprintf(report, "  converged at %d      %ld of %ld (%.2f%%)\n", HIGH_ITERATIONS,
            high_converged_count, n_usable, 100.0 * (double)high_converged_count / (double)n_usable);
    long not_converged_at_base = n_usable - base_converged_count;
    fprintf(report, "  newly converged      %ld", newly_converged);
    if (not_converged_at_base)
        fprintf(report, " (%.2f%% of the %ld that had not)", 
                100.0 * (double)newly_converged / (double)not_converged_at_base, not_converged_at_base);
    fprintf(report, "\n\n");

    fprintf(report, "How far the estimates moved\n\n");
    fprintf(report, "  Distance on the constrained parameters, the scale the paper names them on.\n");
    fprintf(report, "  relative is that distance over the norm of the %d-iteration estimate.\n\n",
            BASE_ITERATIONS);
    fprintf(report, "  nu is reported on its own as well as inside the total, because it is the\n"
                    "  one coordinate the likelihood is nearly flat in once it is large - the t\n"
                    "  is already Gaussian to the data's own precision - so a search wandering\n"
                    "  up it moves the estimate a long way without moving the fit. The column\n"
                    "  without nu is the one that says whether the parameters the impulse\n"
                    "  responses are built from actually changed.\n\n");
    fprintf(report, "  %-14s %16s %16s %16s %16s\n", "quantile", "relative, all",
            "relative, ran on", "relative, no nu", "largest single");
    for (int q = 0; q < 7; q++) {
        fprintf(report, "  %-14s %16.8f", quantile_label[q],
                quantile(sorted_relative, n_all, quantile_points[q]));
        if (n_moved) fprintf(report, " %16.8f", quantile(sorted_relative_moved, n_moved, quantile_points[q]));
        else fprintf(report, " %16s", "-");
        fprintf(report, " %16.8f %16.8f\n", quantile(sorted_relative_no_nu, n_all, quantile_points[q]),
                quantile(sorted_largest, n_all, quantile_points[q]));
    }

    fprintf(report, "\n  Largest absolute change inside each named block, so a distance can be\n"
                    "  attributed rather than left as a number.\n\n");
    fprintf(report, "  %-12s %16s %16s\n", "block", "mean of maxima", "worst single");
    for (int b = 0; b < N_BLOCKS; b++)
        fprintf(report, "  %-12s %16.8f %16.8f\n", block_name[b],
                block_sum[b] / (double)n_usable, block_worst[b]);

    fprintf(report, "\nPer-series rows are in %s.\n", CSV_PATH);
    fclose(report);

    FILE *csv = fopen(CSV_PATH, "w");
    assert(csv && "cannot open the csv path for writing");
    fprintf(csv, "folder,series,log_lik_%d,gain,gradient_%d,gradient_%d,niter_%d,niter_%d,"
                 "converged_%d,converged_%d,nu_%d,nu_%d,divergent,distance,relative,"
                 "relative_no_nu,largest",
            BASE_ITERATIONS, BASE_ITERATIONS, HIGH_ITERATIONS, BASE_ITERATIONS, HIGH_ITERATIONS,
            BASE_ITERATIONS, HIGH_ITERATIONS, BASE_ITERATIONS, HIGH_ITERATIONS);
    for (int b = 0; b < N_BLOCKS; b++) fprintf(csv, ",largest_%s", block_name[b]);
    fputc('\n', csv);
    for (long i = 0; i < total_tasks; i++) {
        if (!usable[i]) continue;
        Task task = tasks[i];
        fprintf(csv, "%s,%d,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g",
                folders[task.folder_index], task.series, base_log_lik[i], gain[i],
                base_gradient[i], high_gradient[i], base_niter[i], high_niter[i],
                base_converged[i], high_converged[i], base_nu[i], high_nu[i], divergent[i],
                distance[i], relative[i], relative_no_nu[i], largest[i]);
        for (int b = 0; b < N_BLOCKS; b++) fprintf(csv, ",%.17g", block_largest[i * N_BLOCKS + b]);
        fputc('\n', csv);
    }
    fclose(csv);

    for (int f = 0; f < n_folders; f++) free(folders[f]);
    free(folders);
    free(n_series);
    free(tasks);
    free(gain); free(base_log_lik); free(base_gradient); free(high_gradient);
    free(distance); free(relative); free(relative_no_nu); free(largest);
    free(base_niter); free(high_niter);
    free(base_nu); free(high_nu);
    free(base_converged); free(high_converged); free(usable); free(divergent);
    free(block_largest);
    free(sorted_gain); free(sorted_relative); free(sorted_largest);
    free(sorted_gain_moved); free(sorted_relative_moved); free(sorted_relative_no_nu);
    return 0;
}
