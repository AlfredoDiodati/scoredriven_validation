/*
How long 500,000 t-QVARMA(1,1,2) fits take on this machine. The auxiliary
model of an indirect-inference loop over an ABM has to be estimated once
per simulated dataset, so what decides whether that loop is feasible is
the wall time of the whole battery, not the cost of one fit. This script
runs the battery and writes the two timestamps that answer the question.

Fits every series applications/abm_system_scale_extract.c wrote under
dataset/abm_system_scale/ - 500 folders of 1000 - with the single spec the
pipeline settled on: t-QVARMA(1,1,2), K_star 3, R 1, shared beta, the same
K_dagger 2 partition and the same build_start convention
applications/abm_system_fit_qvarma.c and
applications/us_qvarma_employment_change.c use. One spec, not the two of
abm_system_fit_qvarma.c, because 500 x 1000 is the size being timed and r
= 4 is not part of it.

Every fitted parameter set is written the moment that one optimization
returns, by the fit's own worker thread, never batched: 500,000 fits are
500,000 independent writes to
out/abm_system_scale_fit_qvarma_i<cap>/<folder>/series_<NNN>_p1q1r2_fit.json,
mirroring the input layout exactly so which series a cached fit belongs
to is never in question. An interrupted run therefore keeps everything it
had finished.

The <cap> in that path is MAX_ITERATIONS, and every output this script
writes carries it. 86.85% of the 500,000 fits at a cap of 2000 stopped at
the cap rather than at the gradient tolerance, so what a larger budget buys
- in wall time and in how far the estimates move - is a question the runs
themselves have to answer, which means two budgets have to coexist on disk
rather than one overwriting the other. Set the budget at compile time
(-DMAX_ITERATIONS=4000); the Makefile generates one target per entry in
SCALE_ITERATION_CAPS and names the binary after the cap too, so make cannot
serve a binary built at one budget for a request at another.
applications/abm_system_scale_iteration_comparison.c reads two such trees
and reports the difference.

out/abm_system_scale_fit_qvarma_i<cap>_timing.txt is the point of the
script. It is opened and the start timestamp flushed before the first fit
begins, a progress line is appended every PROGRESS_EVERY completed fits,
and the end timestamp, the elapsed wall time and the throughput are
appended when the loop returns. Because it is written as the run proceeds
rather than at the end, an interrupted run still leaves a readable record
of how far it got. It is opened with "w", so rerunning a budget destroys
that budget's own previous measurement - the reason the budget is in the
name at all.

The timing file also records how many of the fits were actually estimated
and how many were served from an existing cache, because those two are not
the same measurement. A cache file that loads for the series in question is
reused as it stands, converged or not - this is a throughput test, not the
estimation pipeline, so unlike abm_system_fit_qvarma.c nothing here carries
an unconverged fit on for another budget. A rerun over a complete cache
therefore measures the I/O, not the fits, and the timing file says which of
the two a given run was.

Parallelized with OpenMP over the 500,000 series, schedule(dynamic, 8):
dynamic because a fit that meets the tolerance early and one that runs to
the iteration cap differ by an order of magnitude in cost, chunked at 8 so
that half a million scheduling decisions do not themselves become the
contended resource. Thread count is whatever OMP_NUM_THREADS says, the
machine's own count by default, and is recorded in the timing file - a
throughput number without it means nothing.

openblas_set_num_threads(1) for the reason abm_system_fit_qvarma.c's own
header comment gives: without it OpenBLAS spawns a thread pool inside every
OpenMP worker and the outer loop contends with itself. K = 5 here, so
BLAS-level threading has nothing to work with anyway.

Memory: the per-fit working set is a 5 x 400 series plus a 42-parameter
L-BFGS history, kilobytes, and the CSV and DataFrame behind it are freed
before the next task starts. What scales with the battery is the summary,
four arrays of one entry per fit, about 12 MB at 500,000, and the task
list at 4 MB. Nothing else is held across tasks.

out/abm_system_scale_fit_qvarma_i<cap>_manifest.txt records log-likelihood,
gradient norm, iterations, convergence and whether the fit was estimated
or reused, one line per fit. Written once at the end, a summary of what
the individual JSONs already hold rather than a substitute for them, and
about 45 MB at this size.

Requires dataset/abm_system_scale/ to already exist. Deliberately not made
to depend on app-abm_system_scale_extract in the Makefile: that step writes
21 GB and rerunning the fits must not rewrite it. Nothing printed.
*/

#include "abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <omp.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define START_NU ((mreal)30)
#define P 1
#define Q 1
#define RLAG 2
#define SPEC_LABEL "p1q1r2"

#define PROGRESS_EVERY 10000

/* The solver budget, and the only thing that varies between runs of this
   script. Set it at compile time (-DMAX_ITERATIONS=4000); every path below is
   built from it, so two budgets write to two separate trees and a run at one
   budget can neither overwrite nor silently reuse the other's fits. */
#ifndef MAX_ITERATIONS
#define MAX_ITERATIONS 2000
#endif

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
#define BUDGET_SUFFIX "_i" STRINGIFY(MAX_ITERATIONS)

#define INPUT_DIR "dataset/abm_system_scale"
#define OUTPUT_DIR "out/abm_system_scale_fit_qvarma" BUDGET_SUFFIX
#define TIMING_PATH "out/abm_system_scale_fit_qvarma" BUDGET_SUFFIX "_timing.txt"
#define MANIFEST_PATH "out/abm_system_scale_fit_qvarma" BUDGET_SUFFIX "_manifest.txt"

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_scale_fit_qvarma: mkdir failed");
}

static int compare_names(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static char **list_subdirs(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_scale_fit_qvarma: cannot open dataset/abm_system_scale/ - run abm_system_scale_extract first");

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
    assert(handle && "abm_system_scale_fit_qvarma: cannot open a folder under dataset/abm_system_scale/");
    int n = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL)
        if (strncmp(entry->d_name, "series_", 7) == 0) n++;
    closedir(handle);
    return n;
}

static void timestamp_now(char *buffer, size_t size) {
    time_t now = time(NULL);
    struct tm *local = localtime(&now);
    strftime(buffer, size, "%Y-%m-%d %H:%M:%S", local);
}

static mreal first_difference_sd(Mat y, int row) {
    int periods = y.c;
    Mat difference = mat_new(1, periods - 1);
    for (int t = 1; t < periods; t++) AT(difference, 0, t - 1) = AT(y, row, t) - AT(y, row, t - 1);
    mreal sd = (mreal)sqrt((double)stats_var(difference));
    mat_free(difference);
    return sd;
}

static QvarmaParams build_start(Mat y, int rlag) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, rlag, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    int K_dag = K - K_STAR;

    for (int k = 0; k < K; k++) {
        if (k < K_STAR) {
            Mat row = mat_slice(y, k, k + 1, 0, y.c);
            AT(m.c, k, 0) = stats_mean(row);
        } else {
            AT(m.c, k, 0) = AT(y, k, 0);
        }
    }
    AT(m.Phi_star, 0, 0) = (mreal)0.3;
    for (int a = 0; a < K_STAR; a++) AT(m.Psi_star[0], a, a) = (mreal)0.05;

    for (int a = 0; a < K; a++)
        for (int b = 0; b <= a; b++)
            AT(m.Omega_inv, a, b) = (mreal)(b == a ? first_difference_sd(y, a) : 0);
    m.nu = START_NU;

    for (int l = 0; l < rlag; l++) {
        mreal decay = (mreal)pow(0.6, (double)l);
        for (int i = 0; i < K_dag; i++) AT(m.alpha[l], i, 0) = decay * (mreal)(0.10 - 0.04 * i);
    }
    for (int j = 0; j < K_dag; j++) AT(m.beta[0], 0, j) = (mreal)(j == 0 ? 1 : 0);
    AT(m.beta[0], 0, K_dag - 1) = (mreal)0.3;

    Vec theta = mat_new(qvarma_n_theta(&m), 1);
    _qvarma_unlink(&m, theta);
    qvarma_params_from_theta(theta, &m);
    mat_free(theta);
    return m;
}

/* The cache, if it holds a fit of this shape for this data. The shape has to be
   built before the load: qvarma_load_fit compares every field of it, and
   mu_star_stationary_only and phi_star_bound both change the length of theta it
   expects, so a cache written by this script only loads into params carrying
   this script's own settings. */
static int load_cached(QvarmaFitResult *out, Mat y, const char *path) {
    out->params = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    out->params.phi_star_bound = PHI_STAR_BOUND;
    out->params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    if (qvarma_load_fit(out, y, path)) return 1;
    qvarma_params_free(&out->params);
    return 0;
}

typedef struct { int folder_index, series; } Task;

int main(void) {
    /* See this file's own header comment for why. */
    openblas_set_num_threads(1);

    make_directory(OUTPUT_DIR);

    int n_folders;
    char **folders = list_subdirs(INPUT_DIR, &n_folders);
    assert(n_folders > 0 && "abm_system_scale_fit_qvarma: no folders under dataset/abm_system_scale/");

    int *n_series = malloc((size_t)n_folders * sizeof(int));
    long total_tasks = 0;
    for (int f = 0; f < n_folders; f++) {
        n_series[f] = count_series(folders[f]);
        char out_dir[600];
        snprintf(out_dir, sizeof out_dir, "%s/%s", OUTPUT_DIR, folders[f]);
        make_directory(out_dir);
        total_tasks += n_series[f];
    }
    assert(total_tasks > 0 && "abm_system_scale_fit_qvarma: no series found to fit");

    Task *tasks = malloc((size_t)total_tasks * sizeof(Task));
    long t_idx = 0;
    for (int f = 0; f < n_folders; f++)
        for (int s = 0; s < n_series[f]; s++) tasks[t_idx++] = (Task){ f, s };

    double *log_lik = malloc((size_t)total_tasks * sizeof(double));
    double *gradient_norm = malloc((size_t)total_tasks * sizeof(double));
    int *niter = malloc((size_t)total_tasks * sizeof(int));
    int *converged = malloc((size_t)total_tasks * sizeof(int));
    char *was_cached = malloc((size_t)total_tasks);

    FILE *timing = fopen(TIMING_PATH, "w");
    assert(timing && "cannot open the timing path for writing");
    char started_at[32], finished_at[32];
    timestamp_now(started_at, sizeof started_at);
    int n_threads = omp_get_max_threads();
    fprintf(timing, "t-QVARMA(%d,%d,%d) fits over %s\n", P, Q, RLAG, INPUT_DIR);
    fprintf(timing, "parameters to %s\n\n", OUTPUT_DIR);
    fprintf(timing, "folders          %d\n", n_folders);
    fprintf(timing, "fits             %ld\n", total_tasks);
    fprintf(timing, "threads          %d\n", n_threads);
    fprintf(timing, "max_iterations   %d\n", MAX_ITERATIONS);
    fprintf(timing, "spec             K %d, K_star %d, R %d, shared beta %d, nu start %g\n\n",
            K, K_STAR, R, SHARED_BETA, (double)START_NU);
    fprintf(timing, "start            %s\n", started_at);
    fflush(timing);

    long completed = 0;
    double wall_start = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic, 8)
    for (long i = 0; i < total_tasks; i++) {
        Task task = tasks[i];
        const char *folder = folders[task.folder_index];

        char csv_path[600];
        snprintf(csv_path, sizeof csv_path, "%s/%s/series_%03d.csv", INPUT_DIR, folder, task.series);

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

        char cache_path[640];
        snprintf(cache_path, sizeof cache_path, "%s/%s/series_%03d_%s_fit.json",
                 OUTPUT_DIR, folder, task.series, SPEC_LABEL);

        QvarmaFitResult result;
        if (load_cached(&result, y, cache_path)) {
            was_cached[i] = 1;
        } else {
            QvarmaParams start = build_start(y, RLAG);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = MAX_ITERATIONS;
            result = qvarma_fit(y, &start, options);
            qvarma_save_fit(&result, y, cache_path);
            qvarma_params_free(&start);
            was_cached[i] = 0;
        }

        log_lik[i] = (double)result.log_likelihood;
        gradient_norm[i] = (double)result.gradient_norm;
        niter[i] = result.niter;
        converged[i] = result.is_converged;

        qvarma_fit_result_free(&result);
        mat_free(y);

        #pragma omp critical
        {
            completed++;
            if (completed % PROGRESS_EVERY == 0) {
                double elapsed = omp_get_wtime() - wall_start;
                char now[32];
                timestamp_now(now, sizeof now);
                fprintf(timing, "progress         %s  %ld of %ld fits, %.0f s elapsed, "
                                "%.1f fits/s\n", now, completed, total_tasks, elapsed,
                        (double)completed / elapsed);
                fflush(timing);
            }
        }
    }

    double elapsed = omp_get_wtime() - wall_start;
    timestamp_now(finished_at, sizeof finished_at);

    long n_estimated = 0, n_reused = 0, n_converged = 0;
    for (long i = 0; i < total_tasks; i++) {
        if (was_cached[i]) n_reused++; else n_estimated++;
        if (converged[i]) n_converged++;
    }

    fprintf(timing, "end              %s\n\n", finished_at);
    fprintf(timing, "elapsed          %.1f s (%.2f h)\n", elapsed, elapsed / 3600.0);
    fprintf(timing, "per fit          %.4f s wall, %.4f s of thread time\n",
            elapsed / (double)total_tasks, elapsed * n_threads / (double)total_tasks);
    fprintf(timing, "throughput       %.1f fits/s over %d threads\n\n",
            (double)total_tasks / elapsed, n_threads);
    fprintf(timing, "estimated        %ld fits\n", n_estimated);
    fprintf(timing, "reused           %ld fits already cached from an earlier run\n", n_reused);
    fprintf(timing, "converged        %ld of %ld met the gradient tolerance; the rest stopped "
                    "at the %d iteration cap\n", n_converged, total_tasks, MAX_ITERATIONS);
    if (n_reused > 0)
        fprintf(timing, "\nThe elapsed time above covers %ld estimated fits and %ld cache loads, "
                        "so it is\na throughput number only to the extent the run was not "
                        "already cached.\n", n_estimated, n_reused);
    fclose(timing);

    FILE *manifest = fopen(MANIFEST_PATH, "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "%d folders, %ld fits of t-QVARMA(%d,%d,%d)\n\n", n_folders, total_tasks,
            P, Q, RLAG);
    fprintf(manifest, "origin is what this run did: estimated, fitted from build_start and "
                      "written to its own\nJSON; reused, loaded from a JSON an earlier run "
                      "wrote for the same data.\n\n");
    fprintf(manifest, "%-36s %8s %14s %10s %8s %10s %10s\n", "folder", "series", "log_lik",
            "gradient", "niter", "converged", "origin");
    for (long i = 0; i < total_tasks; i++) {
        Task task = tasks[i];
        fprintf(manifest, "%-36s %8d %14.4f %10.4g %8d %10s %10s\n", folders[task.folder_index],
                task.series, log_lik[i], gradient_norm[i], niter[i],
                converged[i] ? "yes" : "no", was_cached[i] ? "reused" : "estimated");
    }
    fclose(manifest);

    for (int f = 0; f < n_folders; f++) free(folders[f]);
    free(folders);
    free(n_series);
    free(tasks);
    free(log_lik);
    free(gradient_norm);
    free(niter);
    free(converged);
    free(was_cached);
    return 0;
}
