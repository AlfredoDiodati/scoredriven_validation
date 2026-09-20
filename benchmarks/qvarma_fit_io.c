/*
What the file handling around one fit costs, now that the fit itself is fast
enough for it to matter. A benchmark, not a correctness gate.

abm_system_fit_qvarma.c reads one replicate and, through
qvarma_fit_cached, writes one JSON per fit and reads it back on a rerun. At
the taped filter's cost those were rounding errors beside a four second fit.
They are not any more, and a rerun that hits the cache does nothing else at
all, so both directions are timed here.

Writes out/qvarma_fit_io.txt and, while running, one scratch JSON under out/
which it removes before returning.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <time.h>
#include <malloc.h>
#include <stdio.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define START_NU ((mreal)30)
#define MAX_ITERATIONS 2000
#define P 1
#define Q 1
#define SPEC_R 2
#define REPEATS 200
#define SCRATCH "out/qvarma_fit_io_scratch.json"

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
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

static Mat read_series(const char *path) {
    return abm_system_read_replicate(path, 0);
}

static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);
    const char *csv_path = "dataset/abm_system/EstimationSeriesSample1_1";  /* replicate 0 */

    Mat y = read_series(csv_path);
    QvarmaParams start = build_start(y, SPEC_R);
    QvarmaFitOptions options = qvarma_default_fit_options();
    options.max_iterations = MAX_ITERATIONS;

    double t0 = now_seconds();
    QvarmaFitResult result = qvarma_fit(y, &start, options);
    double fit_seconds = now_seconds() - t0;

    t0 = now_seconds();
    for (int i = 0; i < REPEATS; i++) { Mat r = read_series(csv_path); mat_free(r); }
    double csv_read = (now_seconds() - t0) / REPEATS;

    t0 = now_seconds();
    for (int i = 0; i < REPEATS; i++) qvarma_save_fit(&result, y, SCRATCH);
    double json_write = (now_seconds() - t0) / REPEATS;

    t0 = now_seconds();
    for (int i = 0; i < REPEATS; i++) {
        QvarmaFitResult loaded;
        loaded.params = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
        /* qvarma_load_params compares the whole shape, and both of these change
           the length of theta it expects, so a cache written by the pipeline
           only loads into params carrying the pipeline's own settings. */
        loaded.params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
        loaded.params.phi_star_bound = PHI_STAR_BOUND;
        int ok = qvarma_load_fit(&loaded, y, SCRATCH);
        assert(ok && "qvarma_fit_io: the fit just written did not load back");
        qvarma_fit_result_free(&loaded);
    }
    double json_read = (now_seconds() - t0) / REPEATS;

    FILE *report = fopen("out/qvarma_fit_io.txt", "w");
    assert(report && "qvarma_fit_io: cannot open the report path");
    fprintf(report, "series: %s\nK = %d, T = %d, spec r = %d, %d repeats per file "
                    "operation, single threaded\n", csv_path, K, y.c, SPEC_R, REPEATS);
    fprintf(report, "Files are already in the page cache: this is parse and format cost, "
                    "not a cold disk.\n\n");
    fprintf(report, "%-34s %12s %14s\n", "operation", "ms", "share of a fit");
    fprintf(report, "%-34s %12.4f %13.1f%%\n", "fit, 2000 iterations",
            1e3 * fit_seconds, 100.0);
    fprintf(report, "%-34s %12.4f %13.1f%%\n", "read the replicate",
            1e3 * csv_read, 100.0 * csv_read / fit_seconds);
    fprintf(report, "%-34s %12.4f %13.1f%%\n", "write the fit JSON",
            1e3 * json_write, 100.0 * json_write / fit_seconds);
    fprintf(report, "%-34s %12.4f %13.1f%%\n", "read the fit JSON back",
            1e3 * json_read, 100.0 * json_read / fit_seconds);
    fprintf(report, "\nper fit including its own file handling: %.4f ms\n",
            1e3 * (fit_seconds + csv_read + json_write));
    fprintf(report, "per cached fit, nothing recomputed:        %.4f ms\n",
            1e3 * (csv_read + json_read));
    fclose(report);

    remove(SCRATCH);
    qvarma_fit_result_free(&result);
    qvarma_params_free(&start);
    mat_free(y);
    return 0;
}
