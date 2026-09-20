/*
How close the OpenMP loop in applications/abm_system_fit_qvarma.c gets to
linear on this machine's four cores, measured on real fits rather than
assumed. Not part of the pipeline; run explicitly.

Fits the same eight (replicate, spec) pairs at one, two and four threads and
reports wall time and the speedup over one thread. openblas_set_num_threads(1)
is set exactly as the pipeline sets it. Writes out/qvarma_thread_scaling.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>

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
#define N_FITS 8

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

static Mat load_series(const char *path) {
    DataFrame df = df_read_csv(path, csv_read_options_default());
    Mat y = mat_new(K, df.r);
    static const char *row_name[K] = {
        "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
    };
    for (int k = 0; k < K; k++) {
        Mat column = df_col_numeric(&df, row_name[k]);
        for (int t = 0; t < df.r; t++) AT(y, k, t) = AT(column, t, 0);
    }
    df_free(&df);
    return y;
}

/* The tape frees its blocks on every objective evaluation, and handing those
   back to the kernel costs more than the arithmetic does: 22 million minor
   page faults over eight fits, against 2,331 once the arena is kept resident.
   Set here rather than left to MALLOC_TRIM_THRESHOLD_ in the environment, so
   the numbers this writes do not depend on how it was launched. The pipeline
   itself does not set these yet. */
static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);

    Mat series[N_FITS];
    int r_of[N_FITS];
    for (int i = 0; i < N_FITS; i++) {
        char path[256];
        snprintf(path, sizeof path, "dataset/abm_system/EstimationSeriesSample1_1/replicate_%03d.csv",
                 i / 2);
        series[i] = load_series(path);
        r_of[i] = (i % 2) ? 4 : 2;
    }

    FILE *report = fopen("out/qvarma_thread_scaling.txt", "w");
    assert(report);
    fprintf(report, "%d fits (replicates 0 to 3, both specs), T = %d, max_iterations = %d\n\n",
            N_FITS, series[0].c, MAX_ITERATIONS);
    fprintf(report, "%8s %12s %12s %12s\n", "threads", "wall_s", "speedup", "s_per_fit");

    double baseline = 0;
    const int thread_list[3] = { 1, 2, 4 };
    for (int c = 0; c < 3; c++) {
        omp_set_num_threads(thread_list[c]);
        double t0 = now_seconds();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N_FITS; i++) {
            QvarmaParams start = build_start(series[i], r_of[i]);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = MAX_ITERATIONS;
            QvarmaFitResult result = qvarma_fit(series[i], &start, options);
            qvarma_fit_result_free(&result);
            qvarma_params_free(&start);
        }
        double wall = now_seconds() - t0;
        if (c == 0) baseline = wall;
        fprintf(report, "%8d %12.3f %12.2f %12.3f\n", thread_list[c], wall,
                baseline / wall, wall / N_FITS);
        fflush(report);
    }
    fclose(report);
    for (int i = 0; i < N_FITS; i++) mat_free(series[i]);
    return 0;
}
