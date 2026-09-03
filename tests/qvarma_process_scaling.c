/*
Whether independent processes fit faster than the OpenMP threads
tests/qvarma_thread_scaling.c times, which separates contention on something
shared inside one process from a limit of the machine itself. A benchmark, not
a correctness gate.

Each copy fits the same number of replicates single threaded, so the four-copy
run does four times the work of the one-copy run and the speedup column is
copies times the one-copy time over the four-copy time. Copies are separate
processes rather than threads because OpenBLAS's buffer table, which
tests/small_call_scaling.c shows is the contended resource, is one structure
per process.

The spec is t-QVARMA(1,1,2) only, the one the ABM pipeline fits.

Writes out/qvarma_process_scaling.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
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
#define SPEC_R 2
#define FITS_PER_COPY 2

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

static Mat load_series(int replicate) {
    char path[256];
    snprintf(path, sizeof path,
             "dataset/abm_system/EstimationSeriesSample1_1/replicate_%03d.csv", replicate);
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

/* One copy's share: FITS_PER_COPY replicates, offset so no two copies fit the
   same series and any per-series difference in iteration count is spread. */
static void run_share(int copy) {
    openblas_set_num_threads(1);
    for (int f = 0; f < FITS_PER_COPY; f++) {
        Mat y = load_series(copy * FITS_PER_COPY + f);
        QvarmaParams start = build_start(y, SPEC_R);
        QvarmaFitOptions options = qvarma_default_fit_options();
        options.max_iterations = MAX_ITERATIONS;
        QvarmaFitResult result = qvarma_fit(y, &start, options);
        qvarma_fit_result_free(&result);
        qvarma_params_free(&start);
        mat_free(y);
    }
}

static double time_copies(int copies) {
    double t0 = now_seconds();
    for (int c = 0; c < copies; c++) {
        pid_t child = fork();
        assert(child >= 0 && "qvarma_process_scaling: fork failed");
        if (child == 0) { run_share(c); _exit(0); }
    }
    for (int c = 0; c < copies; c++) {
        int status;
        wait(&status);
    }
    return now_seconds() - t0;
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
    FILE *report = fopen("out/qvarma_process_scaling.txt", "w");
    assert(report && "qvarma_process_scaling: cannot open the report path");
    fprintf(report, "spec r = %d, %d fits per copy, max_iterations = %d, "
                    "one BLAS thread per copy\n\n", SPEC_R, FITS_PER_COPY, MAX_ITERATIONS);
    fprintf(report, "%8s %8s %10s %12s %10s\n", "copies", "fits", "wall_s", "s_per_fit", "speedup");

    double baseline = 0;
    const int copy_list[3] = { 1, 2, 4 };
    for (int i = 0; i < 3; i++) {
        int copies = copy_list[i];
        double wall = time_copies(copies);
        if (i == 0) baseline = wall;
        int fits = copies * FITS_PER_COPY;
        fprintf(report, "%8d %8d %10.3f %12.3f %10.2f\n", copies, fits, wall,
                wall / fits, copies * baseline / wall);
        fflush(report);
    }
    fclose(report);
    return 0;
}
