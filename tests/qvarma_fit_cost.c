/*
Where the wall time of one qvarma fit goes, so that a plan for a 500-folder run
can be ranked against measured numbers rather than guesses. Not part of the
pipeline; run explicitly. Writes out/qvarma_fit_cost.txt.

Measures, on one replicate of dataset/abm_system, for both specs the pipeline
fits (r = 2 and r = 4):
  - one objective evaluation with the gradient (forward tape plus backward pass)
  - one objective evaluation without it (what a line search costs)
  - a full fit at the pipeline's own 2000 iteration cap, single threaded
  - the iterations that fit actually used and whether it converged
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
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

int main(int argc, char **argv) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);
    const char *sample = argc > 1 ? argv[1] : "EstimationSeriesSample1_1";
    int replicate = argc > 2 ? atoi(argv[2]) : 0;
    int repeats = argc > 3 ? atoi(argv[3]) : 200;

    char csv_path[560];
    snprintf(csv_path, sizeof csv_path, "dataset/abm_system/%s", sample);
    Mat y = abm_system_read_replicate(csv_path, replicate);

    FILE *report = fopen("out/qvarma_fit_cost.txt", "w");
    assert(report);
    fprintf(report, "series: %s replicate %d\n", sample, replicate);
    fprintf(report, "K = %d, T = %d, evaluation repeats = %d, single threaded\n\n", K, y.c, repeats);
    fprintf(report, "%-6s %8s %14s %14s %12s %8s %10s %12s\n", "spec", "n_theta",
            "value+grad_ms", "value_only_ms", "fit_s", "niter", "converged", "grad_norm");

    const int r_list[2] = { 2, 4 };
    for (int s = 0; s < 2; s++) {
        QvarmaParams start = build_start(y, r_list[s]);
        int n = qvarma_n_theta(&start);
        Vec theta = mat_new(n, 1);
        _qvarma_unlink(&start, theta);
        Vec gradient = mat_new(n, 1);
        /* The workspace qvarma_fit builds and reuses across every evaluation
           of one fit. Leaving it NULL is allowed and makes each call build its
           own, which is not what a fit does and not what this is timing. */
        QvarmaFitContext context = { y, &start, qvarma_analytic_new(&start, y.c) };

        /* warm the allocator and the caches before timing */
        for (int i = 0; i < 5; i++) qvarma_negative_log_likelihood(theta, gradient, &context);

        double t0 = now_seconds();
        for (int i = 0; i < repeats; i++) qvarma_negative_log_likelihood(theta, gradient, &context);
        double with_gradient = (now_seconds() - t0) / repeats;

        Vec no_gradient = { 0 };
        t0 = now_seconds();
        for (int i = 0; i < repeats; i++) qvarma_negative_log_likelihood(theta, no_gradient, &context);
        double value_only = (now_seconds() - t0) / repeats;

        QvarmaFitOptions options = qvarma_default_fit_options();
        options.max_iterations = MAX_ITERATIONS;
        t0 = now_seconds();
        QvarmaFitResult result = qvarma_fit(y, &start, options);
        double fit_seconds = now_seconds() - t0;

        fprintf(report, "%-6d %8d %14.4f %14.4f %12.3f %8d %10s %12.4g\n", r_list[s], n,
                1e3 * with_gradient, 1e3 * value_only, fit_seconds, result.niter,
                result.is_converged ? "yes" : "no", (double)result.gradient_norm);
        fflush(report);

        qvarma_fit_result_free(&result);
        qvarma_analytic_free(context.workspace);
        mat_free(gradient);
        mat_free(theta);
        qvarma_params_free(&start);
    }
    fclose(report);
    mat_free(y);
    return 0;
}
