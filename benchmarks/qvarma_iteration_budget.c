/*
What the 2000 iteration cap in applications/abm_system_fit_qvarma.c buys, per
iteration, so that shortening it can be judged against the log-likelihood it
would give up rather than assumed. Not part of the pipeline; run explicitly.

For each of several replicates and both specs, records the full L-BFGS trace
and reports the iteration at which the log-likelihood first comes within
0.01, 0.1 and 1.0 of the value the fit reaches at iteration 2000.
Writes out/qvarma_iteration_budget.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <string.h>
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

typedef struct { const char *sample; int replicate; } Series;

static Mat load_series(Series series) {
    char path[560];
    snprintf(path, sizeof path, "dataset/abm_system/%s", series.sample);
    return abm_system_read_replicate(path, series.replicate);
}

/* Pull "iter" and "f" out of the solver's own trace lines. */
static int read_trace(const char *path, double *value, int cap) {
    FILE *f = fopen(path, "r");
    assert(f);
    char line[512];
    int n = 0;
    while (fgets(line, sizeof line, f)) {
        int iter; double v;
        if (sscanf(line, "  iter %d  f %lg", &iter, &v) == 2 && iter < cap) {
            value[iter] = v;
            if (iter + 1 > n) n = iter + 1;
        }
    }
    fclose(f);
    return n;
}

static int first_within(const double *nll, int n, double final, double gap) {
    for (int i = 0; i < n; i++)
        if (nll[i] - final <= gap) return i;
    return -1;
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
    static const Series replicate[] = {
        { "EstimationSeriesSample1_1", 0 },
        { "EstimationSeriesSample1_1", 17 },
        { "EstimationSeriesSample1_50", 4 },
        { "EstimationSeriesSample1_100", 31 }
    };
    int n_replicates = (int)(sizeof replicate / sizeof replicate[0]);
    const int r_list[2] = { 2, 4 };

    FILE *report = fopen("out/qvarma_iteration_budget.txt", "w");
    assert(report);
    fprintf(report, "Iterations needed to come within a given log-likelihood of what "
                    "%d iterations reach.\n", MAX_ITERATIONS);
    fprintf(report, "nll is the negative log-likelihood the solver minimises; "
                    "gaps are in the same units.\n\n");
    fprintf(report, "%-46s %4s %14s %14s %8s %8s %8s\n", "series", "r",
            "nll at 2000", "nll at 200", "within1", "within.1", "within.01");

    double *nll = malloc((size_t)MAX_ITERATIONS * sizeof(double));
    for (int i = 0; i < n_replicates; i++) {
        Mat y = load_series(replicate[i]);
        for (int s = 0; s < 2; s++) {
            QvarmaParams start = build_start(y, r_list[s]);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = MAX_ITERATIONS;
            options.trace = fopen("out/qvarma_iteration_budget_trace.txt", "w");
            QvarmaFitResult result = qvarma_fit(y, &start, options);
            fclose(options.trace);

            int n = read_trace("out/qvarma_iteration_budget_trace.txt", nll, MAX_ITERATIONS);
            double final = -(double)result.log_likelihood;
            char label[128];
            snprintf(label, sizeof label, "%s replicate %d",
                     replicate[i].sample, replicate[i].replicate);
            fprintf(report, "%-46s %4d %14.4f %14.4f %8d %8d %8d\n", label, r_list[s],
                    final, n > 200 ? nll[200] : nll[n - 1],
                    first_within(nll, n, final, 1.0),
                    first_within(nll, n, final, 0.1),
                    first_within(nll, n, final, 0.01));
            fflush(report);
            qvarma_fit_result_free(&result);
            qvarma_params_free(&start);
        }
        mat_free(y);
    }
    free(nll);
    fclose(report);
    remove("out/qvarma_iteration_budget_trace.txt");
    return 0;
}
