/*
What et_al.'s analytic filter is worth against the traced one it replaced, and
whether either scales across cores. A benchmark, not a correctness gate,
though it does check that the two compute the same number before timing
either: a fast filter that answers a different question is not a speedup.

Both paths come from et_al. The traced one is _qvarma_link plus
_qvarma_filter plus tape_backward, which is what qvarma_negative_log_likelihood
called before the analytic path existed; the analytic one is
qvarma_analytic_log_likelihood against a workspace built once and reused, which
is what qvarma_fit does now.

Value-only and value-and-gradient are timed separately because a line search
asks for the value alone, and that is most of the objective calls a fit makes.

The spec is t-QVARMA(1,1,2), the one the ABM pipeline fits.
Writes out/qvarma_taped_vs_fused.txt.
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
#define P 1
#define Q 1
#define SPEC_R 2
#define N_THREADS 4

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* The traced evaluation, kept alive here because qvarma_negative_log_likelihood
   no longer takes this path and there is nothing left to compare against
   otherwise. Same sign convention as the analytic one: the log-likelihood, not
   its negation. */
static mreal taped_log_likelihood(Vec theta, const QvarmaParams *shape, Mat y, Vec gradient) {
    Tape *tape = tape_new();
    Node *theta_node = ad_leaf(tape, theta);
    QvarmaLinked linked = _qvarma_link(tape, theta_node, shape);
    Node *objective = _qvarma_filter(tape, &linked, shape, y, NULL, NULL, NULL);
    mreal value = objective->val.d[0];
    if (gradient.d) {
        tape_backward(tape, objective);
        for (int i = 0; i < theta.r; i++) gradient.d[i] = theta_node->grad.d[i];
    }
    qvarma_linked_free(&linked);
    tape_free(tape);
    return value;
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

/* The tape frees its blocks on every evaluation, and handing those back to the
   kernel costs more than the arithmetic does. Set here rather than left to
   MALLOC_TRIM_THRESHOLD_ in the environment, so what this writes does not
   depend on how it was launched. */
static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

/* One worker's share of a timed loop, so the one-thread and four-thread runs
   execute identical code and differ only in how many run at once. */
static void run_evaluations(int fused, int want_gradient, int repeats, Vec theta,
                            const QvarmaParams *shape, Mat y) {
    int n = theta.r;
    Vec local = mat_copy(theta);
    Vec gradient = want_gradient ? mat_new(n, 1) : (Vec){ 0, 0, 0, NULL };
    QvarmaAnalytic *workspace = fused ? qvarma_analytic_new(shape, y.c) : NULL;
    volatile mreal sink = 0;
    for (int i = 0; i < repeats; i++)
        sink += fused ? qvarma_analytic_log_likelihood(workspace, local, y, gradient)
                      : taped_log_likelihood(local, shape, y, gradient);
    (void)sink;
    if (workspace) qvarma_analytic_free(workspace);
    if (want_gradient) mat_free(gradient);
    mat_free(local);
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);

    Mat y = abm_system_read_replicate("dataset/abm_system/EstimationSeriesSample1_1", 0);
    int T = y.c;

    QvarmaParams shape = build_start(y, SPEC_R);
    int n = qvarma_n_theta(&shape);
    Vec theta = mat_new(n, 1);
    _qvarma_unlink(&shape, theta);

    /* Agreement before timing. Both values and both gradients, at the point
       the pipeline starts every fit from. */
    Vec taped_gradient = mat_new(n, 1), fused_gradient = mat_new(n, 1);
    mreal taped_value = taped_log_likelihood(theta, &shape, y, taped_gradient);
    QvarmaAnalytic *check = qvarma_analytic_new(&shape, T);
    mreal fused_value = qvarma_analytic_log_likelihood(check, theta, y, fused_gradient);
    qvarma_analytic_free(check);
    mreal worst_gradient = 0;
    for (int i = 0; i < n; i++) {
        mreal gap = MABS(taped_gradient.d[i] - fused_gradient.d[i]);
        mreal scale = MABS(taped_gradient.d[i]);
        mreal relative = scale > 1 ? gap / scale : gap;
        if (relative > worst_gradient) worst_gradient = relative;
    }

    FILE *report = fopen("out/qvarma_taped_vs_fused.txt", "w");
    assert(report && "qvarma_taped_vs_fused: cannot open the report path");
    fprintf(report, "series: EstimationSeriesSample1_1 replicate 0\n");
    fprintf(report, "K = %d, T = %d, spec r = %d, %d parameters\n\n", K, T, SPEC_R, n);
    fprintf(report, "agreement at the pipeline's own starting theta\n");
    fprintf(report, "  log-likelihood  taped %.10f  analytic %.10f  gap %.3g\n",
            (double)taped_value, (double)fused_value,
            (double)MABS(taped_value - fused_value));
    fprintf(report, "  gradient        worst relative difference over %d coordinates %.3g\n\n",
            n, (double)worst_gradient);

    fprintf(report, "%-24s %8s %8s %12s %10s\n", "evaluation", "path", "threads",
            "per_eval_ms", "scaling");

    struct { const char *name; int gradient; int taped_repeats; int fused_repeats; } job[2] = {
        { "value only", 0, 400, 4000 },
        { "value and gradient", 1, 400, 4000 }
    };

    double per_eval[2][2][2];
    for (int j = 0; j < 2; j++) {
        for (int fused = 0; fused < 2; fused++) {
            int repeats = fused ? job[j].fused_repeats : job[j].taped_repeats;
            double one_thread = 0;
            for (int c = 0; c < 2; c++) {
                int threads = c == 0 ? 1 : N_THREADS;
                /* Warm the allocator and the caches at this configuration. */
                run_evaluations(fused, job[j].gradient, 5, theta, &shape, y);
                double t0 = now_seconds();
                #pragma omp parallel num_threads(threads)
                run_evaluations(fused, job[j].gradient, repeats, theta, &shape, y);
                double wall = now_seconds() - t0;
                double each = wall / repeats;
                if (c == 0) one_thread = each;
                per_eval[j][fused][c] = each;
                fprintf(report, "%-24s %8s %8d %12.4f %10.2f\n", job[j].name,
                        fused ? "fused" : "taped", threads, 1e3 * each,
                        threads * one_thread / each);
            }
        }
    }

    fprintf(report, "\nanalytic against taped, same thread count\n");
    for (int j = 0; j < 2; j++)
        for (int c = 0; c < 2; c++)
            fprintf(report, "  %-20s %2d threads  %6.1fx\n", job[j].name,
                    c == 0 ? 1 : N_THREADS, per_eval[j][0][c] / per_eval[j][1][c]);

    fclose(report);
    mat_free(taped_gradient);
    mat_free(fused_gradient);
    mat_free(theta);
    qvarma_params_free(&shape);
    mat_free(y);
    return 0;
}
