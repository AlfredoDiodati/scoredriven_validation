/*
Whether the fits that never converge are badly scaled or chasing a likelihood
with no maximum. A study, not a speed measurement and not a pass-or-fail check.

out/qvarma_stuck_fits.txt establishes that the resumed fits are not blocked:
a downhill step exists, it just sits ten orders of magnitude closer than the
line search looks, and reaching it gains 3e-4. That is what extreme curvature
looks like from the optimizer's side. This asks where the curvature comes
from, by running the same fit at growing iteration budgets from the same start
and watching four things move: the log-likelihood, the smallest diagonal of
Omega_inv, nu, and the gradient norm.

A residual standard deviation walking toward zero while the log-likelihood
keeps climbing is a likelihood with no interior maximum, and no amount of
optimizer tuning fixes one of those. A residual standard deviation settling
while the log-likelihood flattens is a scaling problem, which preconditioning
or a reparameterization does fix. The two have different remedies, so which
one this is has to be established before either is attempted.

At the largest budget the Hessian's eigenvalues are taken as well, for the
condition number and for which parameter blocks carry the stiff and the flat
directions.

The spec is t-QVARMA(1,1,2), the one the ABM pipeline fits.
Writes out/qvarma_conditioning.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <frame/csv.h>
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
#define P 1
#define Q 1
#define SPEC_R 2

static const int budget[] = { 250, 500, 1000, 2000, 4000, 8000, 16000, 32000 };
#define N_BUDGETS ((int)(sizeof budget / sizeof budget[0]))

static const char *replicate_path[] = {
    "dataset/abm_system/EstimationSeriesSample1_1/replicate_000.csv",
    "dataset/abm_system/EstimationSeriesSample1_1/replicate_017.csv",
    "dataset/abm_system/EstimationSeriesSample1_1/replicate_030.csv",
    "dataset/abm_system/EstimationSeriesSample1_50/replicate_004.csv",
    "dataset/abm_system/EstimationSeriesSample1_100/replicate_031.csv"
};
#define N_REPLICATES ((int)(sizeof replicate_path / sizeof replicate_path[0]))

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

static int smallest_diagonal_row(const QvarmaParams *m, double *value) {
    int at = 0;
    *value = (double)AT(m->Omega_inv, 0, 0);
    for (int k = 1; k < K; k++) {
        double d = (double)AT(m->Omega_inv, k, k);
        if (d < *value) { *value = d; at = k; }
    }
    return at;
}

static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);

    FILE *report = fopen("out/qvarma_conditioning.txt", "w");
    assert(report);
    fprintf(report, "Each fit starts from build_start and runs to the budget in the first "
                    "column.\nBudgets are independent runs, not a continued one.\n");
    fprintf(report, "spec r = %d, K = %d, series named per row below.\n", SPEC_R, K);
    fprintf(report, "min_sd is the smallest diagonal of Omega_inv, which is the residual "
                    "standard\ndeviation of the series named beside it.\n\n");

    static const char *series_name[K] = {
        "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
    };

    for (int rep = 0; rep < N_REPLICATES; rep++) {
        Mat y = load_series(replicate_path[rep]);
        fprintf(report, "%s\n", replicate_path[rep] + 21);
        fprintf(report, "  %8s %14s %12s %12s %14s %12s %10s\n", "budget", "log_lik",
                "min_sd", "at", "nu", "grad_norm", "status");

        Mat last_theta = { 0, 0, 0, NULL };
        QvarmaParams last_shape = build_start(y, SPEC_R);
        for (int b = 0; b < N_BUDGETS; b++) {
            QvarmaParams start = build_start(y, SPEC_R);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = budget[b];
            QvarmaFitResult result = qvarma_fit(y, &start, options);

            double min_sd;
            int at = smallest_diagonal_row(&result.params, &min_sd);
            fprintf(report, "  %8d %14.4f %12.6g %12s %14.6g %12.4g %10s\n", budget[b],
                    (double)result.log_likelihood, min_sd, series_name[at],
                    (double)result.params.nu, (double)result.gradient_norm,
                    result.is_converged ? "converged" : "capped");
            fflush(report);

            if (b == N_BUDGETS - 1) {
                if (last_theta.d) mat_free(last_theta);
                last_theta = mat_new(qvarma_n_theta(&result.params), 1);
                _qvarma_unlink(&result.params, last_theta);
                qvarma_params_free(&last_shape);
                last_shape = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA,
                                               WARMUP_LONGEST);
                last_shape.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
                last_shape.phi_star_bound = PHI_STAR_BOUND;
            }
            qvarma_fit_result_free(&result);
            qvarma_params_free(&start);
        }

        /* The curvature at the last point, which is what the optimizer is
           actually walking on. Eigenvalues of the negative log-likelihood's
           Hessian: at a maximum of the likelihood they are all positive. */
        Mat H = _qvarma_hessian(last_theta, &last_shape, y);
        int n = H.r;
        Vec eigenvalue = mat_new(n, 1);
        Mat copy = mat_copy(H);
        int ok = _syevd(copy.d, n, copy.stride, eigenvalue.d);
        if (ok == 0) {
            double smallest = (double)eigenvalue.d[0], largest = (double)eigenvalue.d[n - 1];
            double smallest_magnitude = fabs(smallest);
            int negative = 0;
            for (int i = 0; i < n; i++) {
                double a = fabs((double)eigenvalue.d[i]);
                if (a < smallest_magnitude) smallest_magnitude = a;
                if ((double)eigenvalue.d[i] < 0) negative++;
            }
            fprintf(report, "  Hessian at the %d iteration point: %d parameters, "
                            "eigenvalues %.4g to %.4g,\n", budget[N_BUDGETS - 1], n,
                    smallest, largest);
            fprintf(report, "    %d negative, condition number |max|/|min| %.4g\n",
                    negative, fabs(largest) / smallest_magnitude);
            fprintf(report, "    curvature per parameter, largest and smallest diagonal:\n");
            int stiff = 0, flat = 0;
            for (int i = 1; i < n; i++) {
                if (fabs((double)AT(H, i, i)) > fabs((double)AT(H, stiff, stiff))) stiff = i;
                if (fabs((double)AT(H, i, i)) < fabs((double)AT(H, flat, flat))) flat = i;
            }
            char stiff_name[64], flat_name[64];
            _qvarma_theta_name(&last_shape, stiff, stiff_name, (int)sizeof stiff_name);
            _qvarma_theta_name(&last_shape, flat, flat_name, (int)sizeof flat_name);
            fprintf(report, "    stiffest %-18s %.4g\n", stiff_name, (double)AT(H, stiff, stiff));
            fprintf(report, "    flattest %-18s %.4g\n", flat_name, (double)AT(H, flat, flat));
        } else {
            fprintf(report, "  Hessian eigenvalues did not converge (_syevd returned %d)\n", ok);
        }
        fprintf(report, "\n");
        mat_free(copy); mat_free(eigenvalue); mat_free(H);
        mat_free(last_theta);
        qvarma_params_free(&last_shape);
        mat_free(y);
    }
    fclose(report);
    return 0;
}
