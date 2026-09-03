/*
The t-QVARMA (no drift) counterpart to us_qvarmad_employment_change.c, same
data, same partition, same two specs that file's own grid settled on -
(p,q,r) = (1,1,2) and (1,1,4) - so the two can be compared directly: same
build_block, same build_start convention, same residual battery, with a
mean squared error check added (docs/MODEL_TEMPLATE.md entry 16 is why this
is its own file rather than qvarma_d.h with the drift term switched off:
qvarma.h and qvarma_d.h define the same names, so a translation unit uses
one or the other, never both).

Data construction, partition, residual checks: identical to
us_qvarmad_employment_change.c's own - see that file's own header comment.
Cached per spec: out/us_qvarma_employment_change_p<p>q<q>r<r>_fit.json.
Output: out/us_qvarma_employment_change.txt,
out/us_qvarma_employment_change_p<p>q<q>r<r>_report.txt,
out/us_qvarma_employment_change_p<p>q<q>r<r>_residuals.csv. In
EXPERIMENT_STEMS. Nothing printed.
*/

#include "us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./inference/unit_root.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <stdlib.h>

#define K 5
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define START_NU ((mreal)30)
#define MAX_ITERATIONS 8000
#define PERIODS (ESTIMATION_PERIODS - 1)

enum { ROW_GDP_GROWTH, ROW_EN_GROWTH, ROW_EMPLOYMENT_CHANGE, ROW_INFLATION, ROW_INTEREST_RATE };
static const char *row_name[K] = {
    "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
};

static Mat build_block(Mat original) {
    Mat y = mat_new(K, PERIODS);
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

static mreal first_difference_sd(Mat y, int row) {
    int periods = y.c;
    Mat difference = mat_new(1, periods - 1);
    for (int t = 1; t < periods; t++) AT(difference, 0, t - 1) = AT(y, row, t) - AT(y, row, t - 1);
    mreal sd = (mreal)sqrt((double)stats_var(difference));
    mat_free(difference);
    return sd;
}

typedef struct { int p, q, r; } Spec;
/* The two specs us_qvarmad_employment_change.c's own 9-point grid settled
   on, not a fresh grid search here - this file exists to compare against
   that choice, not to re-litigate it. */
static const Spec spec_grid[] = { { 1, 1, 2 }, { 1, 1, 4 } };
#define N_SPECS ((int)(sizeof spec_grid / sizeof spec_grid[0]))

static QvarmaParams build_start(Mat y, Spec spec) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, spec.p, spec.q, spec.r, R, SHARED_BETA, WARMUP_LONGEST);
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
    for (int i = 0; i < spec.p; i++) AT(m.Phi_star, i, 0) = (mreal)0.3 / (mreal)(i + 1);
    for (int j = 0; j < spec.q; j++)
        for (int a = 0; a < K_STAR; a++) AT(m.Psi_star[j], a, a) = (mreal)0.05;

    for (int a = 0; a < K; a++)
        for (int b = 0; b <= a; b++)
            AT(m.Omega_inv, a, b) = (mreal)(b == a ? first_difference_sd(y, a) : 0);
    m.nu = START_NU;

    for (int l = 0; l < spec.r; l++) {
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

typedef struct {
    Spec spec;
    QvarmaFitResult result;
    Mat residual;
    int is_maximum;
    mreal smallest_curvature;
} Fitted;

static void write_residual_csv(Mat residual, const char *path) {
    Mat transposed = mat_new(residual.c, residual.r);
    for (int t = 0; t < residual.c; t++)
        for (int a = 0; a < residual.r; a++) AT(transposed, t, a) = AT(residual, a, t);
    DataFrame df = df_from_matrix(transposed, row_name);
    df_write_csv(&df, path, csv_write_options_default());
    df_free(&df);
    mat_free(transposed);
}

static void run_spec(Mat y, Spec spec, Fitted *out) {
    char cache_path[128], report_path[128], residual_path[128];
    snprintf(cache_path, sizeof cache_path, "out/us_qvarma_employment_change_p%dq%dr%d_fit.json",
             spec.p, spec.q, spec.r);
    snprintf(report_path, sizeof report_path,
             "out/us_qvarma_employment_change_p%dq%dr%d_report.txt", spec.p, spec.q, spec.r);
    snprintf(residual_path, sizeof residual_path,
             "out/us_qvarma_employment_change_p%dq%dr%d_residuals.csv", spec.p, spec.q, spec.r);

    QvarmaParams start = build_start(y, spec);
    QvarmaFitOptions options = qvarma_default_fit_options();
    options.max_iterations = MAX_ITERATIONS;
    QvarmaFitResult result = qvarma_fit_cached(y, &start, options, cache_path, 0);
    qvarma_write_report(&result, y, report_path);

    QvarmaStandardErrors errors = qvarma_standard_errors(&result.params, y);
    out->is_maximum = errors.is_maximum;
    out->smallest_curvature = errors.smallest_curvature;
    qvarma_standard_errors_free(&errors);

    int periods = y.c;
    Vec theta = mat_new(qvarma_n_theta(&result.params), 1);
    _qvarma_unlink(&result.params, theta);
    Tape *tape = tape_new();
    Node *theta_node = ad_leaf(tape, theta);
    QvarmaLinked linked = _qvarma_link(tape, theta_node, &result.params);
    Node **v_nodes = (Node**)malloc((size_t)periods * sizeof(Node*));
    _qvarma_filter(tape, &linked, &result.params, y, NULL, NULL, v_nodes);

    Mat residual = mat_new(K, periods);
    for (int t = 0; t < periods; t++)
        for (int a = 0; a < K; a++) AT(residual, a, t) = v_nodes[t]->val.d[a];
    write_residual_csv(residual, residual_path);

    free(v_nodes);
    qvarma_linked_free(&linked);
    tape_free(tape);
    mat_free(theta);
    qvarma_params_free(&start);

    out->spec = spec;
    out->result = result;
    out->residual = residual;
}

static void write_residual_mean_check(FILE *out, const Fitted *f) {
    int periods = f->residual.c;
    int bandwidth = kpss_bandwidth(periods);
    fprintf(out, "  residual mean, null is zero, Newey-West se at bandwidth %d\n", bandwidth);
    for (int a = 0; a < K; a++) {
        mreal mean = 0;
        for (int t = 0; t < periods; t++) mean += AT(f->residual, a, t);
        mean /= (mreal)periods;
        mreal *deviation = (mreal*)malloc((size_t)periods * sizeof(mreal));
        for (int t = 0; t < periods; t++) deviation[t] = AT(f->residual, a, t) - mean;
        mreal long_run_variance = _bartlett_long_run_variance(deviation, periods, bandwidth);
        mreal se = (mreal)sqrt((double)long_run_variance / (double)periods);
        mreal t_statistic = mean / se;
        fprintf(out, "    %-14s mean %8.4f   se %7.4f   t %7.3f   %s\n", row_name[a],
                (double)mean, (double)se, (double)t_statistic,
                MABS(t_statistic) > (mreal)1.96 ? "mean is nonzero at 5%" : "consistent with zero");
        free(deviation);
    }
}

static void write_covariance_check(FILE *out, const Fitted *f) {
    int periods = f->residual.c;
    Mat transposed = mat_new(periods, K);
    for (int t = 0; t < periods; t++)
        for (int a = 0; a < K; a++) AT(transposed, t, a) = AT(f->residual, a, t);
    Mat covariance = stats_autocov(transposed, 0);
    mreal inflation = f->result.params.nu / (f->result.params.nu - 2);
    mreal worst = 0;
    for (int a = 0; a < K; a++)
        for (int b = 0; b < K; b++) {
            mreal expected = inflation * AT(f->result.params.Sigma, a, b);
            mreal difference = MABS(AT(covariance, a, b) - expected);
            if (difference > worst) worst = difference;
        }
    fprintf(out, "  contemporaneous covariance against nu/(nu-2) Sigma: worst |difference| %.4g\n",
            (double)worst);
    mat_free(transposed);
    mat_free(covariance);
}

static void write_autocorrelation_check(FILE *out, const Fitted *f) {
    int periods = f->residual.c;
    mreal band = (mreal)(1.96 / sqrt((double)periods));
    fprintf(out, "  residual autocorrelation, lags 1 to 4, against the approximate 95 per "
                 "cent white-noise band +-%.4f, and the Ljung-Box joint test over those\n"
                 "  same 4 lags, null is no autocorrelation up to lag 4\n", (double)band);
    for (int a = 0; a < K; a++) {
        Mat row = mat_slice(f->residual, a, a + 1, 0, periods);
        fprintf(out, "    %-14s", row_name[a]);
        int flagged = 0;
        for (int lag = 1; lag <= 4; lag++) {
            mreal rho = stats_autocorr(row, lag);
            fprintf(out, " %7.4f%s", (double)rho, MABS(rho) > band ? "*" : " ");
            if (MABS(rho) > band) flagged++;
        }
        StatsLjungBox joint = stats_ljung_box(row, 4);
        fprintf(out, "   %d of 4 outside the band, Ljung-Box Q %7.3f, p %.4g\n", flagged,
                joint.statistic, joint.p_value);
    }
}

static void write_quadratic_form_check(FILE *out, const Fitted *f) {
    int periods = f->residual.c;
    const QvarmaParams *m = &f->result.params;
    mreal *q_over_k = (mreal*)malloc((size_t)periods * sizeof(mreal));
    for (int t = 0; t < periods; t++) {
        Vec residual = mat_new(K, 1);
        for (int a = 0; a < K; a++) residual.d[a] = AT(f->residual, a, t);
        Vec weighted = vec_chol_solve(m->Omega_inv, residual);
        mreal quadratic = 0;
        for (int a = 0; a < K; a++) quadratic += residual.d[a] * weighted.d[a];
        q_over_k[t] = quadratic / (mreal)K;
        mat_free(residual); mat_free(weighted);
    }
    double mean = 0;
    for (int t = 0; t < periods; t++) mean += (double)q_over_k[t];
    mean /= periods;
    double variance = 0;
    for (int t = 0; t < periods; t++) variance += (mean - (double)q_over_k[t]) * (mean - (double)q_over_k[t]);
    variance /= periods;
    double nu = (double)m->nu;
    double expected_mean = nu / (nu - 2);
    double expected_variance = nu > 4
        ? 2 * nu * nu * ((double)K + nu - 2) / ((double)K * (nu - 4) * (nu - 2) * (nu - 2))
        : (double)NAN;
    fprintf(out, "  q_t/K against F(%d, %.2f): sample mean %.4f against %.4f\n", K, nu, mean,
            expected_mean);
    if (nu > 4)
        fprintf(out, "    sample variance %.4f against %.4f\n", variance, expected_variance);
    else
        fprintf(out, "    sample variance %.4f against undefined at this nu (nu <= 4)\n",
                variance);
    free(q_over_k);
}

/* Mean squared residual, et_al.'s own stats_mse against an all-zero target -
   the residual already is actual-minus-fitted, so this is exactly the
   in-sample prediction MSE, pooled over every series and period and also
   broken out per series, in y's own units (not comparable across a series
   whose own scale differs, which is why both are reported). */
static void write_mse_check(FILE *out, const Fitted *f) {
    int periods = f->residual.c;
    Mat zero = mat_new(K, periods);
    fprintf(out, "  mean squared residual (et_al.'s stats_mse against zero)\n");
    fprintf(out, "    pooled, all %d series x %d periods: %.6f\n", K, periods,
            (double)stats_mse(f->residual, zero));
    for (int a = 0; a < K; a++) {
        Mat row = mat_slice(f->residual, a, a + 1, 0, periods);
        Mat zero_row = mat_slice(zero, a, a + 1, 0, periods);
        fprintf(out, "    %-14s %.6f\n", row_name[a], (double)stats_mse(row, zero_row));
    }
    mat_free(zero);
}

int main(void) {
    Mat original = load_us_system();
    Mat y = build_block(original);
    mat_free(original);

    Fitted fitted[N_SPECS];
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N_SPECS; i++) run_spec(y, spec_grid[i], &fitted[i]);

    FILE *out = fopen("out/us_qvarma_employment_change.txt", "w");
    assert(out && "cannot open the output path for writing");
    fprintf(out, "t-QVARMA (no drift), same Fisher-relation partition with Employment "
                 "differenced as us_qvarmad_employment_change.c: GDP_growth, EN_growth, "
                 "Employment_change I(0), Inflation and InterestRate co-integrated "
                 "(K_star 3, K_dagger 2, R 1 - forced) - built to compare against that "
                 "file's own drift-carrying fit, not to re-search the grid\n");
    fprintf(out, "(GDP_growth, EN_growth, Employment_change, Inflation, InterestRate), "
                 "1973Q2 to 2019Q4, %d quarters\n\n", PERIODS);

    fprintf(out, "%-10s %6s %10s %10s %10s %14s %10s %8s %6s %10s\n", "spec", "params",
            "aic", "bic", "hqc", "log_lik", "gradient", "niter", "C_1", "genuine_max");
    for (int i = 0; i < N_SPECS; i++) {
        const QvarmaFitResult *r = &fitted[i].result;
        char label[16];
        snprintf(label, sizeof label, "(%d,%d,%d)", fitted[i].spec.p, fitted[i].spec.q,
                 fitted[i].spec.r);
        mreal c1 = qvarma_max_eigenvalue_modulus(&r->params);
        fprintf(out, "%-10s %6d %10.4f %10.4f %10.4f %14.4f %10.4g %8d %6.3f%s %10s\n", label,
                qvarma_n_theta(&r->params), (double)r->aic, (double)r->bic, (double)r->hannan_quinn,
                (double)r->log_likelihood, (double)r->gradient_norm, r->niter, (double)c1,
                c1 >= 1 ? "  not stationary" : "", fitted[i].is_maximum ? "yes" : "no");
    }

    fprintf(out, "\nresidual analysis, every spec\n");
    for (int i = 0; i < N_SPECS; i++) {
        fprintf(out, "\n(%d,%d,%d)%s, curvature %.4g\n", fitted[i].spec.p, fitted[i].spec.q,
                fitted[i].spec.r,
                fitted[i].is_maximum ? "" : " - NOT A GENUINE MAXIMUM, read with that in mind",
                (double)fitted[i].smallest_curvature);
        write_residual_mean_check(out, &fitted[i]);
        write_covariance_check(out, &fitted[i]);
        write_autocorrelation_check(out, &fitted[i]);
        write_quadratic_form_check(out, &fitted[i]);
        write_mse_check(out, &fitted[i]);
    }

    fclose(out);
    for (int i = 0; i < N_SPECS; i++) {
        qvarma_fit_result_free(&fitted[i].result);
        mat_free(fitted[i].residual);
    }
    mat_free(y);
    return 0;
}
