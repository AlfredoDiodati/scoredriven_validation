/*
How the log-likelihood of each fitted replicate of the winning configuration
moves when nu alone is changed, with every other parameter left where the fit
put it. Nothing is fitted.

A Student-t likelihood in a near-Gaussian series should change little once nu
is large. On some of the winner's fits it does not: moving nu by one per cent
costs thousands of log-likelihood units, and the loss grows roughly in
proportion to nu. That is what the filter does when some periods' quadratic
form q_t = v_t' Sigma^-1 v_t becomes very large, since each period contributes
-((nu+K)/2) log(1 + q_t/nu). The score u_t = nu v_t / (nu + q_t) feeds the
location recursion, so nu changes the path v_t itself and not only the density.

For every replicate of the winner the script evaluates the log-likelihood at the
fitted parameters with nu multiplied by each factor in nu_factor, and with nu
set to each value in nu_fixed. It also records the largest q_t/K over the
periods at the fitted nu and at nu one per cent higher. A fit is called
    smooth       when both one per cent moves change the log-likelihood by less
                 than SMOOTH_LIMIT
    knife-edge   when either changes it by more than KNIFE_EDGE_LIMIT
    between      otherwise
Each fit's loss is read from out/abm_system_mse_qvarma_joint.csv.

Requires out/abm_system_mcs_joint.csv, out/abm_system_mse_qvarma_joint.csv, the
fit cache and dataset/abm_system/. None are Makefile prerequisites, because
rebuilding them reruns a million fits.

Output, none of it printed:
    out/abm_system_winner_nu_likelihood_scan.txt   the summary
    out/abm_system_winner_nu_likelihood_scan.csv   one row per replicate
*/

#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define P 1
#define Q 1
#define RLAG 2
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define SPEC_LABEL "p1q1r2"

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs_joint.csv"
#define LOSS_TABLE_PATH "out/abm_system_mse_qvarma_joint.csv"
#define REPORT_PATH "out/abm_system_winner_nu_likelihood_scan.txt"
#define TABLE_PATH "out/abm_system_winner_nu_likelihood_scan.csv"

static const double nu_factor[] = { 0.5, 0.9, 0.99, 1.01, 1.1, 2 };
#define N_FACTORS ((int)(sizeof nu_factor / sizeof nu_factor[0]))
#define FACTOR_DOWN_ONE_PERCENT 2
#define FACTOR_UP_ONE_PERCENT 3

static const double nu_fixed[] = { 7.198920587, 100, 1e6 };
#define N_FIXED ((int)(sizeof nu_fixed / sizeof nu_fixed[0]))

/* Log-likelihood units. */
#define SMOOTH_LIMIT 1.0
#define KNIFE_EDGE_LIMIT 100.0

enum { SMOOTH, BETWEEN, KNIFE_EDGE, N_CLASSES };
static const char *class_name[N_CLASSES] = { "smooth", "between", "knife-edge" };

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    return m;
}

static void find_winner(char *sample, size_t size) {
    DataFrame df = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    char **name = df_col_string(&df, "model");
    Mat mean_loss = df_col_numeric(&df, "mean_loss");
    Mat in_set = df_col_numeric(&df, "in_set");
    int best = -1;
    for (int i = 0; i < df.r; i++)
        if (AT(in_set, i, 0) > 0 && (best < 0 || AT(mean_loss, i, 0) < AT(mean_loss, best, 0))) best = i;
    assert(best >= 0 && "abm_system_winner_nu_likelihood_scan: no model is in the confidence set");
    const char *marker = strstr(name[best], "_qvarma_" SPEC_LABEL);
    assert(marker && "abm_system_winner_nu_likelihood_scan: the winning model is not the p1q1r2 spec");
    size_t length = (size_t)(marker - name[best]);
    assert(length < size);
    memcpy(sample, name[best], length);
    sample[length] = '\0';
    df_free(&df);
}

/* Log-likelihood at theta with the nu coordinate set to nu, and the largest
   q_t/K along the path that evaluation produced. */
static double likelihood_at_nu(QvarmaAnalytic *workspace, QvarmaParams *working, Vec theta, int nu_offset,
                               double nu, Mat y, double *largest_q_over_k) {
    mreal kept = theta.d[nu_offset];
    theta.d[nu_offset] = (mreal)log(nu - 2.0);
    Vec no_gradient = { 0, 0, 0, NULL };
    double value = (double)qvarma_analytic_log_likelihood(workspace, theta, y, no_gradient);

    if (largest_q_over_k) {
        qvarma_params_from_theta(theta, working);
        Vec residual = mat_new(K, 1);
        double largest = 0;
        for (int t = 0; t < y.c; t++) {
            const mreal *v = qvarma_analytic_v(workspace, t);
            for (int a = 0; a < K; a++) residual.d[a] = v[a];
            Vec whitened = vec_triangular_solve(working->Omega_inv, residual, 'L', 'N', 'N');
            double q = 0;
            for (int a = 0; a < K; a++) q += (double)whitened.d[a] * (double)whitened.d[a];
            mat_free(whitened);
            if (q / K > largest) largest = q / K;
        }
        mat_free(residual);
        *largest_q_over_k = largest;
    }
    theta.d[nu_offset] = kept;
    return value;
}

typedef struct {
    int replicate;
    double nu, log_likelihood, loss;
    double factor_change[N_FACTORS];
    double fixed_change[N_FIXED];
    double largest_q_over_k, largest_q_over_k_up;
    int class_index;
} ScanRow;

static double quantile_of(const double *values, int n, double probability) {
    Mat copy = mat_new(n, 1);
    for (int i = 0; i < n; i++) copy.d[i] = (mreal)values[i];
    double q = (double)stats_quantile(copy, (mreal)probability);
    mat_free(copy);
    return q;
}

int main(void) {
    char sample[64];
    find_winner(sample, sizeof sample);
    char sample_dir[600];
    snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, sample);

    DataFrame losses = df_read_csv(LOSS_TABLE_PATH, csv_read_options_default());
    char column_name[128];
    snprintf(column_name, sizeof column_name, "%s_qvarma_%s", sample, SPEC_LABEL);
    Mat loss_column = df_col_numeric(&losses, column_name);
    Mat replicate_column = df_col_numeric(&losses, "replicate");
    int n_replicates = losses.r;

    QvarmaParams shape = spec_shape();
    int n_theta = qvarma_n_theta(&shape);
    int nu_offset, nu_count;
    qvarma_block_range(&shape, QVARMA_BLOCK_NU, &nu_offset, &nu_count);

    ScanRow *rows = (ScanRow*)calloc((size_t)n_replicates, sizeof(ScanRow));

    #pragma omp parallel
    {
        QvarmaParams working = spec_shape();
        Vec theta = mat_new(n_theta, 1);
        QvarmaAnalytic *workspace = NULL;
        int workspace_periods = -1;

        #pragma omp for schedule(dynamic)
        for (int row = 0; row < n_replicates; row++) {
            ScanRow *s = &rows[row];
            s->replicate = (int)AT(replicate_column, row, 0);
            s->loss = (double)AT(loss_column, row, 0);

            Mat y = abm_system_read_replicate(sample_dir, s->replicate);
            if (y.c != workspace_periods) {
                if (workspace) qvarma_analytic_free(workspace);
                workspace = qvarma_analytic_new(&shape, y.c);
                workspace_periods = y.c;
            }
            char path[700];
            snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, sample, s->replicate,
                     SPEC_LABEL);
            int loaded = qvarma_load_params(&working, path);
            assert(loaded && "abm_system_winner_nu_likelihood_scan: a cached fit of the winner is missing");
            _qvarma_unlink(&working, theta);
            s->nu = (double)working.nu;

            s->log_likelihood = likelihood_at_nu(workspace, &working, theta, nu_offset, s->nu, y,
                                                 &s->largest_q_over_k);
            for (int f = 0; f < N_FACTORS; f++)
                s->factor_change[f] = likelihood_at_nu(workspace, &working, theta, nu_offset,
                                                       2.0 + (s->nu - 2.0) * nu_factor[f], y,
                                                       f == FACTOR_UP_ONE_PERCENT ? &s->largest_q_over_k_up : NULL)
                                      - s->log_likelihood;
            for (int f = 0; f < N_FIXED; f++)
                s->fixed_change[f] = likelihood_at_nu(workspace, &working, theta, nu_offset, nu_fixed[f], y, NULL)
                                     - s->log_likelihood;
            qvarma_params_from_theta(theta, &working);

            double worst_one_percent = fabs(s->factor_change[FACTOR_DOWN_ONE_PERCENT]);
            if (fabs(s->factor_change[FACTOR_UP_ONE_PERCENT]) > worst_one_percent)
                worst_one_percent = fabs(s->factor_change[FACTOR_UP_ONE_PERCENT]);
            if (!isfinite(worst_one_percent) || worst_one_percent > KNIFE_EDGE_LIMIT) s->class_index = KNIFE_EDGE;
            else if (worst_one_percent < SMOOTH_LIMIT) s->class_index = SMOOTH;
            else s->class_index = BETWEEN;
            mat_free(y);
        }

        if (workspace) qvarma_analytic_free(workspace);
        mat_free(theta);
        qvarma_params_free(&working);
    }

    FILE *table = fopen(TABLE_PATH, "w");
    assert(table && "abm_system_winner_nu_likelihood_scan: cannot open the table for writing");
    fprintf(table, "configuration,replicate,nu,log_likelihood,loss,class,largest_q_over_k,largest_q_over_k_nu_up_1pct");
    for (int f = 0; f < N_FACTORS; f++) fprintf(table, ",change_nu_times_%g", nu_factor[f]);
    for (int f = 0; f < N_FIXED; f++) fprintf(table, ",change_nu_at_%g", nu_fixed[f]);
    fprintf(table, "\n");
    for (int row = 0; row < n_replicates; row++) {
        const ScanRow *s = &rows[row];
        fprintf(table, "%s,%d,%.10g,%.10g,%.10g,%s,%.6g,%.6g", sample, s->replicate, s->nu, s->log_likelihood,
                s->loss, class_name[s->class_index], s->largest_q_over_k, s->largest_q_over_k_up);
        for (int f = 0; f < N_FACTORS; f++) fprintf(table, ",%.6g", s->factor_change[f]);
        for (int f = 0; f < N_FIXED; f++) fprintf(table, ",%.6g", s->fixed_change[f]);
        fprintf(table, "\n");
    }
    fclose(table);

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_winner_nu_likelihood_scan: cannot open the report path for writing");
    fprintf(out, "%s: %d fitted replicates, log-likelihood with nu changed and every other parameter at the fit\n",
            sample, n_replicates);
    fprintf(out, "nu times a factor means nu' = 2 + (nu - 2) * factor; change = log-likelihood at nu' minus at the "
                 "fitted nu\n");
    fprintf(out, "smooth: both 1%% moves change it by less than %g; knife-edge: either by more than %g\n\n",
            SMOOTH_LIMIT, KNIFE_EDGE_LIMIT);

    double *values = (double*)malloc((size_t)n_replicates * sizeof(double));
    fprintf(out, "%-12s %6s %10s %10s %10s %10s %14s %14s\n", "class", "fits", "nu p10", "nu p50", "nu p90",
            "mean loss", "q/K max p50", "q/K max +1% p50");
    for (int c = 0; c < N_CLASSES; c++) {
        int members = 0;
        double loss_sum = 0;
        for (int row = 0; row < n_replicates; row++) members += rows[row].class_index == c;
        if (members == 0) { fprintf(out, "%-12s %6d\n", class_name[c], 0); continue; }
        double *nu = (double*)malloc((size_t)members * sizeof(double));
        double *q_own = (double*)malloc((size_t)members * sizeof(double));
        double *q_up = (double*)malloc((size_t)members * sizeof(double));
        int at = 0;
        for (int row = 0; row < n_replicates; row++) {
            if (rows[row].class_index != c) continue;
            nu[at] = rows[row].nu;
            q_own[at] = rows[row].largest_q_over_k;
            q_up[at] = isfinite(rows[row].largest_q_over_k_up) ? rows[row].largest_q_over_k_up : 1e300;
            loss_sum += rows[row].loss;
            at++;
        }
        fprintf(out, "%-12s %6d %10.4g %10.4g %10.4g %10.5f %14.4g %14.4g\n", class_name[c], members,
                quantile_of(nu, at, 0.1), quantile_of(nu, at, 0.5), quantile_of(nu, at, 0.9), loss_sum / members,
                quantile_of(q_own, at, 0.5), quantile_of(q_up, at, 0.5));
        free(nu); free(q_own); free(q_up);
    }

    fprintf(out, "\nlog-likelihood change, quantiles over all %d fits\n", n_replicates);
    fprintf(out, "%-18s %12s %12s %12s %12s %12s\n", "nu moved to", "p10", "p25", "p50", "p75", "p90");
    for (int f = 0; f < N_FACTORS; f++) {
        int n = 0;
        for (int row = 0; row < n_replicates; row++)
            values[n++] = isfinite(rows[row].factor_change[f]) ? rows[row].factor_change[f] : -1e300;
        char label[32];
        snprintf(label, sizeof label, "nu x %g", nu_factor[f]);
        fprintf(out, "%-18s %12.4g %12.4g %12.4g %12.4g %12.4g\n", label, quantile_of(values, n, 0.1),
                quantile_of(values, n, 0.25), quantile_of(values, n, 0.5), quantile_of(values, n, 0.75),
                quantile_of(values, n, 0.9));
    }
    for (int f = 0; f < N_FIXED; f++) {
        int n = 0;
        for (int row = 0; row < n_replicates; row++)
            values[n++] = isfinite(rows[row].fixed_change[f]) ? rows[row].fixed_change[f] : -1e300;
        char label[32];
        snprintf(label, sizeof label, "nu = %g", nu_fixed[f]);
        fprintf(out, "%-18s %12.4g %12.4g %12.4g %12.4g %12.4g\n", label, quantile_of(values, n, 0.1),
                quantile_of(values, n, 0.25), quantile_of(values, n, 0.5), quantile_of(values, n, 0.75),
                quantile_of(values, n, 0.9));
    }
    fclose(out);

    free(values);
    free(rows);
    df_free(&losses);
    qvarma_params_free(&shape);
    return 0;
}
