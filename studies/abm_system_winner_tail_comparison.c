/*
Why the US data give nu near 7 and the fits of the winning configuration give
nu in the thousands: whether the US residuals have tails the simulated ones do
not, and which quarters produce them.

The measure is the quadratic form of each period's residual, q_t = v_t'
Sigma^-1 v_t, divided by K. Under Gaussian innovations it is chi-squared with K
degrees of freedom over K, so its mean is 1, its variance is 2/K and the share
of periods above a threshold c is the chi-squared survival at K c. Under the
fitted Student-t it is F(K, nu), with variance growing without bound as nu
falls towards 4. How much more often q_t/K lands far out than the Gaussian
reference allows is what the likelihood reads as a small nu.

Four sets of residuals, each at the parameters of its own fit:
    US, nu free          out/us_qvarma_spec_choice_p1q1r2_fit.json
    US, nu held          out/us_qvarma_nu_sensitivity_nu<nu>_fit.json, the US
                         data refitted with nu held at the winner's averaged nu,
                         so the US tails are also measured under a near-Gaussian
                         fit
    winner, nu free      every cached replicate fit of the winning configuration
    winner by nu group   the same fits split by their own nu
The largest US values of q_t/K are listed with their quarter and the whitened
residual z_t = Omega_inv^-1 v_t series by series, to show which variable
carries them.

Requires the two US fit caches above, out/us_system.csv,
out/abm_system_mcs.csv, the fit cache and dataset/abm_system/.

Output, none of it printed:
    out/abm_system_winner_tail_comparison.txt
*/

#include "applications/abm_system.h"
#include "applications/us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./special.h>
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
#define BENCHMARK_FIT_PATH "out/us_qvarma_spec_choice_p1q1r2_fit.json"
#define WINNER_PARAMS_PATH "out/abm_system_winner_irf_theta.json"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs.csv"
#define REPORT_PATH "out/abm_system_winner_tail_comparison.txt"

/* The US block starts one quarter after us_system.csv, lost to differencing. */
#define US_FIRST_YEAR 1973
#define US_FIRST_QUARTER 2
#define N_LARGEST_US 12

static const double threshold[] = { 2, 3, 5, 8 };
#define N_THRESHOLDS ((int)(sizeof threshold / sizeof threshold[0]))

static const char *series_name[K] = { "GDP", "energy", "employment", "inflation", "interest" };

#define N_NU_GROUPS 3
static const double nu_group_upper[N_NU_GROUPS] = { 1e3, 1e4, INFINITY };
static const char *nu_group_name[N_NU_GROUPS] = { "nu below 1e3", "nu 1e3 to 1e4", "nu above 1e4" };

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    return m;
}

static Mat build_us_block(void) {
    Mat original = load_us_system();
    int periods = ESTIMATION_PERIODS - 1;
    Mat y = mat_new(K, periods);
    for (int t = 1; t < ESTIMATION_PERIODS; t++) {
        int c = t - 1;
        AT(y, ROW_GDP_GROWTH, c) = AT(original, LOG_GDP, t) - AT(original, LOG_GDP, t - 1);
        AT(y, ROW_EN_GROWTH, c) = AT(original, LOG_ENERGY_DEMAND, t)
                                 - AT(original, LOG_ENERGY_DEMAND, t - 1);
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = AT(original, EMPLOYMENT, t) - AT(original, EMPLOYMENT, t - 1);
        AT(y, ROW_INFLATION, c) = AT(original, LOG_CPI, t) - AT(original, LOG_CPI, t - 1);
        AT(y, ROW_INTEREST_RATE, c) = AT(original, INTEREST_RATE, t);
    }
    mat_free(original);
    return y;
}

static void find_winner(char *sample, size_t size) {
    DataFrame df = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    char **name = df_col_string(&df, "model");
    Mat mean_loss = df_col_numeric(&df, "mean_loss");
    Mat in_set = df_col_numeric(&df, "in_set");
    int best = -1;
    for (int i = 0; i < df.r; i++)
        if (AT(in_set, i, 0) > 0 && (best < 0 || AT(mean_loss, i, 0) < AT(mean_loss, best, 0))) best = i;
    assert(best >= 0 && "abm_system_winner_tail_comparison: no model is in the confidence set");
    const char *marker = strstr(name[best], "_qvarma_" SPEC_LABEL);
    assert(marker && "abm_system_winner_tail_comparison: the winning model is not the p1q1r2 spec");
    size_t length = (size_t)(marker - name[best]);
    assert(length < size);
    memcpy(sample, name[best], length);
    sample[length] = '\0';
    df_free(&df);
}

/* q_t/K for every period, and the whitened residual z_t, at m on y. The caller
   frees both. */
static void residual_path(const QvarmaParams *m, Mat y, Mat *q_over_k, Mat *z) {
    QvarmaAnalytic *workspace = qvarma_analytic_new(m, y.c);
    Vec theta = mat_new(qvarma_n_theta(m), 1);
    _qvarma_unlink(m, theta);
    Vec no_gradient = { 0, 0, 0, NULL };
    qvarma_analytic_log_likelihood(workspace, theta, y, no_gradient);
    *q_over_k = mat_new(y.c, 1);
    *z = mat_new(K, y.c);
    for (int t = 0; t < y.c; t++) {
        const mreal *whitened = workspace->path + (size_t)t * workspace->path_stride
                                + (size_t)QVARMA_ANALYTIC_Z * K;
        mreal q = 0;
        for (int a = 0; a < K; a++) {
            AT(*z, a, t) = whitened[a];
            q += whitened[a] * whitened[a];
        }
        q_over_k->d[t] = q / (mreal)K;
    }
    mat_free(theta);
    qvarma_analytic_free(workspace);
}

typedef struct {
    double mean, variance, largest;
    double share_above[N_THRESHOLDS];
} TailSummary;

static TailSummary summarise(Mat q_over_k) {
    TailSummary s;
    s.mean = (double)stats_mean(q_over_k);
    s.variance = (double)stats_var(q_over_k);
    s.largest = 0;
    for (int i = 0; i < N_THRESHOLDS; i++) s.share_above[i] = 0;
    for (int t = 0; t < q_over_k.r; t++) {
        double value = (double)q_over_k.d[t];
        if (value > s.largest) s.largest = value;
        for (int i = 0; i < N_THRESHOLDS; i++) s.share_above[i] += value > threshold[i];
    }
    for (int i = 0; i < N_THRESHOLDS; i++) s.share_above[i] /= q_over_k.r;
    return s;
}

static void write_header(FILE *out) {
    fprintf(out, "  %-28s %6s %8s %9s %9s", "residuals", "fits", "mean", "variance", "largest");
    for (int i = 0; i < N_THRESHOLDS; i++) fprintf(out, "    >%-4g", threshold[i]);
    fprintf(out, "\n");
}

static void write_row(FILE *out, const char *label, int fits, TailSummary s) {
    fprintf(out, "  %-28s %6d %8.3f %9.3f %9.2f", label, fits, s.mean, s.variance, s.largest);
    for (int i = 0; i < N_THRESHOLDS; i++) fprintf(out, " %8.4f", s.share_above[i]);
    fprintf(out, "\n");
}

/* The median of each summary over a group of fits, as one row. */
static TailSummary median_summary(const TailSummary *summaries, const int *member, int n) {
    int count = 0;
    for (int i = 0; i < n; i++) count += member[i];
    Mat values = mat_new(count, 1);
    TailSummary out;
    #define MEDIAN_OF(field) do { int at = 0; for (int i = 0; i < n; i++) if (member[i]) values.d[at++] = (mreal)summaries[i].field; out.field = (double)stats_median(values); } while (0)
    MEDIAN_OF(mean);
    MEDIAN_OF(variance);
    MEDIAN_OF(largest);
    for (int k = 0; k < N_THRESHOLDS; k++) MEDIAN_OF(share_above[k]);
    #undef MEDIAN_OF
    mat_free(values);
    return out;
}

int main(void) {
    char sample[64];
    find_winner(sample, sizeof sample);

    QvarmaParams winner_average = spec_shape();
    int loaded = qvarma_load_params(&winner_average, WINNER_PARAMS_PATH);
    assert(loaded && "abm_system_winner_tail_comparison: the winner's averaged parameters are missing");
    char held_path[160];
    snprintf(held_path, sizeof held_path, "out/us_qvarma_nu_sensitivity_nu%.4f_fit.json", (double)winner_average.nu);

    Mat us_y = build_us_block();
    QvarmaParams shape = spec_shape();
    QvarmaFitResult us_free = qvarma_fit_result_new(&shape);
    QvarmaFitResult us_held = qvarma_fit_result_new(&shape);
    loaded = qvarma_load_fit(&us_free, us_y, BENCHMARK_FIT_PATH);
    assert(loaded && "abm_system_winner_tail_comparison: the US benchmark fit is missing or stale");
    loaded = qvarma_load_fit(&us_held, us_y, held_path);
    assert(loaded && "abm_system_winner_tail_comparison: the US fit with nu held is missing; run make study-us_qvarma_nu_sensitivity");

    Mat us_free_q, us_free_z, us_held_q, us_held_z;
    residual_path(&us_free.params, us_y, &us_free_q, &us_free_z);
    residual_path(&us_held.params, us_y, &us_held_q, &us_held_z);

    char sample_dir[600];
    snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, sample);
    int n_replicates = 0;
    int *replicate = abm_system_list_replicates(sample_dir, &n_replicates);
    TailSummary *winner = (TailSummary*)calloc((size_t)n_replicates, sizeof(TailSummary));
    double *winner_nu = (double*)calloc((size_t)n_replicates, sizeof(double));

    #pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < n_replicates; row++) {
        Mat y = abm_system_read_replicate(sample_dir, replicate[row]);
        QvarmaParams m = spec_shape();
        char path[700];
        snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, sample, replicate[row], SPEC_LABEL);
        int ok = qvarma_load_params(&m, path);
        assert(ok && "abm_system_winner_tail_comparison: a cached fit of the winner is missing");
        Mat q, z;
        residual_path(&m, y, &q, &z);
        winner[row] = summarise(q);
        winner_nu[row] = (double)m.nu;
        mat_free(q); mat_free(z); mat_free(y);
        qvarma_params_free(&m);
    }

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_winner_tail_comparison: cannot open the report path for writing");
    fprintf(out, "q_t/K, the per-period quadratic form of the residual over K = %d, at each fit's own parameters\n", K);
    fprintf(out, "US: %d quarters, 1973Q2 to 2019Q4. %s: %d replicates of %d periods, one row per group is the median over its fits\n\n",
            us_y.c, sample, n_replicates, 400);

    TailSummary gaussian;
    gaussian.mean = 1;
    gaussian.variance = 2.0 / K;
    gaussian.largest = NAN;
    for (int i = 0; i < N_THRESHOLDS; i++) gaussian.share_above[i] = special_chi_squared_sf(K * threshold[i], K);
    double us_nu = (double)us_free.params.nu;
    TailSummary student = gaussian;
    student.mean = us_nu / (us_nu - 2);
    student.variance = 2 * us_nu * us_nu * (K + us_nu - 2) / (K * (us_nu - 4) * (us_nu - 2) * (us_nu - 2));
    for (int i = 0; i < N_THRESHOLDS; i++) student.share_above[i] = NAN;

    write_header(out);
    fprintf(out, "  reference, Gaussian: chi-squared(%d)/%d                    mean 1, variance %.3f; shares:", K, K, gaussian.variance);
    for (int i = 0; i < N_THRESHOLDS; i++) fprintf(out, " %.4f", gaussian.share_above[i]);
    fprintf(out, "\n");
    fprintf(out, "  reference, F(%d, %.2f), the US fit's nu: mean %.3f, variance %.3f\n", K, us_nu, student.mean, student.variance);

    char label[96];
    snprintf(label, sizeof label, "US, nu free (%.2f)", us_nu);
    write_row(out, label, 1, summarise(us_free_q));
    snprintf(label, sizeof label, "US, nu held at %.0f", (double)us_held.params.nu);
    write_row(out, label, 1, summarise(us_held_q));

    int *member = (int*)malloc((size_t)n_replicates * sizeof(int));
    for (int row = 0; row < n_replicates; row++) member[row] = 1;
    snprintf(label, sizeof label, "%s, all fits", sample);
    write_row(out, label, n_replicates, median_summary(winner, member, n_replicates));
    for (int g = 0; g < N_NU_GROUPS; g++) {
        double lower = g == 0 ? 0 : nu_group_upper[g - 1];
        int count = 0;
        for (int row = 0; row < n_replicates; row++) {
            member[row] = winner_nu[row] > lower && winner_nu[row] <= nu_group_upper[g];
            count += member[row];
        }
        if (count == 0) continue;
        snprintf(label, sizeof label, "%s, %s", sample, nu_group_name[g]);
        write_row(out, label, count, median_summary(winner, member, n_replicates));
    }

    /* How many of the winner's fits reach the US value on each measure. */
    TailSummary us_summary = summarise(us_held_q);
    int variance_reached = 0, largest_reached = 0, share5_reached = 0;
    for (int row = 0; row < n_replicates; row++) {
        variance_reached += winner[row].variance >= us_summary.variance;
        largest_reached += winner[row].largest >= us_summary.largest;
        share5_reached += winner[row].share_above[2] >= us_summary.share_above[2];
    }
    fprintf(out, "\nfits of %s at or above the US value with nu held near Gaussian: variance %d, largest %d, share above 5 %d, of %d\n",
            sample, variance_reached, largest_reached, share5_reached, n_replicates);

    fprintf(out, "\nthe %d largest US values of q_t/K at the nu-free fit, with the whitened residual z_t by series\n", N_LARGEST_US);
    fprintf(out, "  %-8s %8s", "quarter", "q_t/K");
    for (int a = 0; a < K; a++) fprintf(out, " %11s", series_name[a]);
    fprintf(out, "\n");
    int *used = (int*)calloc((size_t)us_y.c, sizeof(int));
    for (int k = 0; k < N_LARGEST_US; k++) {
        int best = -1;
        for (int t = 0; t < us_y.c; t++)
            if (!used[t] && (best < 0 || us_free_q.d[t] > us_free_q.d[best])) best = t;
        used[best] = 1;
        int total = (US_FIRST_QUARTER - 1) + best;
        fprintf(out, "  %4dQ%d   %8.2f", US_FIRST_YEAR + total / 4, 1 + total % 4, (double)us_free_q.d[best]);
        for (int a = 0; a < K; a++) fprintf(out, " %11.2f", (double)AT(us_free_z, a, best));
        fprintf(out, "\n");
    }

    double squared_by_series[K] = { 0 }, squared_total = 0;
    for (int t = 0; t < us_y.c; t++)
        for (int a = 0; a < K; a++) {
            double value = (double)AT(us_free_z, a, t) * (double)AT(us_free_z, a, t);
            squared_by_series[a] += value;
            squared_total += value;
        }
    fprintf(out, "\nshare of the US sum of q_t carried by each whitened series, nu-free fit:");
    for (int a = 0; a < K; a++) fprintf(out, " %s %.3f", series_name[a], squared_by_series[a] / squared_total);
    fprintf(out, "\n");
    fclose(out);

    free(used); free(member); free(winner); free(winner_nu); free(replicate);
    mat_free(us_free_q); mat_free(us_free_z); mat_free(us_held_q); mat_free(us_held_z);
    qvarma_fit_result_free(&us_free); qvarma_fit_result_free(&us_held);
    qvarma_params_free(&shape); qvarma_params_free(&winner_average);
    mat_free(us_y);
    return 0;
}
