/*
Why the configuration the Model Confidence Set keeps wins, given that the
parameter set averaged over its replicates has nu near 2853 while the US fit has
nu near 7.2.

Reads every cached fit under out/abm_system_fit_qvarma/ and recomputes each
replicate's impulse response exactly as applications/abm_system_mse_qvarma.c
does: total component, horizons 0 to 20, stacked into 525 elements, loss the
mean absolute error against the US benchmark's. The recomputed losses are
checked against out/abm_system_mse_qvarma_joint.csv before anything else is
reported. Nothing is fitted.

Four questions, one section of the report each:

1. How the winner's fitted parameters are distributed over its replicates,
   next to the US benchmark and the averaged set of
   applications/abm_system_winner_irf.c. The average is taken on the
   optimizer's unconstrained scale, where nu is log(nu - 2), so it can sit far
   from the typical fit when that scale is spread out.
2. Whether, within the winner, the replicates with smaller nu have smaller or
   larger loss.
3. Whether, across configurations, the typical nu of a configuration's fits
   goes with its mean loss.
4. Whether the loss rewards small responses rather than responses shaped like
   the US ones. A configuration whose responses were all zero would score the
   mean absolute US response, so that number and each configuration's mean
   response magnitude are reported beside the mean losses.

A last section splits the loss of the winner, of the configuration eliminated
last before it, and of the median configuration by response series and by
horizon, to show which elements the win comes from.

The winner is the configuration in the confidence set with the smallest mean
loss, and the runner-up the one with the largest elimination round, both read
from out/abm_system_mcs_joint.csv.

Requires out/us_qvarma_spec_choice_p1q1r2_fit.json, out/us_system.csv,
out/abm_system_mse_qvarma_joint.csv, out/abm_system_mcs_joint.csv, the fit cache
and dataset/abm_system/. None are Makefile prerequisites, because rebuilding
the loss table reruns a million impulse responses.

Output, none of it printed:
    out/abm_system_winner_diagnostics.txt
    out/abm_system_winner_diagnostics_parameters.csv      one row per coordinate
                                                          of theta, winner only,
                                                          with the US and
                                                          averaged values also
                                                          on the unconstrained
                                                          scale and the link
    out/abm_system_winner_diagnostics_configurations.csv  one row per
                                                          configuration
    out/abm_system_winner_diagnostics_winner_fits.csv     one row per replicate
                                                          of the winner: nu and
                                                          loss
    out/abm_system_winner_diagnostics_winner_parameters.csv
                                                          one row per replicate
                                                          of the winner: every
                                                          coordinate of theta
                                                          on the model's scale,
                                                          and whether the fit
                                                          converged
    out/abm_system_winner_diagnostics_winner_theta.csv    the same on the
                                                          unconstrained scale
                                                          the optimizer steps
*/

#include "applications/abm_system.h"
#include "applications/us_data.h"
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
#define HORIZON 20
#define IRF_LENGTH (K * K * (HORIZON + 1))
#define SPEC_LABEL "p1q1r2"
#define COLUMN_SUFFIX "_qvarma_p1q1r2"

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define BENCHMARK_FIT_PATH "out/us_qvarma_spec_choice_p1q1r2_fit.json"
#define LOSS_TABLE_PATH "out/abm_system_mse_qvarma_joint.csv"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs_joint.csv"
#define REPORT_PATH "out/abm_system_winner_diagnostics.txt"
#define PARAMETERS_PATH "out/abm_system_winner_diagnostics_parameters.csv"
#define CONFIGURATIONS_PATH "out/abm_system_winner_diagnostics_configurations.csv"
#define WINNER_FITS_PATH "out/abm_system_winner_diagnostics_winner_fits.csv"
#define WINNER_PARAMETERS_PATH "out/abm_system_winner_diagnostics_winner_parameters.csv"
#define WINNER_THETA_PATH "out/abm_system_winner_diagnostics_winner_theta.csv"

/* Horizon groups for the decomposition: impact, the first year, the rest. */
#define N_HORIZON_GROUPS 3
static const char *horizon_group_name[N_HORIZON_GROUPS] = { "h0", "h1-4", "h5-20" };
static int horizon_group(int h) { return h == 0 ? 0 : (h <= 4 ? 1 : 2); }

static const char *series_name[K] = {
    "GDP growth", "Energy growth", "Employment change", "Inflation", "Interest rate"
};

#define N_NU_BINS 5
static const mreal nu_bin_upper[N_NU_BINS] = { 10, 30, 100, 1000, (mreal)INFINITY };
static const char *nu_bin_name[N_NU_BINS] = {
    "2 to 10", "10 to 30", "30 to 100", "100 to 1000", "above 1000"
};
static int nu_bin(mreal nu) {
    int bin = 0;
    while (bin < N_NU_BINS - 1 && nu > nu_bin_upper[bin]) bin++;
    return bin;
}

/* The link as a formula in theta, for labelling the unconstrained figures. */
static void link_formula(QvarmaLink kind, mreal scale, char *out, size_t size) {
    switch (kind) {
        case QVARMA_LINK_IDENTITY: snprintf(out, size, "theta"); break;
        case QVARMA_LINK_TANH: snprintf(out, size, "%g tanh(theta)", (double)scale); break;
        case QVARMA_LINK_EXP: snprintf(out, size, "exp(theta)"); break;
        case QVARMA_LINK_EXP_PLUS_TWO: snprintf(out, size, "exp(theta) + 2"); break;
    }
}

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    return m;
}

/* The US block applications/us_qvarma_spec_choice.c fits. */
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

/* Returns 0 when nu is not above two or a response is not finite, the two
   cases applications/abm_system_mse_qvarma.c counts as a missing cell. */
static int stacked_total_irf(const QvarmaParams *m, Mat y, mreal *out) {
    if (!(m->nu > 2)) return 0;
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses responses = qvarma_impulse_responses(m, D, options);
    int at = 0, finite = 1;
    for (int h = 0; h <= HORIZON; h++)
        for (int i = 0; i < K * K; i++) {
            mreal value = responses.total[h].d[i];
            if (MISNAN(value) || MISINF(value)) finite = 0;
            out[at++] = value;
        }
    qvarma_impulse_responses_free(&responses);
    mat_free(D);
    return finite;
}

typedef struct {
    char sample[64];
    mreal table_mean_loss;
    mreal mean_loss;
    mreal median_nu;
    mreal share_nu_above_100;
    mreal mean_irf_magnitude;
    mreal distance_by_group[K][N_HORIZON_GROUPS];
    int n_used;
} ConfigurationSummary;

/* Mean of |x_i| over the elements of one response series and one horizon
   group, for a K x K response stacked horizon first, row a the response. */
static void accumulate_by_group(const mreal *difference, mreal by_group[K][N_HORIZON_GROUPS]) {
    for (int h = 0; h <= HORIZON; h++)
        for (int a = 0; a < K; a++)
            for (int b = 0; b < K; b++)
                by_group[a][horizon_group(h)] += MABS(difference[h * K * K + a * K + b]);
}

static int group_element_count(int group) {
    int horizons = group == 0 ? 1 : (group == 1 ? 4 : HORIZON - 4);
    return horizons * K;
}

static mreal quantile_of(const mreal *values, int n, mreal probability) {
    Mat copy = mat_new(n, 1);
    memcpy(copy.d, values, (size_t)n * sizeof(mreal));
    mreal q = stats_quantile(copy, probability);
    mat_free(copy);
    return q;
}

static mreal spearman_of(const mreal *x, const mreal *y, int n) {
    Mat a = mat_new(n, 1), b = mat_new(n, 1);
    memcpy(a.d, x, (size_t)n * sizeof(mreal));
    memcpy(b.d, y, (size_t)n * sizeof(mreal));
    mreal rho = stats_spearman(a, b);
    mat_free(a); mat_free(b);
    return rho;
}

int main(void) {
    Mat us_y = build_us_block();
    QvarmaParams shape = spec_shape();
    QvarmaFitResult benchmark = qvarma_fit_result_new(&shape);
    int loaded = qvarma_load_fit(&benchmark, us_y, BENCHMARK_FIT_PATH);
    assert(loaded && "abm_system_winner_diagnostics: the US benchmark fit is missing or was fitted to other data");
    mreal *us_irf = (mreal*)malloc(IRF_LENGTH * sizeof(mreal));
    int us_ok = stacked_total_irf(&benchmark.params, us_y, us_irf);
    assert(us_ok && "abm_system_winner_diagnostics: the US benchmark response is not finite");
    mreal us_magnitude = 0;
    for (int i = 0; i < IRF_LENGTH; i++) us_magnitude += MABS(us_irf[i]);
    us_magnitude /= (mreal)IRF_LENGTH;

    DataFrame losses = df_read_csv(LOSS_TABLE_PATH, csv_read_options_default());
    DataFrame confidence_set = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    int n_configurations = losses.n_cols - 1;
    int n_replicates = losses.r;
    assert(confidence_set.r == n_configurations
           && "abm_system_winner_diagnostics: the confidence set and the loss table disagree on the configurations");

    Mat replicate_index = df_col_numeric(&losses, "replicate");
    char **model = df_col_string(&confidence_set, "model");
    Mat set_mean_loss = df_col_numeric(&confidence_set, "mean_loss");
    Mat in_set = df_col_numeric(&confidence_set, "in_set");
    Mat elimination_round = df_col_numeric(&confidence_set, "elimination_round");

    ConfigurationSummary *summary = (ConfigurationSummary*)calloc((size_t)n_configurations,
                                                                  sizeof(ConfigurationSummary));
    int winner = -1, runner_up = -1;
    for (int c = 0; c < n_configurations; c++) {
        size_t length = strlen(model[c]) - strlen(COLUMN_SUFFIX);
        assert(strcmp(model[c] + length, COLUMN_SUFFIX) == 0
               && "abm_system_winner_diagnostics: a model name without the expected spec suffix");
        memcpy(summary[c].sample, model[c], length);
        summary[c].sample[length] = '\0';
        summary[c].table_mean_loss = AT(set_mean_loss, c, 0);
        if (AT(in_set, c, 0) > 0 && (winner < 0 || AT(set_mean_loss, c, 0) < AT(set_mean_loss, winner, 0)))
            winner = c;
        if (AT(in_set, c, 0) == 0
            && (runner_up < 0 || AT(elimination_round, c, 0) > AT(elimination_round, runner_up, 0)))
            runner_up = c;
    }
    assert(winner >= 0 && runner_up >= 0);

    int n_theta = qvarma_n_theta(&shape);
    Mat nu_by_configuration = mat_new(n_configurations, n_replicates);
    Mat loss_by_configuration = mat_new(n_configurations, n_replicates);
    Mat winner_theta = mat_new(n_replicates, n_theta);
    int *winner_converged = (int*)calloc((size_t)n_replicates, sizeof(int));
    mreal largest_loss_mismatch = 0;
    int n_missing = 0;

    #pragma omp parallel reduction(max:largest_loss_mismatch) reduction(+:n_missing)
    {
        QvarmaParams working = spec_shape();
        QvarmaFitResult winner_fit = qvarma_fit_result_new(&working);
        Vec theta = mat_new(n_theta, 1);
        mreal *irf = (mreal*)malloc(IRF_LENGTH * sizeof(mreal));
        mreal *difference = (mreal*)malloc(IRF_LENGTH * sizeof(mreal));

        #pragma omp for schedule(dynamic)
        for (int c = 0; c < n_configurations; c++) {
            ConfigurationSummary *s = &summary[c];
            Mat table_column = df_col_numeric(&losses, model[c]);
            char sample_dir[560];
            snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, s->sample);
            mreal magnitude_sum = 0, loss_sum = 0;

            for (int row = 0; row < n_replicates; row++) {
                int replicate = (int)AT(replicate_index, row, 0);
                char cache_path[600];
                snprintf(cache_path, sizeof cache_path, "%s/%s/replicate_%03d_%s_fit.json",
                         FIT_DIR, s->sample, replicate, SPEC_LABEL);
                int ok = qvarma_load_params(&working, cache_path);
                if (ok) {
                    Mat y = abm_system_read_replicate(sample_dir, replicate);
                    ok = stacked_total_irf(&working, y, irf);
                    /* The convergence flag sits with the diagnostics, which only
                       the full load reads, and it checks the fit was made on y. */
                    if (ok && c == winner) {
                        int loaded = qvarma_load_fit(&winner_fit, y, cache_path);
                        assert(loaded && "abm_system_winner_diagnostics: a winner fit has no diagnostics or was fitted to other data");
                        winner_converged[row] = winner_fit.is_converged;
                    }
                    mat_free(y);
                }
                if (!ok) {
                    n_missing++;
                    AT(nu_by_configuration, c, row) = (mreal)NAN;
                    AT(loss_by_configuration, c, row) = (mreal)NAN;
                    continue;
                }

                mreal loss = 0, magnitude = 0;
                for (int i = 0; i < IRF_LENGTH; i++) {
                    difference[i] = irf[i] - us_irf[i];
                    loss += MABS(difference[i]);
                    magnitude += MABS(irf[i]);
                }
                loss /= (mreal)IRF_LENGTH;
                magnitude /= (mreal)IRF_LENGTH;
                accumulate_by_group(difference, s->distance_by_group);

                mreal mismatch = MABS(loss - AT(table_column, row, 0));
                if (mismatch > largest_loss_mismatch) largest_loss_mismatch = mismatch;

                AT(nu_by_configuration, c, row) = working.nu;
                AT(loss_by_configuration, c, row) = loss;
                loss_sum += loss;
                magnitude_sum += magnitude;
                s->n_used++;

                if (c == winner) {
                    _qvarma_unlink(&working, theta);
                    for (int i = 0; i < n_theta; i++) AT(winner_theta, row, i) = theta.d[i];
                }
            }

            s->mean_loss = loss_sum / (mreal)s->n_used;
            s->mean_irf_magnitude = magnitude_sum / (mreal)s->n_used;
            for (int a = 0; a < K; a++)
                for (int g = 0; g < N_HORIZON_GROUPS; g++)
                    s->distance_by_group[a][g] /= (mreal)(s->n_used * group_element_count(g));

            mreal *nu = (mreal*)malloc((size_t)s->n_used * sizeof(mreal));
            int kept = 0, above_100 = 0;
            for (int row = 0; row < n_replicates; row++) {
                mreal value = AT(nu_by_configuration, c, row);
                if (MISNAN(value)) continue;
                nu[kept++] = value;
                above_100 += value > 100;
            }
            s->median_nu = quantile_of(nu, kept, (mreal)0.5);
            s->share_nu_above_100 = (mreal)above_100 / (mreal)kept;
            free(nu);
        }

        free(irf);
        free(difference);
        mat_free(theta);
        qvarma_fit_result_free(&winner_fit);
        qvarma_params_free(&working);
    }
    assert(n_missing == 0 && "abm_system_winner_diagnostics: a cell the loss table holds could not be recomputed");

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_winner_diagnostics: cannot open the report path for writing");
    const ConfigurationSummary *w = &summary[winner];
    const ConfigurationSummary *second = &summary[runner_up];

    fprintf(out, "t-QVARMA(1,1,2) fits of %d configurations x %d replicates, impulse responses recomputed\n",
            n_configurations, n_replicates);
    fprintf(out, "total component, horizons 0 to %d, %d elements, loss = mean absolute error against the US benchmark\n",
            HORIZON, IRF_LENGTH);
    fprintf(out, "winner %s (in the confidence set, smallest mean loss), runner-up %s (eliminated last)\n",
            w->sample, second->sample);
    fprintf(out, "check: largest |recomputed loss - %s| over every cell %.3g\n\n",
            LOSS_TABLE_PATH, (double)largest_loss_mismatch);

    /* Section 1. Quantiles commute with a monotone link, so each is taken on
       the constrained scale directly. */
    QvarmaLink *kind = (QvarmaLink*)malloc((size_t)n_theta * sizeof(QvarmaLink));
    mreal *scale = (mreal*)malloc((size_t)n_theta * sizeof(mreal));
    _qvarma_link_kinds(&shape, kind);
    _qvarma_link_scales(&shape, scale);
    Vec us_theta = mat_new(n_theta, 1);
    _qvarma_unlink(&benchmark.params, us_theta);

    FILE *parameters = fopen(PARAMETERS_PATH, "w");
    assert(parameters && "abm_system_winner_diagnostics: cannot open the parameter table for writing");
    fprintf(parameters, "parameter,us_benchmark,averaged_unconstrained,mean_constrained,p05,p25,p50,p75,p95,"
            "us_benchmark_theta,averaged_theta,link\n");
    fprintf(out, "1. %s: fitted parameters over its %d replicates, on the model's scale\n", w->sample, n_replicates);
    fprintf(out, "   averaged = link of the mean on the unconstrained scale, what abm_system_winner_irf uses\n");
    fprintf(out, "   %-18s %10s %10s %10s %10s %10s %10s %10s %10s\n", "parameter", "US", "averaged",
            "mean", "p05", "p25", "p50", "p75", "p95");
    mreal *column = (mreal*)malloc((size_t)n_replicates * sizeof(mreal));
    for (int i = 0; i < n_theta; i++) {
        mreal theta_mean = 0, constrained_mean = 0;
        for (int row = 0; row < n_replicates; row++) {
            theta_mean += AT(winner_theta, row, i);
            column[row] = qvarma_link_forward(kind[i], AT(winner_theta, row, i), scale[i]);
            constrained_mean += column[row];
        }
        theta_mean /= (mreal)n_replicates;
        constrained_mean /= (mreal)n_replicates;
        char name[64];
        _qvarma_theta_name(&shape, i, name, sizeof name);
        mreal us_value = qvarma_link_forward(kind[i], us_theta.d[i], scale[i]);
        mreal averaged = qvarma_link_forward(kind[i], theta_mean, scale[i]);
        mreal q[5] = { quantile_of(column, n_replicates, (mreal)0.05), quantile_of(column, n_replicates, (mreal)0.25),
                       quantile_of(column, n_replicates, (mreal)0.5), quantile_of(column, n_replicates, (mreal)0.75),
                       quantile_of(column, n_replicates, (mreal)0.95) };
        fprintf(out, "   %-18s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", name, (double)us_value,
                (double)averaged, (double)constrained_mean, (double)q[0], (double)q[1], (double)q[2],
                (double)q[3], (double)q[4]);
        char link[32];
        link_formula(kind[i], scale[i], link, sizeof link);
        fprintf(parameters, "\"%s\",%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%s\n", name,
                (double)us_value, (double)averaged, (double)constrained_mean, (double)q[0], (double)q[1],
                (double)q[2], (double)q[3], (double)q[4], (double)us_theta.d[i], (double)theta_mean, link);
    }
    fclose(parameters);

    /* Section 2. */
    int nu_count[N_NU_BINS] = { 0 };
    mreal nu_loss_sum[N_NU_BINS] = { 0 };
    mreal *winner_nu = (mreal*)malloc((size_t)n_replicates * sizeof(mreal));
    mreal *winner_loss = (mreal*)malloc((size_t)n_replicates * sizeof(mreal));
    for (int row = 0; row < n_replicates; row++) {
        winner_nu[row] = AT(nu_by_configuration, winner, row);
        winner_loss[row] = AT(loss_by_configuration, winner, row);
        int bin = nu_bin(winner_nu[row]);
        nu_count[bin]++;
        nu_loss_sum[bin] += winner_loss[row];
    }
    fprintf(out, "\n2. %s: loss against nu, over its replicates\n", w->sample);
    fprintf(out, "   Spearman rank correlation between a replicate's nu and its loss: %.4f\n",
            (double)spearman_of(winner_nu, winner_loss, n_replicates));
    fprintf(out, "   %-14s %8s %10s\n", "nu", "fits", "mean loss");
    for (int bin = 0; bin < N_NU_BINS; bin++) {
        if (nu_count[bin] == 0) fprintf(out, "   %-14s %8d %10s\n", nu_bin_name[bin], 0, "-");
        else fprintf(out, "   %-14s %8d %10.5f\n", nu_bin_name[bin], nu_count[bin],
                     (double)(nu_loss_sum[bin] / (mreal)nu_count[bin]));
    }

    /* Sections 3 and 4. */
    mreal *mean_loss = (mreal*)malloc((size_t)n_configurations * sizeof(mreal));
    mreal *median_nu = (mreal*)malloc((size_t)n_configurations * sizeof(mreal));
    mreal *magnitude = (mreal*)malloc((size_t)n_configurations * sizeof(mreal));
    mreal *magnitude_gap = (mreal*)malloc((size_t)n_configurations * sizeof(mreal));
    int median_nu_rank = 1, magnitude_rank = 1, below_30 = 0;
    mreal largest_table_mismatch = 0;
    for (int c = 0; c < n_configurations; c++) {
        mean_loss[c] = summary[c].mean_loss;
        median_nu[c] = summary[c].median_nu;
        magnitude[c] = summary[c].mean_irf_magnitude;
        magnitude_gap[c] = MABS(summary[c].mean_irf_magnitude - us_magnitude);
        below_30 += median_nu[c] < 30;
        if (median_nu[c] < w->median_nu) median_nu_rank++;
        if (magnitude[c] < w->mean_irf_magnitude) magnitude_rank++;
        mreal mismatch = MABS(summary[c].mean_loss - summary[c].table_mean_loss);
        if (mismatch > largest_table_mismatch) largest_table_mismatch = mismatch;
    }
    fprintf(out, "\n3. across configurations: typical nu against mean loss\n");
    fprintf(out, "   check: largest |recomputed mean loss - %s| %.3g\n", CONFIDENCE_SET_PATH,
            (double)largest_table_mismatch);
    fprintf(out, "   median nu of a configuration's fits: p05 %.4g, p50 %.4g, p95 %.4g; below 30 in %d of %d\n",
            (double)quantile_of(median_nu, n_configurations, (mreal)0.05),
            (double)quantile_of(median_nu, n_configurations, (mreal)0.5),
            (double)quantile_of(median_nu, n_configurations, (mreal)0.95), below_30, n_configurations);
    fprintf(out, "   %s: median nu %.4g (rank %d of %d, 1 = smallest), share of fits with nu above 100 %.3f\n",
            w->sample, (double)w->median_nu, median_nu_rank, n_configurations, (double)w->share_nu_above_100);
    fprintf(out, "   %s: median nu %.4g, share above 100 %.3f\n", second->sample, (double)second->median_nu,
            (double)second->share_nu_above_100);
    fprintf(out, "   Spearman rank correlation, median nu against mean loss: %.4f\n",
            (double)spearman_of(median_nu, mean_loss, n_configurations));

    fprintf(out, "\n4. does the loss reward small responses\n");
    fprintf(out, "   mean |US response| over the %d elements, the loss of an all-zero response: %.5f\n",
            IRF_LENGTH, (double)us_magnitude);
    fprintf(out, "   mean |response| of a configuration's fits: p05 %.5f, p50 %.5f, p95 %.5f\n",
            (double)quantile_of(magnitude, n_configurations, (mreal)0.05),
            (double)quantile_of(magnitude, n_configurations, (mreal)0.5),
            (double)quantile_of(magnitude, n_configurations, (mreal)0.95));
    fprintf(out, "   %s: mean |response| %.5f (rank %d of %d, 1 = smallest), mean loss %.5f\n", w->sample,
            (double)w->mean_irf_magnitude, magnitude_rank, n_configurations, (double)w->mean_loss);
    fprintf(out, "   %s: mean |response| %.5f, mean loss %.5f\n", second->sample,
            (double)second->mean_irf_magnitude, (double)second->mean_loss);
    fprintf(out, "   Spearman rank correlation, mean |response| against mean loss: %.4f\n",
            (double)spearman_of(magnitude, mean_loss, n_configurations));
    fprintf(out, "   Spearman rank correlation, |mean |response| - US value| against mean loss: %.4f\n",
            (double)spearman_of(magnitude_gap, mean_loss, n_configurations));

    /* Section 5: the median configuration by mean loss. */
    int *order = (int*)malloc((size_t)n_configurations * sizeof(int));
    for (int c = 0; c < n_configurations; c++) order[c] = c;
    for (int i = 1; i < n_configurations; i++) {
        int current = order[i], j = i;
        while (j > 0 && mean_loss[order[j - 1]] > mean_loss[current]) { order[j] = order[j - 1]; j--; }
        order[j] = current;
    }
    const ConfigurationSummary *median_configuration = &summary[order[n_configurations / 2]];
    const ConfigurationSummary *shown[3] = { w, second, median_configuration };
    fprintf(out, "\n5. where the loss comes from: mean |response - US response| by series and horizon\n");
    fprintf(out, "   rows are the response series, averaged over the five shocks and the fits\n");
    for (int k = 0; k < 3; k++) {
        fprintf(out, "   %s, mean loss %.5f%s\n", shown[k]->sample, (double)shown[k]->mean_loss,
                k == 2 ? ", median configuration by mean loss" : "");
        fprintf(out, "   %-18s", "series");
        for (int g = 0; g < N_HORIZON_GROUPS; g++) fprintf(out, " %9s", horizon_group_name[g]);
        fprintf(out, "\n");
        for (int a = 0; a < K; a++) {
            fprintf(out, "   %-18s", series_name[a]);
            for (int g = 0; g < N_HORIZON_GROUPS; g++)
                fprintf(out, " %9.5f", (double)shown[k]->distance_by_group[a][g]);
            fprintf(out, "\n");
        }
    }
    fclose(out);

    FILE *configurations = fopen(CONFIGURATIONS_PATH, "w");
    assert(configurations && "abm_system_winner_diagnostics: cannot open the configuration table for writing");
    fprintf(configurations, "configuration,mean_loss,median_nu,share_nu_above_100,mean_irf_magnitude\n");
    for (int c = 0; c < n_configurations; c++)
        fprintf(configurations, "%s,%.10g,%.10g,%.10g,%.10g\n", summary[c].sample, (double)summary[c].mean_loss,
                (double)summary[c].median_nu, (double)summary[c].share_nu_above_100,
                (double)summary[c].mean_irf_magnitude);
    fclose(configurations);

    FILE *winner_fits = fopen(WINNER_FITS_PATH, "w");
    assert(winner_fits && "abm_system_winner_diagnostics: cannot open the winner's fit table for writing");
    fprintf(winner_fits, "configuration,replicate,nu,loss\n");
    for (int row = 0; row < n_replicates; row++)
        fprintf(winner_fits, "%s,%d,%.10g,%.10g\n", w->sample, (int)AT(replicate_index, row, 0),
                (double)winner_nu[row], (double)winner_loss[row]);
    fclose(winner_fits);

    /* The same table twice, on the model's scale and on the optimizer's. */
    const char *winner_table_path[2] = { WINNER_PARAMETERS_PATH, WINNER_THETA_PATH };
    for (int table = 0; table < 2; table++) {
        FILE *winner_table = fopen(winner_table_path[table], "w");
        assert(winner_table && "abm_system_winner_diagnostics: cannot open a winner parameter table for writing");
        fprintf(winner_table, "configuration,replicate");
        for (int i = 0; i < n_theta; i++) {
            char name[64];
            _qvarma_theta_name(&shape, i, name, sizeof name);
            fprintf(winner_table, ",\"%s\"", name);
        }
        fprintf(winner_table, ",is_converged\n");
        for (int row = 0; row < n_replicates; row++) {
            fprintf(winner_table, "%s,%d", w->sample, (int)AT(replicate_index, row, 0));
            for (int i = 0; i < n_theta; i++) {
                mreal theta_value = AT(winner_theta, row, i);
                mreal value = table == 0 ? qvarma_link_forward(kind[i], theta_value, scale[i]) : theta_value;
                fprintf(winner_table, ",%.10g", (double)value);
            }
            fprintf(winner_table, ",%d\n", winner_converged[row]);
        }
        fclose(winner_table);
    }

    free(order); free(mean_loss); free(median_nu); free(magnitude); free(magnitude_gap);
    free(winner_nu); free(winner_loss); free(column); free(kind); free(scale);
    mat_free(us_theta); mat_free(winner_theta); free(winner_converged);
    mat_free(nu_by_configuration); mat_free(loss_by_configuration);
    free(summary); free(us_irf);
    df_free(&losses); df_free(&confidence_set);
    qvarma_fit_result_free(&benchmark);
    qvarma_params_free(&shape);
    mat_free(us_y);
    return 0;
}
