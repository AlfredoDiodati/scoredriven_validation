/*
Whether the degrees of freedom nu change the impulse response function the Model
Confidence Set compares configurations on.

The loss in applications/abm_system_irf_loss.c is the mean absolute error
between the stacked total impulse responses of two t-QVARMA(1,1,2) fits,
horizons 0 to 20. nu enters those responses through the contemporaneous scale
(nu/(nu-2))^(1/2), through the factor ((nu-2) nu)^(1/2) Omega_inv D, and
through the averaged score Jacobian D. The same fit also re-estimates the
score loadings Psi_star and Psi_dag, which can absorb a change in nu. If they
absorb it, a fit with thin tails and a fit with nu near 7 give nearly the same
responses, and the loss cannot tell them apart.

Three parameter sets on the US data, all at the Fisher-relation partition
K_star 3, K_dagger 2, R 1:

    benchmark   the fit with nu free, out/us_qvarma_spec_choice_p1q1r2_fit.json,
                the response every simulated replicate is compared against
    swapped     the benchmark with nu replaced by the winning configuration's
                nu and nothing else changed, which is the direct effect of nu
                alone
    refitted    every parameter but nu re-estimated with nu held at the winning
                configuration's nu, starting from the benchmark, which is the
                effect once the loadings are allowed to respond

The winning configuration's nu is read from out/abm_system_winner_irf_theta.json,
the parameter set applications/abm_system_winner_irf.c averages over that
configuration's replicates.

The distance of swapped and refitted from benchmark is reported in the loss's
own units, over all 525 elements and split into horizon 0 and horizons 1 to 20,
beside the spread of mean losses across the configurations in
out/abm_system_mcs.csv, which is the scale the confidence set separates
configurations on.

The refitted fit is cached to out/us_qvarma_nu_sensitivity_nu<nu>_fit.json. The
data fingerprint and the parameter shape stored in it invalidate it when the
data or the specification change, and the held nu, to four decimals, is in the
file name. Pass --refit to estimate it again regardless.

Requires out/us_qvarma_spec_choice_p1q1r2_fit.json (make
app-us_qvarma_spec_choice), out/abm_system_mcs.csv (make
app-abm_system_mcs) and out/abm_system_winner_irf_theta.json (make
app-abm_system_winner_irf). None is a Makefile prerequisite of this study,
which only reads what the main pipeline already wrote.

Output, none of it printed:
    out/us_qvarma_nu_sensitivity.txt
    out/us_qvarma_nu_sensitivity_nu<nu>_fit.json
*/

#include "applications/us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <stdlib.h>
#include <string.h>

#define K 5
#define K_STAR 3
#define P 1
#define Q 1
#define RLAG 2
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define PERIODS (ESTIMATION_PERIODS - 1)
#define HORIZON 20
#define MAX_ITERATIONS 8000

#define BENCHMARK_FIT_PATH "out/us_qvarma_spec_choice_p1q1r2_fit.json"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs.csv"
#define WINNER_PARAMS_PATH "out/abm_system_winner_irf_theta.json"
#define REPORT_PATH "out/us_qvarma_nu_sensitivity.txt"

enum { ROW_GDP_GROWTH, ROW_EN_GROWTH, ROW_EMPLOYMENT_CHANGE, ROW_INFLATION, ROW_INTEREST_RATE };

/* The same block applications/us_qvarma_spec_choice.c fits, so the benchmark
   cache's data fingerprint matches. */
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

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    return m;
}

/* A copy of source with nu replaced. The round trip through theta rebuilds
   every quantity the parameter struct derives, so nothing computed at the old
   nu survives. */
static QvarmaParams with_nu(const QvarmaParams *source, mreal nu) {
    QvarmaParams copy = spec_shape();
    Vec theta = mat_new(qvarma_n_theta(source), 1);
    _qvarma_unlink(source, theta);
    qvarma_params_from_theta(theta, &copy);
    copy.nu = nu;
    _qvarma_unlink(&copy, theta);
    qvarma_params_from_theta(theta, &copy);
    mat_free(theta);
    return copy;
}

/* total[0..HORIZON], each K x K, stacked horizon 0 first: the vector the loss
   in applications/abm_system_irf_loss.c is taken over. */
static Vec stacked_total_irf(const QvarmaParams *m, Mat y) {
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses responses = qvarma_impulse_responses(m, D, options);

    Vec stacked = mat_new(K * K * (HORIZON + 1), 1);
    int at = 0;
    for (int h = 0; h <= HORIZON; h++)
        for (int i = 0; i < K * K; i++) stacked.d[at++] = responses.total[h].d[i];

    qvarma_impulse_responses_free(&responses);
    mat_free(D);
    return stacked;
}

typedef struct {
    mreal all_horizons;
    mreal impact;
    mreal after_impact;
} IrfDistance;

static IrfDistance irf_distance(Vec reference, Vec other) {
    int impact_length = K * K;
    IrfDistance distance;
    distance.all_horizons = stats_mae(reference, other);
    distance.impact = stats_mae(mat_slice(reference, 0, impact_length, 0, 1),
                                mat_slice(other, 0, impact_length, 0, 1));
    distance.after_impact = stats_mae(mat_slice(reference, impact_length, reference.r, 0, 1),
                                      mat_slice(other, impact_length, other.r, 0, 1));
    return distance;
}

/* Median over horizons 1 to 20 of the ratio of another response to the
   benchmark's, element by element, over benchmark entries large enough that
   the ratio is not noise around zero. It says whether the responses after
   impact were scaled by a common factor, and by how much. */
static mreal median_ratio_after_impact(Vec reference, Vec other) {
    int impact_length = K * K;
    int n = reference.r - impact_length;
    Mat ratios = mat_new(n, 1);
    int kept = 0;
    for (int i = impact_length; i < reference.r; i++)
        if (MABS(reference.d[i]) > (mreal)1e-3) ratios.d[kept++] = other.d[i] / reference.d[i];
    assert(kept > 0 && "us_qvarma_nu_sensitivity: every benchmark response after impact is near zero");
    mreal median = stats_median(mat_slice(ratios, 0, kept, 0, 1));
    mat_free(ratios);
    return median;
}

typedef struct {
    int n_configurations;
    mreal smallest_mean_loss;
    mreal largest_mean_loss;
    mreal first_gap;
    char winner[64];
    char runner_up[64];
} LossScale;

/* The spread of per-configuration mean losses, and the gap between the two
   configurations the confidence set separated last. */
static LossScale read_loss_scale(void) {
    FILE *probe = fopen(CONFIDENCE_SET_PATH, "r");
    assert(probe && "us_qvarma_nu_sensitivity: out/abm_system_mcs.csv is missing; run make app-abm_system_mcs");
    fclose(probe);

    DataFrame table = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    Mat mean_loss = df_col_numeric(&table, "mean_loss");

    int best = 0, second = -1;
    LossScale scale;
    scale.n_configurations = table.r;
    scale.smallest_mean_loss = AT(mean_loss, 0, 0);
    scale.largest_mean_loss = AT(mean_loss, 0, 0);
    for (int i = 1; i < table.r; i++) {
        mreal value = AT(mean_loss, i, 0);
        if (value > scale.largest_mean_loss) scale.largest_mean_loss = value;
        if (value < AT(mean_loss, best, 0)) { second = best; best = i; }
        else if (second < 0 || value < AT(mean_loss, second, 0)) second = i;
    }
    scale.smallest_mean_loss = AT(mean_loss, best, 0);
    scale.first_gap = AT(mean_loss, second, 0) - AT(mean_loss, best, 0);
    snprintf(scale.winner, sizeof scale.winner, "%s", df_col_string(&table, "model")[best]);
    snprintf(scale.runner_up, sizeof scale.runner_up, "%s", df_col_string(&table, "model")[second]);

    df_free(&table);
    return scale;
}

/* nu of the parameter set averaged over the winning configuration's replicates. */
static mreal read_winner_nu(void) {
    QvarmaParams winner = spec_shape();
    int loaded = qvarma_load_params(&winner, WINNER_PARAMS_PATH);
    assert(loaded && "us_qvarma_nu_sensitivity: out/abm_system_winner_irf_theta.json is missing or has "
                     "another shape; run make app-abm_system_winner_irf");
    mreal nu = winner.nu;
    qvarma_params_free(&winner);
    return nu;
}

static QvarmaFitResult fit_with_nu_held(Mat y, const QvarmaParams *benchmark, mreal held_nu,
                                        int force_refit, const char *cache_path) {
    QvarmaFitResult cached = qvarma_fit_result_new(benchmark);
    if (!force_refit && qvarma_load_fit(&cached, y, cache_path)) return cached;
    qvarma_fit_result_free(&cached);

    QvarmaParams start = with_nu(benchmark, held_nu);
    QvarmaFixedParams fixed = qvarma_fixed_params_new(&start);
    qvarma_fix_block(&fixed, &start, QVARMA_BLOCK_NU);

    QvarmaFitOptions options = qvarma_default_fit_options();
    options.max_iterations = MAX_ITERATIONS;
    QvarmaFitResult result = qvarma_fit_with_fixed(y, &start, &fixed, options);
    qvarma_save_fit(&result, y, cache_path);

    qvarma_fixed_params_free(&fixed);
    qvarma_params_free(&start);
    return result;
}

static void write_distance_row(FILE *out, const char *label, IrfDistance distance, mreal ratio) {
    fprintf(out, "  %-10s all %.5f   horizon 0 %.5f   horizons 1-20 %.5f   median ratio after impact %.3f\n",
            label, (double)distance.all_horizons, (double)distance.impact,
            (double)distance.after_impact, (double)ratio);
}

int main(int argc, char **argv) {
    int force_refit = argc > 1 && strcmp(argv[1], "--refit") == 0;

    Mat original = load_us_system();
    Mat y = build_block(original);
    mat_free(original);

    QvarmaParams shape = spec_shape();
    QvarmaFitResult benchmark = qvarma_fit_result_new(&shape);
    qvarma_params_free(&shape);
    int loaded = qvarma_load_fit(&benchmark, y, BENCHMARK_FIT_PATH);
    assert(loaded && "us_qvarma_nu_sensitivity: the benchmark fit is missing or was fitted to other data; "
                     "run make app-us_qvarma_spec_choice");

    mreal winner_nu = read_winner_nu();
    char cache_path[128];
    snprintf(cache_path, sizeof cache_path, "out/us_qvarma_nu_sensitivity_nu%.4f_fit.json",
             (double)winner_nu);
    QvarmaFitResult refitted = fit_with_nu_held(y, &benchmark.params, winner_nu, force_refit, cache_path);
    QvarmaParams swapped = with_nu(&benchmark.params, winner_nu);

    Vec benchmark_irf = stacked_total_irf(&benchmark.params, y);
    Vec swapped_irf = stacked_total_irf(&swapped, y);
    Vec refitted_irf = stacked_total_irf(&refitted.params, y);

    IrfDistance swapped_distance = irf_distance(benchmark_irf, swapped_irf);
    IrfDistance refitted_distance = irf_distance(benchmark_irf, refitted_irf);
    mreal swapped_ratio = median_ratio_after_impact(benchmark_irf, swapped_irf);
    mreal refitted_ratio = median_ratio_after_impact(benchmark_irf, refitted_irf);
    LossScale scale = read_loss_scale();

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "us_qvarma_nu_sensitivity: cannot open the report path for writing");

    fprintf(out, "t-QVARMA(1,1,2) on the US data, K_star 3, K_dagger 2, R 1, 1973Q2 to 2019Q4, %d quarters\n",
            PERIODS);
    fprintf(out, "impulse responses: total component, horizons 0 to %d, stacked into %d elements\n\n",
            HORIZON, K * K * (HORIZON + 1));

    fprintf(out, "parameter sets\n");
    fprintf(out, "  benchmark  nu %.4f free, log-likelihood %.4f, gradient norm %.4g, converged %s\n",
            (double)benchmark.params.nu, (double)benchmark.log_likelihood,
            (double)benchmark.gradient_norm, benchmark.is_converged ? "yes" : "no");
    fprintf(out, "  winner     nu %.4f, from %s\n", (double)winner_nu, WINNER_PARAMS_PATH);
    fprintf(out, "  swapped    benchmark with nu set to the winner's, nothing re-estimated\n");
    fprintf(out, "  refitted   nu held at the winner's, the other %d coordinates re-estimated from the benchmark\n",
            qvarma_n_theta(&refitted.params) - 1);
    fprintf(out, "             log-likelihood %.4f, gradient norm over the free coordinates %.4g,\n",
            (double)refitted.log_likelihood, (double)refitted.gradient_norm);
    fprintf(out, "             converged %s (%s), %d iterations in the run that wrote it\n",
            refitted.is_converged ? "yes" : "no",
            refitted.status_is_known ? lbfgs_status_text(refitted.status) : "reason not stored",
            refitted.niter);
    fprintf(out, "  likelihood ratio, 2 (benchmark - refitted): %.4f\n\n",
            2 * (double)(benchmark.log_likelihood - refitted.log_likelihood));

    fprintf(out, "distance from the benchmark response, mean absolute error in the loss's units\n");
    write_distance_row(out, "swapped", swapped_distance, swapped_ratio);
    write_distance_row(out, "refitted", refitted_distance, refitted_ratio);
    fprintf(out, "  median ratio: refitted or swapped response over benchmark response, element by element,\n");
    fprintf(out, "  horizons 1 to 20, over the benchmark elements with absolute value above 0.001\n\n");

    fprintf(out, "scale the confidence set works on, from %s\n", CONFIDENCE_SET_PATH);
    fprintf(out, "  mean loss over %d configurations: smallest %.5f (%s), largest %.5f, range %.5f\n",
            scale.n_configurations, (double)scale.smallest_mean_loss, scale.winner,
            (double)scale.largest_mean_loss,
            (double)(scale.largest_mean_loss - scale.smallest_mean_loss));
    fprintf(out, "  gap between the two smallest mean losses (%s, %s): %.5f\n",
            scale.winner, scale.runner_up, (double)scale.first_gap);

    fclose(out);

    mat_free(benchmark_irf);
    mat_free(swapped_irf);
    mat_free(refitted_irf);
    qvarma_params_free(&swapped);
    qvarma_fit_result_free(&refitted);
    qvarma_fit_result_free(&benchmark);
    mat_free(y);
    return 0;
}
