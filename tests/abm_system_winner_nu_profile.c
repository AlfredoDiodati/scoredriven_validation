/*
How flat the likelihood is in nu on the replicates of the configuration the
Model Confidence Set kept, and what moving nu does to their loss.

Its 1000 fits put nu anywhere from about 20 to about 10^12, with most between
10^2 and 10^4. When a series is close to Gaussian the Student-t likelihood
changes little once nu is large, so an optimizer can stop anywhere along that
stretch. This script measures that directly rather than inferring it from the
spread of the estimates.

For every replicate of the winner, and for each nu in nu_grid, the t-QVARMA(1,1,2)
is refitted with nu held at that value and the other 41 coordinates free,
starting from the replicate's own cached fit with nu replaced. For each refit it
records the log-likelihood and the loss the confidence set uses: the mean
absolute error between the stacked total impulse response, horizons 0 to 20,
and the US benchmark's.

Two differences are summarised per grid value:
    log-likelihood gain  refit minus the replicate's own fit. Close to zero
                         across a range of nu means the likelihood is flat
                         there; positive means holding nu there found a higher
                         likelihood than the fit the pipeline kept
    loss change          refit minus the replicate's own loss, in the units the
                         confidence set separates configurations on

nu_grid starts at the US benchmark's value, so the first row also says how far
this configuration's data are from the US tails in likelihood terms.

Each refit is cached under out/abm_system_winner_nu_profile/, named for the
replicate and the held nu. The data fingerprint and parameter shape stored in it
invalidate it when either changes; --refit estimates every one again. An
optional first argument limits the run to the first N replicates, for a pilot.

Requires out/abm_system_mcs_joint.csv, out/us_qvarma_spec_choice_p1q1r2_fit.json,
out/us_system.csv, the fit cache and dataset/abm_system/. None are Makefile
prerequisites, because rebuilding them reruns a million fits.

Output, none of it printed:
    out/abm_system_winner_nu_profile.txt   the summary
    out/abm_system_winner_nu_profile.csv   one row per replicate and grid value
    out/abm_system_winner_nu_profile/      the cached refits
*/

#include "applications/abm_system.h"
#include "applications/us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <sys/stat.h>
#include <errno.h>
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
#define MAX_ITERATIONS 20000

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define BENCHMARK_FIT_PATH "out/us_qvarma_spec_choice_p1q1r2_fit.json"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs_joint.csv"
#define PROFILE_DIR "out/abm_system_winner_nu_profile"
#define REPORT_PATH "out/abm_system_winner_nu_profile.txt"
#define TABLE_PATH "out/abm_system_winner_nu_profile.csv"

/* The US benchmark's nu first, then roughly logarithmic steps through the range
   the winner's fits cover. */
static const double nu_grid[] = { 7.198920587, 10, 30, 100, 1e3, 1e4, 1e5, 1e6 };
#define N_GRID ((int)(sizeof nu_grid / sizeof nu_grid[0]))

/* A gain below this is read as no gain: it is the order of the solver's
   function tolerance times a log-likelihood near a thousand. */
#define LIKELIHOOD_TIE 0.01

/* Replicates grouped by their own fitted nu. */
#define N_OWN_GROUPS 3
static const double own_group_upper[N_OWN_GROUPS] = { 1e3, 1e4, INFINITY };
static const char *own_group_name[N_OWN_GROUPS] = { "own nu below 1e3", "own nu 1e3 to 1e4", "own nu above 1e4" };

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

static void stacked_total_irf(const QvarmaParams *m, Mat y, mreal *out) {
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses responses = qvarma_impulse_responses(m, D, options);
    int at = 0;
    for (int h = 0; h <= HORIZON; h++)
        for (int i = 0; i < K * K; i++) out[at++] = responses.total[h].d[i];
    qvarma_impulse_responses_free(&responses);
    mat_free(D);
}

static mreal irf_loss(const QvarmaParams *m, Mat y, const mreal *us_irf) {
    mreal irf[IRF_LENGTH];
    stacked_total_irf(m, y, irf);
    mreal loss = 0;
    for (int i = 0; i < IRF_LENGTH; i++) loss += MABS(irf[i] - us_irf[i]);
    return loss / (mreal)IRF_LENGTH;
}

/* A copy of source with nu replaced, rebuilt through theta so every derived
   quantity is consistent with the new nu. */
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

static void find_winner(char *sample, size_t size) {
    DataFrame df = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    char **name = df_col_string(&df, "model");
    Mat mean_loss = df_col_numeric(&df, "mean_loss");
    Mat in_set = df_col_numeric(&df, "in_set");
    int best = -1;
    for (int i = 0; i < df.r; i++)
        if (AT(in_set, i, 0) > 0 && (best < 0 || AT(mean_loss, i, 0) < AT(mean_loss, best, 0))) best = i;
    assert(best >= 0 && "abm_system_winner_nu_profile: no model is in the confidence set");
    const char *marker = strstr(name[best], "_qvarma_" SPEC_LABEL);
    assert(marker && "abm_system_winner_nu_profile: the winning model is not the p1q1r2 spec");
    size_t length = (size_t)(marker - name[best]);
    assert(length < size);
    memcpy(sample, name[best], length);
    sample[length] = '\0';
    df_free(&df);
}

typedef struct {
    double own_nu, own_log_likelihood, own_loss;
    int own_converged;
} OwnFit;

typedef struct {
    double log_likelihood, loss, gradient_norm;
    int converged, status_known, status, niter;
} Refit;

static double quantile_of(const double *values, int n, double probability) {
    Mat copy = mat_new(n, 1);
    for (int i = 0; i < n; i++) copy.d[i] = (mreal)values[i];
    double q = (double)stats_quantile(copy, (mreal)probability);
    mat_free(copy);
    return q;
}

int main(int argc, char **argv) {
    int force_refit = 0, replicate_limit = -1;
    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--refit") == 0) force_refit = 1;
        else replicate_limit = atoi(argv[a]);
    }

    char sample[64];
    find_winner(sample, sizeof sample);
    char sample_dir[600];
    snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, sample);

    Mat us_y = build_us_block();
    QvarmaParams shape = spec_shape();
    QvarmaFitResult benchmark = qvarma_fit_result_new(&shape);
    int loaded = qvarma_load_fit(&benchmark, us_y, BENCHMARK_FIT_PATH);
    assert(loaded && "abm_system_winner_nu_profile: the US benchmark fit is missing or was fitted to other data");
    mreal us_irf[IRF_LENGTH];
    stacked_total_irf(&benchmark.params, us_y, us_irf);

    int n_replicates = 0;
    int *replicate = abm_system_list_replicates(sample_dir, &n_replicates);
    if (replicate_limit > 0 && replicate_limit < n_replicates) n_replicates = replicate_limit;

    if (mkdir(PROFILE_DIR, 0755) != 0)
        assert(errno == EEXIST && "abm_system_winner_nu_profile: cannot create the refit directory");

    OwnFit *own = (OwnFit*)calloc((size_t)n_replicates, sizeof(OwnFit));
    Refit *refit = (Refit*)calloc((size_t)(n_replicates * N_GRID), sizeof(Refit));

    #pragma omp parallel for schedule(dynamic)
    for (int task = 0; task < n_replicates * N_GRID; task++) {
        int row = task / N_GRID, g = task % N_GRID;
        Mat y = abm_system_read_replicate(sample_dir, replicate[row]);

        char own_path[700];
        snprintf(own_path, sizeof own_path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, sample,
                 replicate[row], SPEC_LABEL);
        QvarmaFitResult own_fit = qvarma_fit_result_new(&shape);
        int own_loaded = qvarma_load_fit(&own_fit, y, own_path);
        assert(own_loaded && "abm_system_winner_nu_profile: a cached fit of the winner is missing or stale");
        if (g == 0) {
            own[row].own_nu = (double)own_fit.params.nu;
            own[row].own_log_likelihood = (double)own_fit.log_likelihood;
            own[row].own_loss = (double)irf_loss(&own_fit.params, y, us_irf);
            own[row].own_converged = own_fit.is_converged;
        }

        char cache_path[700];
        snprintf(cache_path, sizeof cache_path, "%s/replicate_%03d_nu%g_fit.json", PROFILE_DIR,
                 replicate[row], nu_grid[g]);
        QvarmaFitResult result = qvarma_fit_result_new(&shape);
        if (force_refit || !qvarma_load_fit(&result, y, cache_path)) {
            qvarma_fit_result_free(&result);
            QvarmaParams start = with_nu(&own_fit.params, (mreal)nu_grid[g]);
            QvarmaFixedParams fixed = qvarma_fixed_params_new(&start);
            qvarma_fix_block(&fixed, &start, QVARMA_BLOCK_NU);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = MAX_ITERATIONS;
            result = qvarma_fit_with_fixed(y, &start, &fixed, options);
            qvarma_save_fit(&result, y, cache_path);
            qvarma_fixed_params_free(&fixed);
            qvarma_params_free(&start);
        }

        Refit *out = &refit[task];
        out->log_likelihood = (double)result.log_likelihood;
        out->loss = (double)irf_loss(&result.params, y, us_irf);
        out->gradient_norm = (double)result.gradient_norm;
        out->converged = result.is_converged;
        out->status_known = result.status_is_known;
        out->status = (int)result.status;
        out->niter = result.niter;

        qvarma_fit_result_free(&result);
        qvarma_fit_result_free(&own_fit);
        mat_free(y);
    }

    FILE *table = fopen(TABLE_PATH, "w");
    assert(table && "abm_system_winner_nu_profile: cannot open the table for writing");
    fprintf(table, "configuration,replicate,own_nu,own_log_likelihood,own_loss,own_converged,held_nu,"
                   "log_likelihood,loss,gradient_norm,converged,stop_reason,iterations\n");
    for (int row = 0; row < n_replicates; row++)
        for (int g = 0; g < N_GRID; g++) {
            const Refit *f = &refit[row * N_GRID + g];
            fprintf(table, "%s,%d,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%.10g,%.6g,%d,\"%s\",%d\n", sample,
                    replicate[row], own[row].own_nu, own[row].own_log_likelihood, own[row].own_loss,
                    own[row].own_converged, nu_grid[g], f->log_likelihood, f->loss, f->gradient_norm,
                    f->converged, f->status_known ? lbfgs_status_text((LbfgsStatus)f->status) : "unknown",
                    f->niter);
        }
    fclose(table);

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_winner_nu_profile: cannot open the report path for writing");
    fprintf(out, "%s: %d replicates, t-QVARMA(1,1,2) refitted with nu held at each grid value, other 41 "
                 "coordinates free,\n", sample, n_replicates);
    fprintf(out, "starting from the replicate's own cached fit with nu replaced; solver cap %d iterations\n",
            MAX_ITERATIONS);
    fprintf(out, "gain = refit log-likelihood - own fit log-likelihood; loss = MAE of the stacked total IRF, "
                 "horizons 0-%d, against the US benchmark\n", HORIZON);
    fprintf(out, "a gain above %.2g counts as the held nu finding a higher likelihood than the pipeline's fit\n\n",
            LIKELIHOOD_TIE);

    double *gain = (double*)malloc((size_t)n_replicates * sizeof(double));
    double *loss_change = (double*)malloc((size_t)n_replicates * sizeof(double));
    double *loss_level = (double*)malloc((size_t)n_replicates * sizeof(double));
    double own_loss_mean = 0;
    for (int row = 0; row < n_replicates; row++) own_loss_mean += own[row].own_loss;
    own_loss_mean /= n_replicates;

    fprintf(out, "all replicates, own mean loss %.5f\n", own_loss_mean);
    fprintf(out, "  %10s %10s %10s %10s %9s %10s %10s %10s %9s\n", "held nu", "gain p10", "gain p50",
            "gain p90", "higher", "mean loss", "change p50", "change p90", "converged");
    for (int g = 0; g < N_GRID; g++) {
        int higher = 0, converged = 0;
        double mean_loss = 0;
        for (int row = 0; row < n_replicates; row++) {
            const Refit *f = &refit[row * N_GRID + g];
            gain[row] = f->log_likelihood - own[row].own_log_likelihood;
            loss_change[row] = f->loss - own[row].own_loss;
            loss_level[row] = f->loss;
            higher += gain[row] > LIKELIHOOD_TIE;
            converged += f->converged;
            mean_loss += f->loss;
        }
        mean_loss /= n_replicates;
        fprintf(out, "  %10g %10.3f %10.3f %10.3f %4d/%-4d %10.5f %10.5f %10.5f %4d/%-4d\n", nu_grid[g],
                quantile_of(gain, n_replicates, 0.1), quantile_of(gain, n_replicates, 0.5),
                quantile_of(gain, n_replicates, 0.9), higher, n_replicates, mean_loss,
                quantile_of(loss_change, n_replicates, 0.5), quantile_of(loss_change, n_replicates, 0.9),
                converged, n_replicates);
    }

    fprintf(out, "\nmedian gain by the replicate's own fitted nu\n");
    fprintf(out, "  %-20s %6s", "group", "fits");
    for (int g = 0; g < N_GRID; g++) fprintf(out, " %9g", nu_grid[g]);
    fprintf(out, "\n");
    for (int group = 0; group < N_OWN_GROUPS; group++) {
        double lower = group == 0 ? 0 : own_group_upper[group - 1];
        int members = 0;
        for (int row = 0; row < n_replicates; row++)
            members += own[row].own_nu > lower && own[row].own_nu <= own_group_upper[group];
        fprintf(out, "  %-20s %6d", own_group_name[group], members);
        for (int g = 0; g < N_GRID; g++) {
            if (members == 0) { fprintf(out, " %9s", "-"); continue; }
            int at = 0;
            for (int row = 0; row < n_replicates; row++)
                if (own[row].own_nu > lower && own[row].own_nu <= own_group_upper[group])
                    gain[at++] = refit[row * N_GRID + g].log_likelihood - own[row].own_log_likelihood;
            fprintf(out, " %9.3f", quantile_of(gain, at, 0.5));
        }
        fprintf(out, "\n");
    }
    fclose(out);

    free(gain); free(loss_change); free(loss_level);
    free(own); free(refit); free(replicate);
    qvarma_fit_result_free(&benchmark);
    qvarma_params_free(&shape);
    mat_free(us_y);
    return 0;
}
