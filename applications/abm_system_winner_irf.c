/*
The impulse response function of the model the Model Confidence Set kept,
with the sign-restricted confidence bands of Blazsek, Escribano and Licht
(2023), section 4.3.

applications/abm_system_mcs.c reduces the 200 candidate models (100 ABM
parameterizations times the two driftless t-QVARMA specs p1q1r2 and p1q1r4)
to the handful the data cannot separate, and writes them to
out/abm_system_mcs_joint.csv with their mean IRF loss against the real-data
fit. This file takes the surviving model with the smallest mean loss and
produces the object the whole comparison exists to look at - that model's
own impulse responses, with bands - rather than one more table of losses.

The steps, in the order main() runs them:

  1. Read out/abm_system_mcs_joint.csv, keep the rows with in_set = 1, and
     take the one with the smallest mean_loss. Its name carries both the ABM
     sample directory and the spec, "<sample>_qvarma_p1q1r<N>", the naming
     applications/abm_system_mse_qvarma.c gave its own loss columns.
  2. Read every replicate_<NNN>_p1q1r<N>_fit.json that sample's directory
     under out/abm_system_fit_qvarma/ holds, one per replicate.
  3. Average theta over those replicates, coordinate by coordinate, in the
     unconstrained space the JSON stores - not in the constrained space.
     That is what the averaged object is: the mean of what the optimizer
     stepped. Averaging the constrained parameters instead would give a
     different model, since the link is nonlinear and the mean of a
     transform is not the transform of the mean, and it could also leave a
     mean outside the constrained set (Phi_star has to stay inside
     (-b,b) and nu above two, and an average of admissible values need not
     be admissible in a space with corners). params_from_theta then maps the
     one averaged vector through the link, so every quantity the impulse
     response formulae read - Phi_star, Psi_star, Psi_dag, Omega_inv, nu -
     is the link applied to the average, once.
  4. The impulse responses of (14) to (24) need D, the sample average of the
     score Jacobian of (21), which depends on the data and not only on the
     parameters. The averaged parameter set does not belong to any one
     replicate, so D is averaged the same way the parameters were: computed
     on each of the 108 replicates at the averaged parameters and averaged
     over them. Every replicate has the same number of periods, so this is
     exactly the D of the pooled sample - one average over 108 x 400
     observations rather than 108 separate ones.
  5. Bands by 4.3: ten million random rotations Q, kept when the impact
     matrix carries the signs of the paper's Table 1, and the 10, 50 and 90
     percent percentiles of the responses over the kept rotations.

The paper's Table 1 restricts three series against three shocks:

                              supply   demand   monetary
    US real GDP growth          +        +         -
    US inflation rate           -        +         -
    effective federal funds     .        +         +

This system has five series, three of which are the paper's: GDP growth,
inflation and the interest rate. Energy demand growth and employment change
are not in the paper's system and carry no restriction, and the two shocks
beyond the paper's three are unrestricted in every row - so the identified
set here is the paper's, widened by what this system adds and by nothing
else. A shock is labelled by the column it occupies, and only the first
three columns carry a label at all.

What the fitted models being averaged actually are matters for reading the
result: 105 of this sample's own 108 fits stopped at
applications/abm_system_fit_qvarma.c's own 2000-iteration cap rather than
converging, which the manifest records per replicate. The average is over
what those runs reached, and it is not the average of 108 maximum-likelihood
estimates.

Output:

  out/abm_system_winner_irf.csv        one row per (component, horizon,
                                       shock, response): the unrotated
                                       response, and the lower, median and
                                       upper band, with both the index and a
                                       readable name for the shock and the
                                       response. Long rather than wide so a
                                       plot can filter it directly. The
                                       components are qvarma.h's own five -
                                       contemporaneous, stationary,
                                       cointegrated, their total, and the
                                       total cumulated over the horizon,
                                       which is the response of the level a
                                       differenced series differences.
  out/abm_system_winner_irf_theta.json the averaged parameter set, in
                                       save_params' own format, so it
                                       reloads through load_params.
  out/abm_system_winner_irf_manifest.txt
                                       which model won and by what loss, how
                                       many replicates went into the average
                                       and how many of them converged, the
                                       maximum companion eigenvalue modulus
                                       and nu at the averaged parameters,
                                       and how many rotations satisfied
                                       Table 1.

The unrotated column is the response at Q = I, the one point in the
identified set that Omega_inv's own Cholesky orientation picks out. It is
not the estimate the paper plots and it need not satisfy Table 1 itself; the
median over the accepted rotations is what a figure should show, with the
lower and upper as its band. Both are written because the difference between
them is a statement about how much the sign restrictions actually pin down.

shock_name names the sign-restricted reading of a column and nothing else. A
rotation mixes the reduced-form innovations, so column b is the shock whose
impact signs match column b of Table 1, and it is not the shock to series b -
two of the five columns satisfy no restriction at all and so have no name.
The unrotated column in the same row is a different identification, where
Omega_inv is lower triangular and column b is the part of series b's own
innovation orthogonal to the series before it in this file's row order: there,
and only there, a column does name a series. A consumer reading the unrotated
column should label by the shock index and this row order, not by shock_name.

The cumulative band comes out of impulse_bands with the rest rather than
being summed here afterwards, which would be a different and wrong number: a
quantile of a sum is not the sum of quantiles, so the band of the cumulated
path has to be taken over the cumulated path of each accepted rotation. It
costs nothing extra to get right, since a rotation multiplies every horizon's
response on the right and the sum of (response times Q) is (sum of responses)
times Q.

applications/abm_system_winner_irf_plots.py reads the CSV and writes the
figures; it is the only consumer, and nothing here is written for its
convenience rather than for a reader of the CSV itself.

Requires out/abm_system_mcs_joint.csv (applications/abm_system_mcs.c),
out/abm_system_fit_qvarma/ (applications/abm_system_fit_qvarma.c) and
dataset/abm_system/ (applications/abm_system_extract.c). Nothing printed.
*/

#include "abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define P 1
#define Q 1
#define N_REPLICATES 108
#define HORIZON 20

/* Ten times the paper's own million draws, at its own 10 and 90 percent
   percentiles. The paper's three series and three shocks kept 9561 of a
   million; this system's five series and five shocks put the same eight
   restrictions on a larger orthogonal group and keep a few hundred, which
   is too few for a stable tenth percentile. Ten million draws cost about
   half a minute here and keep a few thousand. The seed is this file's own,
   so a rerun reproduces the same band. */
#define BAND_DRAWS 10000000
#define BAND_SEED 20260822

#define MCS_PATH "out/abm_system_mcs_joint.csv"
#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define IRF_PATH "out/abm_system_winner_irf.csv"
#define THETA_PATH "out/abm_system_winner_irf_theta.json"
#define MANIFEST_PATH "out/abm_system_winner_irf_manifest.txt"

/* Column order of a response matrix, the three the paper's Table 1 names
   followed by the two this system adds. */
enum { SHOCK_SUPPLY, SHOCK_DEMAND, SHOCK_MONETARY };

static const char *shock_name[K] = {
    "Supply shock", "Demand shock", "Monetary policy shock",
    "Unrestricted shock 4", "Unrestricted shock 5"
};

/* Row order of a response matrix, abm_system.h's own ROW_ enum. */
static const char *series_name[K] = {
    "GDP growth", "Energy demand growth", "Employment change",
    "Inflation", "Interest rate"
};

static const char *component_name[QVARMA_N_IMPULSE_COMPONENTS] = {
    "contemporaneous", "stationary", "cointegrated", "total", "cumulative"
};

static QvarmaParams spec_shape(int r) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, r, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

typedef struct {
    char sample[256];
    char spec_label[32];
    int r;
    mreal mean_loss;
    int n_in_set;
    int n_models;
} Winner;

/* The kept model with the smallest mean loss. The name splits at
   "_qvarma_p", so a drift-carrying row would fail the assert rather than be
   read as a driftless one - this file includes qvarma.h and can only
   produce a driftless model's impulse responses. */
static Winner find_winner(void) {
    DataFrame df = df_read_csv(MCS_PATH, csv_read_options_default());
    char **name = df_col_string(&df, "model");
    Mat mean_loss = df_col_numeric(&df, "mean_loss");
    Mat in_set = df_col_numeric(&df, "in_set");

    Winner w;
    w.n_models = df.r;
    w.n_in_set = 0;
    int best = -1;
    for (int i = 0; i < df.r; i++) {
        if (AT(in_set, i, 0) == 0) continue;
        w.n_in_set++;
        if (best < 0 || AT(mean_loss, i, 0) < AT(mean_loss, best, 0)) best = i;
    }
    assert(best >= 0 && "abm_system_winner_irf: no model in out/abm_system_mcs_joint.csv "
                        "is in the confidence set");

    const char *marker = strstr(name[best], "_qvarma_p");
    assert(marker && "abm_system_winner_irf: the winning model is not a driftless t-QVARMA - "
                     "this file includes qvarma.h and cannot read a qvarma_d fit");
    size_t sample_length = (size_t)(marker - name[best]);
    assert(sample_length < sizeof w.sample);
    memcpy(w.sample, name[best], sample_length);
    w.sample[sample_length] = '\0';

    const char *spec = marker + strlen("_qvarma_");
    snprintf(w.spec_label, sizeof w.spec_label, "%s", spec);
    const char *lag_count = strrchr(spec, 'r');
    assert(lag_count && "abm_system_winner_irf: the winning spec label carries no r<N>");
    w.r = atoi(lag_count + 1);
    assert(w.r >= 1);
    w.mean_loss = AT(mean_loss, best, 0);

    df_free(&df);
    return w;
}

/* The stored theta of one cached fit, appended into total, with the
   replicate's own convergence flag reported to the manifest. Returns 0 when
   the file is missing or was written for a different shape. */
static int accumulate_theta(const char *path, int n, Vec total, int *converged) {
    JsonValue *root = json_parse_file(path);
    if (!root) return 0;
    JsonValue *values = json_object_get(root, "theta");
    if (!values || json_array_len(values) != n) {
        json_free(root);
        return 0;
    }
    for (int i = 0; i < n; i++)
        total.d[i] += (mreal)json_as_number(json_array_get(values, i));

    JsonValue *diagnostics = json_object_get(root, "fit");
    JsonValue *flag = diagnostics ? json_object_get(diagnostics, "is_converged") : NULL;
    *converged = flag ? (int)json_as_number(flag) : 0;
    json_free(root);
    return 1;
}

/* One replicate's own K x T series, same convention
   applications/abm_system_mse_qvarma.c's own read_y uses. */
static Mat read_y(const char *sample, int replicate) {
    char csv_path[560];
    snprintf(csv_path, sizeof csv_path, "%s/%s/replicate_%03d.csv", INPUT_DIR, sample, replicate);
    DataFrame df = df_read_csv(csv_path, csv_read_options_default());
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

/* Table 1, on the five rows and five columns this system has. Zero is
   unrestricted, which is what the two series the paper does not have and the
   two shocks beyond its three carry everywhere. */
static Mat paper_sign_restrictions(void) {
    Mat s = mat_new(K, K);
    AT(s, ROW_GDP_GROWTH, SHOCK_SUPPLY) = 1;
    AT(s, ROW_INFLATION, SHOCK_SUPPLY) = -1;
    AT(s, ROW_GDP_GROWTH, SHOCK_DEMAND) = 1;
    AT(s, ROW_INFLATION, SHOCK_DEMAND) = 1;
    AT(s, ROW_INTEREST_RATE, SHOCK_DEMAND) = 1;
    AT(s, ROW_GDP_GROWTH, SHOCK_MONETARY) = -1;
    AT(s, ROW_INFLATION, SHOCK_MONETARY) = -1;
    AT(s, ROW_INTEREST_RATE, SHOCK_MONETARY) = 1;
    return s;
}

static void write_irf_csv(const QvarmaImpulseResponses *point, const QvarmaImpulseBands *bands,
                          const char *path) {
    FILE *out = fopen(path, "w");
    assert(out && "abm_system_winner_irf: cannot open the impulse response path for writing");
    fprintf(out, "component,horizon,shock,response,shock_name,response_name,"
                 "unrotated,lower,median,upper\n");

    const Mat *point_component[QVARMA_N_IMPULSE_COMPONENTS] = {
        point->contemporaneous, point->stationary, point->cointegrated,
        point->total, point->cumulative };
    const Mat *lower_component[QVARMA_N_IMPULSE_COMPONENTS] = {
        bands->lower.contemporaneous, bands->lower.stationary,
        bands->lower.cointegrated, bands->lower.total, bands->lower.cumulative };
    const Mat *median_component[QVARMA_N_IMPULSE_COMPONENTS] = {
        bands->median.contemporaneous, bands->median.stationary,
        bands->median.cointegrated, bands->median.total, bands->median.cumulative };
    const Mat *upper_component[QVARMA_N_IMPULSE_COMPONENTS] = {
        bands->upper.contemporaneous, bands->upper.stationary,
        bands->upper.cointegrated, bands->upper.total, bands->upper.cumulative };

    for (int c = 0; c < QVARMA_N_IMPULSE_COMPONENTS; c++)
        for (int h = 0; h <= point->horizon; h++)
            for (int b = 0; b < K; b++)
                for (int a = 0; a < K; a++)
                    fprintf(out, "%s,%d,%d,%d,%s,%s,%.10g,%.10g,%.10g,%.10g\n",
                            component_name[c], h, b + 1, a + 1, shock_name[b], series_name[a],
                            (double)AT(point_component[c][h], a, b),
                            (double)AT(lower_component[c][h], a, b),
                            (double)AT(median_component[c][h], a, b),
                            (double)AT(upper_component[c][h], a, b));
    fclose(out);
}

int main(void) {
    FILE *manifest = fopen(MANIFEST_PATH, "w");
    assert(manifest && "abm_system_winner_irf: cannot open the manifest path for writing");

    Winner winner = find_winner();
    fprintf(manifest, "winning model  %s_qvarma_%s\n", winner.sample, winner.spec_label);
    fprintf(manifest, "mean IRF loss  %.10g, the smallest of the %d models in the confidence "
                      "set of %d\n", (double)winner.mean_loss, winner.n_in_set, winner.n_models);
    fprintf(manifest, "source         %s\n\n", MCS_PATH);

    QvarmaParams m = spec_shape(winner.r);
    int n = qvarma_n_theta(&m);
    Vec average = mat_new(n, 1);
    int n_loaded = 0, n_converged = 0;
    for (int replicate = 0; replicate < N_REPLICATES; replicate++) {
        char cache_path[560];
        snprintf(cache_path, sizeof cache_path, "%s/%s/replicate_%03d_%s_fit.json",
                 FIT_DIR, winner.sample, replicate, winner.spec_label);
        int converged = 0;
        if (accumulate_theta(cache_path, n, average, &converged)) {
            n_loaded++;
            n_converged += converged;
        } else {
            fprintf(manifest, "missing or wrong shape: replicate %03d\n", replicate);
        }
    }
    assert(n_loaded > 0 && "abm_system_winner_irf: no cached fit of the winning model could "
                           "be read - run abm_system_fit_qvarma first");
    for (int i = 0; i < n; i++) average.d[i] /= (mreal)n_loaded;

    fprintf(manifest, "averaged over  %d of %d replicates, %d parameters each, in the "
                      "unconstrained space\n", n_loaded, N_REPLICATES, n);
    fprintf(manifest, "converged      %d of the %d averaged fits; the rest stopped at the "
                      "fitting script's own iteration cap\n\n", n_converged, n_loaded);

    qvarma_params_from_theta(average, &m);
    qvarma_save_params(&m, THETA_PATH);
    fprintf(manifest, "nu             %.6f\n", (double)m.nu);
    fprintf(manifest, "companion      maximum eigenvalue modulus %.6f\n\n",
            (double)qvarma_max_eigenvalue_modulus(&m));

    /* D on each replicate at the averaged parameters, averaged over them -
       every replicate has the same number of periods, so this is the score
       Jacobian of the pooled sample. */
    Mat D = mat_new(K, K);
    int n_series = 0, n_periods = 0;
    for (int replicate = 0; replicate < N_REPLICATES; replicate++) {
        Mat y = read_y(winner.sample, replicate);
        Mat own = qvarma_mean_score_jacobian(&m, y);
        for (int i = 0; i < K * K; i++) D.d[i] += own.d[i];
        n_periods = y.c;
        n_series++;
        mat_free(own);
        mat_free(y);
    }
    for (int i = 0; i < K * K; i++) D.d[i] /= (mreal)n_series;
    fprintf(manifest, "score Jacobian averaged over %d replicates of %d periods at the "
                      "averaged parameters\n\n", n_series, n_periods);

    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses point = qvarma_impulse_responses(&m, D, options);

    Mat restrictions = paper_sign_restrictions();
    QvarmaImpulseBandOptions band_options = qvarma_default_impulse_band_options();
    band_options.n_draws = BAND_DRAWS;
    Rng rng = rng_new(BAND_SEED, 0);
    QvarmaImpulseBands bands = qvarma_impulse_bands(&rng, &m, D, restrictions, options, band_options);
    assert(bands.n_accepted > 0 && "abm_system_winner_irf: no rotation satisfied Table 1, so "
                                   "every band entry would be not-a-number");

    fprintf(manifest, "bands          %d of %d rotations satisfied Table 1, percentiles "
                      "%.0f%%, 50%%, %.0f%%\n", bands.n_accepted, bands.n_draws,
            (double)band_options.lower_percentile * 100.0,
            (double)band_options.upper_percentile * 100.0);
    fprintf(manifest, "horizon        0 to %d, seed %d\n", HORIZON, BAND_SEED);
    fprintf(manifest, "restrictions   GDP growth (+,+,-), inflation (-,+,-), interest rate "
                      "(.,+,+) on the supply, demand and monetary policy shocks; energy "
                      "demand growth, employment change and shocks 4 and 5 unrestricted\n");
    fprintf(manifest, "shock names    apply to the band alone. The unrotated column is the "
                      "Cholesky orientation, where shock b is series b's own orthogonalized "
                      "innovation under the row order %s, %s, %s, %s, %s\n",
            series_name[0], series_name[1], series_name[2], series_name[3], series_name[4]);

    write_irf_csv(&point, &bands, IRF_PATH);

    qvarma_impulse_bands_free(&bands);
    qvarma_impulse_responses_free(&point);
    mat_free(restrictions);
    mat_free(D);
    mat_free(average);
    qvarma_params_free(&m);
    fclose(manifest);
    return 0;
}
