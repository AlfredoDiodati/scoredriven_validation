/*
Whether the t-QVARMA(1,1,2) estimates of the configuration the Model Confidence
Set keeps are jointly normal across its replicates, on the unconstrained scale
theta the optimizer steps.

Reads out/abm_system_winner_diagnostics_winner_theta.csv, written by
studies/abm_system_winner_diagnostics.c: one row per replicate, the d = 42
coordinates of theta and whether the fit converged. Nothing is fitted and no
fit cache is read. The test is run twice, on every replicate and on the
replicates whose fit converged.

The test is Mardia's (1970). With xbar the sample mean, S the sample covariance
divided by n, and

    g_ij = (x_i - xbar)^T S^-1 (x_j - xbar),

the multivariate skewness and kurtosis are

    b1 = (1/n^2) sum_i sum_j g_ij^3,    b2 = (1/n) sum_i g_ii^2.

Under multivariate normality n b1 / 6 is asymptotically chi-squared on
d(d+1)(d+2)/6 degrees of freedom, and b2 is asymptotically normal with mean
d(d+2) and variance 8 d(d+2) / n. Each is reported in two forms: the asymptotic
one above, and a finite-sample one, the skewness statistic multiplied by
Mardia's (1974) correction k = (d+1)(n+1)(n+3) / (n((n+1)(d+1) - 6)) and the
kurtosis standardized with its exact null mean d(d+2)(n-1)/(n+1) and variance
8d(d+2)(n-3)(n-d-1)(n-d+1) / ((n+1)^2 (n+3)(n+5)). Both tests are invariant to
an affine change of the data, so the scale of each coordinate does not matter.

S^-1 is never formed: with S = L L^T, y_i solves L y_i = x_i - xbar and
g_ij = y_i^T y_j. The diagonal g_ii is the Mahalanobis form inside the
multivariate normal log-density, so it is checked against dist/mv/gauss.h's
mvgauss_logpdf evaluated at N(xbar, S).

How far the chi-squared and normal approximations can be trusted at d = 42 is
measured rather than assumed: the same four tests are run on samples drawn from
a 42-dimensional standard normal with mvgauss_sample, at both sample sizes, and
the share rejected at the 5 per cent level is reported. By affine invariance a
standard normal stands for every normal. Seed and replication count are the
constants below.

Mardia's test is a general statistical tool, not specific to this project. It
is written here against et_al's primitives and will be ported to et_al's
inference/ directory later; the functions take a plain n x d sample for that
reason.

Output, none of it printed:
    out/abm_system_winner_normality.txt
*/

#include <et_al./dist/mv/gauss.h>
#include <et_al./linalg/solver.h>
#include <et_al./random/random.h>
#include <et_al./stats.h>
#include <et_al./special.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define THETA_PATH "out/abm_system_winner_diagnostics_winner_theta.csv"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs.csv"
#define REPORT_PATH "out/abm_system_winner_normality.txt"

#define SIGNIFICANCE_LEVEL 0.05
#define CALIBRATION_REPLICATIONS 1000
#define CALIBRATION_SEED 20260917

typedef struct {
    int n, d;
    double b1, b2;
    double skewness_statistic, skewness_p;
    double skewness_corrected_statistic, skewness_corrected_p;
    double skewness_df;
    double kurtosis_z, kurtosis_p;
    double kurtosis_exact_z, kurtosis_exact_p;
    /* Largest |g_ii - Mahalanobis form from mvgauss_logpdf| / max(1, g_ii). */
    double largest_logpdf_mismatch;
    /* |mean_i g_ii - d|, zero up to rounding when S divides by n. */
    double trace_identity_gap;
} MardiaResult;

/* Mardia's skewness and kurtosis tests on an n x d sample, one observation per
   row. The sample covariance must be positive definite, so n > d and no
   coordinate may be an exact linear combination of the others. */
static MardiaResult mardia_test(Mat x) {
    int n = x.r, d = x.c;
    assert(n > d && d >= 1 && "mardia_test: needs more observations than dimensions");
    assert(mat_all_finite(x) && "mardia_test: the sample holds a NaN or an infinity");

    Mat mean = stats_vec_mean(x);
    Mat covariance = stats_autocov(x, 0);
    Mat cholesky_factor = mat_chol(covariance);

    Mat whitened = mat_new(n, d);
    Vec centered = mat_new(d, 1);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < d; k++) centered.d[k] = AT(x, i, k) - AT(mean, 0, k);
        Vec solved = vec_triangular_solve(cholesky_factor, centered, 'L', 'N', 'N');
        for (int k = 0; k < d; k++) AT(whitened, i, k) = solved.d[k];
        mat_free(solved);
    }
    Mat whitened_transpose = mat_T(whitened);
    Mat g = mat_mul(whitened, whitened_transpose);

    double cube_sum = 0, square_sum = 0, trace_sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double value = AT(g, i, j);
            cube_sum += value * value * value;
        }
        square_sum += AT(g, i, i) * AT(g, i, i);
        trace_sum += AT(g, i, i);
    }

    MardiaResult result;
    result.n = n;
    result.d = d;
    result.b1 = cube_sum / ((double)n * (double)n);
    result.b2 = square_sum / (double)n;

    double dd = (double)d, nn = (double)n;
    result.skewness_df = dd * (dd + 1) * (dd + 2) / 6;
    result.skewness_statistic = nn * result.b1 / 6;
    result.skewness_p = special_chi_squared_sf(result.skewness_statistic, result.skewness_df);
    double correction = (dd + 1) * (nn + 1) * (nn + 3) / (nn * ((nn + 1) * (dd + 1) - 6));
    result.skewness_corrected_statistic = correction * result.skewness_statistic;
    result.skewness_corrected_p = special_chi_squared_sf(result.skewness_corrected_statistic, result.skewness_df);

    double kurtosis_mean = dd * (dd + 2);
    result.kurtosis_z = (result.b2 - kurtosis_mean) / sqrt(8 * dd * (dd + 2) / nn);
    result.kurtosis_p = 2 * special_norm_cdf(-fabs(result.kurtosis_z));
    double exact_mean = dd * (dd + 2) * (nn - 1) / (nn + 1);
    double exact_variance = 8 * dd * (dd + 2) * (nn - 3) * (nn - dd - 1) * (nn - dd + 1)
                          / ((nn + 1) * (nn + 1) * (nn + 3) * (nn + 5));
    result.kurtosis_exact_z = (result.b2 - exact_mean) / sqrt(exact_variance);
    result.kurtosis_exact_p = 2 * special_norm_cdf(-fabs(result.kurtosis_exact_z));

    /* logpdf_i = -g_ii/2 - log det(S)/2 - (d/2) log(2 pi), with log det(S)
       from the factor already in hand. */
    double log_determinant = 0;
    for (int k = 0; k < d; k++) log_determinant += 2 * log(AT(cholesky_factor, k, k));
    Mat log_density = mvgauss_logpdf(x, mean, covariance);
    result.largest_logpdf_mismatch = 0;
    for (int i = 0; i < n; i++) {
        double quadratic_form = -2 * AT(log_density, i, 0) - log_determinant - 2 * dd * MVGAUSS_HALF_LOG_2PI;
        double mismatch = fabs(quadratic_form - AT(g, i, i)) / fmax(1.0, AT(g, i, i));
        if (mismatch > result.largest_logpdf_mismatch) result.largest_logpdf_mismatch = mismatch;
    }
    result.trace_identity_gap = fabs(trace_sum / nn - dd);

    mat_free(log_density);
    mat_free(g);
    mat_free(whitened_transpose);
    mat_free(centered);
    mat_free(whitened);
    mat_free(cholesky_factor);
    mat_free(covariance);
    mat_free(mean);
    return result;
}

/* Largest ratio of eigenvalues of the sample correlation matrix, which is
   scale free, unlike the condition number of the covariance itself. */
static double correlation_condition_number(Mat x) {
    Mat covariance = stats_autocov(x, 0);
    Mat correlation = mat_new(x.c, x.c);
    for (int a = 0; a < x.c; a++)
        for (int b = 0; b < x.c; b++)
            AT(correlation, a, b) = AT(covariance, a, b) / sqrt(AT(covariance, a, a) * AT(covariance, b, b));
    double condition = (double)mat_cond(correlation);
    mat_free(correlation);
    mat_free(covariance);
    return condition;
}

typedef struct {
    int n, replications;
    int skewness_rejections, skewness_corrected_rejections;
    int kurtosis_rejections, kurtosis_exact_rejections;
} Calibration;

/* Rejections at SIGNIFICANCE_LEVEL over samples of n draws from a
   d-dimensional standard normal. Replication r draws from stream r of the
   fixed seed, so the whole table reproduces. */
static Calibration calibrate(int n, int d) {
    Calibration calibration = { n, CALIBRATION_REPLICATIONS, 0, 0, 0, 0 };
    Mat loc = mat_new(1, d);
    Mat identity = mat_eye(d);
    for (int r = 0; r < CALIBRATION_REPLICATIONS; r++) {
        Rng rng = rng_new(CALIBRATION_SEED, (uint64_t)r);
        Mat sample = mvgauss_sample(&rng, loc, identity, n);
        MardiaResult result = mardia_test(sample);
        calibration.skewness_rejections += result.skewness_p < SIGNIFICANCE_LEVEL;
        calibration.skewness_corrected_rejections += result.skewness_corrected_p < SIGNIFICANCE_LEVEL;
        calibration.kurtosis_rejections += result.kurtosis_p < SIGNIFICANCE_LEVEL;
        calibration.kurtosis_exact_rejections += result.kurtosis_exact_p < SIGNIFICANCE_LEVEL;
        mat_free(sample);
    }
    mat_free(identity);
    mat_free(loc);
    return calibration;
}

/* The configuration the confidence set keeps: in the set, smallest mean loss. */
static char *winning_model(void) {
    DataFrame confidence_set = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    char **model = df_col_string(&confidence_set, "model");
    Mat mean_loss = df_col_numeric(&confidence_set, "mean_loss");
    Mat in_set = df_col_numeric(&confidence_set, "in_set");
    int winner = -1;
    for (int c = 0; c < confidence_set.r; c++)
        if (AT(in_set, c, 0) > 0 && (winner < 0 || AT(mean_loss, c, 0) < AT(mean_loss, winner, 0)))
            winner = c;
    assert(winner >= 0 && "abm_system_winner_normality: the confidence set keeps no configuration");
    char *name = frame_strdup(model[winner]);
    df_free(&confidence_set);
    return name;
}

/* A p-value as text. Beyond about z = 38 or a chi-squared statistic far in its
   tail the probability is below the smallest double and evaluates to zero,
   which is not the same claim as zero. */
static void format_p(double p, char *out, size_t size) {
    if (p > 0) snprintf(out, size, "%.4g", p);
    else snprintf(out, size, "below %.1e (underflows double precision)", DBL_MIN);
}

static void write_result(FILE *out, const char *label, const MardiaResult *result, double condition) {
    char skewness_p[64], skewness_corrected_p[64], kurtosis_p[64], kurtosis_exact_p[64];
    format_p(result->skewness_p, skewness_p, sizeof skewness_p);
    format_p(result->skewness_corrected_p, skewness_corrected_p, sizeof skewness_corrected_p);
    format_p(result->kurtosis_p, kurtosis_p, sizeof kurtosis_p);
    format_p(result->kurtosis_exact_p, kurtosis_exact_p, sizeof kurtosis_exact_p);
    fprintf(out, "%s: n = %d fits, d = %d coordinates of theta\n", label, result->n, result->d);
    fprintf(out, "   condition number of the sample correlation matrix %.4g\n", condition);
    fprintf(out, "   check: largest relative gap between g_ii and the Mahalanobis form of mvgauss_logpdf %.3g\n",
            result->largest_logpdf_mismatch);
    fprintf(out, "   check: |mean of g_ii - d| %.3g\n", result->trace_identity_gap);
    fprintf(out, "   skewness b1 = %.6g, degrees of freedom %.0f\n", result->b1, result->skewness_df);
    fprintf(out, "      asymptotic n b1 / 6 = %.6g, p %s\n", result->skewness_statistic, skewness_p);
    fprintf(out, "      small sample k n b1 / 6 = %.6g, p %s\n", result->skewness_corrected_statistic,
            skewness_corrected_p);
    fprintf(out, "   kurtosis b2 = %.6g, d(d+2) = %d\n", result->b2, result->d * (result->d + 2));
    fprintf(out, "      asymptotic moments z = %.6g, p %s\n", result->kurtosis_z, kurtosis_p);
    fprintf(out, "      exact moments z = %.6g, p %s\n", result->kurtosis_exact_z, kurtosis_exact_p);
}

static void write_calibration(FILE *out, const Calibration *calibration) {
    double denominator = (double)calibration->replications;
    fprintf(out, "   n = %d: skewness asymptotic %.3f, small sample %.3f; kurtosis asymptotic %.3f, exact %.3f\n",
            calibration->n, calibration->skewness_rejections / denominator,
            calibration->skewness_corrected_rejections / denominator,
            calibration->kurtosis_rejections / denominator, calibration->kurtosis_exact_rejections / denominator);
}

int main(void) {
    DataFrame table = df_read_csv(THETA_PATH, csv_read_options_default());
    char **configuration = df_col_string(&table, "configuration");
    for (int row = 1; row < table.r; row++)
        assert(strcmp(configuration[row], configuration[0]) == 0
               && "abm_system_winner_normality: the theta table holds more than one configuration");
    char *winner = winning_model();
    size_t name_length = strlen(configuration[0]);
    assert(strncmp(winner, configuration[0], name_length) == 0 && winner[name_length] == '_'
           && "abm_system_winner_normality: the theta table is not for the configuration the confidence set keeps; "
              "run make study-abm_system_winner_diagnostics");

    int n_parameters = 0;
    for (int k = 0; k < table.n_cols; k++) {
        const char *name = table.columns[k].name;
        if (table.columns[k].type == COL_NUMERIC && strcmp(name, "replicate") != 0 && strcmp(name, "is_converged") != 0)
            n_parameters++;
    }
    Mat is_converged = df_col_numeric(&table, "is_converged");
    int n_converged = 0;
    for (int row = 0; row < table.r; row++) n_converged += AT(is_converged, row, 0) > 0;

    Mat all_fits = mat_new(table.r, n_parameters);
    Mat converged_fits = mat_new(n_converged, n_parameters);
    int column = 0;
    for (int k = 0; k < table.n_cols; k++) {
        const char *name = table.columns[k].name;
        if (table.columns[k].type != COL_NUMERIC || strcmp(name, "replicate") == 0 || strcmp(name, "is_converged") == 0)
            continue;
        Mat values = df_col_numeric(&table, name);
        int kept = 0;
        for (int row = 0; row < table.r; row++) {
            AT(all_fits, row, column) = AT(values, row, 0);
            if (AT(is_converged, row, 0) > 0) AT(converged_fits, kept++, column) = AT(values, row, 0);
        }
        column++;
    }

    MardiaResult all_result = mardia_test(all_fits);
    MardiaResult converged_result = mardia_test(converged_fits);
    Calibration all_calibration = calibrate(all_fits.r, n_parameters);
    Calibration converged_calibration = calibrate(converged_fits.r, n_parameters);

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_winner_normality: cannot open the report path for writing");
    fprintf(out, "Mardia's test of joint normality of the t-QVARMA(1,1,2) estimates of %s\n", configuration[0]);
    fprintf(out, "data: %s, one row per replicate, theta on the unconstrained scale\n", THETA_PATH);
    fprintf(out, "null hypothesis: the rows are independent draws from one %d-dimensional normal\n", n_parameters);
    fprintf(out, "rejection: p below %.2f; p-values are two-sided for kurtosis, upper tail for skewness\n\n",
            SIGNIFICANCE_LEVEL);

    write_result(out, "all replicates", &all_result, correlation_condition_number(all_fits));
    fprintf(out, "\n");
    write_result(out, "replicates whose fit converged", &converged_result, correlation_condition_number(converged_fits));

    fprintf(out, "\ncalibration: share of %d samples from a %d-dimensional standard normal rejected at %.2f\n",
            CALIBRATION_REPLICATIONS, n_parameters, SIGNIFICANCE_LEVEL);
    fprintf(out, "   draws by mvgauss_sample, seed %d, stream r for replication r, same n as each test above\n",
            CALIBRATION_SEED);
    fprintf(out, "   a correct test rejects %.2f; the binomial standard error of the share is %.4f\n",
            SIGNIFICANCE_LEVEL,
            sqrt(SIGNIFICANCE_LEVEL * (1 - SIGNIFICANCE_LEVEL) / CALIBRATION_REPLICATIONS));
    write_calibration(out, &all_calibration);
    write_calibration(out, &converged_calibration);
    fclose(out);

    free(winner);
    mat_free(converged_fits);
    mat_free(all_fits);
    df_free(&table);
    return 0;
}
