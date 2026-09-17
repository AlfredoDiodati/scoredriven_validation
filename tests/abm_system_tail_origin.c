/*
Whether the simulations lack the tails of the US data, or whether the tails are
there and the fitted t-QVARMA does not see them.

The US fit gives nu near 7 and the fits to the simulated replicates give nu in
the hundreds to thousands. Two explanations predict different things:
    the simulations have thin tails   their raw series, with no model involved,
                                      have smaller skewness and excess kurtosis
                                      than the US series, and the replicates
                                      with the heaviest raw tails are the ones
                                      fitted with the smallest nu
    the fit misses the tails          the raw series are as heavy-tailed as the
                                      US ones, and a replicate's fitted nu does
                                      not follow its raw kurtosis

Part 1 compares the raw series. Every replicate of every configuration under
dataset/abm_system/ is read as stored, which is already the transformation the
auxiliary model is fitted on: GDP growth, energy demand growth and the change in
employment as 100 times a first difference, inflation as 100 times the first
difference of the log price level, and the interest rate in percentage points,
in levels. The interest rate is a persistent level, so its moments mostly
describe where the level wanders rather than how large its shocks are; its
first difference is reported as a sixth series for that reason. The US series
are built the same way from out/us_system.csv, 187 quarters.

A replicate has 400 periods and the US sample 187. The sample kurtosis of a
shorter series is more variable, and from finite samples of heavy-tailed data
it is also biased downward, so every simulated moment is computed twice: over
all 400 periods, and over the first 187 so the comparison with the US value is
at equal length.

Part 2 sets the same moments against the fitted nu of each replicate, read from
the cached fits under out/abm_system_fit_qvarma/. Configurations are ranked by
the median nu of their fits, and the rank correlation between nu and each
series' excess kurtosis is taken across configurations and across all
replicates.

Moments are the plain moment ratios: skewness m3 / m2^(3/2) and excess
kurtosis m4 / m2^2 - 3, with m_k the k-th central sample moment divided by n.
Under Gaussian data their standard errors are close to (6/n)^(1/2) and
(24/n)^(1/2): 0.18 and 0.36 at n = 187, 0.12 and 0.24 at n = 400.
sample_skewness and sample_excess_kurtosis are general statistics and are meant
to move to et_al's stats.h; they are checked here against hand-computed values
before anything else runs, and the study stops if a check fails.

Reads only. Writes:
    out/abm_system_tail_origin_report.txt
    out/abm_system_tail_origin_configurations.csv   one row per configuration
*/
#include "applications/abm_system.h"
#include "applications/us_data.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./random/random.h>
#include <et_al./frame/npz.h>
#include <dirent.h>
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

#define INPUT_DIR "dataset/abm_system"
#define FIT_DIR "out/abm_system_fit_qvarma"
#define REPORT_PATH "out/abm_system_tail_origin_report.txt"
#define CONFIGURATIONS_PATH "out/abm_system_tail_origin_configurations.csv"

#define US_PERIODS (ESTIMATION_PERIODS - 1)
#define N_RANKED_SHOWN 10

/* The five stored series and the first difference of the interest rate. */
#define N_SERIES (K + 1)
#define SERIES_INTEREST_CHANGE K
static const char *series_name[N_SERIES] = {
    "GDP growth", "energy growth", "employment change", "inflation", "interest rate",
    "interest rate change"
};
static const char *series_column[N_SERIES] = {
    "gdp_growth", "energy_growth", "employment_change", "inflation", "interest_rate",
    "interest_rate_change"
};

enum { MOMENT_VARIANCE, MOMENT_SKEWNESS, MOMENT_EXCESS_KURTOSIS, N_MOMENTS };
static const char *moment_name[N_MOMENTS] = { "variance", "skewness", "excess kurtosis" };
static const char *moment_column[N_MOMENTS] = { "variance", "skewness", "excess_kurtosis" };

enum { LENGTH_FULL, LENGTH_US, N_LENGTHS };
static const char *length_name[N_LENGTHS] = { "400 periods", "first 187 periods" };
static const char *length_column[N_LENGTHS] = { "t400", "t187" };

/* Central sample moments of every element of x, divided by n, walking the
   stride so a view is read correctly. Returns 0 when there are fewer than two
   elements. */
static int central_moments(Mat x, double *m2, double *m3, double *m4) {
    long n = (long)x.r * x.c;
    if (n < 2) return 0;
    double mean = 0;
    for (int i = 0; i < x.r; i++)
        for (int j = 0; j < x.c; j++) mean += (double)AT(x, i, j);
    mean /= (double)n;
    double s2 = 0, s3 = 0, s4 = 0;
    for (int i = 0; i < x.r; i++)
        for (int j = 0; j < x.c; j++) {
            double d = (double)AT(x, i, j) - mean;
            double d2 = d * d;
            s2 += d2;
            s3 += d2 * d;
            s4 += d2 * d2;
        }
    *m2 = s2 / (double)n;
    *m3 = s3 / (double)n;
    *m4 = s4 / (double)n;
    return 1;
}

/* Sample skewness m3 / m2^(3/2). Not a number when x has fewer than two
   elements or no spread, where the ratio is undefined. */
static double sample_skewness(Mat x) {
    double m2, m3, m4;
    if (!central_moments(x, &m2, &m3, &m4) || !(m2 > 0)) return NAN;
    return m3 / pow(m2, 1.5);
}

/* Sample excess kurtosis m4 / m2^2 - 3, zero in expectation for Gaussian data
   as n grows. Not a number under the same conditions as sample_skewness. */
static double sample_excess_kurtosis(Mat x) {
    double m2, m3, m4;
    if (!central_moments(x, &m2, &m3, &m4) || !(m2 > 0)) return NAN;
    return m4 / (m2 * m2) - 3.0;
}

static double sample_variance(Mat x) {
    double m2, m3, m4;
    if (!central_moments(x, &m2, &m3, &m4)) return NAN;
    return m2;
}

/* The checks the two moment functions have to pass before the study uses them.
   Hand-computed cases, a strided view, the degenerate inputs, and large
   samples from distributions whose moments are known. */
static void check_moment_functions(FILE *report) {
    double tolerance = 1e-9;

    double symmetric_values[] = { 1, 2, 3, 4, 5 };
    Mat symmetric = mat_new(5, 1);
    for (int i = 0; i < 5; i++) symmetric.d[i] = (mreal)symmetric_values[i];
    /* m2 = 2, m3 = 0, m4 = 34 / 5 = 6.8, so excess kurtosis 6.8 / 4 - 3. */
    assert(fabs(sample_skewness(symmetric)) < tolerance);
    assert(fabs(sample_excess_kurtosis(symmetric) - (-1.3)) < tolerance);

    /* {0, 0, 0, 1}: m2 = 3/16, m3 = 3/32, m4 = 21/256, so skewness
       (3/32) / (3/16)^(3/2) = 2 / sqrt(3) and excess kurtosis 7/3 - 3. */
    Mat skewed = mat_new(4, 1);
    skewed.d[3] = 1;
    assert(fabs(sample_skewness(skewed) - 2.0 / sqrt(3.0)) < tolerance);
    assert(fabs(sample_excess_kurtosis(skewed) - (7.0 / 3.0 - 3.0)) < tolerance);
    /* Reflecting the data flips the sign of the skewness and leaves the
       kurtosis alone. */
    Mat reflected = mat_new(4, 1);
    reflected.d[3] = -1;
    assert(fabs(sample_skewness(reflected) + 2.0 / sqrt(3.0)) < tolerance);
    assert(fabs(sample_excess_kurtosis(reflected) - (7.0 / 3.0 - 3.0)) < tolerance);
    /* Location and scale leave both unchanged. */
    Mat shifted = mat_new(4, 1);
    for (int i = 0; i < 4; i++) shifted.d[i] = (mreal)(-7.5 + 1e3 * skewed.d[i]);
    assert(fabs(sample_skewness(shifted) - 2.0 / sqrt(3.0)) < 1e-7);
    assert(fabs(sample_excess_kurtosis(shifted) - (7.0 / 3.0 - 3.0)) < 1e-7);

    /* A column view of a wider matrix has stride 2 and must give the same
       answer as the contiguous copy. */
    Mat wide = mat_new(5, 2);
    for (int i = 0; i < 5; i++) { AT(wide, i, 0) = (mreal)(100 * i * i); AT(wide, i, 1) = (mreal)symmetric_values[i]; }
    Mat column_view = mat_slice(wide, 0, 5, 1, 2);
    assert(column_view.stride == 2);
    assert(fabs(sample_excess_kurtosis(column_view) - (-1.3)) < tolerance);
    assert(fabs(sample_skewness(column_view)) < tolerance);

    /* No spread and a single element have no defined ratio. */
    Mat constant = mat_new(6, 1);
    for (int i = 0; i < 6; i++) constant.d[i] = 3;
    Mat single = mat_new(1, 1);
    assert(MISNAN(sample_skewness(constant)) && MISNAN(sample_excess_kurtosis(constant)));
    assert(MISNAN(sample_skewness(single)) && MISNAN(sample_excess_kurtosis(single)));

    /* Large samples: Gaussian (skewness 0, excess kurtosis 0) and Student t
       with 10 degrees of freedom (skewness 0, excess kurtosis 6 / (10 - 4) = 1).
       Tolerances are about six standard errors at this size. */
    int n = 2000000;
    Rng rng = rng_new(20260916, 0);
    Mat gaussian = mat_new(n, 1), student = mat_new(n, 1);
    for (int i = 0; i < n; i++) {
        gaussian.d[i] = (mreal)rng_normal(&rng);
        double chi_squared = 2.0 * rng_gamma(&rng, 5.0);
        student.d[i] = (mreal)(rng_normal(&rng) / sqrt(chi_squared / 10.0));
    }
    double gaussian_skewness = sample_skewness(gaussian), gaussian_kurtosis = sample_excess_kurtosis(gaussian);
    double student_skewness = sample_skewness(student), student_kurtosis = sample_excess_kurtosis(student);
    assert(fabs(gaussian_skewness) < 0.01);
    assert(fabs(gaussian_kurtosis) < 0.02);
    assert(fabs(student_skewness) < 0.03);
    assert(fabs(student_kurtosis - 1.0) < 0.15);

    fprintf(report, "moment functions: every check passed. Large samples, n = %d, seed 20260916:\n", n);
    fprintf(report, "  Gaussian      skewness %8.4f (0), excess kurtosis %8.4f (0)\n", gaussian_skewness, gaussian_kurtosis);
    fprintf(report, "  Student t(10) skewness %8.4f (0), excess kurtosis %8.4f (1)\n\n", student_skewness, student_kurtosis);

    mat_free(symmetric); mat_free(skewed); mat_free(reflected); mat_free(shifted); mat_free(wide);
    mat_free(constant); mat_free(single); mat_free(gaussian); mat_free(student);
}

/* The moments of the six series of one K x T block over its first `periods`
   columns, written as out[series][moment]. */
static void block_moments(Mat y, int periods, double out[N_SERIES][N_MOMENTS]) {
    assert(periods >= 3 && periods <= y.c);
    for (int k = 0; k < K; k++) {
        Mat row = mat_slice(y, k, k + 1, 0, periods);
        out[k][MOMENT_VARIANCE] = sample_variance(row);
        out[k][MOMENT_SKEWNESS] = sample_skewness(row);
        out[k][MOMENT_EXCESS_KURTOSIS] = sample_excess_kurtosis(row);
    }
    Mat change = mat_new(1, periods - 1);
    for (int t = 1; t < periods; t++)
        AT(change, 0, t - 1) = AT(y, ROW_INTEREST_RATE, t) - AT(y, ROW_INTEREST_RATE, t - 1);
    out[SERIES_INTEREST_CHANGE][MOMENT_VARIANCE] = sample_variance(change);
    out[SERIES_INTEREST_CHANGE][MOMENT_SKEWNESS] = sample_skewness(change);
    out[SERIES_INTEREST_CHANGE][MOMENT_EXCESS_KURTOSIS] = sample_excess_kurtosis(change);
    mat_free(change);
}

static Mat build_us_block(void) {
    Mat original = load_us_system();
    Mat y = mat_new(K, US_PERIODS);
    for (int t = 1; t < ESTIMATION_PERIODS; t++) {
        int c = t - 1;
        AT(y, ROW_GDP_GROWTH, c) = AT(original, LOG_GDP, t) - AT(original, LOG_GDP, t - 1);
        AT(y, ROW_EN_GROWTH, c) = AT(original, LOG_ENERGY_DEMAND, t) - AT(original, LOG_ENERGY_DEMAND, t - 1);
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = AT(original, EMPLOYMENT, t) - AT(original, EMPLOYMENT, t - 1);
        AT(y, ROW_INFLATION, c) = AT(original, LOG_CPI, t) - AT(original, LOG_CPI, t - 1);
        AT(y, ROW_INTEREST_RATE, c) = AT(original, INTEREST_RATE, t);
    }
    mat_free(original);
    return y;
}

static char **list_configurations(int *count) {
    DIR *handle = opendir(INPUT_DIR);
    assert(handle && "abm_system_tail_origin: cannot open dataset/abm_system");
    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (strncmp(entry->d_name, "cop_", 4) != 0) continue;
        if (n == cap) { cap = cap ? 2 * cap : 1024; names = (char**)realloc(names, (size_t)cap * sizeof(char*)); }
        names[n++] = frame_strdup(entry->d_name);
    }
    closedir(handle);
    for (int i = 1; i < n; i++) {
        char *current = names[i];
        int j = i;
        while (j > 0 && strcmp(names[j - 1], current) > 0) { names[j] = names[j - 1]; j--; }
        names[j] = current;
    }
    *count = n;
    return names;
}

/* Everything the study keeps about one replicate. */
typedef struct {
    int replicate;
    double nu;
    double moment[N_LENGTHS][N_SERIES][N_MOMENTS];
} ReplicateRecord;

typedef struct {
    char name[32];
    int n_replicates;
    ReplicateRecord *record;
    double median_nu;
    double median_moment[N_LENGTHS][N_SERIES][N_MOMENTS];
} Configuration;

static double median_of(const double *values, int n) {
    Mat copy = mat_new(n, 1);
    int kept = 0;
    for (int i = 0; i < n; i++) if (!MISNAN(values[i])) copy.d[kept++] = (mreal)values[i];
    double median = kept ? (double)stats_median(mat_slice(copy, 0, kept, 0, 1)) : NAN;
    mat_free(copy);
    return median;
}

static double quantile_of(const double *values, int n, double probability) {
    Mat copy = mat_new(n, 1);
    int kept = 0;
    for (int i = 0; i < n; i++) if (!MISNAN(values[i])) copy.d[kept++] = (mreal)values[i];
    double q = kept ? (double)stats_quantile(mat_slice(copy, 0, kept, 0, 1), (mreal)probability) : NAN;
    mat_free(copy);
    return q;
}

static double spearman_of(const double *x, const double *y, int n) {
    Mat a = mat_new(n, 1), b = mat_new(n, 1);
    int kept = 0;
    for (int i = 0; i < n; i++)
        if (!MISNAN(x[i]) && !MISNAN(y[i]) && !MISINF(x[i]) && !MISINF(y[i])) {
            a.d[kept] = (mreal)x[i];
            b.d[kept] = (mreal)y[i];
            kept++;
        }
    double rho = kept >= 3 ? (double)stats_spearman(mat_slice(a, 0, kept, 0, 1), mat_slice(b, 0, kept, 0, 1)) : NAN;
    mat_free(a); mat_free(b);
    return rho;
}

/* One configuration: every archive read once, every replicate's moments at
   both lengths, and every replicate's fitted nu. */
static void read_configuration(Configuration *c) {
    char dir[256];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, c->name);
    int *replicate = abm_system_list_replicates(dir, &c->n_replicates);
    c->record = (ReplicateRecord*)calloc((size_t)c->n_replicates, sizeof(ReplicateRecord));

    QvarmaParams params = qvarma_params_new(K, K_STAR, P, Q, RLAG, R, SHARED_BETA, WARMUP_LONGEST);
    params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    params.phi_star_bound = PHI_STAR_BOUND;

    int row = 0;
    while (row < c->n_replicates) {
        char path[512];
        abm_system_batch_path(path, sizeof path, dir, replicate[row]);
        DataFrame df = df_read_npz(path);
        Mat index = df_col_numeric(&df, abm_system_column_names()[K]);
        int batch = replicate[row] / ABM_SYSTEM_BATCH;

        for (; row < c->n_replicates && replicate[row] / ABM_SYSTEM_BATCH == batch; row++) {
            ReplicateRecord *record = &c->record[row];
            record->replicate = replicate[row];
            int start = -1, periods = 0;
            for (int t = 0; t < df.r; t++)
                if ((int)AT(index, t, 0) == replicate[row]) { if (start < 0) start = t; periods++; }
            assert(start >= 0 && periods >= US_PERIODS);
            Mat y = mat_new(K, periods);
            for (int k = 0; k < K; k++) {
                Mat column = df_col_numeric(&df, abm_system_column_names()[k]);
                for (int t = 0; t < periods; t++) AT(y, k, t) = AT(column, start + t, 0);
            }
            block_moments(y, periods, record->moment[LENGTH_FULL]);
            block_moments(y, US_PERIODS, record->moment[LENGTH_US]);
            mat_free(y);

            char fit_path[512];
            snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, c->name,
                     replicate[row], SPEC_LABEL);
            record->nu = qvarma_load_params(&params, fit_path) ? (double)params.nu : NAN;
        }
        df_free(&df);
    }

    double *values = (double*)malloc((size_t)c->n_replicates * sizeof(double));
    for (int i = 0; i < c->n_replicates; i++) values[i] = c->record[i].nu;
    c->median_nu = median_of(values, c->n_replicates);
    for (int l = 0; l < N_LENGTHS; l++)
        for (int s = 0; s < N_SERIES; s++)
            for (int m = 0; m < N_MOMENTS; m++) {
                for (int i = 0; i < c->n_replicates; i++) values[i] = c->record[i].moment[l][s][m];
                c->median_moment[l][s][m] = median_of(values, c->n_replicates);
            }
    free(values);
    qvarma_params_free(&params);
    free(replicate);
}

int main(void) {
    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_tail_origin: cannot open the report path for writing");
    check_moment_functions(report);

    Mat us_y = build_us_block();
    double us_moment[N_SERIES][N_MOMENTS];
    block_moments(us_y, US_PERIODS, us_moment);
    mat_free(us_y);

    int n_configurations;
    char **names = list_configurations(&n_configurations);
    assert(n_configurations > 0);
    Configuration *configuration = (Configuration*)calloc((size_t)n_configurations, sizeof(Configuration));
    for (int i = 0; i < n_configurations; i++) snprintf(configuration[i].name, sizeof configuration[i].name, "%s", names[i]);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_configurations; i++) read_configuration(&configuration[i]);

    long n_total = 0;
    for (int i = 0; i < n_configurations; i++) n_total += configuration[i].n_replicates;
    int missing_nu = 0;
    for (int i = 0; i < n_configurations; i++)
        for (int j = 0; j < configuration[i].n_replicates; j++) missing_nu += MISNAN(configuration[i].record[j].nu);

    fprintf(report, "data: %d configurations, %ld replicates; US %d quarters, 1973Q2 to 2019Q4\n", n_configurations,
            n_total, US_PERIODS);
    fprintf(report, "fitted nu read for %ld replicates, missing for %d\n\n", n_total - missing_nu, missing_nu);

    /* Part 1. */
    fprintf(report, "PART 1. Raw series, no model\n");
    fprintf(report, "Simulated columns: quantiles over the %d configurations of each configuration's median over its replicates.\n",
            n_configurations);
    fprintf(report, "'reps >= US' is the share of all %ld replicates at or above the US value; 'configs >= US' the number of\n", n_total);
    fprintf(report, "configurations whose median is at or above it. For skewness both use the absolute value.\n");
    double *across = (double*)malloc((size_t)n_configurations * sizeof(double));
    for (int l = 0; l < N_LENGTHS; l++) {
        fprintf(report, "\nsimulated moments over %s\n", length_name[l]);
        fprintf(report, "  %-22s %-16s %10s %10s %10s %10s %11s %13s\n", "series", "moment", "US", "sim p05", "sim p50",
                "sim p95", "reps >= US", "configs >= US");
        for (int s = 0; s < N_SERIES; s++)
            for (int m = 0; m < N_MOMENTS; m++) {
                int absolute = m == MOMENT_SKEWNESS;
                double us_value = absolute ? fabs(us_moment[s][m]) : us_moment[s][m];
                long reached = 0;
                int configurations_reached = 0;
                for (int i = 0; i < n_configurations; i++) {
                    across[i] = configuration[i].median_moment[l][s][m];
                    double median = absolute ? fabs(across[i]) : across[i];
                    configurations_reached += median >= us_value;
                    for (int j = 0; j < configuration[i].n_replicates; j++) {
                        double value = configuration[i].record[j].moment[l][s][m];
                        reached += (absolute ? fabs(value) : value) >= us_value;
                    }
                }
                fprintf(report, "  %-22s %-16s %10.4g %10.4g %10.4g %10.4g %10.4f%% %13d\n", s == 0 || m == 0 ? series_name[s] : "",
                        moment_name[m], us_moment[s][m], quantile_of(across, n_configurations, 0.05),
                        quantile_of(across, n_configurations, 0.5), quantile_of(across, n_configurations, 0.95),
                        100.0 * (double)reached / (double)n_total, configurations_reached);
            }
    }

    /* Part 2. */
    int *order = (int*)malloc((size_t)n_configurations * sizeof(int));
    for (int i = 0; i < n_configurations; i++) order[i] = i;
    for (int i = 1; i < n_configurations; i++) {
        int current = order[i], j = i;
        while (j > 0 && configuration[order[j - 1]].median_nu > configuration[current].median_nu) { order[j] = order[j - 1]; j--; }
        order[j] = current;
    }

    fprintf(report, "\nPART 2. Fitted nu against the raw moments\n");
    fprintf(report, "Rank correlation (Spearman) between nu and each series' excess kurtosis, 400 periods.\n");
    fprintf(report, "Across configurations: median nu against median excess kurtosis. Across replicates: each fit's nu\n");
    fprintf(report, "against its own series. Thin simulated tails fitted correctly predict a negative correlation.\n");
    double *nu_across = (double*)malloc((size_t)n_configurations * sizeof(double));
    double *nu_all = (double*)malloc((size_t)n_total * sizeof(double));
    double *kurtosis_all = (double*)malloc((size_t)n_total * sizeof(double));
    for (int i = 0; i < n_configurations; i++) nu_across[i] = log(configuration[i].median_nu);
    fprintf(report, "  %-22s %16s %16s\n", "series", "configurations", "replicates");
    for (int s = 0; s < N_SERIES; s++) {
        for (int i = 0; i < n_configurations; i++) across[i] = configuration[i].median_moment[LENGTH_FULL][s][MOMENT_EXCESS_KURTOSIS];
        long at = 0;
        for (int i = 0; i < n_configurations; i++)
            for (int j = 0; j < configuration[i].n_replicates; j++) {
                nu_all[at] = configuration[i].record[j].nu;
                kurtosis_all[at] = configuration[i].record[j].moment[LENGTH_FULL][s][MOMENT_EXCESS_KURTOSIS];
                at++;
            }
        fprintf(report, "  %-22s %16.4f %16.4f\n", series_name[s], spearman_of(nu_across, across, n_configurations),
                spearman_of(nu_all, kurtosis_all, (int)n_total));
    }

    fprintf(report, "\nconfigurations ranked by median fitted nu: the %d lowest and the %d highest, with the median excess\n",
            N_RANKED_SHOWN, N_RANKED_SHOWN);
    fprintf(report, "kurtosis of each series over 400 periods. US row: nu of the US fit, 7.20, and the US excess kurtosis.\n");
    fprintf(report, "  %-5s %-10s %10s", "rank", "config", "median nu");
    for (int s = 0; s < N_SERIES; s++) fprintf(report, " %12.12s", series_name[s]);
    fprintf(report, "\n  %-5s %-10s %10.2f", "", "US", 7.20);
    for (int s = 0; s < N_SERIES; s++) fprintf(report, " %12.3f", us_moment[s][MOMENT_EXCESS_KURTOSIS]);
    fprintf(report, "\n");
    for (int k = 0; k < n_configurations; k++) {
        if (k >= N_RANKED_SHOWN && k < n_configurations - N_RANKED_SHOWN) {
            if (k == N_RANKED_SHOWN) fprintf(report, "  ...\n");
            continue;
        }
        const Configuration *c = &configuration[order[k]];
        fprintf(report, "  %-5d %-10s %10.1f", k + 1, c->name, c->median_nu);
        for (int s = 0; s < N_SERIES; s++) fprintf(report, " %12.3f", c->median_moment[LENGTH_FULL][s][MOMENT_EXCESS_KURTOSIS]);
        fprintf(report, "\n");
    }
    fclose(report);

    FILE *table = fopen(CONFIGURATIONS_PATH, "w");
    assert(table && "abm_system_tail_origin: cannot open the configuration table for writing");
    fprintf(table, "configuration,nu_rank,median_nu,n_replicates");
    for (int l = 0; l < N_LENGTHS; l++)
        for (int s = 0; s < N_SERIES; s++)
            for (int m = 0; m < N_MOMENTS; m++) fprintf(table, ",%s_%s_%s", series_column[s], moment_column[m], length_column[l]);
    fprintf(table, "\n");
    int *rank = (int*)malloc((size_t)n_configurations * sizeof(int));
    for (int k = 0; k < n_configurations; k++) rank[order[k]] = k + 1;
    for (int i = 0; i < n_configurations; i++) {
        fprintf(table, "%s,%d,%.10g,%d", configuration[i].name, rank[i], configuration[i].median_nu, configuration[i].n_replicates);
        for (int l = 0; l < N_LENGTHS; l++)
            for (int s = 0; s < N_SERIES; s++)
                for (int m = 0; m < N_MOMENTS; m++) fprintf(table, ",%.8g", configuration[i].median_moment[l][s][m]);
        fprintf(table, "\n");
    }
    fclose(table);

    free(rank); free(order); free(across); free(nu_across); free(nu_all); free(kurtosis_all);
    for (int i = 0; i < n_configurations; i++) { free(configuration[i].record); free(names[i]); }
    free(configuration); free(names);
    return 0;
}
