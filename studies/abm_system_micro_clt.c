/*
Whether the cross-section of firms behind GDP is one a central limit theorem
acts on.

GDP_r is Q1.Sum() * dim_mach + Q2.Sum(), a sum over 20 K-firms and 200
C-firms. applications/abm_system_micro_simulate.c writes those per-firm
quantities for 5 configurations of the design at seeds 1 to 5, which are
replicates 0 to 4 of the same configurations in dataset/abm_system. The two
firm types are kept apart throughout: they differ in number and in the weight
they carry into GDP.

A central limit theorem needs two things of such a sum, and each is measured
here rather than assumed.

    no firm dominates    with shares s_j = q_j / sum_k q_k, the concentration
                         H = sum_j s_j^2 is 1/N when every firm is the same
                         size and 1 when one firm is everything. H * N is
                         reported, so 1 is the even case and N the degenerate
                         one. This is the Lindeberg condition in the form the
                         data can show it.
    the terms are not    firms in this model share a wage, an interest rate and
    moving together      the same demand, so they need not be independent. The
                         average correlation over the N (N-1) / 2 pairs of
                         firm output changes is reported, with the effective
                         count N / (1 + (N - 1) rho) it implies: the number of
                         independent firms whose sum would have the variance
                         this correlated sum has.

What the cross-section would give if it were independent draws is then
arithmetic. A sum of N independent draws from a distribution with skewness
gamma and excess kurtosis kappa has skewness gamma / sqrt(N) and excess
kurtosis kappa / N, so the shape of the cross-section at one date fixes how far
from Gaussian the sum of N such firms can be. Both are reported at N and at the
effective count.

Measured on the cross-section at each period: skewness, excess kurtosis,
concentration, and Anderson-Darling against a Gaussian with the null simulated
at that N from 100,000 Gaussian samples rather than read off an asymptotic
approximation.

Part 1 does that on the cross-section of output levels at a date. That is the
cross-section GDP itself is a sum of, but it is not the one GDP growth is
driven by: a firm's level barely moves from one quarter to the next, so the
spread of levels is mostly firm size and says little about what moves the
aggregate. Part 2 uses the cross-section of output changes, which is the sum
GDP's change is.

Part 2 also settles what Part 1 cannot. Writing x_k for a firm's contribution
to GDP, K-firm output times dim_mach and C-firm output as it stands, GDP is
sum_k x_k over the 220 firms and its change is sum_k of the changes. So the
aggregate change is exactly 220 times the cross-sectional mean change, and a
central limit theorem over the cross-section is a statement about that mean.
If the firms are independent at a date with cross-sectional standard deviation
sigma_t, the aggregate change has standard deviation sqrt(220) sigma_t at that
date, and dividing by it should leave something Gaussian. So GDP growth is
measured three ways:

    as it stands             excess kurtosis of the growth series itself,
                             which is the number dataset/abm_system reports
    standardised             the same after dividing each quarter's growth by
                             the cross-sectional standard deviation of that
                             quarter, which is what a central limit theorem
                             acting at each date separately would leave
    dispersion alone         3 (E[v^2] / E[v]^2 - 1) with v_t the squared
                             conditional standard deviation, the excess
                             kurtosis a Gaussian whose variance moves like this
                             one's would have on its own

If the standardised series is Gaussian and the dispersion term accounts for the
raw excess kurtosis, then the central limit theorem holds at every date and the
shape of the growth series is the variance moving over time, not the
aggregation failing.

dim_mach is not read from the model's parameter file. It is recovered as
(GDP - sum of C-firm output) / (sum of K-firm output) at each period and
required to be the same number throughout, which also checks that the per-firm
files and the aggregate file are the same run.

Periods 201 to 600 of each run, the window dataset/abm_system keeps.

Output, none of it printed:
    out/abm_system_micro_clt_report.txt
    out/abm_system_micro_clt_by_run.csv
    out/abm_system_micro_clt_growth.csv
*/

#include "applications/abm_system.h"
#include <et_al./stats.h>
#include <et_al./special.h>
#include <et_al./random/random.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define INPUT_DIR "dataset/abm_system_micro"
#define REPORT_PATH "out/abm_system_micro_clt_report.txt"
#define BY_RUN_PATH "out/abm_system_micro_clt_by_run.csv"
#define GROWTH_PATH "out/abm_system_micro_clt_growth.csv"

#define FIRST_PERIOD ABM_SYSTEM_BURN_IN
#define SEEDS 5
#define N_TYPES 2
static const char *type_name[N_TYPES] = {"K-firms", "C-firms"};
static const char *type_column[N_TYPES] = {"k_firms", "c_firms"};
static const char *type_file[N_TYPES] = {"Q1all", "Q2all"};

#define NULL_DRAWS 100000
#define NULL_SEED 20260923

static int ascending(const void *a, const void *b) {
    double left = *(const double*)a, right = *(const double*)b;
    return left < right ? -1 : left > right ? 1 : 0;
}

/* Skewness m3 / m2^(3/2), excess kurtosis m4 / m2^2 - 3 and the
   Anderson-Darling statistic of one sample, with m_k the k-th central sample
   moment divided by n. buffer holds at least n doubles. Returns 0 when the
   sample has no spread, where none of the three is defined. */
static int shape_of(const double *x, int n, double *buffer,
                    double *skewness, double *excess_kurtosis, double *anderson_darling) {
    if (n < 4) return 0;

    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;

    double m2 = 0, m3 = 0, m4 = 0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean, d2 = d * d;
        m2 += d2;
        m3 += d2 * d;
        m4 += d2 * d2;
    }
    m2 /= (double)n;
    m3 /= (double)n;
    m4 /= (double)n;
    if (!(m2 > 0)) return 0;

    *skewness = m3 / pow(m2, 1.5);
    *excess_kurtosis = m4 / (m2 * m2) - 3.0;

    double sd = sqrt(m2 * (double)n / (double)(n - 1));
    for (int i = 0; i < n; i++) buffer[i] = x[i];
    qsort(buffer, (size_t)n, sizeof(double), ascending);

    double sum = 0;
    for (int i = 0; i < n; i++) {
        double lower = special_norm_cdf((buffer[i] - mean) / sd);
        double upper = special_norm_cdf((buffer[n - 1 - i] - mean) / sd);
        if (lower < 1e-300) lower = 1e-300;
        if (upper > 1.0 - 1e-16) upper = 1.0 - 1e-16;
        sum += (double)(2 * i + 1) * (log(lower) + log1p(-upper));
    }
    *anderson_darling = -(double)n - sum / (double)n;
    return 1;
}

/* The simulated null of the Anderson-Darling statistic at one sample size,
   ascending. */
static double *null_table_new(int n, uint64_t stream) {
    double *table = (double*)malloc((size_t)NULL_DRAWS * sizeof(double));
    assert(table && "abm_system_micro_clt: out of memory");

    Rng rng = rng_new(NULL_SEED, stream);
    double *sample = (double*)malloc((size_t)n * sizeof(double));
    double *buffer = (double*)malloc((size_t)n * sizeof(double));
    assert(sample && buffer && "abm_system_micro_clt: out of memory");

    for (int draw = 0; draw < NULL_DRAWS; draw++) {
        for (int i = 0; i < n; i++) sample[i] = rng_normal(&rng);
        double skewness, excess_kurtosis, anderson_darling;
        int ok = shape_of(sample, n, buffer, &skewness, &excess_kurtosis, &anderson_darling);
        assert(ok && "abm_system_micro_clt: a Gaussian sample had no spread");
        table[draw] = anderson_darling;
    }
    qsort(table, NULL_DRAWS, sizeof(double), ascending);

    free(sample);
    free(buffer);
    return table;
}

static double p_value_of(const double *sorted, double observed) {
    int low = 0, high = NULL_DRAWS;
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (sorted[mid] < observed) low = mid + 1;
        else high = mid;
    }
    return (double)(NULL_DRAWS - low) / (double)NULL_DRAWS;
}

/* One run's per-firm file as rows x columns of doubles, row-major. Caller must
   free. Returns NULL when the file is not there. */
static double *read_panel(const char *path, int *rows, int *columns) {
    FILE *f = fopen(path, "r");
    if (!f) return NULL;

    char line[1 << 16];
    if (!fgets(line, sizeof line, f)) { fclose(f); return NULL; }
    int n_columns = 0;
    for (const char *p = line; *p; ) {
        while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
        if (!*p) break;
        n_columns++;
        while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
    }
    assert(n_columns > 0 && "abm_system_micro_clt: a per-firm file has no columns");
    rewind(f);

    int capacity = n_columns * 1024, filled = 0;
    double *value = (double*)malloc((size_t)capacity * sizeof(double));
    assert(value && "abm_system_micro_clt: out of memory");
    double here;
    while (fscanf(f, "%lf", &here) == 1) {
        if (filled == capacity) {
            capacity *= 2;
            value = (double*)realloc(value, (size_t)capacity * sizeof(double));
            assert(value && "abm_system_micro_clt: out of memory");
        }
        value[filled++] = here;
    }
    fclose(f);

    assert(filled % n_columns == 0 && "abm_system_micro_clt: a per-firm file is ragged");
    *rows = filled / n_columns;
    *columns = n_columns;
    return value;
}

/* What one run of one firm type comes to. */
typedef struct {
    int configuration;
    int seed;
    int firms;
    int periods;
    double skewness;
    double excess_kurtosis;
    double concentration;         /* H * N, 1 when every firm is the same size */
    double anderson_darling;
    double rejected_05;
    double correlation;           /* average over the pairs of firm output changes */
} Run;

/* The average correlation over the N (N-1) / 2 pairs of columns of first
   differences, taken over the rows given. */
static double average_pair_correlation(const double *panel, int columns,
                                       int first, int last) {
    int n = last - first;
    if (n < 3 || columns < 2) return NAN;

    double *mean = (double*)calloc((size_t)columns, sizeof(double));
    double *sd = (double*)calloc((size_t)columns, sizeof(double));
    double *change = (double*)malloc((size_t)columns * (size_t)n * sizeof(double));
    assert(mean && sd && change && "abm_system_micro_clt: out of memory");

    for (int j = 0; j < columns; j++) {
        for (int i = 0; i < n; i++) {
            int t = first + i;
            change[j * n + i] = panel[t * columns + j] - panel[(t - 1) * columns + j];
            mean[j] += change[j * n + i];
        }
        mean[j] /= (double)n;
        for (int i = 0; i < n; i++) {
            double d = change[j * n + i] - mean[j];
            sd[j] += d * d;
        }
        sd[j] = sqrt(sd[j]);
    }

    double total = 0;
    long pairs = 0;
    for (int a = 0; a < columns; a++) {
        if (!(sd[a] > 0)) continue;
        for (int b = a + 1; b < columns; b++) {
            if (!(sd[b] > 0)) continue;
            double cross = 0;
            for (int i = 0; i < n; i++)
                cross += (change[a * n + i] - mean[a]) * (change[b * n + i] - mean[b]);
            total += cross / (sd[a] * sd[b]);
            pairs++;
        }
    }

    free(mean); free(sd); free(change);
    return pairs ? total / (double)pairs : NAN;
}

/* What one run's GDP growth comes to once the cross-section of changes is
   allowed for. */
typedef struct {
    int configuration;
    int seed;
    int periods;
    double dim_mach;
    double growth_excess_kurtosis;        /* the series as dataset/abm_system reports it */
    double standardised_excess_kurtosis;  /* after dividing by each quarter's own dispersion */
    double dispersion_excess_kurtosis;    /* what the moving variance gives on its own */
    double growth_skewness;
    double standardised_skewness;
    double change_excess_kurtosis[N_TYPES];  /* the cross-section of changes, by type */
    double change_concentration[N_TYPES];
} Growth;

static double excess_kurtosis_of(const double *x, int n) {
    if (n < 4) return NAN;
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;
    double m2 = 0, m4 = 0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean, d2 = d * d;
        m2 += d2;
        m4 += d2 * d2;
    }
    m2 /= (double)n;
    m4 /= (double)n;
    return m2 > 0 ? m4 / (m2 * m2) - 3.0 : NAN;
}

static double skewness_of(const double *x, int n) {
    if (n < 4) return NAN;
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;
    double m2 = 0, m3 = 0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean, d2 = d * d;
        m2 += d2;
        m3 += d2 * d;
    }
    m2 /= (double)n;
    m3 /= (double)n;
    return m2 > 0 ? m3 / pow(m2, 1.5) : NAN;
}

/* Cross-sectional standard deviation, and the share concentration of the
   absolute deviations, of one column of values. */
static double cross_sd(const double *x, int n) {
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;
    double m2 = 0;
    for (int i = 0; i < n; i++) { double d = x[i] - mean; m2 += d * d; }
    return sqrt(m2 / (double)n);
}

static int *list_configurations(int *count) {
    DIR *handle = opendir(INPUT_DIR);
    assert(handle &&
           "abm_system_micro_clt: cannot open dataset/abm_system_micro - run make app-abm_system_micro_simulate");

    int *cop = NULL, n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int here;
        if (sscanf(entry->d_name, "cop_%d", &here) != 1) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 64;
            cop = (int*)realloc(cop, (size_t)cap * sizeof(int));
            assert(cop && "abm_system_micro_clt: out of memory");
        }
        cop[n++] = here;
    }
    closedir(handle);

    for (int i = 1; i < n; i++) {
        int current = cop[i], j = i;
        while (j > 0 && cop[j - 1] > current) { cop[j] = cop[j - 1]; j--; }
        cop[j] = current;
    }
    *count = n;
    return cop;
}

int main(void) {
    int n_configurations;
    int *cop = list_configurations(&n_configurations);
    assert(n_configurations > 0 && "abm_system_micro_clt: no configurations stored");

    Run *run = (Run*)calloc((size_t)n_configurations * SEEDS * N_TYPES, sizeof(Run));
    assert(run && "abm_system_micro_clt: out of memory");
    int n_runs = 0;

    /* One null per firm count, built once the counts are known from the first
       file of each type. */
    double *null_table[N_TYPES] = {NULL, NULL};
    int null_firms[N_TYPES] = {0, 0};

    for (int type = 0; type < N_TYPES; type++) {
        for (int c = 0; c < n_configurations; c++) {
            for (int seed = 1; seed <= SEEDS; seed++) {
                char path[640];
                snprintf(path, sizeof path, "%s/cop_%04d/%s_%d.txt",
                         INPUT_DIR, cop[c], type_file[type], seed);
                int rows, columns;
                double *panel = read_panel(path, &rows, &columns);
                if (!panel) continue;
                assert(rows > FIRST_PERIOD && "abm_system_micro_clt: a run is shorter than the burn-in");

                if (!null_table[type]) {
                    null_firms[type] = columns;
                    null_table[type] = null_table_new(columns, (uint64_t)type + 1);
                }
                assert(columns == null_firms[type] &&
                       "abm_system_micro_clt: two runs of one type have different firm counts");

                double *buffer = (double*)malloc((size_t)columns * sizeof(double));
                double *cross = (double*)malloc((size_t)columns * sizeof(double));
                assert(buffer && cross && "abm_system_micro_clt: out of memory");

                Run here = {0};
                here.configuration = cop[c];
                here.seed = seed;
                here.firms = columns;

                int counted = 0;
                for (int t = FIRST_PERIOD; t < rows; t++) {
                    double total = 0;
                    for (int j = 0; j < columns; j++) {
                        cross[j] = panel[t * columns + j];
                        total += cross[j];
                    }

                    double skewness, excess_kurtosis, anderson_darling;
                    if (!shape_of(cross, columns, buffer, &skewness, &excess_kurtosis,
                                  &anderson_darling)) continue;

                    double concentration = 0;
                    if (total > 0) {
                        for (int j = 0; j < columns; j++) {
                            double share = cross[j] / total;
                            concentration += share * share;
                        }
                        concentration *= (double)columns;
                    } else concentration = NAN;

                    here.skewness += skewness;
                    here.excess_kurtosis += excess_kurtosis;
                    here.anderson_darling += anderson_darling;
                    here.concentration += concentration;
                    here.rejected_05 += p_value_of(null_table[type], anderson_darling) < 0.05;
                    counted++;
                }
                assert(counted > 0 && "abm_system_micro_clt: a run had no usable period");

                here.periods = counted;
                here.skewness /= counted;
                here.excess_kurtosis /= counted;
                here.anderson_darling /= counted;
                here.concentration /= counted;
                here.rejected_05 /= counted;
                here.correlation = average_pair_correlation(panel, columns,
                                                            FIRST_PERIOD, rows);
                run[n_runs++] = here;

                free(buffer); free(cross); free(panel);
            }
        }
    }
    assert(n_runs > 0 && "abm_system_micro_clt: nothing to read");

    /* Part 2. The two per-firm files of a run and its aggregate file together,
       so GDP growth can be set against the cross-section of changes that makes
       it. */
    Growth *growth = (Growth*)calloc((size_t)n_configurations * SEEDS, sizeof(Growth));
    assert(growth && "abm_system_micro_clt: out of memory");
    int n_growth = 0;

    for (int c = 0; c < n_configurations; c++) {
        for (int seed = 1; seed <= SEEDS; seed++) {
            char path[640];
            int rows1, columns1, rows2, columns2, rows_res, columns_res;

            snprintf(path, sizeof path, "%s/cop_%04d/Q1all_%d.txt", INPUT_DIR, cop[c], seed);
            double *q1 = read_panel(path, &rows1, &columns1);
            snprintf(path, sizeof path, "%s/cop_%04d/Q2all_%d.txt", INPUT_DIR, cop[c], seed);
            double *q2 = read_panel(path, &rows2, &columns2);
            snprintf(path, sizeof path, "%s/cop_%04d/results_%d.txt", INPUT_DIR, cop[c], seed);
            double *results = read_panel(path, &rows_res, &columns_res);
            if (!q1 || !q2 || !results) { free(q1); free(q2); free(results); continue; }
            assert(rows1 == rows2 && rows1 == rows_res &&
                   "abm_system_micro_clt: the files of one run disagree on length");

            /* dim_mach from the identity itself, which is also the check that
               these files are the same run. Column 2 of the aggregate file is
               real GDP. */
            double dim_mach = 0;
            for (int t = 0; t < rows1; t++) {
                double total1 = 0, total2 = 0;
                for (int j = 0; j < columns1; j++) total1 += q1[t * columns1 + j];
                for (int j = 0; j < columns2; j++) total2 += q2[t * columns2 + j];
                double gdp = results[t * columns_res + 1];
                if (!(total1 > 0)) continue;
                double here = (gdp - total2) / total1;
                if (dim_mach == 0) dim_mach = here;
                assert(fabs(here - dim_mach) < 1e-6 * dim_mach &&
                       "abm_system_micro_clt: the per-firm files do not add up to the aggregate file");
            }
            assert(dim_mach > 0 && "abm_system_micro_clt: could not recover dim_mach");

            int firms = columns1 + columns2;
            int n = rows1 - FIRST_PERIOD;
            double *contribution_change = (double*)malloc((size_t)firms * sizeof(double));
            double *g = (double*)malloc((size_t)n * sizeof(double));
            double *z = (double*)malloc((size_t)n * sizeof(double));
            double *v = (double*)malloc((size_t)n * sizeof(double));
            assert(contribution_change && g && z && v && "abm_system_micro_clt: out of memory");

            Growth here = {0};
            here.configuration = cop[c];
            here.seed = seed;
            here.dim_mach = dim_mach;

            double change_kurtosis[N_TYPES] = {0, 0};
            double change_concentration[N_TYPES] = {0, 0};
            int counted = 0;

            for (int i = 0; i < n; i++) {
                int t = FIRST_PERIOD + i;
                int at = 0;
                for (int j = 0; j < columns1; j++)
                    contribution_change[at++] = dim_mach *
                        (q1[t * columns1 + j] - q1[(t - 1) * columns1 + j]);
                for (int j = 0; j < columns2; j++)
                    contribution_change[at++] =
                        q2[t * columns2 + j] - q2[(t - 1) * columns2 + j];

                double gdp = results[t * columns_res + 1];
                double previous = results[(t - 1) * columns_res + 1];
                g[i] = 100.0 * (log(gdp) - log(previous));

                /* What a central limit theorem acting at this date alone would
                   give the aggregate change: sqrt(firms) times the
                   cross-sectional spread, carried into the growth rate by the
                   level it is divided by. */
                double sigma = cross_sd(contribution_change, firms);
                double conditional = 100.0 * sqrt((double)firms) * sigma / previous;
                z[i] = conditional > 0 ? g[i] / conditional : NAN;
                v[i] = conditional * conditional;

                change_kurtosis[0] += excess_kurtosis_of(contribution_change, columns1);
                change_kurtosis[1] += excess_kurtosis_of(contribution_change + columns1, columns2);

                for (int type = 0; type < N_TYPES; type++) {
                    const double *piece = type == 0 ? contribution_change
                                                    : contribution_change + columns1;
                    int count = type == 0 ? columns1 : columns2;
                    double total = 0;
                    for (int j = 0; j < count; j++) total += fabs(piece[j]);
                    double concentration = 0;
                    if (total > 0) {
                        for (int j = 0; j < count; j++) {
                            double share = fabs(piece[j]) / total;
                            concentration += share * share;
                        }
                        concentration *= (double)count;
                    }
                    change_concentration[type] += concentration;
                }
                counted++;
            }

            here.periods = counted;
            here.growth_excess_kurtosis = excess_kurtosis_of(g, n);
            here.growth_skewness = skewness_of(g, n);
            here.standardised_excess_kurtosis = excess_kurtosis_of(z, n);
            here.standardised_skewness = skewness_of(z, n);

            /* A Gaussian whose variance moves over time has excess kurtosis
               3 (E[v^2] / E[v]^2 - 1) once the variance path is taken as
               given. */
            double mean_v = 0, mean_v2 = 0;
            for (int i = 0; i < n; i++) { mean_v += v[i]; mean_v2 += v[i] * v[i]; }
            mean_v /= (double)n;
            mean_v2 /= (double)n;
            here.dispersion_excess_kurtosis = mean_v > 0 ? 3.0 * (mean_v2 / (mean_v * mean_v) - 1.0) : NAN;

            for (int type = 0; type < N_TYPES; type++) {
                here.change_excess_kurtosis[type] = change_kurtosis[type] / counted;
                here.change_concentration[type] = change_concentration[type] / counted;
            }
            growth[n_growth++] = here;

            free(contribution_change); free(g); free(z); free(v);
            free(q1); free(q2); free(results);
        }
    }

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_micro_clt: cannot open the report path");
    fprintf(report, "Whether the cross-section of firms behind GDP is one a central limit theorem\n"
                    "acts on.\n\n");
    fprintf(report, "%d configurations of the design, seeds 1 to %d, periods %d to the end of\n"
                    "each run. The two firm types are kept apart.\n",
            n_configurations, SEEDS, FIRST_PERIOD + 1);
    fprintf(report, "Anderson-Darling null simulated from %d Gaussian samples at each firm\n"
                    "count, seed %d.\n\n", NULL_DRAWS, NULL_SEED);

    FILE *by_run = fopen(BY_RUN_PATH, "w");
    assert(by_run && "abm_system_micro_clt: cannot open the run table path");
    fprintf(by_run, "firm_type,configuration,seed,firms,periods,cross_skewness,"
                    "cross_excess_kurtosis,concentration_times_n,mean_anderson_darling,"
                    "share_ad_below_05,mean_pair_correlation,effective_firms,"
                    "implied_sum_skewness,implied_sum_excess_kurtosis,"
                    "implied_sum_excess_kurtosis_effective\n");
    for (int i = 0; i < n_runs; i++) {
        const Run *r = &run[i];
        int type = r->firms == null_firms[0] ? 0 : 1;
        double effective = (double)r->firms / (1.0 + ((double)r->firms - 1.0) * r->correlation);
        fprintf(by_run, "%s,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                type_column[type], r->configuration, r->seed, r->firms, r->periods,
                r->skewness, r->excess_kurtosis, r->concentration, r->anderson_darling,
                r->rejected_05, r->correlation, effective,
                r->skewness / sqrt((double)r->firms), r->excess_kurtosis / (double)r->firms,
                effective > 0 ? r->excess_kurtosis / effective : NAN);
    }
    fclose(by_run);

    fprintf(report, "PART 1. The cross-section of output levels at a date, averaged over the\n"
                    "periods of a run and then\n"
                    "over the runs. Concentration is H * N: 1 when every firm is the same size,\n"
                    "N when one firm is everything. AD<.05 is 0.05 under a Gaussian\n"
                    "cross-section.\n\n");
    fprintf(report, "  firm type      N   runs   skewness  exc kurt  concentr   AD<.05\n");
    for (int type = 0; type < N_TYPES; type++) {
        double skewness = 0, excess_kurtosis = 0, concentration = 0, rejected = 0;
        int counted = 0;
        for (int i = 0; i < n_runs; i++) {
            if (run[i].firms != null_firms[type]) continue;
            skewness += run[i].skewness;
            excess_kurtosis += run[i].excess_kurtosis;
            concentration += run[i].concentration;
            rejected += run[i].rejected_05;
            counted++;
        }
        if (!counted) continue;
        fprintf(report, "  %-10s %5d %6d %10.4f %9.4f %9.3f %8.4f\n",
                type_name[type], null_firms[type], counted, skewness / counted,
                excess_kurtosis / counted, concentration / counted, rejected / counted);
    }

    fprintf(report, "\nWhether the terms move together, and what the sum inherits. The average\n"
                    "correlation is over the pairs of firm output changes within a run. The\n"
                    "effective count is N / (1 + (N - 1) rho). A sum of m independent draws\n"
                    "from the cross-section would have skewness gamma / sqrt(m) and excess\n"
                    "kurtosis kappa / m.\n\n");
    fprintf(report, "  firm type      N   mean rho  effective N   sum skew at N  sum kurt at N"
                    "  sum kurt at eff\n");
    for (int type = 0; type < N_TYPES; type++) {
        double correlation = 0, skewness = 0, excess_kurtosis = 0;
        int counted = 0;
        for (int i = 0; i < n_runs; i++) {
            if (run[i].firms != null_firms[type]) continue;
            correlation += run[i].correlation;
            skewness += run[i].skewness;
            excess_kurtosis += run[i].excess_kurtosis;
            counted++;
        }
        if (!counted) continue;
        correlation /= counted;
        skewness /= counted;
        excess_kurtosis /= counted;
        double n = (double)null_firms[type];
        double effective = n / (1.0 + (n - 1.0) * correlation);
        fprintf(report, "  %-10s %5d %10.4f %12.2f %15.4f %14.4f %16.4f\n",
                type_name[type], null_firms[type], correlation, effective,
                skewness / sqrt(n), excess_kurtosis / n,
                effective > 0 ? excess_kurtosis / effective : NAN);
    }

    fprintf(report, "\nPART 2. The cross-section of changes, which is what GDP's change is a sum\n"
                    "of. Averaged over the periods of a run and then over the runs.\n\n");
    fprintf(report, "  firm type      N   runs   exc kurt of change  concentr of |change|\n");
    for (int type = 0; type < N_TYPES; type++) {
        double excess_kurtosis = 0, concentration = 0;
        for (int i = 0; i < n_growth; i++) {
            excess_kurtosis += growth[i].change_excess_kurtosis[type];
            concentration += growth[i].change_concentration[type];
        }
        fprintf(report, "  %-10s %5d %6d %20.4f %21.3f\n",
                type_name[type], null_firms[type], n_growth,
                excess_kurtosis / n_growth, concentration / n_growth);
    }

    fprintf(report, "\nGDP growth three ways, over the %d quarters of each run, averaged over the\n"
                    "runs. 'as it stands' is the series dataset/abm_system reports.\n"
                    "'standardised' divides each quarter by sqrt(N) times that quarter's own\n"
                    "cross-sectional spread, which is what a central limit theorem acting at\n"
                    "each date separately would leave. 'dispersion alone' is the excess\n"
                    "kurtosis a Gaussian whose variance moved like this one's would have.\n\n",
            n_growth ? growth[0].periods : 0);
    fprintf(report, "  quantity                    skewness   exc kurt\n");
    {
        double raw_k = 0, raw_s = 0, std_k = 0, std_s = 0, dispersion = 0;
        for (int i = 0; i < n_growth; i++) {
            raw_k += growth[i].growth_excess_kurtosis;
            raw_s += growth[i].growth_skewness;
            std_k += growth[i].standardised_excess_kurtosis;
            std_s += growth[i].standardised_skewness;
            dispersion += growth[i].dispersion_excess_kurtosis;
        }
        double inverse = n_growth ? 1.0 / (double)n_growth : NAN;
        fprintf(report, "  GDP growth as it stands   %10.4f %10.4f\n", raw_s * inverse, raw_k * inverse);
        fprintf(report, "  standardised              %10.4f %10.4f\n", std_s * inverse, std_k * inverse);
        fprintf(report, "  dispersion alone                     %10.4f\n", dispersion * inverse);
    }

    FILE *by_growth = fopen(GROWTH_PATH, "w");
    assert(by_growth && "abm_system_micro_clt: cannot open the growth table path");
    fprintf(by_growth, "configuration,seed,periods,dim_mach,growth_skewness,"
                       "growth_excess_kurtosis,standardised_skewness,"
                       "standardised_excess_kurtosis,dispersion_excess_kurtosis,"
                       "k_firm_change_excess_kurtosis,c_firm_change_excess_kurtosis,"
                       "k_firm_change_concentration,c_firm_change_concentration\n");
    for (int i = 0; i < n_growth; i++) {
        const Growth *r = &growth[i];
        fprintf(by_growth, "%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                r->configuration, r->seed, r->periods, r->dim_mach,
                r->growth_skewness, r->growth_excess_kurtosis,
                r->standardised_skewness, r->standardised_excess_kurtosis,
                r->dispersion_excess_kurtosis,
                r->change_excess_kurtosis[0], r->change_excess_kurtosis[1],
                r->change_concentration[0], r->change_concentration[1]);
    }
    fclose(by_growth);

    fprintf(report, "\nEvery run on its own is in %s and %s.\n", BY_RUN_PATH, GROWTH_PATH);
    fclose(report);

    for (int type = 0; type < N_TYPES; type++) free(null_table[type]);
    free(run);
    free(growth);
    free(cop);
    return EXIT_SUCCESS;
}
