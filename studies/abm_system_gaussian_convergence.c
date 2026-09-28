/*
Whether the departure from a Gaussian in a run's own time series shrinks as the
sample grows.

Each replicate is one economy's history. This reads down that history and
nothing else: no sample here is built by putting different replicates side by
side. The system is not ergodic, so a sample taken across replicates at a fixed
date is a mixture over economies in different states and says nothing about a
central limit theorem.

For every replicate and every variable the first n quarters are taken, n = 100,
150, 200, 250, 300, 350, 400, and the same four quantities are measured at each
n. What the sweep over n shows is whether a departure holds up as the sample
grows or thins out with it.

The variables are the five the experiment stores: GDP growth, energy growth,
employment change, inflation and the interest rate. The first four are already
first differences. The interest rate is stored as a level, so its first
difference is measured as a sixth variable; it has 399 values rather than 400
and so reports up to n = 350.

What is measured, for every sample:

    skewness             m3 / m2^(3/2), with m_k the k-th central sample moment
                         divided by n
    excess kurtosis      m4 / m2^2 - 3, zero for a Gaussian
    Anderson-Darling     A2 = -n - (1/n) sum_i (2i-1) [ln z_i + ln(1 - z_(n+1-i))]
                         on z_i = Phi((x_(i) - xbar) / s), the sample ordered
                         and standardised by its own mean and standard
                         deviation. It weights the tails.
    Jarque-Bera          (n/6) (S^2 + (K - 3)^2 / 4), which reads the two
                         moments above and nothing else.

Two statistics rather than one because they fail in different places: a sample
can match the first four moments of a Gaussian and have the wrong shape in the
tails, and it can have the right shape everywhere but the tails.

Neither is read off its asymptotic distribution. At every n the null
distribution of both is simulated from Gaussian samples of that size, and an
observed statistic becomes a p-value by where it falls in that simulated null.
The mean of each statistic under that null is reported beside the observed one,
because sample excess kurtosis is biased downward and a Gaussian sample does
not give exactly zero.

Two things this does not control for. The quarters of one run are dependent
while the simulated null draws them independently, so a rejection can come from
the dependence rather than from a non-Gaussian shape. And a longer n reaches
later into a run, so the sample grows and the stretch of history it covers
changes together.

Under the null a p-value is uniform on 0 to 1, so the share of tests below 0.05
should be 0.05 and the share below 0.01 should be 0.01.

    ./bin/abm_system_gaussian_convergence [CONFIGURATIONS]

Reads dataset/abm_system only. Writes, none of it printed:
    out/abm_system_gaussian_convergence_report.txt
    out/abm_system_gaussian_convergence_by_sample_size.csv
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

#define K ABM_SYSTEM_K
#define INPUT_DIR "dataset/abm_system"
#define REPORT_PATH "out/abm_system_gaussian_convergence_report.txt"
#define BY_SAMPLE_SIZE_PATH "out/abm_system_gaussian_convergence_by_sample_size.csv"

#define PERIODS 400
#define REPLICATES 1000

#define N_SAMPLE_SIZES 7
static const int sample_size[N_SAMPLE_SIZES] = {100, 150, 200, 250, 300, 350, 400};

/* The five stored series and the first difference of the interest rate. */
#define N_SERIES (K + 1)
#define SERIES_INTEREST_CHANGE K

#define NULL_DRAWS 100000
#define NULL_SEED 20260923

static const char *series_name[N_SERIES] = {
    "GDP growth", "energy growth", "employment change", "inflation", "interest rate",
    "interest rate change"
};
static const char *series_column[N_SERIES] = {
    "gdp_growth", "energy_growth", "employment_change", "inflation", "interest_rate",
    "interest_rate_change"
};

typedef struct {
    double skewness;
    double excess_kurtosis;
    double anderson_darling;
    double jarque_bera;
} Departure;

static int ascending(const void *a, const void *b) {
    double left = *(const double*)a, right = *(const double*)b;
    return left < right ? -1 : left > right ? 1 : 0;
}

/* The four quantities of one sample. buffer holds at least n doubles and takes
   the sort Anderson-Darling needs. Returns 0 when the sample has no spread,
   where none of the four is defined. */
static int departure_of(const double *x, int n, double *buffer, Departure *out) {
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

    out->skewness = m3 / pow(m2, 1.5);
    out->excess_kurtosis = m4 / (m2 * m2) - 3.0;
    out->jarque_bera = (double)n / 6.0 *
                       (out->skewness * out->skewness +
                        out->excess_kurtosis * out->excess_kurtosis / 4.0);

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
    out->anderson_darling = -(double)n - sum / (double)n;
    return 1;
}

/* The simulated null of both statistics at one sample size, ascending, and
   what excess kurtosis and Anderson-Darling average to there. */
typedef struct {
    int n;
    double *anderson_darling;
    double *jarque_bera;
    double mean_excess_kurtosis;
    double mean_anderson_darling;
} NullTable;

static NullTable null_table_new(int n, uint64_t stream) {
    NullTable table;
    table.n = n;
    table.anderson_darling = (double*)malloc((size_t)NULL_DRAWS * sizeof(double));
    table.jarque_bera = (double*)malloc((size_t)NULL_DRAWS * sizeof(double));
    assert(table.anderson_darling && table.jarque_bera &&
           "abm_system_gaussian_convergence: out of memory");

    Rng rng = rng_new(NULL_SEED, stream);
    double *sample = (double*)malloc((size_t)n * sizeof(double));
    double *buffer = (double*)malloc((size_t)n * sizeof(double));
    assert(sample && buffer && "abm_system_gaussian_convergence: out of memory");

    table.mean_excess_kurtosis = 0;
    table.mean_anderson_darling = 0;
    int kept = 0;
    for (int draw = 0; draw < NULL_DRAWS; draw++) {
        for (int i = 0; i < n; i++) sample[i] = rng_normal(&rng);
        Departure d;
        if (!departure_of(sample, n, buffer, &d)) continue;
        table.anderson_darling[kept] = d.anderson_darling;
        table.jarque_bera[kept] = d.jarque_bera;
        table.mean_excess_kurtosis += d.excess_kurtosis;
        table.mean_anderson_darling += d.anderson_darling;
        kept++;
    }
    assert(kept == NULL_DRAWS && "abm_system_gaussian_convergence: a Gaussian sample had no spread");
    table.mean_excess_kurtosis /= (double)kept;
    table.mean_anderson_darling /= (double)kept;

    qsort(table.anderson_darling, NULL_DRAWS, sizeof(double), ascending);
    qsort(table.jarque_bera, NULL_DRAWS, sizeof(double), ascending);

    free(sample);
    free(buffer);
    return table;
}

static void null_table_free(NullTable *table) {
    free(table->anderson_darling);
    free(table->jarque_bera);
}

/* The share of the simulated null at or above an observed statistic, which is
   the p-value both take, since both reject in the upper tail. */
static double p_value_of(const double *sorted, double observed) {
    int low = 0, high = NULL_DRAWS;
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (sorted[mid] < observed) low = mid + 1;
        else high = mid;
    }
    return (double)(NULL_DRAWS - low) / (double)NULL_DRAWS;
}

/* Running totals over the replicates in one (variable, n) cell. */
typedef struct {
    long tests;
    double skewness;
    double excess_kurtosis;
    double anderson_darling;
    long anderson_darling_below_05;
    long anderson_darling_below_01;
    long jarque_bera_below_05;
    long jarque_bera_below_01;
} Cell;

static void cell_add(Cell *cell, const Departure *d, const NullTable *table) {
    cell->tests++;
    cell->skewness += d->skewness;
    cell->excess_kurtosis += d->excess_kurtosis;
    cell->anderson_darling += d->anderson_darling;
    double ad = p_value_of(table->anderson_darling, d->anderson_darling);
    double jb = p_value_of(table->jarque_bera, d->jarque_bera);
    cell->anderson_darling_below_05 += ad < 0.05;
    cell->anderson_darling_below_01 += ad < 0.01;
    cell->jarque_bera_below_05 += jb < 0.05;
    cell->jarque_bera_below_01 += jb < 0.01;
}

static void cell_merge(Cell *into, const Cell *from) {
    into->tests += from->tests;
    into->skewness += from->skewness;
    into->excess_kurtosis += from->excess_kurtosis;
    into->anderson_darling += from->anderson_darling;
    into->anderson_darling_below_05 += from->anderson_darling_below_05;
    into->anderson_darling_below_01 += from->anderson_darling_below_01;
    into->jarque_bera_below_05 += from->jarque_bera_below_05;
    into->jarque_bera_below_01 += from->jarque_bera_below_01;
}

static int *list_configurations(int *count) {
    DIR *handle = opendir(INPUT_DIR);
    assert(handle && "abm_system_gaussian_convergence: cannot open dataset/abm_system");

    int *cop = NULL, n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int here;
        if (sscanf(entry->d_name, "cop_%d", &here) != 1) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 1024;
            cop = (int*)realloc(cop, (size_t)cap * sizeof(int));
            assert(cop && "abm_system_gaussian_convergence: out of memory");
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

int main(int argc, char **argv) {
    int n_configurations;
    int *cop = list_configurations(&n_configurations);
    assert(n_configurations > 0 && "abm_system_gaussian_convergence: no configurations stored");
    if (argc > 1) {
        int wanted = atoi(argv[1]);
        assert(wanted >= 1 && "abm_system_gaussian_convergence: nothing to read");
        if (wanted < n_configurations) n_configurations = wanted;
    }

    const double null_started = omp_get_wtime();
    NullTable *null_table = (NullTable*)malloc((size_t)N_SAMPLE_SIZES * sizeof(NullTable));
    assert(null_table && "abm_system_gaussian_convergence: out of memory");
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N_SAMPLE_SIZES; i++)
        null_table[i] = null_table_new(sample_size[i], (uint64_t)i + 1);
    const double null_elapsed = omp_get_wtime() - null_started;

    Cell *cell = (Cell*)calloc((size_t)N_SERIES * N_SAMPLE_SIZES, sizeof(Cell));
    assert(cell && "abm_system_gaussian_convergence: out of memory");

    const double sweep_started = omp_get_wtime();
    long replicates_read = 0;

    #pragma omp parallel
    {
        Cell *mine = (Cell*)calloc((size_t)N_SERIES * N_SAMPLE_SIZES, sizeof(Cell));
        double *series = (double*)malloc((size_t)PERIODS * sizeof(double));
        double *buffer = (double*)malloc((size_t)PERIODS * sizeof(double));
        assert(mine && series && buffer && "abm_system_gaussian_convergence: out of memory");
        long my_replicates = 0;

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n_configurations; i++) {
            char dir[512];
            snprintf(dir, sizeof dir, "%s/cop_%04d", INPUT_DIR, cop[i]);

            for (int batch = 0; batch < REPLICATES / ABM_SYSTEM_BATCH; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);

                for (int b = 0; b < count; b++) {
                    Mat y = block[b];
                    if (y.c < PERIODS) { mat_free(y); continue; }
                    my_replicates++;

                    for (int s = 0; s < N_SERIES; s++) {
                        int available;
                        if (s < K) {
                            available = PERIODS;
                            for (int t = 0; t < PERIODS; t++) series[t] = (double)AT(y, s, t);
                        } else {
                            available = PERIODS - 1;
                            for (int t = 1; t < PERIODS; t++)
                                series[t - 1] = (double)AT(y, ROW_INTEREST_RATE, t) -
                                                (double)AT(y, ROW_INTEREST_RATE, t - 1);
                        }

                        for (int q = 0; q < N_SAMPLE_SIZES; q++) {
                            int n = sample_size[q];
                            if (n > available) continue;
                            Departure d;
                            if (!departure_of(series, n, buffer, &d)) continue;
                            cell_add(&mine[s * N_SAMPLE_SIZES + q], &d, &null_table[q]);
                        }
                    }
                    mat_free(y);
                }
            }
        }

        #pragma omp critical(merge)
        {
            for (int i = 0; i < N_SERIES * N_SAMPLE_SIZES; i++) cell_merge(&cell[i], &mine[i]);
            replicates_read += my_replicates;
        }

        free(mine); free(series); free(buffer);
    }

    const double sweep_elapsed = omp_get_wtime() - sweep_started;

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_gaussian_convergence: cannot open the report path");

    fprintf(report, "Whether the departure from a Gaussian in a run's own time series shrinks as\n"
                    "the sample grows.\n\n");
    fprintf(report, "%d configurations, %ld replicates, the first n of each run's 400 quarters.\n",
            n_configurations, replicates_read);
    fprintf(report, "Nothing is pooled across replicates: every sample is one run's own history.\n");
    fprintf(report, "Null distributions simulated from %d Gaussian samples at each n, seed %d,\n"
                    "%.1f minutes. Sweep %.1f minutes on %d threads.\n\n",
            NULL_DRAWS, NULL_SEED, null_elapsed / 60.0, sweep_elapsed / 60.0,
            omp_get_max_threads());
    fprintf(report, "'null' is what that statistic averages to on Gaussian samples of the same n.\n"
                    "The 0.05 and 0.01 columns read against 0.05 and 0.01.\n");

    FILE *table = fopen(BY_SAMPLE_SIZE_PATH, "w");
    assert(table && "abm_system_gaussian_convergence: cannot open the table path");
    fprintf(table, "series,sample_size,replicates,mean_skewness,mean_excess_kurtosis,"
                   "null_excess_kurtosis,mean_anderson_darling,null_anderson_darling,"
                   "share_ad_below_05,share_ad_below_01,share_jb_below_05,share_jb_below_01\n");

    fprintf(report, "\n  series                     n     tests   skewness  exc kurt      null"
                    "     AD    null   AD<.05  AD<.01   JB<.05  JB<.01\n");
    for (int s = 0; s < N_SERIES; s++) {
        int shown = 0;
        for (int q = 0; q < N_SAMPLE_SIZES; q++) {
            const Cell *c = &cell[s * N_SAMPLE_SIZES + q];
            if (!c->tests) continue;
            double inverse = 1.0 / (double)c->tests;
            fprintf(report, "  %-20s %5d %9ld %10.4f %9.4f %9.4f %6.2f %7.2f %8.4f %7.4f %8.4f %7.4f\n",
                    shown++ == 0 ? series_name[s] : "", sample_size[q], c->tests,
                    c->skewness * inverse, c->excess_kurtosis * inverse,
                    null_table[q].mean_excess_kurtosis,
                    c->anderson_darling * inverse, null_table[q].mean_anderson_darling,
                    (double)c->anderson_darling_below_05 * inverse,
                    (double)c->anderson_darling_below_01 * inverse,
                    (double)c->jarque_bera_below_05 * inverse,
                    (double)c->jarque_bera_below_01 * inverse);
            fprintf(table, "%s,%d,%ld,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    series_column[s], sample_size[q], c->tests,
                    c->skewness * inverse, c->excess_kurtosis * inverse,
                    null_table[q].mean_excess_kurtosis,
                    c->anderson_darling * inverse, null_table[q].mean_anderson_darling,
                    (double)c->anderson_darling_below_05 * inverse,
                    (double)c->anderson_darling_below_01 * inverse,
                    (double)c->jarque_bera_below_05 * inverse,
                    (double)c->jarque_bera_below_01 * inverse);
        }
    }

    fclose(table);
    fclose(report);

    for (int i = 0; i < N_SAMPLE_SIZES; i++) null_table_free(&null_table[i]);
    free(null_table);
    free(cell);
    free(cop);
    return EXIT_SUCCESS;
}
