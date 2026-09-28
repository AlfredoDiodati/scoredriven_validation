/*
Whether the macro series are stable with index below 2, or have a finite
variance.

studies/abm_system_gaussian_convergence.c settles that they are not Gaussian
and that the departure grows with the sample rather than shrinking. A growing
sample excess kurtosis is what an infinite fourth moment looks like. It is also
what a fixed leptokurtic distribution looks like, because sample excess
kurtosis is bounded above by roughly n and biased downward, so that measurement
does not separate the two.

Stable laws are indexed by alpha in (0, 2]. Alpha = 2 is the Gaussian. Every
alpha below 2 has P(|X| > x) ~ c x^-alpha and an infinite variance, so within
the family finite variance and Gaussian are the same condition. Three
measurements separate the cases, and none of them reads a normality test.

    variance against n     for a finite variance the sample variance settles as
                           n grows; for a stable alpha below 2 it grows like
                           n^(2/alpha - 1). The slope of log variance on log n
                           is 0 in the first case and 2/alpha - 1 in the second.
    the tail index         Hill's estimator on the k largest of |x|,
                           alpha_hat = k / sum_i (ln |x|_(i) - ln |x|_(k+1)),
                           with |x|_(1) the largest. A power-law tail with
                           exponent alpha drives it to alpha; a Gaussian, which
                           has no power-law tail, drives it upward with n and
                           with k.
    scale against block    summing k consecutive quarters multiplies the scale
                           of a stable law by k^(1/alpha) and leaves its shape
                           alone. With a finite variance and independent terms
                           the factor is k^(1/2). The slope of log scale on
                           log k is 1/alpha in the first case and 1/2 in the
                           second. The scale is the interquartile range as well
                           as the standard deviation, because the standard
                           deviation of a sample from an infinite-variance law
                           estimates nothing while the interquartile range is
                           defined for every stable law.

None of the three is read off a formula for what it should give. Part 1 runs
all three on samples of the same length drawn from laws whose answer is known:
a Gaussian, and symmetric stable laws at alpha = 1.8 and alpha = 1.5, drawn by
the Chambers-Mallows-Stuck method. Whatever those return is what the model's
numbers are read against.

Three things none of this controls for. A run's quarters are dependent, and
dependence moves the block-scale slope on its own, upward for positive
dependence, so a slope above 1/2 is not by itself a tail index below 2. A
series that is not stationary makes the sample variance grow with n whatever
its tails are. And Hill's estimator assumes the tail is already power-law over
the range of order statistics used, which at n = 400 and k = 40 it need not be.

Part 5 runs Hill on the cross-section of firm output changes from
dataset/abm_system_micro, which is the sum GDP's change is, rather than on the
aggregate.

    ./bin/abm_system_stable_tails [CONFIGURATIONS]

Reads dataset/abm_system and dataset/abm_system_micro. Writes, none of it
printed:
    out/abm_system_stable_tails_report.txt
    out/abm_system_stable_tails_by_series.csv
*/

#include "applications/abm_system.h"
#include <et_al./stats.h>
#include <et_al./random/random.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define K ABM_SYSTEM_K
#define INPUT_DIR "dataset/abm_system"
#define MICRO_DIR "dataset/abm_system_micro"
#define REPORT_PATH "out/abm_system_stable_tails_report.txt"
#define BY_SERIES_PATH "out/abm_system_stable_tails_by_series.csv"

#define PERIODS 400
#define REPLICATES 1000
#define MICRO_SEEDS 5

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

/* Sample sizes the variance is taken over. */
#define N_LENGTHS 8
static const int length[N_LENGTHS] = {50, 100, 150, 200, 250, 300, 350, 400};

/* Block lengths the scale is taken at. The last incomplete block is dropped,
   so k = 32 leaves 12 values of a 400-quarter run. */
#define N_BLOCKS 6
static const int block_length[N_BLOCKS] = {1, 2, 4, 8, 16, 32};

/* Fractions of the sample Hill's estimator uses. */
#define N_HILL 3
static const double hill_fraction[N_HILL] = {0.10, 0.05, 0.025};

#define CALIBRATION_DRAWS 20000
#define CALIBRATION_SEED 20260923

static int ascending(const void *a, const void *b) {
    double left = *(const double*)a, right = *(const double*)b;
    return left < right ? -1 : left > right ? 1 : 0;
}

static int descending(const void *a, const void *b) { return ascending(b, a); }

/* Hill's estimator of the tail index of |x| from its k largest values. buffer
   holds at least n doubles. NAN when a value is zero, where the logarithm it
   needs does not exist, or when the k largest are not distinct enough to give
   a positive denominator. */
static double hill_of(const double *x, int n, int k, double *buffer) {
    if (k < 2 || k + 1 > n) return NAN;
    for (int i = 0; i < n; i++) buffer[i] = fabs(x[i]);
    qsort(buffer, (size_t)n, sizeof(double), descending);
    if (!(buffer[k] > 0)) return NAN;

    double sum = 0;
    for (int i = 0; i < k; i++) sum += log(buffer[i]) - log(buffer[k]);
    return sum > 0 ? (double)k / sum : NAN;
}

static double variance_of(const double *x, int n) {
    if (n < 2) return NAN;
    double mean = 0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= (double)n;
    double m2 = 0;
    for (int i = 0; i < n; i++) { double d = x[i] - mean; m2 += d * d; }
    return m2 / (double)(n - 1);
}

/* The interquartile range, which is defined whatever the tails are. buffer
   holds at least n doubles. */
static double iqr_of(const double *x, int n, double *buffer) {
    if (n < 4) return NAN;
    for (int i = 0; i < n; i++) buffer[i] = x[i];
    qsort(buffer, (size_t)n, sizeof(double), ascending);
    double lower = buffer[(int)(0.25 * (n - 1))];
    double upper = buffer[(int)(0.75 * (n - 1))];
    return upper - lower;
}

/* Slope of log y on log x over the points where both are positive. NAN when
   fewer than two such points. */
static double log_slope(const double *x, const double *y, int n) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int used = 0;
    for (int i = 0; i < n; i++) {
        if (!(x[i] > 0) || !(y[i] > 0)) continue;
        double lx = log(x[i]), ly = log(y[i]);
        sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
        used++;
    }
    if (used < 2) return NAN;
    double denominator = (double)used * sxx - sx * sx;
    return denominator != 0 ? ((double)used * sxy - sx * sy) / denominator : NAN;
}

/* One symmetric stable draw by Chambers, Mallows and Stuck. alpha = 2 is the
   Gaussian, up to the scale convention, and is drawn directly instead. */
static double stable_draw(Rng *rng, double alpha) {
    if (alpha >= 2.0) return rng_normal(rng);
    const double pi = 3.14159265358979323846;
    double u = pi * (rng_uniform(rng) - 0.5);
    double w = -log(rng_uniform(rng) + 1e-300);
    double term = sin(alpha * u) / pow(cos(u), 1.0 / alpha);
    return term * pow(cos(u - alpha * u) / w, (1.0 - alpha) / alpha);
}

/* The three measurements of one sample of PERIODS values. */
typedef struct {
    double variance_slope;
    double hill[N_HILL];
    double iqr_slope;
    double sd_slope;
} Measurement;

static void measure(const double *x, double *buffer, double *work, Measurement *out) {
    double n_points[N_LENGTHS], variance[N_LENGTHS];
    for (int i = 0; i < N_LENGTHS; i++) {
        n_points[i] = (double)length[i];
        variance[i] = variance_of(x, length[i]);
    }
    out->variance_slope = log_slope(n_points, variance, N_LENGTHS);

    for (int i = 0; i < N_HILL; i++)
        out->hill[i] = hill_of(x, PERIODS, (int)(hill_fraction[i] * PERIODS), buffer);

    double k_points[N_BLOCKS], iqr[N_BLOCKS], sd[N_BLOCKS];
    for (int b = 0; b < N_BLOCKS; b++) {
        int k = block_length[b], blocks = PERIODS / k;
        for (int j = 0; j < blocks; j++) {
            double total = 0;
            for (int i = 0; i < k; i++) total += x[j * k + i];
            work[j] = total;
        }
        k_points[b] = (double)k;
        iqr[b] = iqr_of(work, blocks, buffer);
        double v = variance_of(work, blocks);
        sd[b] = v > 0 ? sqrt(v) : NAN;
    }
    out->iqr_slope = log_slope(k_points, iqr, N_BLOCKS);
    out->sd_slope = log_slope(k_points, sd, N_BLOCKS);
}

/* Running totals over the samples of one group. */
typedef struct {
    long samples;
    double variance_slope;
    double hill[N_HILL];
    long hill_counted[N_HILL];
    double iqr_slope;
    double sd_slope;
} Group;

static void group_add(Group *group, const Measurement *m) {
    group->samples++;
    group->variance_slope += m->variance_slope;
    group->iqr_slope += m->iqr_slope;
    group->sd_slope += m->sd_slope;
    for (int i = 0; i < N_HILL; i++)
        if (!isnan(m->hill[i])) { group->hill[i] += m->hill[i]; group->hill_counted[i]++; }
}

static void group_merge(Group *into, const Group *from) {
    into->samples += from->samples;
    into->variance_slope += from->variance_slope;
    into->iqr_slope += from->iqr_slope;
    into->sd_slope += from->sd_slope;
    for (int i = 0; i < N_HILL; i++) {
        into->hill[i] += from->hill[i];
        into->hill_counted[i] += from->hill_counted[i];
    }
}

static void group_print(FILE *report, const char *label, const Group *g) {
    if (!g->samples) return;
    double inverse = 1.0 / (double)g->samples;
    fprintf(report, "  %-24s %9ld %10.3f", label, g->samples, g->variance_slope * inverse);
    for (int i = 0; i < N_HILL; i++)
        fprintf(report, " %8.2f", g->hill_counted[i] ? g->hill[i] / (double)g->hill_counted[i] : NAN);
    fprintf(report, " %9.3f %9.3f\n", g->iqr_slope * inverse, g->sd_slope * inverse);
}

static int *list_configurations(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    if (!handle) { *count = 0; return NULL; }

    int *cop = NULL, n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int here;
        if (sscanf(entry->d_name, "cop_%d", &here) != 1) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 1024;
            cop = (int*)realloc(cop, (size_t)cap * sizeof(int));
            assert(cop && "abm_system_stable_tails: out of memory");
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
    if (n_columns <= 0) { fclose(f); return NULL; }
    rewind(f);

    int capacity = n_columns * 1024, filled = 0;
    double *value = (double*)malloc((size_t)capacity * sizeof(double));
    assert(value && "abm_system_stable_tails: out of memory");
    double here;
    while (fscanf(f, "%lf", &here) == 1) {
        if (filled == capacity) {
            capacity *= 2;
            value = (double*)realloc(value, (size_t)capacity * sizeof(double));
            assert(value && "abm_system_stable_tails: out of memory");
        }
        value[filled++] = here;
    }
    fclose(f);
    *rows = filled / n_columns;
    *columns = n_columns;
    return value;
}

int main(int argc, char **argv) {
    int n_configurations;
    int *cop = list_configurations(INPUT_DIR, &n_configurations);
    assert(n_configurations > 0 && "abm_system_stable_tails: no configurations stored");
    if (argc > 1) {
        int wanted = atoi(argv[1]);
        assert(wanted >= 1 && "abm_system_stable_tails: nothing to read");
        if (wanted < n_configurations) n_configurations = wanted;
    }

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_stable_tails: cannot open the report path");
    fprintf(report, "Whether the macro series are stable with index below 2, or have a finite\n"
                    "variance.\n\n");
    fprintf(report, "Columns. 'var slope' is the slope of log sample variance on log n over\n"
                    "n = %d to %d: 0 for a finite variance, 2/alpha - 1 for a stable alpha.\n",
            length[0], length[N_LENGTHS - 1]);
    fprintf(report, "'Hill' is the tail index from the largest 10, 5 and 2.5 per cent of |x|.\n");
    fprintf(report, "'iqr slope' and 'sd slope' are the slope of log scale on log block length\n"
                    "over k = %d to %d: 1/2 for a finite variance with independent terms,\n"
                    "1/alpha for a stable alpha.\n\n",
            block_length[0], block_length[N_BLOCKS - 1]);

    /* Part 1. The same three measurements on laws whose answer is known. */
    fprintf(report, "PART 1. Calibration. %d samples of %d values from each law, seed %d.\n",
            CALIBRATION_DRAWS, PERIODS, CALIBRATION_SEED);
    fprintf(report, "A stable alpha should give var slope 2/alpha - 1, Hill alpha, and\n"
                    "iqr slope 1/alpha. The Gaussian row is alpha = 2: 0, no power-law tail,\n"
                    "and 0.5.\n\n");
    fprintf(report, "  law                        samples  var slope  Hill 10%%  Hill 5%%  Hill2.5%%"
                    "  iqr slope   sd slope\n");

    static const double calibration_alpha[3] = {2.0, 1.8, 1.5};
    static const char *calibration_name[3] = {
        "Gaussian (alpha = 2)", "stable alpha = 1.8", "stable alpha = 1.5"
    };
    for (int a = 0; a < 3; a++) {
        Group group = {0};
        #pragma omp parallel
        {
            Group mine = {0};
            Rng rng = rng_new(CALIBRATION_SEED, (uint64_t)(a * 64 + omp_get_thread_num() + 1));
            double *sample = (double*)malloc((size_t)PERIODS * sizeof(double));
            double *buffer = (double*)malloc((size_t)PERIODS * sizeof(double));
            double *work = (double*)malloc((size_t)PERIODS * sizeof(double));
            assert(sample && buffer && work && "abm_system_stable_tails: out of memory");

            #pragma omp for schedule(static)
            for (int draw = 0; draw < CALIBRATION_DRAWS; draw++) {
                for (int i = 0; i < PERIODS; i++) sample[i] = stable_draw(&rng, calibration_alpha[a]);
                Measurement m;
                measure(sample, buffer, work, &m);
                group_add(&mine, &m);
            }
            #pragma omp critical(calibrate)
            group_merge(&group, &mine);
            free(sample); free(buffer); free(work);
        }
        group_print(report, calibration_name[a], &group);
    }

    /* Parts 2 to 4. The same three on every replicate of the experiment. */
    Group *series_group = (Group*)calloc((size_t)N_SERIES, sizeof(Group));
    assert(series_group && "abm_system_stable_tails: out of memory");
    long replicates_read = 0;
    const double sweep_started = omp_get_wtime();

    #pragma omp parallel
    {
        Group *mine = (Group*)calloc((size_t)N_SERIES, sizeof(Group));
        double *sample = (double*)malloc((size_t)PERIODS * sizeof(double));
        double *buffer = (double*)malloc((size_t)PERIODS * sizeof(double));
        double *work = (double*)malloc((size_t)PERIODS * sizeof(double));
        assert(mine && sample && buffer && work && "abm_system_stable_tails: out of memory");
        long my_replicates = 0;

        #pragma omp for schedule(dynamic)
        for (int c = 0; c < n_configurations; c++) {
            char dir[512];
            snprintf(dir, sizeof dir, "%s/cop_%04d", INPUT_DIR, cop[c]);

            for (int batch = 0; batch < REPLICATES / ABM_SYSTEM_BATCH; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);

                for (int b = 0; b < count; b++) {
                    Mat y = block[b];
                    if (y.c < PERIODS) { mat_free(y); continue; }
                    my_replicates++;

                    for (int s = 0; s < N_SERIES; s++) {
                        if (s < K) {
                            for (int t = 0; t < PERIODS; t++) sample[t] = (double)AT(y, s, t);
                        } else {
                            /* One value short, so the last is repeated rather
                               than the block grid changing for one series. */
                            for (int t = 1; t < PERIODS; t++)
                                sample[t - 1] = (double)AT(y, ROW_INTEREST_RATE, t) -
                                                (double)AT(y, ROW_INTEREST_RATE, t - 1);
                            sample[PERIODS - 1] = sample[PERIODS - 2];
                        }
                        Measurement m;
                        measure(sample, buffer, work, &m);
                        group_add(&mine[s], &m);
                    }
                    mat_free(y);
                }
            }
        }

        #pragma omp critical(merge)
        {
            for (int s = 0; s < N_SERIES; s++) group_merge(&series_group[s], &mine[s]);
            replicates_read += my_replicates;
        }
        free(mine); free(sample); free(buffer); free(work);
    }
    const double sweep_elapsed = omp_get_wtime() - sweep_started;

    fprintf(report, "\nPARTS 2 TO 4. The experiment's own series, %d configurations, %ld\n"
                    "replicates, 400 quarters each, one sample per replicate. %.1f minutes.\n\n",
            n_configurations, replicates_read, sweep_elapsed / 60.0);
    fprintf(report, "  series                     samples  var slope  Hill 10%%  Hill 5%%  Hill2.5%%"
                    "  iqr slope   sd slope\n");
    for (int s = 0; s < N_SERIES; s++) group_print(report, series_name[s], &series_group[s]);

    /* Part 5. Hill on the cross-section of firm output changes. */
    int n_micro;
    int *micro_cop = list_configurations(MICRO_DIR, &n_micro);
    if (n_micro > 0) {
        double hill_total[2][N_HILL] = {{0}};
        long hill_counted[2][N_HILL] = {{0}};
        long cross_sections = 0;
        int firms[2] = {0, 0};
        static const char *micro_file[2] = {"Q1all", "Q2all"};
        static const char *micro_name[2] = {"K-firm changes", "C-firm changes"};

        for (int type = 0; type < 2; type++) {
            for (int c = 0; c < n_micro; c++) {
                for (int seed = 1; seed <= MICRO_SEEDS; seed++) {
                    char path[640];
                    snprintf(path, sizeof path, "%s/cop_%04d/%s_%d.txt",
                             MICRO_DIR, micro_cop[c], micro_file[type], seed);
                    int rows, columns;
                    double *panel = read_panel(path, &rows, &columns);
                    if (!panel) continue;
                    firms[type] = columns;

                    double *change = (double*)malloc((size_t)columns * sizeof(double));
                    double *buffer = (double*)malloc((size_t)columns * sizeof(double));
                    assert(change && buffer && "abm_system_stable_tails: out of memory");

                    for (int t = ABM_SYSTEM_BURN_IN; t < rows; t++) {
                        for (int j = 0; j < columns; j++)
                            change[j] = panel[t * columns + j] - panel[(t - 1) * columns + j];
                        for (int i = 0; i < N_HILL; i++) {
                            int k = (int)(hill_fraction[i] * columns);
                            double alpha = hill_of(change, columns, k, buffer);
                            if (!isnan(alpha)) { hill_total[type][i] += alpha; hill_counted[type][i]++; }
                        }
                        if (type == 0) cross_sections++;
                    }
                    free(change); free(buffer); free(panel);
                }
            }
        }

        fprintf(report, "\nPART 5. Hill on the cross-section of firm output changes at each date,\n"
                        "%ld cross-sections from dataset/abm_system_micro. The 10, 5 and 2.5 per\n"
                        "cent of a cross-section this small is a handful of firms, so these are\n"
                        "the noisiest numbers here.\n\n", cross_sections);
        fprintf(report, "  cross-section              firms  Hill 10%%  Hill 5%%  Hill2.5%%\n");
        for (int type = 0; type < 2; type++) {
            if (!firms[type]) continue;
            fprintf(report, "  %-24s %7d", micro_name[type], firms[type]);
            for (int i = 0; i < N_HILL; i++)
                fprintf(report, " %9.2f", hill_counted[type][i]
                        ? hill_total[type][i] / (double)hill_counted[type][i] : NAN);
            fprintf(report, "\n");
        }
        free(micro_cop);
    }

    FILE *table = fopen(BY_SERIES_PATH, "w");
    assert(table && "abm_system_stable_tails: cannot open the table path");
    fprintf(table, "group,samples,variance_slope,hill_10,hill_5,hill_2p5,iqr_slope,sd_slope\n");
    for (int s = 0; s < N_SERIES; s++) {
        const Group *g = &series_group[s];
        if (!g->samples) continue;
        double inverse = 1.0 / (double)g->samples;
        fprintf(table, "%s,%ld,%.9g", series_column[s], g->samples, g->variance_slope * inverse);
        for (int i = 0; i < N_HILL; i++)
            fprintf(table, ",%.9g", g->hill_counted[i] ? g->hill[i] / (double)g->hill_counted[i] : NAN);
        fprintf(table, ",%.9g,%.9g\n", g->iqr_slope * inverse, g->sd_slope * inverse);
    }
    fclose(table);

    fprintf(report, "\nEvery series' own row is in %s.\n", BY_SERIES_PATH);
    fclose(report);

    free(series_group);
    free(cop);
    return EXIT_SUCCESS;
}
