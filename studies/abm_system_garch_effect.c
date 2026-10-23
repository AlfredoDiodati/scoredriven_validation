/*
Whether the simulated macro series have a GARCH effect, meaning that their
shocks cluster: a large shock this quarter makes a large one next quarter more
likely.

The series are the five the experiment stores, as stored: GDP growth, energy
growth, employment change and inflation in per cent changes, and the interest
rate in per cent levels.

The test is McLeod and Li's. Fit an autoregression with an intercept and
P_MEAN lags to the series, square its residuals, and run the Ljung-Box test on
the squares. With no conditional heteroskedasticity the squares are
uncorrelated and the statistic is chi-squared with as many degrees of freedom
as lags, with no adjustment for the autoregression's own coefficients. It is
run at 4 and 8 lags. Rejection is a p-value below 0.05.

The simulator sets the interest rate by a smoothed Taylor rule,

    r_t = taylor r_(t-1) + (1 - taylor) ((1 + r_a)^(1/4) - 1),

with r_a a linear function of annual inflation and the unemployment gap, and
replaces any value at or below zero by 1e-6, which is 1e-4 in the per cent
units stored here. Nothing in the rule makes the variance of the shocks depend
on past shocks, but the floor can make the squared residuals correlated on its
own: while the rate sits on the floor its residuals are near zero quarter after
quarter, and away from it they are not. The interest rate's rejection rates are
therefore also reported by the share of a replicate's quarters on the floor, so
that the replicates that never touch it give the answer free of that mechanism.

Part 1 runs the same procedure, the same autoregression and the same lags, on
samples of PERIODS values from laws whose answer is known. All four are
x_t = (1 - phi) mu + phi x_(t-1) + e_t with phi = 0.9:
    Gaussian        e_t standard normal; no GARCH effect, so the rejection
                    rate should be 0.05.
    Student t4      e_t Student t with 4 degrees of freedom scaled to unit
                    variance; no GARCH effect, heavy tails, so whatever it
                    rejects beyond 0.05 is what heavy tails alone produce.
    GARCH(1,1)      e_t = sqrt(h_t) z_t, h_t = 0.05 + 0.10 e_(t-1)^2 +
                    0.85 h_(t-1); a GARCH effect, so its rejection rate is
                    the power of the test at this sample size.
    floored AR(1)   Gaussian e_t, mu = 2, and x_t replaced by 0 whenever it
                    falls below 0, which keeps it on the floor about a fifth
                    of the time; no GARCH effect, so whatever it rejects is
                    what a floor alone produces.

Part 4 runs the same test on the US series of applications/us_data.h, 187
quarters from 1973Q2 to 2019Q4 in the same layout as the simulated ones. That is
less than half the simulated length, so Part 1 is also run at 187 values, and
the US p-values are read against those rates rather than the 400-value ones.
The US rate sat between 0.07 and 0.24 per cent from 2008Q3 to 2015Q3, its own
floor, so it is tested again on the quarters before the first one below
US_FLOOR_THRESHOLD.

    ./bin/abm_system_garch_effect [CONFIGURATIONS]

Reads out/dataset/abm_system and out/applications/us_system.csv. Writes, none of
it printed:
    out/studies/abm_system_garch_effect_report.txt
    out/studies/abm_system_garch_effect_by_configuration.csv
    out/studies/abm_system_garch_effect_us.csv
*/

#include "applications/abm_system.h"
#include "applications/us_data.h"
#include <et_al./regression.h>
#include <et_al./stats.h>
#include <et_al./random/random.h>
#include <dirent.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define K ABM_SYSTEM_K
#define INPUT_DIR "out/dataset/abm_system"
#define REPORT_PATH "out/studies/abm_system_garch_effect_report.txt"
#define BY_CONFIGURATION_PATH "out/studies/abm_system_garch_effect_by_configuration.csv"
#define US_PATH "out/studies/abm_system_garch_effect_us.csv"

#define PERIODS 400
#define REPLICATES 1000

static const char *series_name[K] = {
    "GDP growth", "energy growth", "employment change", "inflation", "interest rate"
};
static const char *series_column[K] = {
    "gdp_growth", "energy_growth", "employment_change", "inflation", "interest_rate"
};

/* Lags of the autoregression the residuals come from: one year of quarters,
   which is also how far back the Taylor rule reads inflation. */
#define P_MEAN 4

#define N_LAGS 2
static const int test_lag[N_LAGS] = {4, 8};
#define MAX_TEST_LAG 8
#define SIGNIFICANCE 0.05

/* A stored rate at or below this, in per cent, counts as on the floor. The
   floor itself is 1e-4; the margin is half of it. */
#define FLOOR_THRESHOLD 1.5e-4

/* A US rate below this, in per cent, counts as on its floor. */
#define US_FLOOR_THRESHOLD 0.25

/* Groups of replicates by the share of their PERIODS quarters on the floor:
   none, up to 5 per cent, up to 25 per cent, more. */
#define N_FLOOR_GROUPS 4
static const double floor_group_upper[N_FLOOR_GROUPS - 1] = {0.0, 0.05, 0.25};
static const char *floor_group_name[N_FLOOR_GROUPS] = {
    "never on the floor", "floor share (0, 0.05]", "floor share (0.05, 0.25]",
    "floor share above 0.25"
};

#define CALIBRATION_DRAWS 10000
#define CALIBRATION_SEED 20261016
#define CALIBRATION_BURN_IN 200
#define CALIBRATION_PHI 0.9

enum { LAW_GAUSSIAN, LAW_STUDENT, LAW_GARCH, LAW_FLOORED, N_LAWS };
static const char *law_name[N_LAWS] = {
    "Gaussian", "Student t4", "GARCH(1,1)", "floored AR(1)"
};

typedef struct {
    int degenerate;
    double p_value[N_LAGS];
    /* First-order autocorrelation of the squared residuals, the size of the
       effect the test is looking for. */
    double squared_autocorrelation;
} ArchTest;

/* stats_autocorr asserts on a constant window, so every pair of windows the
   Ljung-Box sum will correlate is checked first. */
static int lag_windows_vary(const mreal *value, int n, int max_lag) {
    for (int lag = 1; lag <= max_lag; lag++) {
        double low_min = (double)value[0], low_max = low_min;
        double high_min = (double)value[lag], high_max = high_min;
        for (int i = 1; i < n - lag; i++) {
            low_min = fmin(low_min, value[i]);
            low_max = fmax(low_max, value[i]);
        }
        for (int i = lag + 1; i < n; i++) {
            high_min = fmin(high_min, value[i]);
            high_max = fmax(high_max, value[i]);
        }
        if (!(low_max > low_min) || !(high_max > high_min)) return 0;
    }
    return 1;
}

/* McLeod-Li on x_0..x_(n-1): Ljung-Box on the squared residuals of an
   autoregression with an intercept and P_MEAN lags. */
static ArchTest arch_test(const double *x, int n) {
    ArchTest result = {0};
    int rows = n - P_MEAN;
    Mat design = mat_new(rows, P_MEAN + 1);
    Mat response = mat_new(rows, 1);
    for (int i = 0; i < rows; i++) {
        int t = i + P_MEAN;
        AT(design, i, 0) = (mreal)1;
        for (int lag = 1; lag <= P_MEAN; lag++) AT(design, i, lag) = (mreal)x[t - lag];
        AT(response, i, 0) = (mreal)x[t];
    }

    OlsFit fit = ols(design, response);
    mat_free(design);
    mat_free(response);
    if (fit.status < 0) {
        result.degenerate = 1;
        return result;
    }

    Mat squared = mat_new(rows, 1);
    for (int i = 0; i < rows; i++) {
        double residual = (double)AT(fit.residuals, i, 0);
        AT(squared, i, 0) = (mreal)(residual * residual);
    }
    ols_free(&fit);

    if (!lag_windows_vary(squared.d, rows, MAX_TEST_LAG)) {
        result.degenerate = 1;
        mat_free(squared);
        return result;
    }
    for (int l = 0; l < N_LAGS; l++)
        result.p_value[l] = stats_ljung_box(squared, test_lag[l]).p_value;
    result.squared_autocorrelation = (double)stats_autocorr(squared, 1);
    mat_free(squared);
    return result;
}

typedef struct {
    long samples;
    long degenerate;
    long rejected[N_LAGS];
    double squared_autocorrelation;
} Group;

static void group_add(Group *group, const ArchTest *test) {
    group->samples++;
    if (test->degenerate) { group->degenerate++; return; }
    for (int l = 0; l < N_LAGS; l++)
        if (test->p_value[l] < SIGNIFICANCE) group->rejected[l]++;
    group->squared_autocorrelation += test->squared_autocorrelation;
}

static void group_merge(Group *into, const Group *from) {
    into->samples += from->samples;
    into->degenerate += from->degenerate;
    for (int l = 0; l < N_LAGS; l++) into->rejected[l] += from->rejected[l];
    into->squared_autocorrelation += from->squared_autocorrelation;
}

static double group_rejection(const Group *g, int l) {
    long tested = g->samples - g->degenerate;
    return tested ? (double)g->rejected[l] / (double)tested : NAN;
}

static void group_print(FILE *report, const char *label, const Group *g) {
    long tested = g->samples - g->degenerate;
    fprintf(report, "  %-26s %9ld %10ld", label, g->samples, g->degenerate);
    for (int l = 0; l < N_LAGS; l++) fprintf(report, " %10.4f", group_rejection(g, l));
    fprintf(report, " %12.4f\n", tested ? g->squared_autocorrelation / (double)tested : NAN);
}

static void simulate_law(Rng *rng, int law, double *x, int n) {
    const double mu = law == LAW_FLOORED ? 2.0 : 0.0;
    const double garch_omega = 0.05, garch_alpha = 0.10, garch_beta = 0.85;
    double previous = mu;
    double h = garch_omega / (1.0 - garch_alpha - garch_beta), previous_shock = 0;

    for (int t = -CALIBRATION_BURN_IN; t < n; t++) {
        double shock;
        if (law == LAW_STUDENT) {
            /* t4 has variance 2, hence the sqrt(2). */
            double chi_squared = 2.0 * rng_gamma(rng, 2.0);
            shock = rng_normal(rng) / sqrt(chi_squared / 4.0) / sqrt(2.0);
        } else if (law == LAW_GARCH) {
            h = garch_omega + garch_alpha * previous_shock * previous_shock + garch_beta * h;
            shock = sqrt(h) * rng_normal(rng);
            previous_shock = shock;
        } else {
            shock = rng_normal(rng);
        }
        double value = (1.0 - CALIBRATION_PHI) * mu + CALIBRATION_PHI * previous + shock;
        if (law == LAW_FLOORED && value < 0) value = 0;
        previous = value;
        if (t >= 0) x[t] = value;
    }
}

/* Part 1 at one sample length. The seed does not depend on n, so the shorter
   samples are not the leading values of the longer ones but fresh draws. */
static void calibrate(FILE *report, int n) {
    fprintf(report, "\n  %d values per sample\n", n);
    fprintf(report, "  law                          samples degenerate   reject %d   reject %d"
                    "      acf1 sq\n", test_lag[0], test_lag[1]);

    for (int law = 0; law < N_LAWS; law++) {
        Group group = {0};
        #pragma omp parallel
        {
            Group mine = {0};
            Rng rng = rng_new(CALIBRATION_SEED, (uint64_t)(law * 64 + omp_get_thread_num() + 1));
            double *sample = (double*)malloc((size_t)n * sizeof(double));
            assert(sample && "abm_system_garch_effect: out of memory");

            #pragma omp for schedule(static)
            for (int draw = 0; draw < CALIBRATION_DRAWS; draw++) {
                simulate_law(&rng, law, sample, n);
                ArchTest test = arch_test(sample, n);
                group_add(&mine, &test);
            }
            #pragma omp critical(calibrate)
            group_merge(&group, &mine);
            free(sample);
        }
        group_print(report, law_name[law], &group);
    }
}

static int floor_group_of(double floor_share) {
    for (int g = 0; g < N_FLOOR_GROUPS - 1; g++)
        if (floor_share <= floor_group_upper[g]) return g;
    return N_FLOOR_GROUPS - 1;
}

static int ascending(const void *a, const void *b) {
    double left = *(const double*)a, right = *(const double*)b;
    return left < right ? -1 : left > right ? 1 : 0;
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
            assert(cop && "abm_system_garch_effect: out of memory");
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
    int *cop = list_configurations(INPUT_DIR, &n_configurations);
    assert(n_configurations > 0 && "abm_system_garch_effect: no configurations stored");
    if (argc > 1) {
        int wanted = atoi(argv[1]);
        assert(wanted >= 1 && "abm_system_garch_effect: nothing to read");
        if (wanted < n_configurations) n_configurations = wanted;
    }

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_garch_effect: cannot open the report path");
    fprintf(report, "Whether the simulated macro series have a GARCH effect.\n\n");
    fprintf(report, "Test: McLeod-Li, the Ljung-Box test on the squared residuals of an\n"
                    "autoregression of the series as stored, with an intercept and %d lags, at\n"
                    "%d and %d lags of the squares. Rejection is p < %.2f.\n\n",
            P_MEAN, test_lag[0], test_lag[1], SIGNIFICANCE);
    fprintf(report, "Columns. 'samples' is how many series were tested. 'degenerate' is how\n"
                    "many could not be, because the squared residuals were constant over a\n"
                    "window the test correlates, as they are for a rate on the floor\n"
                    "throughout. 'reject L' is the share of the others rejected at L lags.\n"
                    "'acf1 sq' is the mean first-order autocorrelation of the squared\n"
                    "residuals.\n\n");

    Mat us_original = load_us_system();
    Mat us = us_system_growth_block(us_original);
    mat_free(us_original);
    const int us_periods = us.c;

    fprintf(report, "PART 1. Calibration. %d samples from each law, seed %d,\n"
                    "x_t = (1 - phi) mu + phi x_(t-1) + e_t with phi = %.1f. The Gaussian, the\n"
                    "Student t4 and the floored AR(1) have no GARCH effect; the GARCH(1,1) has\n"
                    "omega = 0.05, alpha = 0.10, beta = 0.85. Run at the simulated length and\n"
                    "at the US length.\n",
            CALIBRATION_DRAWS, CALIBRATION_SEED, CALIBRATION_PHI);
    calibrate(report, PERIODS);
    calibrate(report, us_periods);

    /* Parts 2 and 3. Every replicate of the experiment, every series. The
       configuration groups are indexed c * K + s. */
    Group *configuration_group = (Group*)calloc((size_t)n_configurations * K, sizeof(Group));
    double *configuration_floor_share = (double*)calloc((size_t)n_configurations, sizeof(double));
    Group series_group[K] = {{0}};
    Group floor_group[N_FLOOR_GROUPS] = {{0}};
    assert(configuration_group && configuration_floor_share &&
           "abm_system_garch_effect: out of memory");
    const double sweep_started = omp_get_wtime();

    #pragma omp parallel
    {
        Group mine_by_series[K] = {{0}};
        Group mine_by_floor[N_FLOOR_GROUPS] = {{0}};
        double *sample = (double*)malloc((size_t)PERIODS * sizeof(double));
        assert(sample && "abm_system_garch_effect: out of memory");

        #pragma omp for schedule(dynamic)
        for (int c = 0; c < n_configurations; c++) {
            char dir[512];
            snprintf(dir, sizeof dir, "%s/cop_%04d", INPUT_DIR, cop[c]);
            double floor_share_total = 0;
            long replicates_here = 0;

            for (int batch = 0; batch < REPLICATES / ABM_SYSTEM_BATCH; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);

                for (int b = 0; b < count; b++) {
                    Mat y = block[b];
                    if (y.c < PERIODS) { mat_free(y); continue; }
                    replicates_here++;

                    for (int s = 0; s < K; s++) {
                        int on_floor = 0;
                        for (int t = 0; t < PERIODS; t++) {
                            sample[t] = (double)AT(y, s, t);
                            if (sample[t] <= FLOOR_THRESHOLD) on_floor++;
                        }
                        ArchTest test = arch_test(sample, PERIODS);
                        group_add(&configuration_group[c * K + s], &test);
                        group_add(&mine_by_series[s], &test);

                        if (s == ROW_INTEREST_RATE) {
                            double floor_share = (double)on_floor / (double)PERIODS;
                            floor_share_total += floor_share;
                            group_add(&mine_by_floor[floor_group_of(floor_share)], &test);
                        }
                    }
                    mat_free(y);
                }
            }
            if (replicates_here)
                configuration_floor_share[c] = floor_share_total / (double)replicates_here;
        }

        #pragma omp critical(merge)
        {
            for (int s = 0; s < K; s++) group_merge(&series_group[s], &mine_by_series[s]);
            for (int g = 0; g < N_FLOOR_GROUPS; g++) group_merge(&floor_group[g], &mine_by_floor[g]);
        }
        free(sample);
    }
    const double sweep_elapsed = omp_get_wtime() - sweep_started;

    fprintf(report, "\nPART 2. The experiment's own series, %d configurations, %ld replicates,\n"
                    "the %d quarters after burn-in of each, one test per replicate and series.\n"
                    "%.1f minutes.\n\n",
            n_configurations, series_group[0].samples, PERIODS, sweep_elapsed / 60.0);
    fprintf(report, "  series                       samples degenerate   reject %d   reject %d"
                    "      acf1 sq\n", test_lag[0], test_lag[1]);
    for (int s = 0; s < K; s++) group_print(report, series_name[s], &series_group[s]);

    fprintf(report, "\nThe interest rate by the share of a replicate's quarters on the floor. A\n"
                    "quarter is on the floor when the stored rate is at most %.1e per cent.\n\n",
            FLOOR_THRESHOLD);
    fprintf(report, "  interest rate replicates     samples degenerate   reject %d   reject %d"
                    "      acf1 sq\n", test_lag[0], test_lag[1]);
    for (int g = 0; g < N_FLOOR_GROUPS; g++) group_print(report, floor_group_name[g], &floor_group[g]);

    /* Whether rejection is spread over the design or concentrated in a few
       configurations. */
    fprintf(report, "\nPART 3. Rejection share at %d lags across the %d configurations, each\n"
                    "configuration's share over its own replicates: quartiles, and how many\n"
                    "configurations reject in fewer than 10 per cent and in more than 90 per\n"
                    "cent of their replicates.\n\n", test_lag[0], n_configurations);
    fprintf(report, "  series                   lower q   median  upper q   below 0.10   above 0.90\n");
    double *share = (double*)malloc((size_t)n_configurations * sizeof(double));
    assert(share && "abm_system_garch_effect: out of memory");
    for (int s = 0; s < K; s++) {
        int counted = 0, below_tenth = 0, above_ninety = 0;
        for (int c = 0; c < n_configurations; c++) {
            double here = group_rejection(&configuration_group[c * K + s], 0);
            if (here != here) continue;
            share[counted++] = here;
            if (here < 0.1) below_tenth++;
            if (here > 0.9) above_ninety++;
        }
        if (!counted) continue;
        qsort(share, (size_t)counted, sizeof(double), ascending);
        fprintf(report, "  %-22s %9.3f %8.3f %8.3f %12d %12d\n", series_name[s],
                share[(int)(0.25 * (counted - 1))], share[(int)(0.5 * (counted - 1))],
                share[(int)(0.75 * (counted - 1))], below_tenth, above_ninety);
    }
    free(share);

    FILE *table = fopen(BY_CONFIGURATION_PATH, "w");
    assert(table && "abm_system_garch_effect: cannot open the table path");
    fprintf(table, "configuration,series,replicates,degenerate,mean_floor_share,reject_lag_%d,"
                   "reject_lag_%d,mean_squared_autocorrelation\n", test_lag[0], test_lag[1]);
    for (int c = 0; c < n_configurations; c++)
        for (int s = 0; s < K; s++) {
            const Group *g = &configuration_group[c * K + s];
            long tested = g->samples - g->degenerate;
            fprintf(table, "%d,%s,%ld,%ld,", cop[c], series_column[s], g->samples, g->degenerate);
            if (s == ROW_INTEREST_RATE) fprintf(table, "%.9g", configuration_floor_share[c]);
            for (int l = 0; l < N_LAGS; l++) fprintf(table, ",%.9g", group_rejection(g, l));
            fprintf(table, ",%.9g\n", tested ? g->squared_autocorrelation / (double)tested : NAN);
        }
    fclose(table);

    /* Part 4. The US series. Column t of the growth block is quarter t + 1 of
       a sample that starts in 1973Q1. */
    const int us_first_year = 1973;
    int us_before_floor = us_periods;
    for (int t = 0; t < us_periods; t++)
        if ((double)AT(us, ROW_INTEREST_RATE, t) < US_FLOOR_THRESHOLD) { us_before_floor = t; break; }
    int last_year = us_first_year + us_before_floor / 4, last_quarter = us_before_floor % 4 + 1;

    fprintf(report, "\nPART 4. The US series, %d quarters from %dQ2. Each row is one test, so\n"
                    "the columns are p-values rather than rejection shares; read them against\n"
                    "Part 1 at %d values. The last row is the interest rate on the %d quarters\n"
                    "up to %dQ%d, before the first quarter below %.2f per cent.\n\n",
            us_periods, us_first_year, us_periods, us_before_floor, last_year, last_quarter,
            US_FLOOR_THRESHOLD);
    fprintf(report, "  series                       values   p at %d   p at %d      acf1 sq\n",
            test_lag[0], test_lag[1]);

    FILE *us_table = fopen(US_PATH, "w");
    assert(us_table && "abm_system_garch_effect: cannot open the US table path");
    fprintf(us_table, "series,values,p_value_lag_%d,p_value_lag_%d,squared_autocorrelation\n",
            test_lag[0], test_lag[1]);

    double *us_sample = (double*)malloc((size_t)us_periods * sizeof(double));
    assert(us_sample && "abm_system_garch_effect: out of memory");
    for (int s = 0; s <= K; s++) {
        int row = s < K ? s : ROW_INTEREST_RATE;
        int n = s < K ? us_periods : us_before_floor;
        for (int t = 0; t < n; t++) us_sample[t] = (double)AT(us, row, t);
        ArchTest test = arch_test(us_sample, n);
        char label[64];
        snprintf(label, sizeof label, "%s", s < K ? series_name[s] : "interest rate before floor");
        const char *column = s < K ? series_column[s] : "interest_rate_before_floor";
        if (test.degenerate) {
            fprintf(report, "  %-26s %8d   degenerate\n", label, n);
            fprintf(us_table, "%s,%d,nan,nan,nan\n", column, n);
            continue;
        }
        fprintf(report, "  %-26s %8d %9.4f %9.4f %12.4f\n", label, n, test.p_value[0], test.p_value[1],
                test.squared_autocorrelation);
        fprintf(us_table, "%s,%d,%.9g,%.9g,%.9g\n", column, n, test.p_value[0], test.p_value[1],
                test.squared_autocorrelation);
    }
    free(us_sample);
    fclose(us_table);
    mat_free(us);

    fprintf(report, "\nEvery configuration's own row is in %s, the US rows in %s.\n",
            BY_CONFIGURATION_PATH, US_PATH);
    fclose(report);

    free(configuration_floor_share);
    free(configuration_group);
    free(cop);
    return EXIT_SUCCESS;
}
