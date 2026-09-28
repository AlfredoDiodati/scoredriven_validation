/*
Why the local projections fail to recover the configuration a simulated run
came from, when the t-QVARMA succeeds.

docs/MONTECARLO_LP_VALIDATION.md records the failure: with each of cop_0191's
1000 runs standing in for the real data in turn, the confidence set holds
cop_0191 in 5, 72 and 42 of 1000 repetitions under the linear model and the two
states, against 876 of 1000 for the t-QVARMA. This study checks, one at a time,
the things in the local-projection protocol that could be wrong.

1. Is the stored result what the code computes? Every run is fitted again from
   the dataset and its loss against cop_0191 run 706 recomputed, and compared
   with montecarlo/out/lp_irf_loss_<model>.csv.
2. The shock size. The local projections answer a one-unit shock, as R's
   pipeline asks lpirfs to (shock_type 1); the t-QVARMA answers a shock of one
   standard deviation. A unit shock divides every response to shock j by that
   variable's own residual standard deviation, so configurations differ by how
   volatile their variables are, whatever their dynamics. Refitted with shocks
   of one standard deviation.
3. The transformation. R's transform_data builds log GDP as
   100 log(100 + cumsum(GDP_growth)), which is not the log level when the
   growth rates are already 100 log differences, as they are here. Refitted
   with the log levels themselves, the cumulated growth rates.
4. The loss. R uses the mean squared difference, the Monte Carlo the mean
   absolute difference. The linear model's table is recomputed as squared.
5. Noise against signal. A configuration's expected loss is the distance from
   its average response to the stand-in's plus how widely its own runs scatter.
   Both are measured for every configuration, together with how far a vector
   of zeros is from the stand-in, and which (response, shock) pairs carry the
   loss.
6. The volatility of each variable per configuration, which is what a unit
   shock divides by.

Each variant (items 2 and 3, and both) is run the way the Monte Carlo runs:
cop_0191 run 706 as the stand-in, run 706 of every configuration left out, the
mean absolute difference over the 400 response values, MCS_TR at level 0.05
with 10000 bootstrap resamples of single runs and the bootstrap variance, the
three models. Then, for the linear model, runs 0 to REPETITIONS - 1 of cop_0191
each take the stand-in's place in turn, with the same settings, and the count
of sets holding cop_0191 is reported.

Nothing here is stored beyond the report: every fit is recomputed from
dataset/abm_system/. Requires montecarlo/out/benchmark.env and, for item 1,
montecarlo/out/lp_irf_loss_<model>.csv. Writes
out/lp_recovery_diagnostics_report.txt. Progress to stderr.
*/

#include "applications/abm_system_lp.h"
#include "montecarlo/benchmark.h"
#include <et_al./inference/mcs.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <math.h>
#include "studies/recovery_noise.h"

#define REPORT_PATH "out/lp_recovery_diagnostics_report.txt"
#define REPETITIONS 100
#define K LP_SYSTEM_K
#define DIM LP_RESPONSE_DIM

static const char *const r_order[K] = { "GDP", "Employment", "Inflation", "InterestRate", "Energy" };
static const char *const stored_order[K] = { "GDP_growth", "EN_growth", "Employment_change", "Inflation",
                                             "InterestRate" };

typedef enum { TRANSFORM_R, TRANSFORM_LOG_LEVEL, TRANSFORM_GROWTH } TransformChoice;

typedef struct {
    const char *label;
    const char *description;
    TransformChoice transform;
    LpShockType shock;
    const char *const *names;
} Variant;

static const Variant variants[] = {
    { "r_unit", "R's transform_data, unit shocks: the protocol as run", TRANSFORM_R, LP_SHOCK_UNIT, r_order },
    { "r_sd", "R's transform_data, shocks of one standard deviation", TRANSFORM_R,
      LP_SHOCK_STANDARD_DEVIATION, r_order },
    { "loglevel_unit", "log levels as cumulated growth, unit shocks", TRANSFORM_LOG_LEVEL, LP_SHOCK_UNIT,
      r_order },
    { "loglevel_sd", "log levels as cumulated growth, shocks of one standard deviation",
      TRANSFORM_LOG_LEVEL, LP_SHOCK_STANDARD_DEVIATION, r_order },
    { "growth_unit", "the stored series the t-QVARMA is fitted on, in its order (GDP growth, energy growth, "
      "employment change, inflation, interest rate), unit shocks; state weight as in R", TRANSFORM_GROWTH,
      LP_SHOCK_UNIT, stored_order },
};
#define N_VARIANTS ((int)(sizeof variants / sizeof variants[0]))

static int n_samples, n_replicates, stand_in_sample, stand_in_replicate;
static LpSample *samples;
static const char *INPUT_DIR = LP_INPUT_DIR_DEFAULT;

/* R's layout; the same layout with GDP and energy as the cumulated growth
   rates, which are 100 log levels up to a constant, and the switching series
   their moving average; or the stored series themselves, the ones the
   t-QVARMA is fitted on, in its order, from the period the moving average
   first exists, with R's switching series. */
static void transform(Mat stored, TransformChoice choice, Mat *y, Mat *switching) {
    if (choice == TRANSFORM_GROWTH) {
        lp_system_from_stored(stored, y, switching);
        return;
    }
    lp_system_r_transform(stored, y, switching);
    if (choice == TRANSFORM_R) return;
    int first = LP_SYSTEM_SWITCHING_WINDOW - 1;
    double gdp = 0, energy = 0;
    double *level = malloc((size_t)stored.c * sizeof(double));
    for (int t = 0; t < stored.c; t++) {
        gdp += (double)AT(stored, STORED_GDP_GROWTH, t);
        energy += (double)AT(stored, STORED_EN_GROWTH, t);
        level[t] = gdp;
        if (t < first) continue;
        AT(*y, LP_R_ROW_GDP, t - first) = (mreal)gdp;
        AT(*y, LP_R_ROW_ENERGY, t - first) = (mreal)energy;
        double sum = 0;
        for (int w = 0; w < LP_SYSTEM_SWITCHING_WINDOW; w++) sum += level[t - w];
        AT(*switching, t - first, 0) = (mreal)(sum / LP_SYSTEM_SWITCHING_WINDOW);
    }
    free(level);
}

static int fit(Mat stored, const Variant *v, mreal out[LP_N_LOSSES][DIM]) {
    Mat y, switching;
    transform(stored, v->transform, &y, &switching);
    LpNlSpec spec = lp_system_spec();
    spec.lin.shock_type = v->shock;
    LpLinNlFit f = lp_lin_and_nl(y, switching, spec, (LpBands){0});
    int ok = f.lin.d.d && f.nl.d.d;
    if (ok) {
        memcpy(out[LP_LOSS_LIN], f.lin.irf_lin_mean.d, DIM * sizeof(mreal));
        memcpy(out[LP_LOSS_S1], f.nl.irf_s1_mean.d, DIM * sizeof(mreal));
        memcpy(out[LP_LOSS_S2], f.nl.irf_s2_mean.d, DIM * sizeof(mreal));
        for (int loss = 0; loss < LP_N_LOSSES && ok; loss++)
            for (int i = 0; i < DIM; i++)
                if (MISNAN(out[loss][i]) || MISINF(out[loss][i])) { ok = 0; break; }
    }
    lp_lin_nl_fit_free(&f);
    mat_free(y);
    mat_free(switching);
    return ok;
}

static double mean_absolute(const mreal *a, const mreal *b) {
    double s = 0;
    for (int i = 0; i < DIM; i++) s += fabs((double)a[i] - (double)b[i]);
    return s / DIM;
}

static double mean_squared(const mreal *a, const mreal *b) {
    double s = 0;
    for (int i = 0; i < DIM; i++) s += ((double)a[i] - (double)b[i]) * ((double)a[i] - (double)b[i]);
    return s / DIM;
}

static float *cache_slot(float *cache, int sample, int replicate) {
    return cache + ((size_t)sample * n_replicates + replicate) * DIM;
}

/* A confidence set over one loss table, rows replicate and columns
   configuration, with the rows in skip_row and any row holding a NaN left out. */
typedef struct {
    int in_set, rank, set_size, n_rows_used;
    int winner;
    double stand_in_config_mean, winner_mean;
    char members[256];
} SetOutcome;

static SetOutcome confidence_set(Mat values, int skip_row) {
    int keep = 0;
    int *usable = malloc((size_t)values.r * sizeof(int));
    for (int r = 0; r < values.r; r++) {
        usable[r] = r != skip_row;
        for (int c = 0; c < values.c && usable[r]; c++)
            if (MISNAN(AT(values, r, c))) usable[r] = 0;
        keep += usable[r];
    }
    Mat kept = mat_new(keep, values.c);
    for (int r = 0, at = 0; r < values.r; r++) {
        if (!usable[r]) continue;
        for (int c = 0; c < values.c; c++) AT(kept, at, c) = AT(values, r, c);
        at++;
    }
    free(usable);

    const char **names = malloc((size_t)n_samples * sizeof(char *));
    for (int c = 0; c < n_samples; c++) names[c] = samples[c].name;
    DataFrame losses = df_from_matrix(kept, names);
    mat_free(kept);

    MCSOptions opt = mcs_options_default();
    opt.bootstrap = 10000;
    opt.block_length = 1;
    opt.variance = MCS_VARIANCE_BOOTSTRAP;
    opt.stat = MCS_TR;
    MCSResult res = mcs(&losses, opt);

    SetOutcome o;
    memset(&o, 0, sizeof o);
    o.n_rows_used = keep;
    o.in_set = mcs_in_set(&res, stand_in_sample);
    o.set_size = res.n_surviving;
    double best = INFINITY;
    o.stand_in_config_mean = (double)stats_mean(df_col_numeric(&losses, names[stand_in_sample]));
    o.rank = 1;
    for (int c = 0; c < n_samples; c++) {
        double m = (double)stats_mean(df_col_numeric(&losses, names[c]));
        if (c != stand_in_sample && m < o.stand_in_config_mean) o.rank++;
        if (m < best) { best = m; o.winner = c; }
    }
    o.winner_mean = best;
    int written = 0, listed = 0;
    for (int c = 0; c < n_samples && listed < 8; c++)
        if (mcs_in_set(&res, c)) {
            written += snprintf(o.members + written, sizeof o.members - (size_t)written, "%s%s",
                                listed ? " " : "", names[c]);
            listed++;
        }
    if (o.set_size > listed) snprintf(o.members + written, sizeof o.members - (size_t)written, " ...");

    mcs_free(&res);
    df_free(&losses);
    free(names);
    return o;
}

static void report_set(FILE *f, const char *what, SetOutcome o) {
    fprintf(f, "  %-28s set of %d: %s\n", what, o.set_size, o.members);
    fprintf(f, "  %-28s cop_0191 %s the set, rank %d of %d by mean loss (%.4g); smallest mean loss %s %.4g\n",
            "", o.in_set ? "in" : "NOT in", o.rank, n_samples, o.stand_in_config_mean, samples[o.winner].name,
            o.winner_mean);
}

/* Section 1: the protocol as run must reproduce the stored loss tables. */
static void compare_with_stored(FILE *f, Mat values[LP_N_LOSSES]) {
    fprintf(f, "\n1. Recomputed losses against montecarlo/out/lp_irf_loss_<model>.csv\n");
    for (int loss = 0; loss < LP_N_LOSSES; loss++) {
        char path[256];
        snprintf(path, sizeof path, "montecarlo/out/lp_irf_loss_%s.csv", lp_loss_name(loss));
        DataFrame stored = df_read_csv(path, csv_read_options_default());
        Mat replicate = df_col_numeric(&stored, "replicate");
        double worst = 0;
        long compared = 0;
        for (int c = 0; c < n_samples; c++) {
            char name[128];
            snprintf(name, sizeof name, "%s_lp_%s", samples[c].name, lp_loss_name(loss));
            Mat column = df_col_numeric(&stored, name);
            for (int r = 0; r < stored.r; r++) {
                int row = (int)AT(replicate, r, 0);
                double a = (double)AT(column, r, 0), b = (double)AT(values[loss], row, c);
                double gap = fabs(a - b) / fmax(1.0, fabs(a));
                if (gap > worst) worst = gap;
                compared++;
            }
        }
        fprintf(f, "  %s: %ld cells compared, largest difference relative to max(1, |stored|) %.3g\n",
                lp_loss_name(loss), compared, worst);
        df_free(&stored);
    }
}

/* How precise the linear model's responses are (studies/recovery_noise.h),
   and which (response, shock) pairs the stand-in configuration's and the
   winner's losses come from. */
static int noise_and_signal(FILE *f, float *cache, const mreal *reference, const Variant *variant) {
    double reference_double[DIM];
    for (int i = 0; i < DIM; i++) reference_double[i] = (double)reference[i];
    const char **names = malloc((size_t)n_samples * sizeof(char *));
    for (int c = 0; c < n_samples; c++) names[c] = samples[c].name;
    RecoveryResponses responses = { cache, DIM, n_samples, n_replicates, stand_in_sample, stand_in_replicate,
                                    names };
    int winner = recovery_noise_report(f, &responses, reference_double);

    int tracked[2] = { stand_in_sample, winner };
    int n_tracked = winner == stand_in_sample ? 1 : 2;
    for (int t = 0; t < n_tracked; t++) {
        double pair_loss[K][K], total = 0;
        memset(pair_loss, 0, sizeof pair_loss);
        for (int r = 0; r < n_replicates; r++) {
            if (r == stand_in_replicate) continue;
            const float *x = cache_slot(cache, tracked[t], r);
            if (MISNAN(x[0])) continue;
            for (int k = 0; k < K; k++)
                for (int h = 0; h <= LP_SYSTEM_HOR; h++)
                    for (int j = 0; j < K; j++) {
                        int i = (k * (LP_SYSTEM_HOR + 1) + h) * K + j;
                        double e = fabs((double)x[i] - (double)reference[i]);
                        pair_loss[k][j] += e;
                        total += e;
                    }
        }
        fprintf(f, "  share of %s's loss by response (row) and shock (column), percent:\n", samples[tracked[t]].name);
        fprintf(f, "    %-18s", "");
        for (int j = 0; j < K; j++) fprintf(f, " %18s", variant->names[j]);
        fprintf(f, "\n");
        for (int k = 0; k < K; k++) {
            fprintf(f, "    %-18s", variant->names[k]);
            for (int j = 0; j < K; j++) fprintf(f, " %18.1f", 100.0 * pair_loss[k][j] / total);
            fprintf(f, "\n");
        }
    }
    free(names);
    return winner;
}

/* Section 6: the standard deviation of each series' period-to-period change,
   averaged over a configuration's runs, in the layout the protocol fits. */
static void volatility(FILE *f, int winner) {
    double *sd = calloc((size_t)n_samples * K, sizeof(double));
    #pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < n_samples; c++) {
        char dir[560];
        snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[c].name);
        int count = 0;
        for (int batch = 0; batch * ABM_SYSTEM_BATCH < n_replicates; batch++) {
            Mat block[ABM_SYSTEM_BATCH];
            int replicate[ABM_SYSTEM_BATCH];
            int n = abm_system_read_batch(dir, batch, block, replicate);
            for (int b = 0; b < n; b++) {
                Mat y, switching;
                lp_system_r_transform(block[b], &y, &switching);
                for (int k = 0; k < K; k++) {
                    Mat change = mat_new(1, y.c - 1);
                    for (int t = 1; t < y.c; t++) AT(change, 0, t - 1) = AT(y, k, t) - AT(y, k, t - 1);
                    sd[(size_t)c * K + k] += sqrt((double)stats_var(change));
                    mat_free(change);
                }
                count++;
                mat_free(y); mat_free(switching); mat_free(block[b]);
            }
        }
        for (int k = 0; k < K; k++) sd[(size_t)c * K + k] /= count;
    }
    fprintf(f, "\n6. Standard deviation of each series' one-period change, averaged over a configuration's runs, "
               "R's layout\n");
    fprintf(f, "  %-13s %11s %11s %11s %11s %11s %11s\n", "series", "min", "median", "max", "max/min",
            "cop_0191", samples[winner].name);
    Vec across = vec_new(n_samples);
    for (int k = 0; k < K; k++) {
        for (int c = 0; c < n_samples; c++) AT(across, c, 0) = (mreal)sd[(size_t)c * K + k];
        double lo = (double)stats_quantile(across, 0), hi = (double)stats_quantile(across, 1);
        fprintf(f, "  %-13s %11.4g %11.4g %11.4g %11.1f %11.4g %11.4g\n", r_order[k], lo,
                (double)stats_median(across), hi, hi / lo, sd[(size_t)stand_in_sample * K + k],
                sd[(size_t)winner * K + k]);
    }
    mat_free(across);
    free(sd);
}

int main(int argc, char **argv) {
    openblas_set_num_threads(1);

    /* No argument runs every variant into REPORT_PATH; one argument runs the
       variant of that label alone into a report named after it. */
    int only = -1;
    if (argc > 1) {
        for (int v = 0; v < N_VARIANTS; v++) if (strcmp(argv[1], variants[v].label) == 0) only = v;
        assert(only >= 0 && "lp_recovery_diagnostics: unknown variant label");
    }
    char report_path[256];
    if (only < 0) snprintf(report_path, sizeof report_path, "%s", REPORT_PATH);
    else snprintf(report_path, sizeof report_path, "out/lp_recovery_diagnostics_%s_report.txt", variants[only].label);
    Benchmark stand_in = benchmark_read();
    samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);
    stand_in_sample = -1;
    for (int c = 0; c < n_samples; c++) if (strcmp(samples[c].name, stand_in.sample) == 0) stand_in_sample = c;
    assert(stand_in_sample >= 0 && "lp_recovery_diagnostics: the stand-in's configuration is missing");
    stand_in_replicate = stand_in.replicate;

    FILE *f = fopen(report_path, "w");
    assert(f && "lp_recovery_diagnostics: cannot open the report");
    fprintf(f, "Why the local projections do not recover the stand-in's configuration\n\n");
    fprintf(f, "Data: dataset/abm_system, %d configurations x %d runs, 397 periods each after R's transform.\n",
            n_samples, n_replicates);
    fprintf(f, "Stand-in for the real data: %s run %d; run %d of every configuration left out.\n",
            stand_in.sample, stand_in_replicate, stand_in_replicate);
    fprintf(f, "Local projections: 4 lags, horizons 0 to 15, recursive ordering GDP, Employment, Inflation, "
               "InterestRate, Energy; state weight from the HP cycle (lambda 1600) of the 4-period moving average "
               "of the GDP row, logistic slope 2, lagged one period.\n");
    fprintf(f, "Loss: mean absolute difference over the 400 response values (5 responses x 16 horizons x 5 "
               "shocks), unless stated.\n");
    fprintf(f, "Confidence set: MCS_TR, level 0.05, 10000 bootstrap resamples of single runs, bootstrap "
               "variance, et_al seed 123.\n");
    fprintf(f, "Repetitions: runs 0 to %d of %s each take the stand-in's place in turn, linear model only, "
               "that run left out of every configuration.\n", REPETITIONS - 1, stand_in.sample);

    float *cache = malloc((size_t)n_samples * n_replicates * DIM * sizeof(float));
    assert(cache && "lp_recovery_diagnostics: out of memory for the linear responses");

    for (int v = 0; v < N_VARIANTS; v++) {
        if (only >= 0 && v != only) continue;
        const Variant *variant = &variants[v];
        fprintf(stderr, "variant %s\n", variant->label);

        mreal reference[LP_N_LOSSES][DIM];
        {
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, stand_in.sample);
            Mat stored = abm_system_read_replicate(dir, stand_in_replicate);
            int ok = fit(stored, variant, reference);
            assert(ok && "lp_recovery_diagnostics: the stand-in's own fit failed");
            mat_free(stored);
        }

        Mat values[LP_N_LOSSES], squared = mat_new(n_replicates, n_samples);
        for (int loss = 0; loss < LP_N_LOSSES; loss++) values[loss] = mat_new(n_replicates, n_samples);
        long n_failed = 0;

        #pragma omp parallel for schedule(dynamic) reduction(+:n_failed)
        for (int c = 0; c < n_samples; c++) {
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[c].name);
            mreal response[LP_N_LOSSES][DIM];
            for (int batch = 0; batch * ABM_SYSTEM_BATCH < n_replicates; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int n = abm_system_read_batch(dir, batch, block, replicate);
                for (int b = 0; b < n; b++) {
                    int r = replicate[b];
                    float *slot = cache_slot(cache, c, r);
                    if (fit(block[b], variant, response)) {
                        for (int i = 0; i < DIM; i++) slot[i] = (float)response[LP_LOSS_LIN][i];
                        for (int loss = 0; loss < LP_N_LOSSES; loss++)
                            AT(values[loss], r, c) = (mreal)mean_absolute(reference[loss], response[loss]);
                        AT(squared, r, c) = (mreal)mean_squared(reference[LP_LOSS_LIN], response[LP_LOSS_LIN]);
                    } else {
                        slot[0] = (float)NAN;
                        for (int loss = 0; loss < LP_N_LOSSES; loss++) AT(values[loss], r, c) = (mreal)NAN;
                        AT(squared, r, c) = (mreal)NAN;
                        n_failed++;
                    }
                    mat_free(block[b]);
                }
            }
        }

        fprintf(f, "\n\nVariant %s: %s\n", variant->label, variant->description);
        fprintf(f, "  runs whose fit had no responses or a non-finite one: %ld\n", n_failed);
        if (v == 0) compare_with_stored(f, values);

        fprintf(f, "\n  stand-in %s run %d\n", stand_in.sample, stand_in_replicate);
        SetOutcome linear = confidence_set(values[LP_LOSS_LIN], stand_in_replicate);
        report_set(f, "linear", linear);
        report_set(f, "state 1", confidence_set(values[LP_LOSS_S1], stand_in_replicate));
        report_set(f, "state 2", confidence_set(values[LP_LOSS_S2], stand_in_replicate));
        if (v == 0) {
            fprintf(f, "\n4. The linear model with the mean squared difference, R's loss\n");
            report_set(f, "linear, squared", confidence_set(squared, stand_in_replicate));
        }
        fprintf(f, "\n  5. Noise against signal, linear model, stand-in %s run %d\n", stand_in.sample,
                stand_in_replicate);
        noise_and_signal(f, cache, reference[LP_LOSS_LIN], variant);
        fflush(f);

        int found = 0, first = 0;
        Mat table = mat_new(n_replicates, n_samples);
        for (int b = 0; b < REPETITIONS; b++) {
            const float *ref = cache_slot(cache, stand_in_sample, b);
            if (MISNAN(ref[0])) continue;
            #pragma omp parallel for schedule(static)
            for (int c = 0; c < n_samples; c++)
                for (int r = 0; r < n_replicates; r++) {
                    const float *x = cache_slot(cache, c, r);
                    if (MISNAN(x[0])) { AT(table, r, c) = (mreal)NAN; continue; }
                    double s = 0;
                    for (int i = 0; i < DIM; i++) s += fabs((double)x[i] - (double)ref[i]);
                    AT(table, r, c) = (mreal)(s / DIM);
                }
            SetOutcome o = confidence_set(table, b);
            found += o.in_set;
            first += o.rank == 1;
            fprintf(stderr, "  %s repetition %d: in set %d, rank %d\n", variant->label, b, o.in_set, o.rank);
        }
        mat_free(table);
        fprintf(f, "\n  linear model, runs 0 to %d of %s as the stand-in: cop_0191 in the set %d times of %d, "
                   "smallest mean loss %d times\n", REPETITIONS - 1, stand_in.sample, found, REPETITIONS, first);
        fflush(f);

        if (v == 0) volatility(f, linear.winner);
        for (int loss = 0; loss < LP_N_LOSSES; loss++) mat_free(values[loss]);
        mat_free(squared);
    }

    fclose(f);
    free(cache);
    lp_free_samples(samples, n_samples);
    return 0;
}
