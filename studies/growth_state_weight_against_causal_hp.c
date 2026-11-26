/*
How close two recession weights are, on every simulated run and on the US
data. A recession weight F_t lies in (0, 1), near 1 when GDP is low and near
0 when it is high, and splits each period between the two states.

The NLP weight is the one the state-dependent local projection uses:

    s_t = (1/4) sum_{j=0..3} 100 log(100 + GDP_growth_1 + ... + GDP_growth_(t-j))
    c_t = causal HP cycle of s at t, from s_1..s_t only, at lambda
    z_t = (c_t - m_t) / sd_t, m_t and sd_t the mean and standard deviation of c_1..c_t
    F_t = exp(-2 z_(t-1)) / (1 + exp(-2 z_(t-1)))

where lambda is the maximum likelihood estimate on that run's own s. The
estimate is read back from the local projection fits, not estimated again,
so the weight is exactly the one those fits used.

The growth weight is the candidate for the state-conditional t-QVARMA
responses, the same logistic and lag without the HP filter:

    g_t = (GDP_growth_(t-3) + GDP_growth_(t-2) + GDP_growth_(t-1) + GDP_growth_t) / 4
    z_t = (g_t - mean(g)) / sd(g), mean and sd (divided by n - 1) over all periods
    F_t = exp(-2 z_(t-1)) / (1 + exp(-2 z_(t-1)))

Both start at the fourth period, the first where the 4-period averages exist:
397 of a simulated run's 400 periods, and 184 US quarters from 1974Q1. Lagged,
both are defined from the period after, so the comparison covers 396 periods of
a run and 183 US quarters. Per run:

  correlation    Pearson correlation of the two weights over the periods
  mean_abs_diff  mean over the periods of |growth weight - NLP weight|
  agreement      share of periods both weights put on the same side of 0.5,
                 both calling it recession or both expansion
  mean_growth, mean_nlp
                 the mean of each weight, the share of the sample it puts in
                 recession

A run is left out, and counted, when its local projection fit is missing or
was computed on other data, when s or lambda is not finite, or when either
weight is constant, which leaves the correlation undefined.

Writes the report (quantiles of each statistic and of lambda over all runs,
and the US values), a table with the median of each statistic per
configuration, and the US weights quarter by quarter.
*/

#include "applications/abm_system_lp.h"
#include "applications/us_data.h"
#include <et_al./stats.h>
#include <omp.h>
#include <time.h>

#define NLP_FIT_DIR "out/applications/abm_system_fit_lp_causal_hp_estimated_lambda"
#define US_NLP_FIT "out/applications/us_lp_fit_causal_hp_estimated_lambda.npz"
#define INPUT_DIR "out/dataset/abm_system"
#define REPORT_PATH "out/studies/growth_state_weight_against_causal_hp_report.txt"
#define BY_CONFIGURATION_PATH "out/studies/growth_state_weight_against_causal_hp_by_configuration.csv"
#define US_PATH "out/studies/growth_state_weight_against_causal_hp_us.csv"

enum { CORRELATION, MEAN_ABS_DIFF, AGREEMENT, MEAN_GROWTH, MEAN_NLP, N_STATISTICS };
static const char *const statistic_name[N_STATISTICS] = {
    "correlation", "mean_abs_diff", "agreement", "mean_growth", "mean_nlp"
};

enum { EXCLUDED_NO_FIT_ROW, EXCLUDED_NOT_FINITE, EXCLUDED_CONSTANT_WEIGHT, N_EXCLUSIONS };
static const char *const exclusion_name[N_EXCLUSIONS] = {
    "no fit row matching the series", "switching series or lambda not finite",
    "a weight constant over the periods"
};

/* What a local projection fit row holds that this study reads: the run, the
   fingerprints of the series it was fitted on, and its estimated lambda. */
typedef struct {
    int replicate;
    double data_fingerprint, switching_fingerprint, lambda;
} LambdaRow;

/* The rows of the fit archive at path, or NULL when it is absent or was not
   fitted with the NLP weight: causal HP cycle, estimated lambda, logistic of
   slope spec.gamma lagged spec.lag_switching periods. Only the columns used are
   read. Caller frees. */
static LambdaRow *read_lambda_rows(const char *path, LpNlSpec spec, int *count) {
    *count = 0;
    FILE *probe = fopen(path, "rb");
    if (!probe) return NULL;
    fclose(probe);
    static const char *const names[] = {
        "replicate", "data_fingerprint", "switching_fingerprint", "spec_lambda",
        "spec_use_logistic", "spec_use_hp", "spec_gamma", "spec_lag_switching", "spec_causal_hp", "spec_estimate_lambda"
    };
    int n_names = (int)(sizeof names / sizeof names[0]);
    DataFrame frame = df_read_npz_columns(path, names, n_names);
    LambdaRow *rows = NULL;
    if (frame.n_cols == n_names && frame.r >= 1) {
        Mat column[sizeof names / sizeof names[0]];
        for (int j = 0; j < n_names; j++) column[j] = df_col_numeric(&frame, names[j]);
        int same_weight = 1;
        for (int i = 0; i < frame.r && same_weight; i++)
            same_weight = AT(column[4], i, 0) == 1 && AT(column[5], i, 0) == 1 && AT(column[6], i, 0) == spec.gamma
                       && AT(column[7], i, 0) == spec.lag_switching && AT(column[8], i, 0) == 1
                       && AT(column[9], i, 0) == 1;
        if (same_weight) {
            rows = malloc((size_t)frame.r * sizeof(LambdaRow));
            assert(rows && "growth_state_weight_against_causal_hp: out of memory");
            for (int i = 0; i < frame.r; i++)
                rows[i] = (LambdaRow){ (int)AT(column[0], i, 0), (double)AT(column[1], i, 0),
                                       (double)AT(column[2], i, 0), (double)AT(column[3], i, 0) };
            *count = frame.r;
        }
    }
    df_free(&frame);
    return rows;
}

static int row_matches(const LambdaRow *row, Mat y, Mat switching) {
    return row->data_fingerprint == lp_data_fingerprint(y) && row->switching_fingerprint == lp_data_fingerprint(switching);
}

/* z_t of the growth weight, n x 1, from the fourth period on. Caller must
   mat_free. */
static Mat standardised_growth_ma4(Mat stored) {
    int first = LP_SYSTEM_SWITCHING_WINDOW - 1;
    int n = lp_system_n_periods(stored.c);
    Mat growth = mat_slice(stored, STORED_GDP_GROWTH, STORED_GDP_GROWTH + 1, 0, stored.c);
    Mat moving_average = mat_rolling_mean(growth, LP_SYSTEM_SWITCHING_WINDOW, 1);
    Mat z = mat_new(n, 1);
    for (int t = 0; t < n; t++) z.d[t] = AT(moving_average, 0, t + first);
    double mean = (double)stats_mean(z);
    double sd = sqrt((double)stats_var(z) * n / (n - 1));
    for (int t = 0; t < n; t++) z.d[t] = (mreal)(((double)z.d[t] - mean) / sd);
    mat_free(moving_average);
    return z;
}

/* The growth weight from its z_t. With the HP filter off, the local
   projection's own weight function applies only the logistic and the lag. */
static Mat growth_weight(Mat z, LpNlSpec spec) {
    spec.use_hp = 0;
    int unused;
    return _lp_weights(z, &spec, &unused);
}

/* The NLP weight at the run's estimated lambda. */
static Mat nlp_weight(Mat switching, LpNlSpec spec, double lambda) {
    spec.lambda = lambda;
    int unused;
    return _lp_causal_weights(switching, &spec, &unused);
}

static int all_finite(Mat x) {
    for (int i = 0; i < x.r * x.c; i++)
        if (MISNAN(x.d[i]) || MISINF(x.d[i])) return 0;
    return 1;
}

/* The statistics over the periods where both lagged weights are defined.
   Returns 0 when either weight is constant there. */
static int compare_weights(Mat growth, Mat nlp, double *statistic) {
    int first = 1, n = growth.r - first;
    Mat a = mat_slice(growth, first, growth.r, 0, 1);
    Mat b = mat_slice(nlp, first, nlp.r, 0, 1);
    if (stats_var(a) == 0 || stats_var(b) == 0) return 0;
    double abs_diff = 0, agree = 0;
    for (int t = 0; t < n; t++) {
        double x = (double)AT(a, t, 0), y = (double)AT(b, t, 0);
        abs_diff += fabs(x - y);
        agree += (x > 0.5) == (y > 0.5);
    }
    statistic[CORRELATION] = (double)stats_corr(a, b);
    statistic[MEAN_ABS_DIFF] = abs_diff / n;
    statistic[AGREEMENT] = agree / n;
    statistic[MEAN_GROWTH] = (double)stats_mean(a);
    statistic[MEAN_NLP] = (double)stats_mean(b);
    return 1;
}

/* Median of the finite entries of values, which is reordered; NAN when none. */
static double finite_median(mreal *values, int n) {
    int kept = 0;
    for (int i = 0; i < n; i++)
        if (!MISNAN(values[i])) values[kept++] = values[i];
    return kept ? (double)stats_quantile_inplace(values, kept, (mreal)0.5) : NAN;
}

int main(void) {
    openblas_set_num_threads(1);
    double started = omp_get_wtime();
    LpNlSpec spec = lp_system_spec();

    int n_configurations;
    LpSample *samples = lp_list_samples(NLP_FIT_DIR, &n_configurations);
    assert(n_configurations > 0 && "growth_state_weight_against_causal_hp: no fits of the NLP weight");
    int n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);

    char first_dir[560];
    snprintf(first_dir, sizeof first_dir, "%s/%s", INPUT_DIR, samples[0].name);
    Mat first_replicate = abm_system_read_replicate(first_dir, 0);
    int n_periods = lp_system_n_periods(first_replicate.c);
    mat_free(first_replicate);

    size_t n_slots = (size_t)n_configurations * n_replicates;
    /* statistic[s * n_slots + slot], NAN where the replicate is left out. */
    mreal *statistic = malloc((size_t)N_STATISTICS * n_slots * sizeof(mreal));
    mreal *lambda = malloc(n_slots * sizeof(mreal));
    int *excluded = calloc((size_t)n_configurations * N_EXCLUSIONS, sizeof(int));
    assert(statistic && lambda && excluded && "growth_state_weight_against_causal_hp: out of memory");
    for (size_t i = 0; i < (size_t)N_STATISTICS * n_slots; i++) statistic[i] = (mreal)NAN;
    for (size_t i = 0; i < n_slots; i++) lambda[i] = (mreal)NAN;

    #pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < n_configurations; c++) {
        char path[640], dir[560];
        snprintf(path, sizeof path, "%s/%s.npz", NLP_FIT_DIR, samples[c].name);
        snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[c].name);
        int n_records = 0;
        LambdaRow *records = read_lambda_rows(path, spec, &n_records);
        unsigned char seen[n_replicates];
        memset(seen, 0, sizeof seen);
        int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
        for (int batch = 0; batch < n_batches; batch++) {
            Mat block[ABM_SYSTEM_BATCH];
            int replicate[ABM_SYSTEM_BATCH];
            int count = abm_system_read_batch(dir, batch, block, replicate);
            for (int b = 0; b < count; b++) {
                int index = replicate[b];
                assert(index >= 0 && index < n_replicates && "growth_state_weight_against_causal_hp: a replicate index past the count");
                Mat y, switching;
                lp_system_from_stored(block[b], &y, &switching);
                const LambdaRow *found = NULL;
                for (int i = 0; i < n_records && !found; i++) {
                    int at = (index + i) % n_records;
                    if (records[at].replicate == index && row_matches(&records[at], y, switching))
                        found = &records[at];
                }
                seen[index] = 1;
                size_t slot = (size_t)c * n_replicates + index;
                int *out = &excluded[c * N_EXCLUSIONS];
                if (!found) {
                    out[EXCLUDED_NO_FIT_ROW]++;
                } else if (!all_finite(switching) || MISNAN(found->lambda) || MISINF(found->lambda)) {
                    out[EXCLUDED_NOT_FINITE]++;
                } else {
                    lambda[slot] = (mreal)found->lambda;
                    Mat z = standardised_growth_ma4(block[b]);
                    Mat growth = growth_weight(z, spec);
                    Mat nlp = nlp_weight(switching, spec, found->lambda);
                    double value[N_STATISTICS];
                    if (compare_weights(growth, nlp, value))
                        for (int s = 0; s < N_STATISTICS; s++) statistic[(size_t)s * n_slots + slot] = (mreal)value[s];
                    else
                        out[EXCLUDED_CONSTANT_WEIGHT]++;
                    mat_free(nlp);
                    mat_free(growth);
                    mat_free(z);
                }
                mat_free(y);
                mat_free(switching);
                mat_free(block[b]);
            }
        }
        for (int i = 0; i < n_replicates; i++) excluded[c * N_EXCLUSIONS + EXCLUDED_NO_FIT_ROW] += !seen[i];
        free(records);
    }

    FILE *by_configuration = fopen(BY_CONFIGURATION_PATH, "w");
    assert(by_configuration && "growth_state_weight_against_causal_hp: cannot open the configuration table");
    fprintf(by_configuration, "configuration,replicates_compared");
    for (int s = 0; s < N_STATISTICS; s++) fprintf(by_configuration, ",median_%s", statistic_name[s]);
    fprintf(by_configuration, ",median_lambda\n");
    mreal *column = malloc((size_t)n_replicates * sizeof(mreal));
    assert(column && "growth_state_weight_against_causal_hp: out of memory");
    for (int c = 0; c < n_configurations; c++) {
        int compared = 0;
        for (int i = 0; i < n_replicates; i++) compared += !MISNAN(statistic[(size_t)c * n_replicates + i]);
        fprintf(by_configuration, "%s,%d", samples[c].name, compared);
        for (int s = 0; s < N_STATISTICS; s++) {
            memcpy(column, statistic + (size_t)s * n_slots + (size_t)c * n_replicates, (size_t)n_replicates * sizeof(mreal));
            fprintf(by_configuration, ",%.6g", finite_median(column, n_replicates));
        }
        memcpy(column, lambda + (size_t)c * n_replicates, (size_t)n_replicates * sizeof(mreal));
        fprintf(by_configuration, ",%.6g\n", finite_median(column, n_replicates));
    }
    fclose(by_configuration);
    free(column);

    /* The US data, with the same weights. */
    Mat original = load_us_system();
    Mat us_stored = us_system_growth_block(original);
    Mat us_y, us_switching;
    lp_system_from_stored(us_stored, &us_y, &us_switching);
    int us_count = 0;
    LambdaRow *us_record = read_lambda_rows(US_NLP_FIT, spec, &us_count);
    assert(us_record && us_count == 1 && row_matches(&us_record[0], us_y, us_switching)
           && "growth_state_weight_against_causal_hp: the US fit of the NLP weight is missing or was fitted on other data");
    double us_lambda = us_record[0].lambda;
    Mat us_z = standardised_growth_ma4(us_stored);
    Mat us_growth = growth_weight(us_z, spec);
    Mat us_nlp = nlp_weight(us_switching, spec, us_lambda);
    Mat us_cycle = mat_hp_cycle_causal(us_switching, us_lambda, 0);
    double us_statistic[N_STATISTICS];
    int us_compared = compare_weights(us_growth, us_nlp, us_statistic);

    FILE *us = fopen(US_PATH, "w");
    assert(us && "growth_state_weight_against_causal_hp: cannot open the US table");
    fprintf(us, "quarter,standardised_growth_ma4,growth_weight,causal_hp_cycle,nlp_weight\n");
    for (int t = 0; t < us_switching.r; t++) {
        int quarter = 4 * 1974 + t;
        fprintf(us, "%dQ%d,%.10g,%.10g,%.10g,%.10g\n", quarter / 4, quarter % 4 + 1, (double)us_z.d[t],
                (double)us_growth.d[t], (double)us_cycle.d[t], (double)us_nlp.d[t]);
    }
    fclose(us);

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "growth_state_weight_against_causal_hp: cannot open the report");
    fprintf(report, "Recession weight of the 4-period moving average of GDP growth against the NLP's\n\n");
    fprintf(report, "growth weight: exp(-2 z) / (1 + exp(-2 z)) of the previous period, z the 4-period moving average "
                    "of GDP growth standardised with its mean and standard deviation over the run\n");
    fprintf(report, "NLP weight: the same logistic of the previous period, z the causal HP cycle of the 4-period "
                    "moving average of 100 log(100 + cumsum(GDP_growth)), at the lambda estimated by maximum "
                    "likelihood on that run, standardised with the mean and standard deviation up to each period\n");
    fprintf(report, "periods compared: %d per simulated run (periods 5 to %d), %d US quarters (1974Q2 to 2019Q4)\n",
            n_periods - 1, n_periods + LP_SYSTEM_SWITCHING_WINDOW - 1, us_switching.r - 1);
    fprintf(report, "estimated lambdas read from the local projection fits in %s and %s\n\n", NLP_FIT_DIR, US_NLP_FIT);

    int total_excluded[N_EXCLUSIONS] = { 0 };
    for (int c = 0; c < n_configurations; c++)
        for (int e = 0; e < N_EXCLUSIONS; e++) total_excluded[e] += excluded[c * N_EXCLUSIONS + e];
    size_t compared = 0;
    for (size_t i = 0; i < n_slots; i++) compared += !MISNAN(statistic[i]);
    fprintf(report, "simulated: %d configurations x %d replicates, compared %zu\n", n_configurations, n_replicates,
            compared);
    for (int e = 0; e < N_EXCLUSIONS; e++) fprintf(report, "  left out, %s: %d\n", exclusion_name[e], total_excluded[e]);
    fprintf(report, "\n");

    static const double probability[] = { 0, 0.05, 0.25, 0.5, 0.75, 0.95, 1 };
    int n_probabilities = (int)(sizeof probability / sizeof probability[0]);
    fprintf(report, "quantiles over the compared replicates:");
    for (int p = 0; p < n_probabilities; p++) fprintf(report, " %g", probability[p]);
    fprintf(report, ", then the mean\n");
    mreal *pooled = malloc(n_slots * sizeof(mreal));
    assert(pooled && "growth_state_weight_against_causal_hp: out of memory");
    for (int s = 0; s <= N_STATISTICS; s++) {
        const mreal *source = s < N_STATISTICS ? statistic + (size_t)s * n_slots : lambda;
        int kept = 0;
        double sum = 0;
        for (size_t i = 0; i < n_slots; i++)
            if (!MISNAN(source[i]) && !MISNAN(statistic[i])) {
                pooled[kept++] = source[i];
                sum += (double)source[i];
            }
        fprintf(report, "  %s:", s < N_STATISTICS ? statistic_name[s] : "lambda");
        for (int p = 0; p < n_probabilities && kept; p++)
            fprintf(report, " %.4g", (double)stats_quantile_inplace(pooled, kept, (mreal)probability[p]));
        fprintf(report, "; mean %.4g\n", kept ? sum / kept : NAN);
    }
    free(pooled);

    fprintf(report, "\nUS data, estimated lambda %.6g:\n", us_lambda);
    if (us_compared)
        for (int s = 0; s < N_STATISTICS; s++) fprintf(report, "  %s: %.4f\n", statistic_name[s], us_statistic[s]);
    else
        fprintf(report, "  not compared, a weight is constant over the quarters\n");
    fprintf(report, "\nwall time %.0f s on %d threads\n", omp_get_wtime() - started, omp_get_max_threads());
    fclose(report);

    mat_free(us_cycle);
    mat_free(us_nlp);
    mat_free(us_growth);
    mat_free(us_z);
    free(us_record);
    mat_free(us_y);
    mat_free(us_switching);
    mat_free(us_stored);
    mat_free(original);
    free(excluded);
    free(lambda);
    free(statistic);
    lp_free_samples(samples, n_configurations);
    return 0;
}
