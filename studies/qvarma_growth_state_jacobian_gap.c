/*
How far apart D_recession and D_expansion are, the two conditional means of the
t-QVARMA's score Jacobian D_t that the growth-state responses are computed at,
and whether that gap is large enough to show in the responses the confidence
set compares.

The states are those of the growth-state responses: period t is in recession
when the mean of GDP growth over periods t-4..t-1 is at or below its 20th
percentile over the run, in expansion otherwise, and periods 1 to 4 have no
state. Every run of every configuration, its cached fit, its stored series.

Per run, with D_pooled the mean of D_t over the periods with a state:

  relative_gap     ||D_recession - D_expansion|| / ||D_pooled||, Frobenius norms
  trace_ratio      trace(D_recession) / trace(D_expansion) - 1
  trace_t          Welch t statistic of the difference between the mean of
                   trace(D_t) / K over recession periods and over expansion
                   periods; near 0 when the states do not differ, beyond 1.96
                   in about 5% of runs when they do not differ at all
  off_diagonal     ||off-diagonal part of D_pooled|| / ||D_pooled||
  nu               the fitted degrees of freedom

Per configuration, from the stored responses (expansion followed by recession
for the growth states, the responses at the mean of D_t for the plain
t-QVARMA), all as means over the 525 response entries:

  state_gap        mean over runs of (response at D_recession - response at D_expansion)^2
  within_variance  variance over the configuration's runs of the plain t-QVARMA response
  ratio            state_gap / within_variance: how large the change from one
                   state to the other is against the run-to-run noise the loss sees

and, over all configurations, between_variance, the variance across
configurations of their mean plain response, mean over entries.

A run is left out, and counted, when its fit is missing, nu <= 2, or its growth
state is undefined.

Writes the report (quantiles of each statistic) and a table with the median of
each statistic per configuration.
*/

#include "montecarlo/response_cache.h"
#include <et_al./stats.h>
#include <omp.h>

#define REPORT_PATH "out/studies/qvarma_growth_state_jacobian_gap_report.txt"
#define BY_CONFIGURATION_PATH "out/studies/qvarma_growth_state_jacobian_gap_by_configuration.csv"

enum { RELATIVE_GAP, TRACE_RATIO, TRACE_T, OFF_DIAGONAL, NU, N_STATISTICS };
static const char *const statistic_name[N_STATISTICS] = {
    "relative_gap", "trace_ratio", "trace_t", "off_diagonal", "nu"
};

static double frobenius(const mreal *a, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) s += (double)a[i] * a[i];
    return sqrt(s);
}

/* The statistics of one run; returns 0 when the run has no growth state or a
   state has fewer than two periods. */
static int run_statistics(const QvarmaParams *m, Mat y, double *statistic) {
    const int K = m->K, T = y.c, KK = K * K;
    signed char *state = malloc((size_t)T);
    assert(state && "qvarma_growth_state_jacobian_gap: out of memory");
    if (!response_qvarma_growth_states(y, state)) {
        free(state);
        return 0;
    }
    Mat *D_path = qvarma_score_jacobian(m, y);
    double sum[2][64] = { { 0 } }, trace_sum[2] = { 0, 0 }, trace_square[2] = { 0, 0 };
    int count[2] = { 0, 0 };
    assert(KK <= 64 && "qvarma_growth_state_jacobian_gap: K larger than the buffers");
    for (int t = 0; t < T; t++) {
        if (state[t] < 0) continue;
        int s = state[t];
        double trace = 0;
        for (int i = 0; i < KK; i++) sum[s][i] += (double)D_path[t].d[i];
        for (int k = 0; k < K; k++) trace += (double)AT(D_path[t], k, k);
        trace /= K;
        trace_sum[s] += trace;
        trace_square[s] += trace * trace;
        count[s]++;
    }
    qvarma_free_mat_array(D_path, T);
    free(state);
    if (count[0] < 2 || count[1] < 2) return 0;

    mreal recession[64], expansion[64], pooled[64], gap[64], off[64];
    for (int i = 0; i < KK; i++) {
        recession[i] = (mreal)(sum[1][i] / count[1]);
        expansion[i] = (mreal)(sum[0][i] / count[0]);
        pooled[i] = (mreal)((sum[0][i] + sum[1][i]) / (count[0] + count[1]));
        gap[i] = recession[i] - expansion[i];
        off[i] = i % (K + 1) == 0 ? 0 : pooled[i];
    }
    double mean[2], variance[2], trace_recession = 0, trace_expansion = 0;
    for (int s = 0; s < 2; s++) {
        mean[s] = trace_sum[s] / count[s];
        variance[s] = (trace_square[s] - count[s] * mean[s] * mean[s]) / (count[s] - 1);
    }
    for (int k = 0; k < K; k++) {
        trace_recession += (double)recession[k * (K + 1)];
        trace_expansion += (double)expansion[k * (K + 1)];
    }
    statistic[RELATIVE_GAP] = frobenius(gap, KK) / frobenius(pooled, KK);
    statistic[TRACE_RATIO] = trace_recession / trace_expansion - 1;
    statistic[TRACE_T] = (mean[1] - mean[0]) / sqrt(variance[1] / count[1] + variance[0] / count[0]);
    statistic[OFF_DIAGONAL] = frobenius(off, KK) / frobenius(pooled, KK);
    statistic[NU] = (double)m->nu;
    return 1;
}

static void write_quantiles(FILE *out, const char *name, mreal *values, size_t n) {
    static const double probability[] = { 0, 0.05, 0.25, 0.5, 0.75, 0.95, 1 };
    size_t kept = 0;
    double sum = 0;
    for (size_t i = 0; i < n; i++)
        if (!MISNAN(values[i])) {
            values[kept++] = values[i];
            sum += (double)values[i];
        }
    fprintf(out, "  %s:", name);
    for (size_t p = 0; p < sizeof probability / sizeof probability[0] && kept; p++)
        fprintf(out, " %.4g", (double)stats_quantile_inplace(values, (int)kept, (mreal)probability[p]));
    fprintf(out, "; mean %.4g (n %zu)\n", kept ? sum / kept : NAN, kept);
}

int main(void) {
    openblas_set_num_threads(1);
    omp_set_max_active_levels(1);
    double started = omp_get_wtime();
    int n_samples = 0;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    size_t n_slots = (size_t)n_samples * n_replicates;

    /* statistic[s * n_slots + slot], NAN where the run is left out. */
    mreal *statistic = malloc((size_t)N_STATISTICS * n_slots * sizeof(mreal));
    assert(statistic && "qvarma_growth_state_jacobian_gap: out of memory");
    for (size_t i = 0; i < (size_t)N_STATISTICS * n_slots; i++) statistic[i] = (mreal)NAN;
    long n_left_out = 0;

    #pragma omp parallel reduction(+:n_left_out)
    {
        QvarmaParams working = _response_qvarma_shape();
        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < n_samples; sample++) {
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", LP_INPUT_DIR_DEFAULT, samples[sample].name);
            int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
            int seen = 0;
            for (int batch = 0; batch < n_batches; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);
                for (int b = 0; b < count; b++) {
                    if (replicate[b] < n_replicates) {
                        char fit_path[640];
                        snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json",
                                 RESPONSE_QVARMA_FIT_DIR, samples[sample].name, replicate[b],
                                 RESPONSE_QVARMA_SPEC_LABEL);
                        double value[N_STATISTICS];
                        int ok = qvarma_load_params(&working, fit_path) && working.nu > 2
                              && run_statistics(&working, block[b], value);
                        size_t slot = (size_t)sample * n_replicates + replicate[b];
                        if (ok)
                            for (int s = 0; s < N_STATISTICS; s++)
                                statistic[(size_t)s * n_slots + slot] = (mreal)value[s];
                        else
                            n_left_out++;
                        seen++;
                    }
                    mat_free(block[b]);
                }
            }
            n_left_out += n_replicates - seen;
        }
        qvarma_params_free(&working);
    }

    /* The responses the confidence set compares, from the two caches. */
    ResponseCache plain = response_cache_load_or_build(MODEL_QVARMA, samples, n_samples, n_replicates,
                                                       LP_FIT_DIR_DEFAULT, LP_INPUT_DIR_DEFAULT,
                                                       RESPONSE_CACHE_DIR_DEFAULT);
    ResponseCache states = response_cache_load_or_build(MODEL_QVARMA_GROWTH_STATES, samples, n_samples,
                                                        n_replicates, LP_FIT_DIR_DEFAULT, LP_INPUT_DIR_DEFAULT,
                                                        RESPONSE_CACHE_DIR_DEFAULT);
    const int dim = RESPONSE_QVARMA_DIM;
    double *state_gap = malloc((size_t)n_samples * sizeof(double));
    double *within_variance = malloc((size_t)n_samples * sizeof(double));
    double *configuration_mean = calloc((size_t)n_samples * dim, sizeof(double));
    assert(state_gap && within_variance && configuration_mean && "qvarma_growth_state_jacobian_gap: out of memory");
    #pragma omp parallel for schedule(dynamic)
    for (int sample = 0; sample < n_samples; sample++) {
        double gap = 0, *mean = configuration_mean + (size_t)sample * dim;
        double *square = calloc((size_t)dim, sizeof(double));
        assert(square && "qvarma_growth_state_jacobian_gap: out of memory");
        int n_plain = 0, n_states = 0;
        for (int run = 0; run < n_replicates; run++) {
            const float *p = response_slot(&plain, sample, run), *s = response_slot(&states, sample, run);
            if (!MISNAN(p[0])) {
                for (int i = 0; i < dim; i++) {
                    mean[i] += p[i];
                    square[i] += (double)p[i] * p[i];
                }
                n_plain++;
            }
            if (!MISNAN(s[0])) {
                for (int i = 0; i < dim; i++) {
                    double d = (double)s[dim + i] - s[i];
                    gap += d * d;
                }
                n_states++;
            }
        }
        double variance = 0;
        for (int i = 0; i < dim; i++) {
            mean[i] /= n_plain;
            variance += (square[i] - n_plain * mean[i] * mean[i]) / (n_plain - 1);
        }
        state_gap[sample] = gap / ((double)n_states * dim);
        within_variance[sample] = variance / dim;
        free(square);
    }
    double between_variance = 0;
    for (int i = 0; i < dim; i++) {
        double mean = 0, square = 0;
        for (int sample = 0; sample < n_samples; sample++) {
            double v = configuration_mean[(size_t)sample * dim + i];
            mean += v;
            square += v * v;
        }
        mean /= n_samples;
        between_variance += (square - n_samples * mean * mean) / (n_samples - 1);
    }
    between_variance /= dim;
    response_cache_free(&plain);
    response_cache_free(&states);

    FILE *by_configuration = fopen(BY_CONFIGURATION_PATH, "w");
    assert(by_configuration && "qvarma_growth_state_jacobian_gap: cannot open the configuration table");
    fprintf(by_configuration, "configuration");
    for (int s = 0; s < N_STATISTICS; s++) fprintf(by_configuration, ",median_%s", statistic_name[s]);
    fprintf(by_configuration, ",state_gap,within_variance,ratio\n");
    mreal *column = malloc((size_t)n_replicates * sizeof(mreal));
    assert(column && "qvarma_growth_state_jacobian_gap: out of memory");
    for (int sample = 0; sample < n_samples; sample++) {
        fprintf(by_configuration, "%s", samples[sample].name);
        for (int s = 0; s < N_STATISTICS; s++) {
            int kept = 0;
            for (int run = 0; run < n_replicates; run++) {
                mreal v = statistic[(size_t)s * n_slots + (size_t)sample * n_replicates + run];
                if (!MISNAN(v)) column[kept++] = v;
            }
            fprintf(by_configuration, ",%.6g", kept ? (double)stats_quantile_inplace(column, kept, (mreal)0.5) : NAN);
        }
        fprintf(by_configuration, ",%.6g,%.6g,%.6g\n", state_gap[sample], within_variance[sample],
                state_gap[sample] / within_variance[sample]);
    }
    fclose(by_configuration);
    free(column);

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "qvarma_growth_state_jacobian_gap: cannot open the report");
    fprintf(report, "Gap between D_recession and D_expansion of the t-QVARMA growth-state responses\n\n");
    fprintf(report, "recession: mean GDP growth over periods t-4..t-1 at or below its 20th percentile over the run; "
                    "expansion otherwise; periods 1 to 4 in neither\n");
    fprintf(report, "%d configurations x %d runs, left out %ld\n\n", n_samples, n_replicates, n_left_out);
    fprintf(report, "per run, quantiles 0 0.05 0.25 0.5 0.75 0.95 1, then the mean:\n");
    mreal *pooled = malloc(n_slots * sizeof(mreal));
    assert(pooled && "qvarma_growth_state_jacobian_gap: out of memory");
    for (int s = 0; s < N_STATISTICS; s++) {
        memcpy(pooled, statistic + (size_t)s * n_slots, n_slots * sizeof(mreal));
        write_quantiles(report, statistic_name[s], pooled, n_slots);
    }
    long beyond = 0, compared = 0;
    for (size_t i = 0; i < n_slots; i++) {
        mreal t = statistic[(size_t)TRACE_T * n_slots + i];
        if (MISNAN(t)) continue;
        compared++;
        beyond += fabs((double)t) > 1.96;
    }
    fprintf(report, "  share of runs with |trace_t| > 1.96: %.4f (about 0.05 if the states do not differ)\n\n",
            (double)beyond / compared);

    fprintf(report, "per configuration, means over the 525 response entries, quantiles over configurations:\n");
    mreal *values = malloc((size_t)n_samples * sizeof(mreal));
    assert(values && "qvarma_growth_state_jacobian_gap: out of memory");
    for (int sample = 0; sample < n_samples; sample++) values[sample] = (mreal)state_gap[sample];
    write_quantiles(report, "state_gap", values, (size_t)n_samples);
    for (int sample = 0; sample < n_samples; sample++) values[sample] = (mreal)within_variance[sample];
    write_quantiles(report, "within_variance", values, (size_t)n_samples);
    for (int sample = 0; sample < n_samples; sample++) values[sample] = (mreal)(state_gap[sample] / within_variance[sample]);
    write_quantiles(report, "ratio", values, (size_t)n_samples);
    fprintf(report, "between_variance, across configurations of their mean plain response: %.4g\n", between_variance);
    fprintf(report, "\nwall time %.0f s on %d threads\n", omp_get_wtime() - started, omp_get_max_threads());
    fclose(report);

    free(values);
    free(pooled);
    free(configuration_mean);
    free(within_variance);
    free(state_gap);
    free(statistic);
    lp_free_samples(samples, n_samples);
    return 0;
}
