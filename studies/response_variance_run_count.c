/*
How many runs of a configuration are enough to estimate its response variance,
taking the runs from the last one backwards.

For configuration i and m runs, the last m of its n = 1000 (runs n-m .. n-1),

    V_i(m) = sum_h s2_ih(m),

with s2_ih(m) the sample variance of entry h of the response vector over
those m runs, divisor m - 1, and the sum over every entry of the vector
(every horizon, response and shock). V_i(n) uses all runs. The distance from
the full-sample value is

    d_i(m) = |V_i(m) - V_i(n)| / V_i(n).

The runs are taken from the last because the Seri grid with a held-out mean
(montecarlo/seri_split_true_configuration_recovery.c) estimates the mean from
the last runs and forms the losses on the first ones.

Per model and m from 2 to n, across the configurations: the median, the 95th
percentile and the largest d_i(m), and the share of configurations with
d_i(m) below TOLERANCE. A model's run count is the smallest m from which on,
for that m and every larger one, at least COVERAGE of the configurations are
below TOLERANCE. The common run count is the largest of the models' counts.

Models: those of the Monte Carlo grids (model_in_main_grid in
montecarlo/response_cache.h), the linear projection, both states under the
causal HP cycle at the estimated lambda, and the t-QVARMA with its two growth
states. Reads their response caches in out/montecarlo/sweep_grid_response_cache/,
and asserts that no run lacks a response.

Writes out/studies/response_variance_run_count_by_size.csv, one row per model
and m, out/studies/response_variance_run_count_by_configuration.csv, one row
per model and configuration holding V_i(m) for m = 2 .. n, which
studies/response_variance_run_count_plots.py draws, and
out/studies/response_variance_run_count_report.txt.
*/

#include "montecarlo/response_cache.h"
#include <et_al./stats.h>
#include <omp.h>

#define BY_SIZE_PATH "out/studies/response_variance_run_count_by_size.csv"
#define BY_CONFIGURATION_PATH "out/studies/response_variance_run_count_by_configuration.csv"
#define REPORT_PATH "out/studies/response_variance_run_count_report.txt"
#define TOLERANCE 0.05
#define COVERAGE 0.95
/* The most runs per configuration the grid may give to the variance. */
#define RUN_COUNT_BOUND 100

/* V_i(m) for every configuration i and m = 1..n, at summed_variance[i * n + m - 1];
   m = 1 has no sample variance and is left at 0. */
static void summed_variance_by_run_count(const ResponseCache *cache, double *summed_variance, int n_threads) {
    const int dim = cache->dim, n = cache->n_replicates;
    #pragma omp parallel num_threads(n_threads)
    {
        /* Row h holds entry h of runs n-1, n-2, .., 0, so its prefix of length
           m is the last m runs. */
        Mat entries = mat_new(dim, n);
        double *mean = malloc((size_t)n * sizeof(double)), *var = malloc((size_t)n * sizeof(double));
        assert(mean && var);
        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < cache->n_samples; sample++) {
            for (int position = 0; position < n; position++) {
                const float *cell = response_slot(cache, sample, n - 1 - position);
                for (int h = 0; h < dim; h++) AT(entries, h, position) = (mreal)cell[h];
            }
            double *out = summed_variance + (size_t)sample * n;
            for (int position = 0; position < n; position++) out[position] = 0;
            for (int h = 0; h < dim; h++) {
                stats_running_mean_var(mat_slice(entries, h, h + 1, 0, n), mean, var);
                /* var is the population variance of the prefix, rescaled to divisor m - 1. */
                for (int position = 1; position < n; position++)
                    out[position] += var[position] * (position + 1) / position;
            }
        }
        mat_free(entries);
        free(mean);
        free(var);
    }
}

int main(void) {
    int n_samples;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    int n_threads = omp_get_num_procs();
    mkdir("out/studies", 0755);

    double *summed_variance = malloc((size_t)n_samples * n * sizeof(double));
    double *distance = malloc((size_t)n_samples * sizeof(double));
    double *share_below = malloc((size_t)n * sizeof(double));
    assert(summed_variance && distance && share_below);

    FILE *by_size = fopen(BY_SIZE_PATH, "w");
    FILE *by_configuration = fopen(BY_CONFIGURATION_PATH, "w");
    FILE *report = fopen(REPORT_PATH, "w");
    assert(by_size && by_configuration && report && "response_variance_run_count: cannot open an output file");
    fprintf(by_size, "model,runs,median_distance,p95_distance,largest_distance,share_below_tolerance\n");
    fprintf(by_configuration, "model,cop");
    for (int m = 2; m <= n; m++) fprintf(by_configuration, ",V_%d", m);
    fprintf(by_configuration, "\n");
    fprintf(report,
            "Setup: %d configurations, %d runs each, the cached responses of every run (no refit). For "
            "configuration i, V_i(m) is the sum over every entry of the response vector of the entry's sample "
            "variance (divisor m - 1) over the last m runs, runs %d-m .. %d. d_i(m) = |V_i(m) - V_i(%d)| / "
            "V_i(%d), the relative distance from the value over all runs. No randomness, no seed.\n",
            n_samples, n, n, n - 1, n, n);
    fprintf(report,
            "Criterion: a model's run count is the smallest m such that, for m and every larger m, d_i(m) < %.2f "
            "in at least %.0f%% of the configurations. The common count is the largest over the models, capped "
            "at %d.\n\n",
            TOLERANCE, 100 * COVERAGE, RUN_COUNT_BOUND);

    int common = 0;
    for (int model = 0; model < N_MODELS; model++) {
        if (!model_in_main_grid(model)) continue;
        double started = omp_get_wtime();
        ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n, LP_FIT_DIR_DEFAULT,
                                                           LP_INPUT_DIR_DEFAULT, RESPONSE_CACHE_DIR_DEFAULT);
        assert(cache.n_missing == 0 && "response_variance_run_count: a run has no response");
        summed_variance_by_run_count(&cache, summed_variance, n_threads);
        const int dim = cache.dim;
        response_cache_free(&cache);
        for (int sample = 0; sample < n_samples; sample++) {
            fprintf(by_configuration, "%s,%d", model_label[model], samples[sample].index);
            for (int m = 2; m <= n; m++) fprintf(by_configuration, ",%.9g", summed_variance[(size_t)sample * n + m - 1]);
            fprintf(by_configuration, "\n");
        }

        double quoted[3][4];
        const int quoted_runs[3] = { 50, 100, 200 };
        for (int m = 2; m <= n; m++) {
            int below = 0;
            for (int sample = 0; sample < n_samples; sample++) {
                const double *v = summed_variance + (size_t)sample * n;
                distance[sample] = fabs(v[m - 1] - v[n - 1]) / v[n - 1];
                below += distance[sample] < TOLERANCE;
            }
            share_below[m - 1] = (double)below / n_samples;
            Mat values = mat_new(n_samples, 1);
            for (int sample = 0; sample < n_samples; sample++) AT(values, sample, 0) = distance[sample];
            double median = stats_quantile(values, 0.5), p95 = stats_quantile(values, 0.95), largest = 0;
            for (int sample = 0; sample < n_samples; sample++) if (distance[sample] > largest) largest = distance[sample];
            mat_free(values);
            fprintf(by_size, "%s,%d,%.6g,%.6g,%.6g,%.4f\n", model_label[model], m, median, p95, largest,
                    share_below[m - 1]);
            for (int q = 0; q < 3; q++)
                if (m == quoted_runs[q]) {
                    quoted[q][0] = median;
                    quoted[q][1] = p95;
                    quoted[q][2] = largest;
                    quoted[q][3] = share_below[m - 1];
                }
        }

        int run_count = n;
        while (run_count > 2 && share_below[run_count - 2] >= COVERAGE) run_count--;
        if (run_count > common) common = run_count;

        fprintf(report, "%s (%d entries per response vector): run count %d", model_label[model], dim, run_count);
        fprintf(report, "; share of configurations below %.2f at that count %.3f\n", TOLERANCE,
                share_below[run_count - 1]);
        for (int q = 0; q < 3; q++)
            fprintf(report, "  m = %d: median d %.4f, 95th percentile %.4f, largest %.4f, share below %.2f %.3f\n",
                    quoted_runs[q], quoted[q][0], quoted[q][1], quoted[q][2], TOLERANCE, quoted[q][3]);
        fprintf(report, "  %.0f s\n", omp_get_wtime() - started);
        fflush(report);
        fflush(by_size);
    }

    int chosen = common < RUN_COUNT_BOUND ? common : RUN_COUNT_BOUND;
    fprintf(report, "\nCommon run count over the models: %d. Runs given to the variance in the grid: %d%s.\n", common,
            chosen, common > RUN_COUNT_BOUND ? " (the bound, below the common count)" : "");
    fclose(by_size);
    fclose(by_configuration);
    fclose(report);
    free(summed_variance);
    free(distance);
    free(share_below);
    lp_free_samples(samples, n_samples);
    return 0;
}
