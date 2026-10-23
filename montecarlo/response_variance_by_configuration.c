/*
How noisy each configuration's impulse responses are, and whether that noise
decides how often the confidence set holds the configuration, under the Seri
grid and the naive-loss grid.

A configuration's noise is its mean response variance, (1/K) sum_h sigma2_ih,
the term montecarlo/recovery_loss.h subtracts: sigma2_ih is the variance of
entry h across all 1000 runs of configuration i, divisor n - 1, and K the
entries per run. No run is held out, since this describes the configuration
rather than scoring a benchmark.

Reads the response caches of montecarlo/response_cache.h and
out/montecarlo/seri_against_naive_by_configuration_identifiability.csv, which
bin/seri_against_naive_by_configuration writes. Writes
out/montecarlo/response_variance_by_configuration.csv, one row per model and
configuration, and out/montecarlo/response_variance_by_configuration_report.txt:
per model, where cop_0191 ranks by noise, and per grid the Spearman correlation
across configurations between noise and the own rate (share of its own sets
holding it) and the off-target rate (share of the other configurations' sets
holding it).
*/

/* getline is POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/recovery_loss.h"
#include <stdio.h>
#include <string.h>

#define IDENTIFIABILITY_PATH "out/montecarlo/seri_against_naive_by_configuration_identifiability.csv"
#define VARIANCE_PATH "out/montecarlo/response_variance_by_configuration.csv"
#define REPORT_PATH "out/montecarlo/response_variance_by_configuration_report.txt"
#define MAX_COP 1001
#define FOCUS_COP 191
#define N_GRIDS 2

static const char *const grid_label[N_GRIDS] = { "Seri", "naive" };

/* own and off-target percentages, indexed [grid][model][cop], NaN when absent. */
static double own_percent[N_GRIDS][N_MODELS][MAX_COP], off_percent[N_GRIDS][N_MODELS][MAX_COP];

static void read_rates(void) {
    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++)
            for (int cop = 0; cop < MAX_COP; cop++) own_percent[grid][model][cop] = off_percent[grid][model][cop] = NAN;
    FILE *f = fopen(IDENTIFIABILITY_PATH, "r");
    assert(f && "response_variance_by_configuration: run bin/seri_against_naive_by_configuration first");
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, f) > 0) {
        if (header) { header = 0; continue; }
        char grid_name[16], model_name[16];
        int cop;
        long runs;
        double own;
        if (sscanf(line, "%15[^,],%15[^,],%d,%ld,%lf,", grid_name, model_name, &cop, &runs, &own) != 5) continue;
        const char *last_comma = strrchr(line, ',');
        int grid = -1, model = -1;
        for (int g = 0; g < N_GRIDS; g++) if (strcmp(grid_name, grid_label[g]) == 0) grid = g;
        for (int m = 0; m < N_MODELS; m++) if (strcmp(model_name, model_label[m]) == 0) model = m;
        if (grid < 0 || model < 0 || cop < 0 || cop >= MAX_COP || !last_comma) continue;
        own_percent[grid][model][cop] = own;
        off_percent[grid][model][cop] = strtod(last_comma + 1, NULL);
    }
    free(line);
    fclose(f);
}

int main(void) {
    int n_samples;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    int n_threads = omp_get_num_procs();
    read_rates();

    int *kept = malloc((size_t)n_replicates * sizeof(int));
    double *mean_variance = malloc((size_t)n_samples * sizeof(double));
    assert(kept && mean_variance);
    for (int run = 0; run < n_replicates; run++) kept[run] = run;

    FILE *table = fopen(VARIANCE_PATH, "w");
    FILE *report = fopen(REPORT_PATH, "w");
    assert(table && report && "response_variance_by_configuration: cannot open an output file");
    fprintf(table, "model,cop,mean_response_variance\n");
    fprintf(report, "Noise: a configuration's mean response variance over its %d runs, averaged over the entries "
                    "of the response vector. Rank 1 is the least noisy of the %d configurations.\n", n_replicates,
            n_samples);
    fprintf(report, "Own rate: share of a configuration's own sets holding it. Off-target rate: share of the other "
                    "configurations' sets holding it. Correlations are Spearman, across configurations.\n\n");
    fprintf(report, "model, cop_%04d noise, median noise, cop_%04d noise rank, noise against own rate Seri, naive, "
                    "noise against off-target rate Seri, naive\n", FOCUS_COP, FOCUS_COP);

    for (int model = 0; model < N_MODELS; model++) {
        ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates,
                                                           LP_FIT_DIR_DEFAULT,
                                                           LP_INPUT_DIR_DEFAULT,
                                                           RESPONSE_CACHE_DIR_DEFAULT);
        assert(cache.n_missing == 0 && "response_variance_by_configuration: a run has no response");
        recovery_loss_mean_variance(&cache, kept, n_replicates, mean_variance, n_threads);
        response_cache_free(&cache);

        int focus = -1, focus_rank = 1;
        for (int sample = 0; sample < n_samples; sample++) {
            fprintf(table, "%s,%d,%.9g\n", model_label[model], samples[sample].index, mean_variance[sample]);
            if (samples[sample].index == FOCUS_COP) focus = sample;
        }
        assert(focus >= 0);
        for (int sample = 0; sample < n_samples; sample++) focus_rank += mean_variance[sample] < mean_variance[focus];
        Mat noise = mat_new(n_samples, 1);
        for (int sample = 0; sample < n_samples; sample++) AT(noise, sample, 0) = mean_variance[sample];
        mreal median = stats_median(noise);

        double rho[2][N_GRIDS];
        for (int grid = 0; grid < N_GRIDS; grid++) {
            Mat own = mat_new(n_samples, 1), off = mat_new(n_samples, 1);
            for (int sample = 0; sample < n_samples; sample++) {
                AT(own, sample, 0) = own_percent[grid][model][samples[sample].index];
                AT(off, sample, 0) = off_percent[grid][model][samples[sample].index];
            }
            rho[0][grid] = stats_spearman(noise, own);
            rho[1][grid] = stats_spearman(noise, off);
            mat_free(own);
            mat_free(off);
        }
        mat_free(noise);
        fprintf(report, "%s, %.4g, %.4g, %d of %d, %.3f, %.3f, %.3f, %.3f\n", model_label[model],
                mean_variance[focus], (double)median, focus_rank, n_samples, rho[0][0], rho[0][1], rho[1][0],
                rho[1][1]);
        fflush(report);
    }
    fclose(table);
    fclose(report);
    free(kept);
    free(mean_variance);
    lp_free_samples(samples, n_samples);
    return 0;
}
