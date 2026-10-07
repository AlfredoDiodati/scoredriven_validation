/*
Whether montecarlo/fit_metric.h computes the fit metric it describes.

1. sigma2 on a hand-computed case: three vectors of two entries, (1, 2), (3, 6)
   and (5, 10), have variances 4 and 16 with divisor n - 1, so sigma2 is 10;
   and s and v at d = 2, sigma2 = 4 are 1 and 1/2, the same when d is
   multiplied by 100 and sigma2 by 100 squared, as rescaling the responses does.
2. The oracle shortcut, each run's moments taken out of the whole
   configuration's, against sigma2 computed directly over the other runs, on
   responses badly scaled on purpose (an offset of 1000 against a spread of
   0.01, where a difference of large sums would lose the digits), for the first,
   a middle and the last run.
3. d leaves the benchmark's run out: changing the best configuration's
   responses at that run to 1e6 does not move it, changing them at any other run
   does, and it equals a plain double loop over the other runs.
4. The responses of one real series, cop_0191 run 706 (the single benchmark of
   docs/MONTECARLO_VALIDATION.md): the four local-projection vectors
   fit_metric_responses builds against the ones stored in
   out/sweep_grid_response_cache/, and the t-QVARMA vector from the cached fit's
   own parameters against the stored one. The t-QVARMA refitted from those
   parameters is reported, not checked, since the solver may still move.
5. A bootstrap whose blocks are the whole series resamples the series itself
   every time, so sigma2_star is exactly 0 under every method; and two bootstraps
   with the same seed give the same sigma2_star.

Reads dataset/abm_system/, out/abm_system_fit_qvarma/ and the response caches.
Writes out/fit_metric_correctness_report.txt. Exits 1 on any failure.
*/

#include "montecarlo/fit_metric.h"
#include <stdio.h>
#include <math.h>

#define REPORT "out/fit_metric_correctness_report.txt"
#define REAL_COP 191
#define REAL_RUN 706

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) failures++;
}

static void hand_computed(FILE *report) {
    double vectors[] = { 1, 2, 3, 6, 5, 10 };
    check(report, fabs(fit_metric_sigma2(vectors, 3, 2) - 10) < 1e-12, "sigma2 of (1,2), (3,6), (5,10) is 10");
    double s = fit_metric_s(2, 4);
    check(report, fabs(s - 1) < 1e-15 && fabs(fit_metric_v(s) - 0.5) < 1e-15, "s = 1 and v = 1/2 at d = 2, sigma2 = 4");
    check(report, fabs(fit_metric_s(200, 40000) - s) < 1e-15,
          "s is unchanged when the responses are multiplied by 100");
}

/* A cache of n_samples configurations by n runs of dim entries in memory,
   offset plus spread times a uniform draw. */
static ResponseCache synthetic_cache(int n_samples, int n, int dim, double offset, double spread, Rng *rng) {
    ResponseCache cache = { MODEL_LP_LIN, n_samples, n, dim, NULL, 0, NULL, 0 };
    cache.data = malloc((size_t)n_samples * n * dim * sizeof(float));
    for (size_t i = 0; i < (size_t)n_samples * n * dim; i++)
        cache.data[i] = (float)(offset + spread * rng_uniform(rng));
    return cache;
}

static void oracle_shortcut(FILE *report) {
    Rng rng = rng_new(42, 0);
    const int n = 1000, dim = 37;
    ResponseCache cache = synthetic_cache(2, n, dim, 1000, 0.01, &rng);
    double *mean = malloc(dim * sizeof(double)), *squares = malloc(dim * sizeof(double));
    double *held_out = malloc(dim * sizeof(double)), *others = malloc((size_t)(n - 1) * dim * sizeof(double));
    fit_metric_moments(&cache, 1, mean, squares);

    int runs[] = { 0, 517, n - 1 };
    double worst = 0;
    for (int r = 0; r < 3; r++) {
        int at = 0;
        for (int run = 0; run < n; run++) {
            const float *x = response_slot(&cache, 1, run);
            for (int k = 0; k < dim; k++) {
                if (run == runs[r]) held_out[k] = (double)x[k];
                else others[(size_t)at * dim + k] = (double)x[k];
            }
            at += run != runs[r];
        }
        double direct = fit_metric_sigma2(others, n - 1, dim);
        double shortcut = fit_metric_oracle_sigma2(mean, squares, held_out, n, dim);
        double gap = fabs(shortcut - direct) / direct;
        if (gap > worst) worst = gap;
    }
    fprintf(report, "      oracle shortcut against the direct variance, largest relative gap %.2e\n", worst);
    check(report, worst < 1e-9, "the oracle shortcut equals the direct variance over the other runs");
    free(mean); free(squares); free(held_out); free(others); free(cache.data);
}

static double plain_d(const ResponseCache *cache, int best, int held_out_run, const double *benchmark) {
    double total = 0;
    for (int run = 0; run < cache->n_replicates; run++) {
        if (run == held_out_run) continue;
        double sum = 0;
        for (int k = 0; k < cache->dim; k++) sum += fabs((double)response_slot(cache, best, run)[k] - benchmark[k]);
        total += sum / cache->dim;
    }
    return total / (cache->n_replicates - 1);
}

static void d_holds_out(FILE *report) {
    Rng rng = rng_new(43, 0);
    const int n = 50, dim = 11, held = 7;
    ResponseCache cache = synthetic_cache(3, n, dim, 0, 1, &rng);
    double benchmark[11];
    for (int k = 0; k < dim; k++) benchmark[k] = (double)response_slot(&cache, 0, held)[k];
    double before = fit_metric_d(&cache, 2, held, benchmark);
    check(report, fabs(before - plain_d(&cache, 2, held, benchmark)) < 1e-12, "d equals a plain loop over the other runs");
    for (int k = 0; k < dim; k++) response_slot(&cache, 2, held)[k] = 1e6f;
    check(report, fit_metric_d(&cache, 2, held, benchmark) == before, "d does not read the benchmark's run");
    response_slot(&cache, 2, held + 1)[0] = 1e6f;
    check(report, fit_metric_d(&cache, 2, held, benchmark) > before + 1, "d reads every other run");
    free(cache.data);
}

/* The stored response of one run, read from its cache file directly. */
static int stored_response(int model, int sample, int n_samples, int n_replicates, int run, float *out) {
    char path[256];
    snprintf(path, sizeof path, "%s/%s.f32", RESPONSE_CACHE_DIR_DEFAULT, model_label[model]);
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    int header[4];
    int ok = fread(header, sizeof(int), 4, f) == 4 && header[0] == model && header[1] == n_samples
          && header[2] == n_replicates && header[3] == response_dim(model);
    long offset = (long)(4 * sizeof(int) + ((size_t)sample * n_replicates + run) * response_dim(model) * sizeof(float));
    ok = ok && fseek(f, offset, SEEK_SET) == 0
         && fread(out, sizeof(float), (size_t)response_dim(model), f) == (size_t)response_dim(model);
    fclose(f);
    return ok;
}

static double largest_relative_gap(const double *got, const float *stored, int dim) {
    double scale = 0, gap = 0;
    for (int k = 0; k < dim; k++) if (fabs((double)stored[k]) > scale) scale = fabs((double)stored[k]);
    for (int k = 0; k < dim; k++) {
        double g = fabs(got[k] - (double)stored[k]);
        if (g > gap) gap = g;
    }
    return gap / scale;
}

static void real_series(FILE *report) {
    int n_samples;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    int sample = -1;
    for (int i = 0; i < n_samples; i++) if (samples[i].index == REAL_COP) sample = i;
    check(report, sample >= 0, "cop_0191 is in the design");
    if (sample < 0) return;

    char dir[600], fit_path[700];
    snprintf(dir, sizeof dir, "%s/%s", LP_INPUT_DIR_DEFAULT, samples[sample].name);
    snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json", RESPONSE_QVARMA_FIT_DIR,
             samples[sample].name, REAL_RUN, RESPONSE_QVARMA_SPEC_LABEL);
    Mat series = abm_system_read_replicate(dir, REAL_RUN);
    QvarmaParams start = _response_qvarma_shape();
    check(report, qvarma_load_params(&start, fit_path), "the cached t-QVARMA fit of cop_0191 run 706 loads");

    double *out[N_MODELS];
    float *stored = malloc(2 * LP_RESPONSE_DIM * sizeof(float));
    for (int model = 0; model < N_MODELS; model++) out[model] = malloc((size_t)response_dim(model) * sizeof(double));
    int usable[N_MODELS], converged, iterations;
    fit_metric_responses(series, &start, qvarma_default_fit_options(), out, usable, &converged, &iterations);

    for (int model = MODEL_LP_LIN; model < N_MODELS; model++) {
        int read = stored_response(model, sample, n_samples, n_replicates, REAL_RUN, stored);
        double gap = read && usable[model] ? largest_relative_gap(out[model], stored, response_dim(model)) : INFINITY;
        char what[160];
        snprintf(what, sizeof what, "%s response of cop_0191 run 706 equals the stored one to float32 "
                 "(largest gap %.1e of the largest entry)", model_label[model], gap);
        check(report, gap < 1e-6, what);
    }

    float qvarma_from_cache[RESPONSE_QVARMA_DIM];
    double from_cache[RESPONSE_QVARMA_DIM];
    int read = stored_response(MODEL_QVARMA, sample, n_samples, n_replicates, REAL_RUN, stored = realloc(stored,
                               RESPONSE_QVARMA_DIM * sizeof(float)));
    int computed = _response_qvarma_cell(&start, series, qvarma_from_cache);
    for (int k = 0; k < RESPONSE_QVARMA_DIM; k++) from_cache[k] = (double)qvarma_from_cache[k];
    double gap = read && computed ? largest_relative_gap(from_cache, stored, RESPONSE_QVARMA_DIM) : INFINITY;
    char what[160];
    snprintf(what, sizeof what, "qvarma response from the cached parameters equals the stored one "
             "(largest gap %.1e of the largest entry)", gap);
    check(report, gap < 1e-6, what);
    fprintf(report, "      qvarma refitted from the cached parameters: converged %d after %d iterations, "
            "largest gap to the stored response %.1e of the largest entry\n", converged, iterations,
            usable[MODEL_QVARMA] ? largest_relative_gap(out[MODEL_QVARMA], stored, RESPONSE_QVARMA_DIM) : INFINITY);

    FitMetricBootstrapOptions options = fit_metric_bootstrap_options_default();
    options.n_bootstrap = 3;
    options.block_length = series.c;
    Rng rng = rng_new(123, 1);
    FitMetricBootstrap whole = fit_metric_bootstrap(series, &start, options, &rng);
    /* Not exactly 0: the mean of three equal doubles can differ from them in
       the last bit, so a deviation of order 1e-17 survives. */
    int all_zero = 1;
    double largest = 0;
    for (int model = 0; model < N_MODELS; model++) {
        all_zero = all_zero && whole.n_usable[model] == 3 && whole.sigma2_star[model] < 1e-25;
        if (whole.sigma2_star[model] > largest) largest = whole.sigma2_star[model];
    }
    char zero_what[160];
    snprintf(zero_what, sizeof zero_what, "blocks the length of the series give sigma2_star 0 under every method "
             "(largest %.1e)", largest);
    check(report, all_zero, zero_what);

    options.block_length = 8;
    Rng first_rng = rng_new(123, 2), second_rng = rng_new(123, 2);
    FitMetricBootstrap first = fit_metric_bootstrap(series, &start, options, &first_rng);
    FitMetricBootstrap second = fit_metric_bootstrap(series, &start, options, &second_rng);
    int same = 1;
    for (int model = 0; model < N_MODELS; model++) {
        same = same && first.sigma2_star[model] == second.sigma2_star[model];
        fprintf(report, "      %s sigma2_star over 3 resamples, blocks of 8: %.6g (%d usable)\n", model_label[model],
                first.sigma2_star[model], first.n_usable[model]);
    }
    check(report, same, "the same seed gives the same sigma2_star");

    for (int model = 0; model < N_MODELS; model++) free(out[model]);
    free(stored);
    qvarma_params_free(&start);
    mat_free(series);
    lp_free_samples(samples, n_samples);
}

int main(void) {
    openblas_set_num_threads(1);
    FILE *report = fopen(REPORT, "w");
    assert(report && "fit_metric_correctness: cannot open the report");
    hand_computed(report);
    oracle_shortcut(report);
    d_holds_out(report);
    real_series(report);
    fprintf(report, "\n%d failures\n", failures);
    fclose(report);
    return failures ? 1 : 0;
}
