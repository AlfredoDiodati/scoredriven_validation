#ifndef MONTECARLO_FIT_METRIC_H
#define MONTECARLO_FIT_METRIC_H

/*
The fit metric of docs/MONTECARLO_FIT_METRIC.md: how far the best configuration
of a confidence set is from the data, measured against how much the impulse
responses move from one sample to another.

For a benchmark whose true configuration is j and a method m (one of the five
response vectors of montecarlo/response_cache.h), with d the mean loss of the
configuration with the smallest mean loss in the final confidence set and sigma2
the sampling variance of one entry of the response vector, averaged over its K
entries,

    s = d / sigma,    sigma = sqrt(sigma2),    v = 1 / (1 + s).

The loss is a mean absolute difference, in the units of the responses, so it
is divided by a standard deviation rather than a variance: s then does not
change when every response is multiplied by a constant.

sigma2 comes in two versions. The oracle one is the variance across the true
configuration's other runs, which only a simulation has; montecarlo/fit_metric_oracle.c
computes it for every benchmark of the grid from the response caches. The
feasible one, sigma2_star, is the variance across B bootstrap resamples of the
benchmark's own series, re-estimated by every method, which can be computed on
the US data too; fit_metric_bootstrap below does it for one series.

The resamples are moving-block bootstrap draws of the stored series, the
K x T matrix the t-QVARMA is fitted on, whole periods at a time so the five
series keep their contemporaneous relation, through et_al's mcs_block_indices.
Every method is re-estimated on the same resamples: the t-QVARMA from the
estimate on the original series as its starting point, the local projections
through applications/lp_system.h exactly as applications/abm_system_fit_lp.c
fits them, so a resample's four local-projection vectors come from one fit.
*/

#include "montecarlo/response_cache.h"
#include <et_al./inference/mcs.h>
#include <et_al./random/random.h>

static inline double fit_metric_s(double d, double sigma2) { return d / sqrt(sigma2); }
static inline double fit_metric_v(double s) { return 1 / (1 + s); }

/*
sigma2 over n vectors of dim entries stored one after another: per entry the
mean over the n vectors and the variance with divisor n - 1, then the average
of the dim variances. Two passes, so the variance is not formed as a difference
of large sums. n must be at least 2.
*/
static inline double fit_metric_sigma2(const double *vectors, int n, int dim) {
    assert(n >= 2 && dim >= 1);
    double total = 0;
    for (int k = 0; k < dim; k++) {
        double mean = 0;
        for (int b = 0; b < n; b++) mean += vectors[(size_t)b * dim + k];
        mean /= n;
        double squares = 0;
        for (int b = 0; b < n; b++) {
            double deviation = vectors[(size_t)b * dim + k] - mean;
            squares += deviation * deviation;
        }
        total += squares / (n - 1);
    }
    return total / dim;
}

/*
The oracle sigma2 of every benchmark of one configuration comes from two
arrays of dim entries over all n of its runs: the mean, and the sum of squared
deviations from it. Taking run x out of both gives the mean over the other
n - 1 runs, mean' = (n mean - x) / (n - 1), and their sum of squares,
squares - (x - mean)(x - mean'), which divided by n - 2 is the variance with
divisor (n - 1) - 1. That costs dim operations per benchmark instead of n dim.
*/
static inline void fit_metric_moments(const ResponseCache *cache, int sample, double *mean, double *squares) {
    const int n = cache->n_replicates, dim = cache->dim;
    for (int k = 0; k < dim; k++) mean[k] = squares[k] = 0;
    for (int run = 0; run < n; run++) {
        const float *x = response_slot(cache, sample, run);
        for (int k = 0; k < dim; k++) mean[k] += (double)x[k];
    }
    for (int k = 0; k < dim; k++) mean[k] /= n;
    for (int run = 0; run < n; run++) {
        const float *x = response_slot(cache, sample, run);
        for (int k = 0; k < dim; k++) {
            double deviation = (double)x[k] - mean[k];
            squares[k] += deviation * deviation;
        }
    }
}

static inline double fit_metric_oracle_sigma2(const double *mean, const double *squares, const double *held_out,
                                              int n, int dim) {
    assert(n >= 3);
    double total = 0;
    for (int k = 0; k < dim; k++) {
        double mean_without = (n * mean[k] - held_out[k]) / (n - 1);
        total += (squares[k] - (held_out[k] - mean[k]) * (held_out[k] - mean_without)) / (n - 2);
    }
    return total / dim;
}

/*
d for one benchmark: the mean, over every run of configuration best_sample but
held_out_run, of the mean absolute difference between that run's responses and
the benchmark's. The same sums in the same order as montecarlo/sweep_grid.c
forms its loss matrix, whose column for best_sample this averages.
*/
static inline double fit_metric_d(const ResponseCache *cache, int best_sample, int held_out_run,
                                  const double *benchmark) {
    const int dim = cache->dim;
    double total = 0;
    for (int run = 0; run < cache->n_replicates; run++) {
        if (run == held_out_run) continue;
        const float *cell = response_slot(cache, best_sample, run);
        double sum = 0;
        for (int k = 0; k < dim; k++) sum += fabs((double)cell[k] - benchmark[k]);
        total += (double)(mreal)(sum / dim);
    }
    return total / (cache->n_replicates - 1);
}

typedef struct {
    int n_bootstrap;            /* B, resamples per series */
    int block_length;           /* periods per block of the moving-block bootstrap */
    int qvarma_max_iterations;  /* the t-QVARMA solver's budget per resample */
} FitMetricBootstrapOptions;

/* B = 200, blocks of 8 periods (the cube root of 400 periods, rounded up),
   4000 iterations, et_al's default budget for the t-QVARMA. */
static inline FitMetricBootstrapOptions fit_metric_bootstrap_options_default(void) {
    return (FitMetricBootstrapOptions){ 200, 8, 4000 };
}

typedef struct {
    double sigma2_star[N_MODELS];
    int n_usable[N_MODELS];      /* resamples whose response under the method is finite */
    int n_qvarma_converged;      /* resamples whose t-QVARMA fit met the tolerance */
    double qvarma_mean_iterations;
} FitMetricBootstrap;

/*
The five response vectors of one series, written to out[model], which holds
response_dim(model) entries. usable[model] is 0 where that response is missing
or not finite, under lp_nl when either state's is. series is ABM_SYSTEM_K x T,
the layout the t-QVARMA is fitted on; the t-QVARMA is fitted on it from
qvarma_start, and its convergence and iterations go to fit_converged and
fit_iterations. The local projections are fitted through lp_system_series and
lp_lin_and_nl as applications/abm_system_fit_lp.c fits every simulated run.
*/
static inline void fit_metric_responses(Mat series, const QvarmaParams *qvarma_start, QvarmaFitOptions fit_options,
                                        double *out[N_MODELS], int usable[N_MODELS], int *fit_converged,
                                        int *fit_iterations) {
    float qvarma_response[RESPONSE_QVARMA_DIM];
    QvarmaFitResult fit = qvarma_fit(series, qvarma_start, fit_options);
    *fit_converged = fit.is_converged;
    *fit_iterations = fit.niter;
    usable[MODEL_QVARMA] = _response_qvarma_cell(&fit.params, series, qvarma_response);
    for (int i = 0; i < RESPONSE_QVARMA_DIM; i++) out[MODEL_QVARMA][i] = (double)qvarma_response[i];
    qvarma_fit_result_free(&fit);

    Mat y, switching;
    lp_system_series(series, LP_LAYOUT_GROWTH, &y, &switching);
    LpLinNlFit lp = lp_lin_and_nl(y, switching, lp_system_spec(), (LpBands){0});
    int has_responses = lp.lin.d.d && lp.nl.d.d;
    const mreal *state[LP_N_LOSSES] = { NULL, NULL, NULL };
    int finite[LP_N_LOSSES] = { 0, 0, 0 };
    if (has_responses) {
        state[LP_LOSS_LIN] = lp.lin.irf_lin_mean.d;
        state[LP_LOSS_S1] = lp.nl.irf_s1_mean.d;
        state[LP_LOSS_S2] = lp.nl.irf_s2_mean.d;
        for (int loss = 0; loss < LP_N_LOSSES; loss++) finite[loss] = _lp_all_finite(state[loss]);
    }
    for (int model = MODEL_LP_LIN; model < N_MODELS; model++) {
        int first = model == MODEL_LP_NL ? LP_LOSS_S1 : model - MODEL_LP_LIN;
        int last = model == MODEL_LP_NL ? LP_LOSS_S2 : first;
        usable[model] = has_responses;
        for (int loss = first; loss <= last; loss++) usable[model] = usable[model] && finite[loss];
        if (!usable[model]) continue;
        for (int loss = first; loss <= last; loss++)
            for (int i = 0; i < LP_RESPONSE_DIM; i++)
                out[model][(size_t)(loss - first) * LP_RESPONSE_DIM + i] = (double)state[loss][i];
    }
    lp_lin_nl_fit_free(&lp);
    mat_free(y);
    mat_free(switching);
}

/*
sigma2_star under every method for one series. stored is ABM_SYSTEM_K x T, the
series the t-QVARMA is fitted on; qvarma_start is the t-QVARMA estimate on it,
where each resample's fit starts. rng supplies the block starts and nothing
else, so a series and a seed give the same resamples whatever else runs.

A resample whose response under a method is missing or not finite is left out
of that method's sigma2_star and of nothing else; n_usable counts the rest, and
a method with fewer than two gets NaN.
*/
static inline FitMetricBootstrap fit_metric_bootstrap(Mat stored, const QvarmaParams *qvarma_start,
                                                      FitMetricBootstrapOptions options, Rng *rng) {
    assert(options.n_bootstrap >= 2 && options.block_length >= 1 && options.block_length <= stored.c);
    const int B = options.n_bootstrap, T = stored.c;
    FitMetricBootstrap result = { { 0 }, { 0 }, 0, 0 };

    double *responses[N_MODELS], *slot[N_MODELS];
    for (int model = 0; model < N_MODELS; model++) {
        responses[model] = malloc((size_t)B * response_dim(model) * sizeof(double));
        assert(responses[model] && "fit_metric_bootstrap: out of memory");
    }
    int *period = malloc((size_t)T * sizeof(int));
    assert(period && "fit_metric_bootstrap: out of memory");
    Mat resample = mat_new(stored.r, T);
    QvarmaFitOptions fit_options = qvarma_default_fit_options();
    fit_options.max_iterations = options.qvarma_max_iterations;
    long iterations = 0;

    for (int b = 0; b < B; b++) {
        mcs_block_indices(rng, T, options.block_length, period);
        for (int k = 0; k < stored.r; k++)
            for (int t = 0; t < T; t++) AT(resample, k, t) = AT(stored, k, period[t]);

        /* Each method writes into its next free row, which is kept only when
           the response is usable. */
        for (int model = 0; model < N_MODELS; model++)
            slot[model] = responses[model] + (size_t)result.n_usable[model] * response_dim(model);
        int usable[N_MODELS], converged, fit_iterations;
        fit_metric_responses(resample, qvarma_start, fit_options, slot, usable, &converged, &fit_iterations);
        for (int model = 0; model < N_MODELS; model++) result.n_usable[model] += usable[model];
        result.n_qvarma_converged += converged;
        iterations += fit_iterations;
    }

    for (int model = 0; model < N_MODELS; model++) {
        int n = result.n_usable[model];
        result.sigma2_star[model] = n >= 2 ? fit_metric_sigma2(responses[model], n, response_dim(model)) : NAN;
        free(responses[model]);
    }
    result.qvarma_mean_iterations = (double)iterations / B;
    mat_free(resample);
    free(period);
    return result;
}

#endif /* MONTECARLO_FIT_METRIC_H */
