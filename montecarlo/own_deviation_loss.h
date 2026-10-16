#ifndef MONTECARLO_OWN_DEVIATION_LOSS_H
#define MONTECARLO_OWN_DEVIATION_LOSS_H

/*
The corrected loss of docs/MCS_ALTERNATIVE_CALIBRATION.md, step 2, over the
response caches of montecarlo/response_cache.h: for replicate s of
configuration i against benchmark d,

    L*_is = (1/K) sum_h (IR_ish - IR_dh)^2
            - n/(n-1) (1/K) sum_h (IR_ish - IR_bar_ih)^2,

with IR_bar_ih the mean of entry h over the configuration's n kept runs.
Every configuration's kept runs are the same rows, the rows of the loss
matrix, so in the Monte Carlo the benchmark's run is never among them.

The subtracted term, the replicate's own deviation, depends on which rows are
kept and not on the benchmark, so it is computed once per run index by
own_deviation_loss_deviations. own_deviation_loss_fill then forms the squared
errors of a batch of benchmarks through montecarlo/recovery_loss.h's
recovery_loss_fill and subtracts it.
*/

#include "montecarlo/recovery_loss.h"

/* n/(n-1) (1/K) sum_h (IR_ish - IR_bar_ih)^2 for every configuration i and
   every kept row s, written to own_deviation[i * n_kept + s], over the n_kept
   runs listed in kept. The entry means are sums in double over the kept rows,
   in order. Runs on n_threads threads, one configuration at a time per
   thread. */
static inline void own_deviation_loss_deviations(const ResponseCache *cache, const int *kept, int n_kept,
                                                 double *own_deviation, int n_threads) {
    assert(n_kept >= 2 && "own_deviation_loss: a deviation from the mean needs at least two runs");
    const int dim = cache->dim;
    const double scale = (double)n_kept / (n_kept - 1) / dim;

    #pragma omp parallel num_threads(n_threads)
    {
        double *entry_mean = malloc((size_t)dim * sizeof(double));
        assert(entry_mean && "own_deviation_loss: out of memory");
        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < cache->n_samples; sample++) {
            for (int h = 0; h < dim; h++) entry_mean[h] = 0;
            for (int row = 0; row < n_kept; row++) {
                const float *cell = response_slot(cache, sample, kept[row]);
                for (int h = 0; h < dim; h++) entry_mean[h] += (double)cell[h];
            }
            for (int h = 0; h < dim; h++) entry_mean[h] /= n_kept;

            for (int row = 0; row < n_kept; row++) {
                const float *cell = response_slot(cache, sample, kept[row]);
                double sum = 0;
                for (int h = 0; h < dim; h++) {
                    double deviation = (double)cell[h] - entry_mean[h];
                    sum += deviation * deviation;
                }
                own_deviation[(size_t)sample * n_kept + row] = scale * sum;
            }
        }
        free(entry_mean);
    }
}

/* The loss matrices of n_batch benchmarks that share their kept runs.
   reference holds the benchmarks' responses entry by entry,
   reference[h * batch_stride + k] for benchmark k; losses[k] is
   n_kept x n_samples and receives L*_is at (row, sample). no_correction is
   n_samples zeros, handed to recovery_loss_fill so that it writes the plain
   squared error. */
static inline void own_deviation_loss_fill(const ResponseCache *cache, const int *kept, int n_kept,
                                           const double *reference, int n_batch, int batch_stride,
                                           const double *own_deviation, const double *no_correction, Mat *losses,
                                           int n_threads) {
    recovery_loss_fill(cache, kept, n_kept, reference, n_batch, batch_stride, no_correction, losses, n_threads);
    const int n_samples = cache->n_samples;

    #pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int row = 0; row < n_kept; row++)
        for (int k = 0; k < n_batch; k++)
            for (int sample = 0; sample < n_samples; sample++)
                AT(losses[k], row, sample) -= (mreal)own_deviation[(size_t)sample * n_kept + row];
}

#endif /* MONTECARLO_OWN_DEVIATION_LOSS_H */
