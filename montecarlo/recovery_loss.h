#ifndef MONTECARLO_RECOVERY_LOSS_H
#define MONTECARLO_RECOVERY_LOSS_H

/*
The loss of docs/MCS_VALIDATION_PROCEDURE.md, "The steps", over the response
caches of montecarlo/response_cache.h: for replicate s of configuration i
against benchmark d,

    L*_is = (1/K) sum_h (IR_ish - IR_dh)^2 - (1/K) sum_h sigma2_ih,

with sigma2_ih the variance of entry h across the configuration's kept runs,
divisor n - 1. Every configuration's kept runs are the same rows, the rows of
the loss matrix, so in the Monte Carlo the benchmark's run is never among them.

Two steps, because the variance depends only on which rows are kept and every
benchmark sharing a run index shares those rows: recovery_loss_mean_variance
once per run index, then recovery_loss_fill once per batch of benchmarks.
*/

#include "montecarlo/response_cache.h"
#include <et_al./stats.h>
#include <omp.h>

/* (1/K) sum_h sigma2_ih for every configuration i, over the n_kept runs listed
   in kept, written to mean_variance[i]. Each entry's variance is stats_var on
   that entry's n_kept values, rescaled from divisor n to n - 1. Runs on
   n_threads threads, one configuration at a time per thread. */
static inline void recovery_loss_mean_variance(const ResponseCache *cache, const int *kept, int n_kept,
                                               double *mean_variance, int n_threads) {
    assert(n_kept >= 2 && "recovery_loss: a variance needs at least two runs");
    const int dim = cache->dim;
    const double unbiased = (double)n_kept / (n_kept - 1);

    #pragma omp parallel num_threads(n_threads)
    {
        /* One row per entry, so each entry's values are contiguous. */
        Mat entries = mat_new(dim, n_kept);
        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < cache->n_samples; sample++) {
            for (int row = 0; row < n_kept; row++) {
                const float *cell = response_slot(cache, sample, kept[row]);
                for (int h = 0; h < dim; h++) AT(entries, h, row) = (mreal)cell[h];
            }
            double total = 0;
            for (int h = 0; h < dim; h++) total += (double)stats_var(mat_slice(entries, h, h + 1, 0, n_kept));
            mean_variance[sample] = unbiased * total / dim;
        }
        mat_free(entries);
    }
}

/* The loss matrices of n_batch benchmarks that share their kept runs, one read
   of each cell serving every benchmark. reference holds the benchmarks'
   responses entry by entry, reference[h * batch_stride + k] for benchmark k;
   losses[k] is n_kept x n_samples and receives L*_is at (row, sample).
   Each loss is the sum over h in order, divided by dim, minus
   mean_variance[sample]. */
static inline void recovery_loss_fill(const ResponseCache *cache, const int *kept, int n_kept,
                                      const double *reference, int n_batch, int batch_stride,
                                      const double *mean_variance, Mat *losses, int n_threads) {
    assert(n_batch >= 1 && n_batch <= batch_stride);
    const int dim = cache->dim, n_samples = cache->n_samples;

    #pragma omp parallel num_threads(n_threads)
    {
        double *sum = malloc((size_t)batch_stride * sizeof(double));
        assert(sum && "recovery_loss: out of memory");
        #pragma omp for schedule(static)
        for (int row = 0; row < n_kept; row++) {
            for (int sample = 0; sample < n_samples; sample++) {
                const float *cell = response_slot(cache, sample, kept[row]);
                for (int k = 0; k < n_batch; k++) sum[k] = 0;
                for (int h = 0; h < dim; h++) {
                    double x = (double)cell[h];
                    const double *ref = reference + (size_t)h * batch_stride;
                    for (int k = 0; k < n_batch; k++) {
                        double difference = x - ref[k];
                        sum[k] += difference * difference;
                    }
                }
                for (int k = 0; k < n_batch; k++)
                    AT(losses[k], row, sample) = (mreal)(sum[k] / dim - mean_variance[sample]);
            }
        }
        free(sum);
    }
}

#endif /* MONTECARLO_RECOVERY_LOSS_H */
