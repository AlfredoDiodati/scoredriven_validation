/*
Whether montecarlo/own_deviation_loss.h computes the corrected loss of
docs/MCS_ALTERNATIVE_CALIBRATION.md, step 2, and whether it has the two
properties the document derives.

1. Against a reference written here from the definition, as plain loops in
   double: the squared error averaged over the entries, minus n/(n-1) times
   the replicate's own squared deviation from its configuration's entry means
   averaged over the entries. On synthetic caches with configurations at
   magnitudes from 1e-6 to 1e6, a held-out run and a dropped run missing from
   the kept rows, a batch smaller than its stride, and at one thread and at
   every hardware thread.
2. The smallest legal case, two kept runs, where n/(n-1) is 2, and a
   configuration whose runs are all equal, whose loss is the plain squared
   error.
3. Each configuration's own deviations average, over the kept rows, to its
   mean response variance as montecarlo/recovery_loss.h computes it, so the
   two losses have the same column means and differ only row by row.
4. The bias and the spread. Configuration A is the truth: entry means mu,
   standard deviation 1. Configuration B is wrong: means mu + 0.6, standard
   deviation 0.2. The benchmark is a fresh run of A held out of every column,
   200 benchmarks, 400 kept runs each, K = 50 entries. The corrected loss
   must prefer A on average, and every configuration's mean corrected loss
   must be within 0.1 of (1/K) sum_h (mu_ih - IR_dh)^2. The spread across
   rows: in A's column the loss of montecarlo/recovery_loss.h, which
   subtracts one common variance, varies across rows by about
   4 D / K + 2 / K, with D the benchmark's squared distance from mu, about 1,
   while the corrected loss varies by about 4 D / K alone. The ratio of the
   two variances, averaged over benchmarks, must be within 0.05 of 2/3.

Writes out/tests/own_deviation_loss_correctness_report.txt. Exits 1 on any
failure.
*/

#include "montecarlo/own_deviation_loss.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>

#define REPORT "out/tests/own_deviation_loss_correctness_report.txt"

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) failures++;
}

static ResponseCache cache_new(int n_samples, int n_replicates, int dim) {
    ResponseCache cache = { MODEL_QVARMA, n_samples, n_replicates, dim, NULL, 0, NULL, 0 };
    cache.data = malloc((size_t)n_samples * n_replicates * dim * sizeof(float));
    assert(cache.data);
    return cache;
}

/* The loss of configuration sample on kept row row against benchmark, from
   the definition. */
static double reference_loss(const ResponseCache *cache, const int *kept, int n_kept, int sample, int row,
                             const float *benchmark) {
    const int dim = cache->dim;
    double squared = 0, own = 0;
    for (int h = 0; h < dim; h++) {
        double mean = 0;
        for (int r = 0; r < n_kept; r++) mean += (double)response_slot(cache, sample, kept[r])[h];
        mean /= n_kept;
        double value = (double)response_slot(cache, sample, kept[row])[h];
        own += (value - mean) * (value - mean);
        squared += (value - (double)benchmark[h]) * (value - (double)benchmark[h]);
    }
    return squared / dim - (double)n_kept / (n_kept - 1) * own / dim;
}

/* The header's loss matrices for the given benchmarks, one Mat each. */
static Mat *header_losses(const ResponseCache *cache, const int *kept, int n_kept, const float *const *benchmark,
                          int n_batch, int batch_stride, int n_threads) {
    double *own_deviation = malloc((size_t)cache->n_samples * n_kept * sizeof(double));
    double *no_correction = calloc((size_t)cache->n_samples, sizeof(double));
    double *reference = calloc((size_t)cache->dim * batch_stride, sizeof(double));
    Mat *losses = malloc((size_t)n_batch * sizeof(Mat));
    assert(own_deviation && no_correction && reference && losses);
    for (int k = 0; k < n_batch; k++) {
        losses[k] = mat_new(n_kept, cache->n_samples);
        for (int h = 0; h < cache->dim; h++) reference[(size_t)h * batch_stride + k] = (double)benchmark[k][h];
    }
    own_deviation_loss_deviations(cache, kept, n_kept, own_deviation, n_threads);
    own_deviation_loss_fill(cache, kept, n_kept, reference, n_batch, batch_stride, own_deviation, no_correction,
                            losses, n_threads);
    free(own_deviation);
    free(no_correction);
    free(reference);
    return losses;
}

/* The largest difference between the header and the reference, relative to
   the larger of the reference's absolute value and 1e-300. */
static double largest_relative_error(const ResponseCache *cache, const int *kept, int n_kept,
                                     const float *const *benchmark, int n_batch, const Mat *losses) {
    double worst = 0;
    for (int k = 0; k < n_batch; k++)
        for (int row = 0; row < n_kept; row++)
            for (int sample = 0; sample < cache->n_samples; sample++) {
                double expected = reference_loss(cache, kept, n_kept, sample, row, benchmark[k]);
                double scale = fabs(expected) > 1e-300 ? fabs(expected) : 1e-300;
                double error = fabs((double)AT(losses[k], row, sample) - expected) / scale;
                if (error > worst) worst = error;
            }
    return worst;
}

static void free_losses(Mat *losses, int n_batch) {
    for (int k = 0; k < n_batch; k++) mat_free(losses[k]);
    free(losses);
}

static void against_reference(FILE *report) {
    const int n_samples = 7, n_replicates = 11, dim = 13, n_batch = 5, batch_stride = 8;
    const double magnitude[7] = { 1, 1e6, 1e-6, 3.5, 1e3, 0.25, 40 };
    ResponseCache cache = cache_new(n_samples, n_replicates, dim);
    Rng rng = rng_new(20261013, 0);
    for (int sample = 0; sample < n_samples; sample++)
        for (int replicate = 0; replicate < n_replicates; replicate++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, replicate)[h] =
                    (float)(magnitude[sample] * (h % 3 + rng_normal(&rng)));

    /* Run 4 is the benchmark's run index and run 7 a dropped one. */
    int kept[11], n_kept = 0;
    for (int r = 0; r < n_replicates; r++) if (r != 4 && r != 7) kept[n_kept++] = r;
    const float *benchmark[5];
    for (int k = 0; k < n_batch; k++) benchmark[k] = response_slot(&cache, k + 1, 4);

    int thread_counts[2] = { 1, omp_get_num_procs() };
    Mat *at_count[2];
    for (int t = 0; t < 2; t++) {
        at_count[t] = header_losses(&cache, kept, n_kept, benchmark, n_batch, batch_stride, thread_counts[t]);
        double worst = largest_relative_error(&cache, kept, n_kept, benchmark, n_batch, at_count[t]);
        fprintf(report, "      %d thread(s): largest relative error against the reference %.3e\n",
                thread_counts[t], worst);
        char what[160];
        snprintf(what, sizeof what, "loss matches the reference to 1e-9 relative at %d thread(s)",
                 thread_counts[t]);
        check(report, worst < 1e-9, what);
    }
    int identical = 1;
    for (int k = 0; k < n_batch; k++)
        for (int row = 0; row < n_kept; row++)
            for (int sample = 0; sample < n_samples; sample++)
                if (fabs((double)AT(at_count[0][k], row, sample) - (double)AT(at_count[1][k], row, sample))
                    > 1e-12 * fabs((double)AT(at_count[0][k], row, sample)))
                    identical = 0;
    check(report, identical, "one thread and every thread give the same losses");
    for (int t = 0; t < 2; t++) free_losses(at_count[t], n_batch);
    free(cache.data);
}

static void smallest_and_constant(FILE *report) {
    const int n_samples = 3, n_replicates = 3, dim = 4;
    ResponseCache cache = cache_new(n_samples, n_replicates, dim);
    Rng rng = rng_new(20261013, 1);
    for (int sample = 0; sample < n_samples; sample++)
        for (int replicate = 0; replicate < n_replicates; replicate++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, replicate)[h] = sample == 2 ? (float)(h - 1.5) : (float)rng_normal(&rng);

    int kept[2] = { 0, 2 };
    const float *benchmark[1] = { response_slot(&cache, 0, 1) };
    Mat *losses = header_losses(&cache, kept, 2, benchmark, 1, 1, omp_get_num_procs());
    double worst = largest_relative_error(&cache, kept, 2, benchmark, 1, losses);
    check(report, worst < 1e-9, "two kept runs, the smallest legal case, match the reference");

    int plain = 1;
    for (int row = 0; row < 2; row++) {
        double squared = 0;
        for (int h = 0; h < dim; h++) {
            double difference = (double)response_slot(&cache, 2, kept[row])[h] - (double)benchmark[0][h];
            squared += difference * difference;
        }
        if (fabs((double)AT(losses[0], row, 2) - squared / dim) > 1e-12 * (squared / dim + 1)) plain = 0;
    }
    check(report, plain, "a configuration with equal runs has no deviation, its loss is the squared error");
    free_losses(losses, 1);
    free(cache.data);
}

static void same_column_means(FILE *report) {
    const int n_samples = 5, n_replicates = 30, dim = 21;
    const double magnitude[5] = { 1, 1e4, 1e-4, 7, 0.5 };
    ResponseCache cache = cache_new(n_samples, n_replicates, dim);
    Rng rng = rng_new(20261013, 2);
    for (int sample = 0; sample < n_samples; sample++)
        for (int replicate = 0; replicate < n_replicates; replicate++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, replicate)[h] = (float)(magnitude[sample] * (h + rng_normal(&rng)));

    int kept[30], n_kept = 0;
    for (int r = 0; r < n_replicates; r++) if (r != 3) kept[n_kept++] = r;
    double own_deviation[5 * 30], mean_variance[5];
    own_deviation_loss_deviations(&cache, kept, n_kept, own_deviation, omp_get_num_procs());
    recovery_loss_mean_variance(&cache, kept, n_kept, mean_variance, omp_get_num_procs());
    double worst = 0;
    for (int sample = 0; sample < n_samples; sample++) {
        double average = 0;
        for (int row = 0; row < n_kept; row++) average += own_deviation[sample * n_kept + row];
        average /= n_kept;
        double error = fabs(average - mean_variance[sample]) / mean_variance[sample];
        if (error > worst) worst = error;
    }
    fprintf(report, "      largest relative difference from recovery_loss_mean_variance %.3e\n", worst);
    check(report, worst < 1e-9, "the own deviations average to the mean response variance, to 1e-9 relative");
    free(cache.data);
}

static double column_variance(const Mat losses, int sample) {
    double mean = 0, sum = 0;
    for (int row = 0; row < losses.r; row++) mean += (double)AT(losses, row, sample);
    mean /= losses.r;
    for (int row = 0; row < losses.r; row++) {
        double deviation = (double)AT(losses, row, sample) - mean;
        sum += deviation * deviation;
    }
    return sum / (losses.r - 1);
}

static void bias_and_spread(FILE *report) {
    const int dim = 50, n_runs = 401, n_benchmarks = 200;
    const double sd[2] = { 1.0, 0.2 }, shift[2] = { 0.0, 0.6 };
    ResponseCache cache = cache_new(2, n_runs, dim);
    double mu[2][50];
    Rng rng = rng_new(20261013, 3);
    for (int h = 0; h < dim; h++) mu[0][h] = mu[1][h] = 0.5 * rng_normal(&rng);
    for (int h = 0; h < dim; h++) mu[1][h] += shift[1];
    for (int sample = 0; sample < 2; sample++)
        for (int run = 0; run < n_runs; run++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, run)[h] = (float)(mu[sample][h] + sd[sample] * rng_normal(&rng));

    double corrected_average[2] = { 0, 0 }, worst_miss = 0;
    double corrected_spread = 0, common_spread = 0;
    int *kept = malloc((size_t)n_runs * sizeof(int));
    float *benchmark_values = malloc((size_t)dim * sizeof(float));
    double *reference = malloc((size_t)dim * sizeof(double));
    assert(kept && benchmark_values && reference);
    int n_kept = 0;
    for (int r = 1; r < n_runs; r++) kept[n_kept++] = r;
    double mean_variance[2];
    recovery_loss_mean_variance(&cache, kept, n_kept, mean_variance, omp_get_num_procs());

    for (int b = 0; b < n_benchmarks; b++) {
        /* The benchmark is a fresh run of the true configuration, written
           into run 0 and held out of both columns, as in the grid. */
        for (int h = 0; h < dim; h++) {
            benchmark_values[h] = (float)(mu[0][h] + sd[0] * rng_normal(&rng));
            response_slot(&cache, 0, 0)[h] = benchmark_values[h];
            reference[h] = (double)benchmark_values[h];
        }
        const float *benchmark[1] = { benchmark_values };
        Mat *losses = header_losses(&cache, kept, n_kept, benchmark, 1, 1, omp_get_num_procs());
        Mat common = mat_new(n_kept, 2);
        recovery_loss_fill(&cache, kept, n_kept, reference, 1, 1, mean_variance, &common, omp_get_num_procs());

        for (int sample = 0; sample < 2; sample++) {
            double mean_loss = 0;
            for (int row = 0; row < n_kept; row++) mean_loss += (double)AT(losses[0], row, sample);
            mean_loss /= n_kept;
            double target = 0;
            for (int h = 0; h < dim; h++) {
                double difference = mu[sample][h] - (double)benchmark_values[h];
                target += difference * difference;
            }
            target /= dim;
            corrected_average[sample] += mean_loss / n_benchmarks;
            if (fabs(mean_loss - target) > worst_miss) worst_miss = fabs(mean_loss - target);
        }
        corrected_spread += column_variance(losses[0], 0) / n_benchmarks;
        common_spread += column_variance(common, 0) / n_benchmarks;
        mat_free(common);
        free_losses(losses, 1);
    }
    double ratio = corrected_spread / common_spread;
    fprintf(report, "      average over %d benchmarks, true A and wrong B:\n", n_benchmarks);
    fprintf(report, "      corrected mean loss  A %.4f  B %.4f   (expected 1.00 and 1.36)\n", corrected_average[0],
            corrected_average[1]);
    fprintf(report, "      largest miss of the corrected mean loss from its target %.4f\n", worst_miss);
    fprintf(report, "      variance across rows in A's column: own deviation %.4f, common variance %.4f\n",
            corrected_spread, common_spread);
    fprintf(report, "      (expected about 0.08 and 0.12), ratio %.4f (expected 2/3)\n", ratio);
    check(report, corrected_average[0] < corrected_average[1], "the corrected loss prefers the true configuration");
    check(report, worst_miss < 0.1, "every corrected mean loss is within 0.1 of its squared distance");
    check(report, fabs(ratio - 2.0 / 3.0) < 0.05,
          "subtracting the own deviation cuts the spread across rows to 2/3 of the common variance's, within 0.05");
    free(kept);
    free(benchmark_values);
    free(reference);
    free(cache.data);
}

int main(void) {
    FILE *report = fopen(REPORT, "w");
    assert(report && "own_deviation_loss_correctness: cannot open the report");
    fprintf(report, "montecarlo/own_deviation_loss.h against its definition, %d hardware threads\n\n",
            omp_get_num_procs());
    fprintf(report, "1. against a reference from the definition\n");
    against_reference(report);
    fprintf(report, "\n2. two kept runs, and a configuration with equal runs\n");
    smallest_and_constant(report);
    fprintf(report, "\n3. the same column means as montecarlo/recovery_loss.h\n");
    same_column_means(report);
    fprintf(report, "\n4. the bias, and the spread across rows\n");
    bias_and_spread(report);
    fprintf(report, "\n%s, %d failure(s)\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    return failures ? 1 : 0;
}
