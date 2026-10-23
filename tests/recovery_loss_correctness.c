/*
Whether montecarlo/recovery_loss.h computes the loss of
docs/MCS_VALIDATION_PROCEDURE.md, and whether that loss removes the bias the
document derives.

1. Against a reference written here from the definitions, as plain loops in
   double: the squared error averaged over the entries, minus the average over
   the entries of each entry's variance across the kept runs with divisor
   n - 1. On synthetic caches with configurations at magnitudes from 1e-6 to
   1e6, a held-out run and a dropped run missing from the kept rows, a batch
   smaller than its stride, and at one thread and at every hardware thread.
2. The smallest legal case, two kept runs, and a configuration whose runs are
   all equal, whose variance is zero, so its loss is the plain squared error.
3. The bias. Configuration A is the truth: entry means mu, standard deviation
   1. Configuration B is wrong: means mu + 0.6, standard deviation 0.2. The
   benchmark is a run of A held out of every column, 200 benchmarks, 400 kept
   runs each, 50 entries. In expectation over benchmarks the uncorrected
   loss of A is 1 + 1 = 2 and of B 1 + 0.36 + 0.04 = 1.40, so it prefers the
   wrong configuration; the corrected loss is 1 for A and 1.36 for B. Checked:
   the uncorrected loss prefers B on average, the corrected one prefers A, and
   for every benchmark each configuration's mean corrected loss is within 0.1
   of the squared distance between its entry means and the benchmark,
   (1/K) sum_h (mu_ih - IR_dh)^2. The uncorrected loss of A misses that target
   by about its variance, 1, which is the negative control that the check can
   fail.

Writes out/tests/recovery_loss_correctness_report.txt. Exits 1 on any failure.
*/

#include "montecarlo/recovery_loss.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>

#define REPORT "out/tests/recovery_loss_correctness_report.txt"

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
   the definitions. */
static double reference_loss(const ResponseCache *cache, const int *kept, int n_kept, int sample, int row,
                             const float *benchmark) {
    const int dim = cache->dim;
    double squared = 0, variance = 0;
    for (int h = 0; h < dim; h++) {
        double mean = 0;
        for (int r = 0; r < n_kept; r++) mean += (double)response_slot(cache, sample, kept[r])[h];
        mean /= n_kept;
        double sum = 0;
        for (int r = 0; r < n_kept; r++) {
            double deviation = (double)response_slot(cache, sample, kept[r])[h] - mean;
            sum += deviation * deviation;
        }
        variance += sum / (n_kept - 1);
        double difference = (double)response_slot(cache, sample, kept[row])[h] - (double)benchmark[h];
        squared += difference * difference;
    }
    return squared / dim - variance / dim;
}

/* The header's loss matrices for the given benchmarks, one Mat each. */
static Mat *header_losses(const ResponseCache *cache, const int *kept, int n_kept, const float *const *benchmark,
                          int n_batch, int batch_stride, int n_threads) {
    double *mean_variance = malloc((size_t)cache->n_samples * sizeof(double));
    double *reference = calloc((size_t)cache->dim * batch_stride, sizeof(double));
    Mat *losses = malloc((size_t)n_batch * sizeof(Mat));
    assert(mean_variance && reference && losses);
    for (int k = 0; k < n_batch; k++) {
        losses[k] = mat_new(n_kept, cache->n_samples);
        for (int h = 0; h < cache->dim; h++) reference[(size_t)h * batch_stride + k] = (double)benchmark[k][h];
    }
    recovery_loss_mean_variance(cache, kept, n_kept, mean_variance, n_threads);
    recovery_loss_fill(cache, kept, n_kept, reference, n_batch, batch_stride, mean_variance, losses, n_threads);
    free(mean_variance);
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
    Rng rng = rng_new(20261012, 0);
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
    Rng rng = rng_new(20261012, 1);
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
    check(report, plain, "a configuration with equal runs has zero variance, its loss is the squared error");
    free_losses(losses, 1);
    free(cache.data);
}

static void removes_the_bias(FILE *report) {
    const int dim = 50, n_runs = 401, n_benchmarks = 200;
    const double sd[2] = { 1.0, 0.2 }, shift[2] = { 0.0, 0.6 };
    ResponseCache cache = cache_new(2, n_runs, dim);
    double mu[2][50];
    Rng rng = rng_new(20261012, 2);
    for (int h = 0; h < dim; h++) mu[0][h] = mu[1][h] = 0.5 * rng_normal(&rng);
    for (int h = 0; h < dim; h++) mu[1][h] += shift[1];
    for (int sample = 0; sample < 2; sample++)
        for (int run = 0; run < n_runs; run++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, run)[h] = (float)(mu[sample][h] + sd[sample] * rng_normal(&rng));

    double corrected_average[2] = { 0, 0 }, uncorrected_average[2] = { 0, 0 };
    double worst_miss = 0, smallest_uncorrected_miss_of_truth = INFINITY;
    int *kept = malloc((size_t)n_runs * sizeof(int));
    float *benchmark_values = malloc((size_t)dim * sizeof(float));
    assert(kept && benchmark_values);
    for (int b = 0; b < n_benchmarks; b++) {
        /* The benchmark is a fresh run of the true configuration, written
           into run 0 and held out of both columns, as in the grid. */
        for (int h = 0; h < dim; h++) {
            benchmark_values[h] = (float)(mu[0][h] + sd[0] * rng_normal(&rng));
            response_slot(&cache, 0, 0)[h] = benchmark_values[h];
        }
        int n_kept = 0;
        for (int r = 1; r < n_runs; r++) kept[n_kept++] = r;
        const float *benchmark[1] = { benchmark_values };
        Mat *losses = header_losses(&cache, kept, n_kept, benchmark, 1, 1, omp_get_num_procs());
        double mean_variance[2];
        recovery_loss_mean_variance(&cache, kept, n_kept, mean_variance, omp_get_num_procs());

        for (int sample = 0; sample < 2; sample++) {
            double mean_loss = 0;
            for (int row = 0; row < n_kept; row++) mean_loss += (double)AT(losses[0], row, sample);
            mean_loss /= n_kept;
            double uncorrected = mean_loss + mean_variance[sample];
            double target = 0;
            for (int h = 0; h < dim; h++) {
                double difference = mu[sample][h] - (double)benchmark_values[h];
                target += difference * difference;
            }
            target /= dim;
            corrected_average[sample] += mean_loss / n_benchmarks;
            uncorrected_average[sample] += uncorrected / n_benchmarks;
            if (fabs(mean_loss - target) > worst_miss) worst_miss = fabs(mean_loss - target);
            if (sample == 0 && fabs(uncorrected - target) < smallest_uncorrected_miss_of_truth)
                smallest_uncorrected_miss_of_truth = fabs(uncorrected - target);
        }
        free_losses(losses, 1);
    }
    fprintf(report, "      average over %d benchmarks, true A and wrong B:\n", n_benchmarks);
    fprintf(report, "      uncorrected  A %.4f  B %.4f   (expected 2.00 and 1.40)\n", uncorrected_average[0],
            uncorrected_average[1]);
    fprintf(report, "      corrected    A %.4f  B %.4f   (expected 1.00 and 1.36)\n", corrected_average[0],
            corrected_average[1]);
    fprintf(report, "      largest miss of the corrected mean loss from its target %.4f\n", worst_miss);
    fprintf(report, "      smallest miss of A's uncorrected mean loss from its target %.4f\n",
            smallest_uncorrected_miss_of_truth);
    check(report, uncorrected_average[1] < uncorrected_average[0],
          "the uncorrected loss prefers the wrong, less noisy configuration");
    check(report, corrected_average[0] < corrected_average[1], "the corrected loss prefers the true configuration");
    check(report, worst_miss < 0.1, "every corrected mean loss is within 0.1 of its squared distance");
    check(report, smallest_uncorrected_miss_of_truth > 0.5,
          "negative control: the uncorrected loss of A misses its target by more than 0.5");
    free(kept);
    free(benchmark_values);
    free(cache.data);
}

int main(void) {
    FILE *report = fopen(REPORT, "w");
    assert(report && "recovery_loss_correctness: cannot open the report");
    fprintf(report, "montecarlo/recovery_loss.h against its definition, %d hardware threads\n\n",
            omp_get_num_procs());
    fprintf(report, "1. against a reference from the definitions\n");
    against_reference(report);
    fprintf(report, "\n2. two kept runs, and a configuration with zero variance\n");
    smallest_and_constant(report);
    fprintf(report, "\n3. the bias of the uncorrected loss, and its removal\n");
    removes_the_bias(report);
    fprintf(report, "\n%s, %d failure(s)\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    return failures ? 1 : 0;
}
