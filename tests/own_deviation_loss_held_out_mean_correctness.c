/*
Whether own_deviation_loss_held_out_mean_deviations of
montecarlo/own_deviation_loss.h computes, for every configuration i and every
deviation run s,

    S/(S+1) (1/K) sum_h (IR_ish - IR_bar_ih)^2,

IR_bar_ih the mean of entry h over the S mean runs, and whether that term is
unbiased for the mean response variance (1/K) sum_h sigma2_ih when s is not
among the mean runs.

1. Against a reference written here from the definition, as plain loops in
   double, on a synthetic cache with configurations at magnitudes from 1e-6 to
   1e6, the mean runs the last 4 runs, the deviation runs the others with one
   left out, at one thread and at every hardware thread.
2. One mean run: the mean is that run, the scale 1/2, so a configuration whose
   runs are all equal has a zero term.
3. The bias. One configuration, K = 20 entries, entry h normal with mean h and
   standard deviation 0.5 + h / 20, S = 5 mean runs and 200 deviation runs,
   2000 independent draws of the whole cache (seed 20261125). The average term
   over runs and draws must be within 1% of (1/K) sum_h sigma2_ih, while the
   term without the factor S/(S+1) must be above it by 1/S = 20%, within 1%.

Writes out/tests/own_deviation_loss_held_out_mean_correctness_report.txt.
Exits 1 on any failure.
*/

#include "montecarlo/own_deviation_loss.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>

#define REPORT "out/tests/own_deviation_loss_held_out_mean_correctness_report.txt"

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

static double reference_term(const ResponseCache *cache, const int *mean_runs, int n_mean, int sample, int run) {
    const int dim = cache->dim;
    double own = 0;
    for (int h = 0; h < dim; h++) {
        double mean = 0;
        for (int r = 0; r < n_mean; r++) mean += (double)response_slot(cache, sample, mean_runs[r])[h];
        mean /= n_mean;
        double deviation = (double)response_slot(cache, sample, run)[h] - mean;
        own += deviation * deviation;
    }
    return (double)n_mean / (n_mean + 1) * own / dim;
}

static void against_reference(FILE *report) {
    const int n_samples = 7, n_replicates = 15, dim = 13;
    const double magnitude[7] = { 1, 1e6, 1e-6, 3.5, 1e3, 0.25, 40 };
    ResponseCache cache = cache_new(n_samples, n_replicates, dim);
    Rng rng = rng_new(20261125, 0);
    for (int sample = 0; sample < n_samples; sample++)
        for (int replicate = 0; replicate < n_replicates; replicate++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, sample, replicate)[h] =
                    (float)(magnitude[sample] * (h % 3 + rng_normal(&rng)));

    /* Runs 11 to 14 give the mean; run 5 is left out of the deviation runs. */
    int mean_runs[4] = { 11, 12, 13, 14 }, deviation_runs[10], n_deviation = 0;
    for (int r = 0; r < 11; r++) if (r != 5) deviation_runs[n_deviation++] = r;

    double at_count[2][7 * 10];
    int thread_counts[2] = { 1, omp_get_num_procs() };
    for (int t = 0; t < 2; t++) {
        own_deviation_loss_held_out_mean_deviations(&cache, mean_runs, 4, deviation_runs, n_deviation, at_count[t],
                                                    thread_counts[t]);
        double worst = 0;
        for (int sample = 0; sample < n_samples; sample++)
            for (int row = 0; row < n_deviation; row++) {
                double expected = reference_term(&cache, mean_runs, 4, sample, deviation_runs[row]);
                double error = fabs(at_count[t][sample * n_deviation + row] - expected) / fabs(expected);
                if (error > worst) worst = error;
            }
        fprintf(report, "      %d thread(s): largest relative error against the reference %.3e\n",
                thread_counts[t], worst);
        char what[160];
        snprintf(what, sizeof what, "term matches the reference to 1e-9 relative at %d thread(s)", thread_counts[t]);
        check(report, worst < 1e-9, what);
    }
    int identical = 1;
    for (int at = 0; at < n_samples * n_deviation; at++)
        if (fabs(at_count[0][at] - at_count[1][at]) > 1e-12 * fabs(at_count[0][at])) identical = 0;
    check(report, identical, "one thread and every thread give the same terms");
    free(cache.data);
}

static void one_mean_run(FILE *report) {
    const int n_samples = 2, n_replicates = 3, dim = 4;
    ResponseCache cache = cache_new(n_samples, n_replicates, dim);
    Rng rng = rng_new(20261125, 1);
    for (int replicate = 0; replicate < n_replicates; replicate++)
        for (int h = 0; h < dim; h++) {
            response_slot(&cache, 0, replicate)[h] = (float)rng_normal(&rng);
            response_slot(&cache, 1, replicate)[h] = (float)(h - 1.5);
        }
    int mean_runs[1] = { 2 }, deviation_runs[2] = { 0, 1 };
    double term[4];
    own_deviation_loss_held_out_mean_deviations(&cache, mean_runs, 1, deviation_runs, 2, term, omp_get_num_procs());
    int ok = 1;
    for (int row = 0; row < 2; row++) {
        double half_squared = 0;
        for (int h = 0; h < dim; h++) {
            double difference = (double)response_slot(&cache, 0, deviation_runs[row])[h]
                              - (double)response_slot(&cache, 0, 2)[h];
            half_squared += 0.5 * difference * difference / dim;
        }
        if (fabs(term[row] - half_squared) > 1e-12 * half_squared) ok = 0;
    }
    check(report, ok, "one mean run: half the squared distance from that run");
    check(report, term[2] == 0 && term[3] == 0, "a configuration with equal runs has a zero term");
    free(cache.data);
}

static void bias(FILE *report) {
    const int dim = 20, n_mean = 5, n_deviation = 200, n_draws = 2000;
    const int n_replicates = n_mean + n_deviation;
    ResponseCache cache = cache_new(1, n_replicates, dim);
    Rng rng = rng_new(20261125, 2);
    int mean_runs[5], deviation_runs[200];
    for (int r = 0; r < n_deviation; r++) deviation_runs[r] = r;
    for (int r = 0; r < n_mean; r++) mean_runs[r] = n_deviation + r;

    double mean_variance = 0;
    for (int h = 0; h < dim; h++) mean_variance += (0.5 + h / 20.0) * (0.5 + h / 20.0) / dim;

    double *term = malloc((size_t)n_deviation * sizeof(double));
    assert(term);
    double total = 0;
    for (int draw = 0; draw < n_draws; draw++) {
        for (int replicate = 0; replicate < n_replicates; replicate++)
            for (int h = 0; h < dim; h++)
                response_slot(&cache, 0, replicate)[h] = (float)(h + (0.5 + h / 20.0) * rng_normal(&rng));
        own_deviation_loss_held_out_mean_deviations(&cache, mean_runs, n_mean, deviation_runs, n_deviation, term,
                                                    omp_get_num_procs());
        for (int row = 0; row < n_deviation; row++) total += term[row];
    }
    double average = total / ((double)n_draws * n_deviation);
    double unscaled = average * (n_mean + 1) / n_mean;
    fprintf(report, "      mean response variance %.6f, average term %.6f (ratio %.4f), without S/(S+1) %.6f "
                    "(ratio %.4f)\n", mean_variance, average, average / mean_variance, unscaled,
            unscaled / mean_variance);
    check(report, fabs(average / mean_variance - 1) < 0.01, "the term is unbiased for the mean response variance");
    check(report, fabs(unscaled / mean_variance - (1 + 1.0 / n_mean)) < 0.01,
          "without S/(S+1) the term is biased upwards by 1/S");
    free(term);
    free(cache.data);
}

int main(void) {
    mkdir("out/tests", 0755);
    FILE *report = fopen(REPORT, "w");
    assert(report && "own_deviation_loss_held_out_mean_correctness: cannot open the report");
    fprintf(report, "1. against the reference\n");
    against_reference(report);
    fprintf(report, "2. one mean run\n");
    one_mean_run(report);
    fprintf(report, "3. bias\n");
    bias(report);
    fprintf(report, "%s, %d failure(s)\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    printf("own_deviation_loss_held_out_mean_correctness: %s, %d failure(s), report in %s\n",
           failures ? "FAILED" : "passed", failures, REPORT);
    return failures ? 1 : 0;
}
