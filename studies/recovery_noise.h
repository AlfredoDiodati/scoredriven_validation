#ifndef STUDIES_RECOVERY_NOISE_H
#define STUDIES_RECOVERY_NOISE_H

#include <et_al./linalg/mat.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
Whether an auxiliary model's impulse responses are precise enough for the
Monte Carlo protocol to tell configurations apart, measured the same way for
any model, so the local projections (studies/lp_recovery_diagnostics.c) and
the t-QVARMA (studies/qvarma_recovery_diagnostics.c) are compared on the same
arithmetic.

The protocol compares every run's response vector with one stand-in run's by
the mean absolute difference over the vector's entries. A configuration's mean
loss then mixes two things: how far its average response is from the
stand-in's (distance), and how widely its own runs scatter around that average
(scatter). Both are measured per configuration, as mean absolute differences
over the entries, together with:

- the Spearman correlation of the mean loss with each, over configurations;
- per entry, the variance of the configurations' average responses over the
  average variance of a configuration's runs around its own average, which is
  below one when a single run says less about its configuration than the
  configurations differ from each other;
- the distance of an all-zero response from the stand-in, against every
  configuration's mean loss.

responses holds n_samples x n_replicates vectors of dim floats, configuration
major, a NaN in the first entry marking a run with no response. The stand-in's
own run index is left out of every configuration, as the protocol leaves it out.
*/

typedef struct {
    const float *responses;
    int dim, n_samples, n_replicates;
    int stand_in_sample, stand_in_replicate;
    const char *const *names;
} RecoveryResponses;

static inline const float *recovery_slot(const RecoveryResponses *r, int sample, int replicate) {
    return r->responses + ((size_t)sample * r->n_replicates + replicate) * r->dim;
}

static inline int _recovery_used(const RecoveryResponses *r, int sample, int replicate) {
    return replicate != r->stand_in_replicate && !MISNAN(recovery_slot(r, sample, replicate)[0]);
}

/* Writes the measurements to f and returns the configuration with the
   smallest mean loss. reference is the stand-in's response vector. */
static inline int recovery_noise_report(FILE *f, const RecoveryResponses *r, const double *reference) {
    int n = r->n_samples, dim = r->dim;
    Vec mean_loss = vec_new(n), distance = vec_new(n), scatter = vec_new(n);
    double *center_all = calloc((size_t)n * dim, sizeof(double));
    double *within_by_config = calloc((size_t)n * dim, sizeof(double));
    long *runs = calloc((size_t)n, sizeof(long));

    #pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < n; c++) {
        double *center = center_all + (size_t)c * dim;
        double *within = within_by_config + (size_t)c * dim;
        double loss = 0;
        for (int rep = 0; rep < r->n_replicates; rep++) {
            if (!_recovery_used(r, c, rep)) continue;
            const float *x = recovery_slot(r, c, rep);
            double l = 0;
            for (int i = 0; i < dim; i++) {
                center[i] += x[i];
                l += fabs((double)x[i] - reference[i]);
            }
            loss += l / dim;
            runs[c]++;
        }
        for (int i = 0; i < dim; i++) center[i] /= runs[c];
        double s = 0, d = 0;
        for (int rep = 0; rep < r->n_replicates; rep++) {
            if (!_recovery_used(r, c, rep)) continue;
            const float *x = recovery_slot(r, c, rep);
            for (int i = 0; i < dim; i++) {
                double e = (double)x[i] - center[i];
                s += fabs(e);
                within[i] += e * e;
            }
        }
        for (int i = 0; i < dim; i++) d += fabs(center[i] - reference[i]);
        AT(mean_loss, c, 0) = (mreal)(loss / runs[c]);
        AT(distance, c, 0) = (mreal)(d / dim);
        AT(scatter, c, 0) = (mreal)(s / dim / runs[c]);
    }

    long total_runs = 0;
    for (int c = 0; c < n; c++) total_runs += runs[c];
    int n_varying = 0;
    Vec ratio = vec_new(dim);
    for (int i = 0; i < dim; i++) {
        double grand = 0, between = 0, within = 0;
        for (int c = 0; c < n; c++) grand += center_all[(size_t)c * dim + i];
        grand /= n;
        for (int c = 0; c < n; c++) {
            double e = center_all[(size_t)c * dim + i] - grand;
            between += e * e;
            within += within_by_config[(size_t)c * dim + i];
        }
        between /= n - 1;
        within /= total_runs;
        /* Entries fixed by construction, such as the impact of a shock on
           itself under unit shocks, vary in no fit and carry no information. */
        if (within > 0) AT(ratio, n_varying++, 0) = (mreal)(between / within);
    }
    Mat varying = mat_slice(ratio, 0, n_varying, 0, 1);

    double zero = 0;
    for (int i = 0; i < dim; i++) zero += fabs(reference[i]);
    zero /= dim;
    int beat_zero = 0, winner = 0;
    for (int c = 0; c < n; c++) {
        beat_zero += (double)AT(mean_loss, c, 0) < zero;
        if (AT(mean_loss, c, 0) < AT(mean_loss, winner, 0)) winner = c;
    }

    int tracked[2] = { r->stand_in_sample, winner };
    for (int t = 0; t < (winner == r->stand_in_sample ? 1 : 2); t++) {
        int c = tracked[t];
        fprintf(f, "  %s: mean loss %.4g, distance %.4g, scatter %.4g%s\n", r->names[c],
                (double)AT(mean_loss, c, 0), (double)AT(distance, c, 0), (double)AT(scatter, c, 0),
                t == 1 ? " (smallest mean loss)"
                       : winner == c ? " (the stand-in's configuration, and the smallest mean loss)"
                                     : " (the stand-in's configuration)");
    }
    fprintf(f, "  over the %d configurations: distance median %.4g (min %.4g, max %.4g), scatter median %.4g "
               "(min %.4g, max %.4g)\n", n,
            (double)stats_median(distance), (double)stats_quantile(distance, 0), (double)stats_quantile(distance, 1),
            (double)stats_median(scatter), (double)stats_quantile(scatter, 0), (double)stats_quantile(scatter, 1));
    fprintf(f, "  Spearman correlation over configurations: mean loss with distance %.3f, with scatter %.3f\n",
            (double)stats_spearman(mean_loss, distance), (double)stats_spearman(mean_loss, scatter));
    fprintf(f, "  variance between configurations over variance within, per entry (%d of %d entries vary): "
               "median %.3g, 10th percentile %.3g, 90th percentile %.3g\n", n_varying, dim,
            (double)stats_median(varying), (double)stats_quantile(varying, 0.1),
            (double)stats_quantile(varying, 0.9));
    fprintf(f, "  an all-zero response is at %.4g from the stand-in; %d of %d configurations have a smaller mean loss\n",
            zero, beat_zero, n);

    mat_free(mean_loss); mat_free(distance); mat_free(scatter); mat_free(ratio);
    free(center_all); free(within_by_config); free(runs);
    return winner;
}

#endif /* STUDIES_RECOVERY_NOISE_H */
