/*
Whether the t-QVARMA's impulse responses are precise enough, compared with the
differences between configurations, for the Monte Carlo protocol to tell
configurations apart: the same measurements studies/lp_recovery_diagnostics.c
makes for the local projections, through the same code
(studies/recovery_noise.h), so the two auxiliary models are compared on one
arithmetic.

The local projections fail to recover the stand-in's configuration because
their estimates scatter far more within a configuration than the
configurations' averages differ, and the protocol then ranks configurations by
their scatter (docs/MONTECARLO_LP_VALIDATION.md). If the t-QVARMA passes only
because its scatter is small, the failure belongs to the protocol rather than
to the local projections.

The responses are the ones montecarlo/irf_loss.c and montecarlo/sweep_irf.c
compare: the total response, horizons 0 to 20, 5 x 5 x 21 = 525 values, from
each run's stored p1q1r2 fit in out/abm_system_fit_qvarma/ and its own series.
A fit with nu <= 2 or a non-finite response counts as missing. The stand-in is
montecarlo/out/benchmark.env's run, cop_0191 run 706, left out of every
configuration. As a check that these are the protocol's own numbers, every
run's mean absolute difference from the stand-in is compared with
montecarlo/out/irf_loss.csv.

Requires out/abm_system_fit_qvarma/, dataset/abm_system/,
montecarlo/out/benchmark.env and montecarlo/out/irf_loss.csv, and rebuilds
none of them. Writes out/qvarma_recovery_diagnostics_report.txt. Progress to
stderr.
*/

#include "applications/abm_system.h"
#include "montecarlo/benchmark.h"
#include "studies/recovery_noise.h"
#include <et_al./sd/qvarma.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define P 1
#define Q 1
#define SPEC_R 2
#define SPEC_LABEL "p1q1r2"
#define HORIZON 20
#define DIM (K * K * (HORIZON + 1))

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define LOSS_TABLE "montecarlo/out/irf_loss.csv"
#define REPORT_PATH "out/qvarma_recovery_diagnostics_report.txt"

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = 1;
    return m;
}

static int compare_names(const void *a, const void *b) {
    return strcmp(*(char *const *)a, *(char *const *)b);
}

static char **list_samples(int *count) {
    DIR *handle = opendir(FIT_DIR);
    assert(handle && "qvarma_recovery_diagnostics: cannot open the t-QVARMA fits");
    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (strncmp(entry->d_name, "cop_", 4) != 0) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 64;
            names = realloc(names, (size_t)cap * sizeof(char *));
            assert(names && "qvarma_recovery_diagnostics: out of memory listing configurations");
        }
        size_t length = strlen(entry->d_name);
        names[n] = malloc(length + 1);
        memcpy(names[n], entry->d_name, length + 1);
        n++;
    }
    closedir(handle);
    qsort(names, (size_t)n, sizeof(char *), compare_names);
    *count = n;
    return names;
}

/* The stacked total response, horizon 0 first, as montecarlo/irf_loss.c
   builds it. Returns 0 for nu <= 2 or a non-finite value. */
static int total_response(const QvarmaParams *m, Mat y, double *out) {
    if (!(m->nu > 2)) return 0;
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = HORIZON;
    QvarmaImpulseResponses r = qvarma_impulse_responses(m, D, options);
    int at = 0, ok = 1;
    for (int h = 0; h <= r.horizon; h++)
        for (int i = 0; i < K * K; i++) {
            double value = (double)r.total[h].d[i];
            if (MISNAN(value) || MISINF(value)) ok = 0;
            out[at++] = value;
        }
    qvarma_impulse_responses_free(&r);
    mat_free(D);
    return ok;
}

int main(void) {
    openblas_set_num_threads(1);
    Benchmark stand_in = benchmark_read();

    int n_samples;
    char **samples = list_samples(&n_samples);
    int stand_in_sample = -1;
    for (int c = 0; c < n_samples; c++) if (strcmp(samples[c], stand_in.sample) == 0) stand_in_sample = c;
    assert(stand_in_sample >= 0 && "qvarma_recovery_diagnostics: the stand-in's configuration is missing");

    /* 1000 runs a configuration, the design's count, confirmed below against
       the loss table's rows. */
    int n_replicates = 1000;
    float *cache = malloc((size_t)n_samples * n_replicates * DIM * sizeof(float));
    assert(cache && "qvarma_recovery_diagnostics: out of memory for the responses");
    long n_missing = 0;

    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams working = spec_shape();
        double response[DIM];
        #pragma omp for schedule(dynamic)
        for (int c = 0; c < n_samples; c++) {
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, samples[c]);
            for (int r = 0; r < n_replicates; r++) cache[((size_t)c * n_replicates + r) * DIM] = (float)NAN;
            for (int batch = 0; batch * ABM_SYSTEM_BATCH < n_replicates; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);
                for (int b = 0; b < count; b++) {
                    int r = replicate[b];
                    char path[640];
                    snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, samples[c], r,
                             SPEC_LABEL);
                    if (r < n_replicates && qvarma_load_params(&working, path)
                        && total_response(&working, block[b], response)) {
                        float *slot = cache + ((size_t)c * n_replicates + r) * DIM;
                        for (int i = 0; i < DIM; i++) slot[i] = (float)response[i];
                    } else {
                        n_missing++;
                    }
                    mat_free(block[b]);
                }
            }
            if (c % 100 == 0) fprintf(stderr, "  responses: configuration %d of %d\n", c, n_samples);
        }
        qvarma_params_free(&working);
    }

    double reference[DIM];
    {
        const float *slot = cache + ((size_t)stand_in_sample * n_replicates + stand_in.replicate) * DIM;
        assert(!MISNAN(slot[0]) && "qvarma_recovery_diagnostics: the stand-in has no response");
        /* The stand-in in double, as montecarlo/irf_loss.c holds it, recomputed
           rather than taken from the float copy. */
        QvarmaParams m = spec_shape();
        char path[640], dir[560];
        snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, stand_in.sample,
                 stand_in.replicate, SPEC_LABEL);
        snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, stand_in.sample);
        Mat y = abm_system_read_replicate(dir, stand_in.replicate);
        int ok = qvarma_load_params(&m, path) && total_response(&m, y, reference);
        assert(ok && "qvarma_recovery_diagnostics: the stand-in's own response failed");
        mat_free(y);
        qvarma_params_free(&m);
    }

    FILE *f = fopen(REPORT_PATH, "w");
    assert(f && "qvarma_recovery_diagnostics: cannot open the report");
    fprintf(f, "How precise the t-QVARMA's impulse responses are for the Monte Carlo protocol\n\n");
    fprintf(f, "Data: dataset/abm_system, %d configurations x %d runs; each run's stored p1q1r2 fit in %s.\n",
            n_samples, n_replicates, FIT_DIR);
    fprintf(f, "Responses: total response, horizons 0 to %d, %d values per run, as montecarlo/irf_loss.c builds them.\n",
            HORIZON, DIM);
    fprintf(f, "Stand-in: %s run %d, left out of every configuration. Runs with no response: %ld.\n\n",
            stand_in.sample, stand_in.replicate, n_missing);

    /* The protocol's own loss table, against the mean absolute differences
       recomputed here from the float copies. */
    DataFrame table = df_read_csv(LOSS_TABLE, csv_read_options_default());
    Mat replicate = df_col_numeric(&table, "replicate");
    double worst = 0;
    for (int c = 0; c < n_samples; c++) {
        char name[128];
        snprintf(name, sizeof name, "%s_qvarma_%s", samples[c], SPEC_LABEL);
        Mat column = df_col_numeric(&table, name);
        for (int row = 0; row < table.r; row++) {
            int r = (int)AT(replicate, row, 0);
            const float *x = cache + ((size_t)c * n_replicates + r) * DIM;
            double s = 0;
            for (int i = 0; i < DIM; i++) s += fabs((double)x[i] - reference[i]);
            double gap = fabs(s / DIM - (double)AT(column, row, 0)) / fmax(1e-12, fabs((double)AT(column, row, 0)));
            if (gap > worst) worst = gap;
        }
    }
    fprintf(f, "Check against %s: %d rows x %d configurations, largest relative difference %.3g "
               "(the responses here are held as float).\n\n", LOSS_TABLE, table.r, n_samples, worst);
    df_free(&table);

    fprintf(f, "Noise against signal, stand-in %s run %d\n", stand_in.sample, stand_in.replicate);
    RecoveryResponses responses = { cache, DIM, n_samples, n_replicates, stand_in_sample, stand_in.replicate,
                                    (const char *const *)samples };
    recovery_noise_report(f, &responses, reference);
    fclose(f);

    for (int c = 0; c < n_samples; c++) free(samples[c]);
    free(samples);
    free(cache);
    return 0;
}
