#ifndef MONTECARLO_RESPONSE_CACHE_H
#define MONTECARLO_RESPONSE_CACHE_H

/*
Every simulated run's impulse responses under one auxiliary model, held in
memory as float32 so a loss matrix against any benchmark is one pass of
absolute differences over it. Shared by montecarlo/sweep_cops.c and
montecarlo/sweep_grid.c.

Four models: the t-QVARMA (p1q1r2, horizon 20, the stacked total responses,
525 entries per run, computed from the cached fits in
out/abm_system_fit_qvarma/) and the three local projections (linear, state 1,
state 2, 400 entries per run, read from out/abm_system_fit_lp/). A run with no
usable response is marked by a NaN in its first entry.

float32 for storage and double for the arithmetic, as in montecarlo/sweep_irf.c
and montecarlo/lp_sweep.c: 2.1 GB for the t-QVARMA, 1.6 GB for each local
projection.
*/

#include "applications/abm_system_lp.h"
#include <et_al./sd/qvarma.h>
#include <math.h>
#include <omp.h>

#define RESPONSE_QVARMA_K ABM_SYSTEM_K
#define RESPONSE_QVARMA_K_STAR 3
#define RESPONSE_QVARMA_R 1
#define RESPONSE_QVARMA_SHARED_BETA 1
#define RESPONSE_QVARMA_WARMUP_LONGEST 0
#define RESPONSE_QVARMA_MU_STAR_STATIONARY_ONLY 1
#define RESPONSE_QVARMA_P 1
#define RESPONSE_QVARMA_Q 1
#define RESPONSE_QVARMA_SPEC_R 2
#define RESPONSE_QVARMA_SPEC_LABEL "p1q1r2"
#define RESPONSE_QVARMA_HORIZON 20
#define RESPONSE_QVARMA_DIM (RESPONSE_QVARMA_K * RESPONSE_QVARMA_K * (RESPONSE_QVARMA_HORIZON + 1))
#define RESPONSE_QVARMA_FIT_DIR "out/abm_system_fit_qvarma"

enum { MODEL_QVARMA, MODEL_LP_LIN, MODEL_LP_S1, MODEL_LP_S2, N_MODELS };

static const char *const model_label[N_MODELS] = { "qvarma", "lp_lin", "lp_s1", "lp_s2" };

typedef struct {
    int model, n_samples, n_replicates, dim;
    float *data;
    long n_missing;
} ResponseCache;

static inline int response_dim(int model) {
    return model == MODEL_QVARMA ? RESPONSE_QVARMA_DIM : LP_RESPONSE_DIM;
}

static inline float *response_slot(const ResponseCache *cache, int sample, int replicate) {
    return cache->data + ((size_t)sample * cache->n_replicates + replicate) * cache->dim;
}

static inline QvarmaParams _response_qvarma_shape(void) {
    QvarmaParams m = qvarma_params_new(RESPONSE_QVARMA_K, RESPONSE_QVARMA_K_STAR, RESPONSE_QVARMA_P,
                                       RESPONSE_QVARMA_Q, RESPONSE_QVARMA_SPEC_R, RESPONSE_QVARMA_R,
                                       RESPONSE_QVARMA_SHARED_BETA, RESPONSE_QVARMA_WARMUP_LONGEST);
    m.mu_star_stationary_only = RESPONSE_QVARMA_MU_STAR_STATIONARY_ONLY;
    return m;
}

/* Fills out with the stacked total responses, horizon 0 first, and returns 1.
   Returns 0 when nu <= 2, which impulse_responses would abort on, or when a
   response comes back non-finite. The same as montecarlo/sweep_irf.c. */
static inline int _response_qvarma_cell(const QvarmaParams *m, Mat y, float *out) {
    if (!(m->nu > 2)) return 0;
    Mat D = qvarma_mean_score_jacobian(m, y);
    QvarmaImpulseOptions options = qvarma_default_impulse_options();
    options.horizon = RESPONSE_QVARMA_HORIZON;
    QvarmaImpulseResponses r = qvarma_impulse_responses(m, D, options);

    int at = 0, ok = 1;
    for (int h = 0; h <= r.horizon; h++)
        for (int i = 0; i < RESPONSE_QVARMA_K * RESPONSE_QVARMA_K; i++) {
            mreal value = r.total[h].d[i];
            if (MISNAN(value) || MISINF(value)) ok = 0;
            out[at++] = (float)value;
        }

    qvarma_impulse_responses_free(&r);
    mat_free(D);
    return ok;
}

/* Every run's response vector under one model, samples in the order given. */
static inline ResponseCache response_cache_build(int model, const LpSample *samples, int n_samples,
                                                 int n_replicates, const char *lp_fit_dir,
                                                 const char *input_dir) {
    ResponseCache cache = { model, n_samples, n_replicates, response_dim(model), NULL, 0 };
    size_t cells = (size_t)n_samples * n_replicates;
    cache.data = malloc(cells * cache.dim * sizeof(float));
    assert(cache.data && "response_cache: out of memory");
    for (size_t cell = 0; cell < cells; cell++) cache.data[cell * cache.dim] = (float)NAN;

    long n_missing = 0;
    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams working = _response_qvarma_shape();

        #pragma omp for schedule(dynamic)
        for (int sample = 0; sample < n_samples; sample++) {
            if (model == MODEL_QVARMA) {
                char dir[560];
                snprintf(dir, sizeof dir, "%s/%s", input_dir, samples[sample].name);
                int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
                for (int batch = 0; batch < n_batches; batch++) {
                    Mat block[ABM_SYSTEM_BATCH];
                    int replicate[ABM_SYSTEM_BATCH];
                    int count = abm_system_read_batch(dir, batch, block, replicate);
                    for (int b = 0; b < count; b++) {
                        if (replicate[b] < n_replicates) {
                            float *slot = response_slot(&cache, sample, replicate[b]);
                            char fit_path[640];
                            snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json",
                                     RESPONSE_QVARMA_FIT_DIR, samples[sample].name, replicate[b],
                                     RESPONSE_QVARMA_SPEC_LABEL);
                            int ok = qvarma_load_params(&working, fit_path)
                                  && _response_qvarma_cell(&working, block[b], slot);
                            if (!ok) slot[0] = (float)NAN;
                        }
                        mat_free(block[b]);
                    }
                }
            } else {
                int loss = model - MODEL_LP_LIN;
                LpConfiguration c = lp_configuration_load(lp_fit_dir, input_dir, samples[sample].name,
                                                          n_replicates, LP_LAYOUT_GROWTH);
                for (int replicate = 0; replicate < n_replicates; replicate++) {
                    if (!lp_configuration_ok(&c, replicate, loss)) continue;
                    const mreal *response = lp_response(c.row[replicate], loss);
                    float *slot = response_slot(&cache, sample, replicate);
                    for (int i = 0; i < LP_RESPONSE_DIM; i++) slot[i] = (float)response[i];
                }
                lp_configuration_free(&c);
            }

            for (int replicate = 0; replicate < n_replicates; replicate++)
                if (MISNAN(response_slot(&cache, sample, replicate)[0])) n_missing++;
            if (sample % 100 == 0)
                fprintf(stderr, "  %s cache: configuration %d of %d\n", model_label[model], sample, n_samples);
        }

        qvarma_params_free(&working);
    }

    cache.n_missing = n_missing;
    return cache;
}

/* The same cache, kept on disk under dir as <model>.f32 after it is first
   built, so a restarted run reads it back instead of rereading the fits and
   the dataset (about 20 minutes for a local projection, 10 for the t-QVARMA).
   The file starts with the model, configuration count, run count and entries
   per run as four ints, and is rebuilt when any of them differs. Delete the
   file to force a rebuild after the fits or the dataset change. */
static inline ResponseCache response_cache_load_or_build(int model, const LpSample *samples, int n_samples,
                                                         int n_replicates, const char *lp_fit_dir,
                                                         const char *input_dir, const char *dir) {
    char path[640];
    snprintf(path, sizeof path, "%s/%s.f32", dir, model_label[model]);
    int expected[4] = { model, n_samples, n_replicates, response_dim(model) };
    size_t values = (size_t)n_samples * n_replicates * expected[3];

    FILE *f = fopen(path, "rb");
    if (f) {
        int header[4];
        ResponseCache cache = { model, n_samples, n_replicates, expected[3], NULL, 0 };
        if (fread(header, sizeof(int), 4, f) == 4 && memcmp(header, expected, sizeof header) == 0) {
            cache.data = malloc(values * sizeof(float));
            assert(cache.data && "response_cache: out of memory");
            if (fread(cache.data, sizeof(float), values, f) == values) {
                fclose(f);
                for (int sample = 0; sample < n_samples; sample++)
                    for (int replicate = 0; replicate < n_replicates; replicate++)
                        if (MISNAN(response_slot(&cache, sample, replicate)[0])) cache.n_missing++;
                fprintf(stderr, "  %s cache read from %s\n", model_label[model], path);
                return cache;
            }
            free(cache.data);
        }
        fclose(f);
    }

    ResponseCache cache = response_cache_build(model, samples, n_samples, n_replicates, lp_fit_dir, input_dir);
    mkdir(dir, 0755);
    f = fopen(path, "wb");
    if (f) {
        size_t written = fwrite(expected, sizeof(int), 4, f);
        written += fwrite(cache.data, sizeof(float), values, f);
        fclose(f);
        if (written != values + 4) remove(path);
    }
    return cache;
}

static inline void response_cache_free(ResponseCache *cache) {
    free(cache->data);
    cache->data = NULL;
}

#endif /* MONTECARLO_RESPONSE_CACHE_H */
