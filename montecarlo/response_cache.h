#ifndef MONTECARLO_RESPONSE_CACHE_H
#define MONTECARLO_RESPONSE_CACHE_H

/*
Every simulated run's impulse responses under one auxiliary model, held in
memory as float32 so a loss matrix against any benchmark is one pass of
absolute differences over it. Shared by montecarlo/sweep_cops.c and
montecarlo/sweep_grid.c.

Five models: the t-QVARMA (p1q1r2, horizon 20, the stacked total responses,
525 entries per run, computed from the cached fits in
out/abm_system_fit_qvarma/), the three local projections (linear, state 1,
state 2, 400 entries per run, read from out/abm_system_fit_lp/), and lp_nl, the
state-dependent local projection with its two states in one vector, state 1's
400 entries followed by state 2's (_temp/Note on non-lin LP.pdf, equation 5).
The mean absolute difference over lp_nl's 800 entries is the average of the
two states' mean absolute differences. A run with no usable response, under
lp_nl a run missing either state, is marked by a NaN in its first entry.

float32 for storage and double for the arithmetic, as in montecarlo/sweep_irf.c:
2.1 GB for the t-QVARMA, 1.6 GB for each local projection, 3.2 GB for lp_nl.
The cache is a file mapped into memory rather than read into it, so the kernel
can drop pages it needs elsewhere and read them back, instead of the machine
running out of memory with a cache larger than the memory left.
*/

#include "applications/abm_system_lp.h"
#include <et_al./sd/qvarma.h>
#include <math.h>
#include <omp.h>
#include <sys/mman.h>

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
#define RESPONSE_CACHE_DIR_DEFAULT "out/sweep_grid_response_cache"

/* New models go at the end: a cache file's header records the model's index. */
enum { MODEL_QVARMA, MODEL_LP_LIN, MODEL_LP_S1, MODEL_LP_S2, MODEL_LP_NL, N_MODELS };

static const char *const model_label[N_MODELS] = { "qvarma", "lp_lin", "lp_s1", "lp_s2", "lp_nl" };

/* A local-projection model's label without its "lp_" prefix, as the
   local-projection programs name their files: lin, s1, s2, nl. */
static inline const char *lp_model_name(int model) {
    assert(model >= MODEL_LP_LIN && model < N_MODELS);
    return model_label[model] + 3;
}

typedef struct {
    int model, n_samples, n_replicates, dim;
    float *data;
    long n_missing;
    void *mapping;
    size_t mapping_length;
} ResponseCache;

static inline int response_dim(int model) {
    return model == MODEL_QVARMA ? RESPONSE_QVARMA_DIM : model == MODEL_LP_NL ? 2 * LP_RESPONSE_DIM : LP_RESPONSE_DIM;
}

/* One replicate's response under a local-projection model, response_dim(model)
   entries written to out. Returns 0 when the response, or under lp_nl either
   state's response, is missing or not finite. */
static inline int lp_model_response(const LpConfiguration *c, int replicate, int model, mreal *out) {
    int first = model == MODEL_LP_NL ? LP_LOSS_S1 : model - MODEL_LP_LIN;
    int last = model == MODEL_LP_NL ? LP_LOSS_S2 : first;
    for (int loss = first; loss <= last; loss++) {
        if (!lp_configuration_ok(c, replicate, loss)) return 0;
        memcpy(out + (size_t)(loss - first) * LP_RESPONSE_DIM, lp_response(c->row[replicate], loss),
               LP_RESPONSE_DIM * sizeof(mreal));
    }
    return 1;
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

/* Every run's response vector under one model, samples in the order given,
   written into cache->data, which holds n_samples x n_replicates x dim floats.
   Counts the runs with no usable response into cache->n_missing. */
static inline void _response_cache_fill(ResponseCache *cache, const LpSample *samples, const char *lp_fit_dir,
                                        const char *input_dir) {
    const int model = cache->model, n_samples = cache->n_samples, n_replicates = cache->n_replicates;
    size_t cells = (size_t)n_samples * n_replicates;
    for (size_t cell = 0; cell < cells; cell++) cache->data[cell * cache->dim] = (float)NAN;

    long n_missing = 0;
    #pragma omp parallel reduction(+:n_missing)
    {
        QvarmaParams working = _response_qvarma_shape();
        mreal *response = malloc((size_t)cache->dim * sizeof(mreal));
        assert(response && "response_cache: out of memory");

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
                            float *slot = response_slot(cache, sample, replicate[b]);
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
                LpConfiguration c = lp_configuration_load(lp_fit_dir, input_dir, samples[sample].name,
                                                          n_replicates, LP_LAYOUT_GROWTH);
                for (int replicate = 0; replicate < n_replicates; replicate++) {
                    if (!lp_model_response(&c, replicate, model, response)) continue;
                    float *slot = response_slot(cache, sample, replicate);
                    for (int i = 0; i < cache->dim; i++) slot[i] = (float)response[i];
                }
                lp_configuration_free(&c);
            }

            for (int replicate = 0; replicate < n_replicates; replicate++)
                if (MISNAN(response_slot(cache, sample, replicate)[0])) n_missing++;
            if (sample % 100 == 0)
                fprintf(stderr, "  %s cache: configuration %d of %d\n", model_label[model], sample, n_samples);
        }

        free(response);
        qvarma_params_free(&working);
    }
    cache->n_missing = n_missing;
}

/* The file mapped whole, or NULL when it cannot be opened or mapped. */
static inline void *_response_cache_map(const char *path, size_t length, int writable) {
    FILE *f = fopen(path, writable ? "r+b" : "rb");
    if (!f) return NULL;
    void *mapping = mmap(NULL, length, writable ? PROT_READ | PROT_WRITE : PROT_READ, MAP_SHARED, fileno(f), 0);
    fclose(f);
    return mapping == MAP_FAILED ? NULL : mapping;
}

static inline void response_cache_free(ResponseCache *cache) {
    munmap(cache->mapping, cache->mapping_length);
    cache->mapping = NULL;
    cache->data = NULL;
}

/* lp_nl's cache from the two state caches, each opened or built as any other,
   so the archives are not read a second time for the same responses. A run
   missing either state is missing. */
static inline ResponseCache response_cache_load_or_build(int model, const LpSample *samples, int n_samples,
                                                         int n_replicates, const char *lp_fit_dir,
                                                         const char *input_dir, const char *dir);

static inline void _response_cache_stack_states(ResponseCache *cache, const LpSample *samples,
                                                const char *lp_fit_dir, const char *input_dir, const char *dir) {
    ResponseCache state[2];
    for (int k = 0; k < 2; k++)
        state[k] = response_cache_load_or_build(MODEL_LP_S1 + k, samples, cache->n_samples, cache->n_replicates,
                                                lp_fit_dir, input_dir, dir);
    long n_missing = 0;
    for (int sample = 0; sample < cache->n_samples; sample++)
        for (int replicate = 0; replicate < cache->n_replicates; replicate++) {
            float *slot = response_slot(cache, sample, replicate);
            for (int k = 0; k < 2; k++)
                memcpy(slot + k * LP_RESPONSE_DIM, response_slot(&state[k], sample, replicate),
                       LP_RESPONSE_DIM * sizeof(float));
            if (MISNAN(slot[0]) || MISNAN(slot[LP_RESPONSE_DIM])) {
                slot[0] = (float)NAN;
                n_missing++;
            }
        }
    for (int k = 0; k < 2; k++) response_cache_free(&state[k]);
    cache->n_missing = n_missing;
}

/* Every run's response vector under one model, kept on disk under dir as
   <model>.f32 and mapped into memory. The file starts with the model,
   configuration count, run count and entries per run as four ints, and is
   rebuilt when any of them or the file's length differs; a restarted run then
   maps it in seconds instead of rereading the fits and the dataset (about 20
   minutes for a local projection, 10 for the t-QVARMA; lp_nl is put together
   from the two state caches instead). A build goes to
   <model>.f32.partial, renamed once complete, so an interrupted build leaves no
   file that looks finished. Delete the file to force a rebuild after the fits
   or the dataset change. */
static inline ResponseCache response_cache_load_or_build(int model, const LpSample *samples, int n_samples,
                                                         int n_replicates, const char *lp_fit_dir,
                                                         const char *input_dir, const char *dir) {
    char path[640], partial[660];
    snprintf(path, sizeof path, "%s/%s.f32", dir, model_label[model]);
    snprintf(partial, sizeof partial, "%s.partial", path);
    int expected[4] = { model, n_samples, n_replicates, response_dim(model) };
    size_t length = sizeof expected + (size_t)n_samples * n_replicates * expected[3] * sizeof(float);
    ResponseCache cache = { model, n_samples, n_replicates, expected[3], NULL, 0, NULL, length };

    struct stat st;
    if (stat(path, &st) == 0 && (size_t)st.st_size == length) {
        cache.mapping = _response_cache_map(path, length, 0);
        if (cache.mapping && memcmp(cache.mapping, expected, sizeof expected) == 0) {
            cache.data = (float *)((char *)cache.mapping + sizeof expected);
            for (int sample = 0; sample < n_samples; sample++)
                for (int replicate = 0; replicate < n_replicates; replicate++)
                    if (MISNAN(response_slot(&cache, sample, replicate)[0])) cache.n_missing++;
            fprintf(stderr, "  %s cache read from %s\n", model_label[model], path);
            return cache;
        }
        if (cache.mapping) munmap(cache.mapping, length);
    }

    mkdir(dir, 0755);
    FILE *f = fopen(partial, "wb");
    assert(f && "response_cache: cannot create the cache file");
    int sized = fseek(f, (long)length - 1, SEEK_SET) == 0 && fputc(0, f) == 0;
    sized = fclose(f) == 0 && sized;
    assert(sized && "response_cache: cannot size the cache file");
    (void)sized;
    cache.mapping = _response_cache_map(partial, length, 1);
    assert(cache.mapping && "response_cache: cannot map the cache file");
    cache.data = (float *)((char *)cache.mapping + sizeof expected);
    if (model == MODEL_LP_NL) _response_cache_stack_states(&cache, samples, lp_fit_dir, input_dir, dir);
    else _response_cache_fill(&cache, samples, lp_fit_dir, input_dir);
    memcpy(cache.mapping, expected, sizeof expected);
    int synced = msync(cache.mapping, length, MS_SYNC) == 0 && rename(partial, path) == 0;
    assert(synced && "response_cache: cannot write the cache file");
    (void)synced;
    return cache;
}

#endif /* MONTECARLO_RESPONSE_CACHE_H */
