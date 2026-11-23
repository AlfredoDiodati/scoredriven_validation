#ifndef MONTECARLO_RESPONSE_CACHE_H
#define MONTECARLO_RESPONSE_CACHE_H

/*
Every simulated run's impulse responses under one auxiliary model, held in
memory as float32 so a loss matrix against any benchmark is one pass over it.
Read by montecarlo/true_configuration_recovery.c through
montecarlo/recovery_loss.h, and by the discarded naive-loss programs under
old/naive_mae_loss/montecarlo/.

Fourteen models. The first five: the t-QVARMA (p1q1r2, horizon 20, the stacked
total responses, 525 entries per run, computed from the cached fits in
out/applications/abm_system_fit_qvarma/), the three local projections (linear,
state 1, state 2, 400 entries per run, read from
out/applications/abm_system_fit_lp/), and lp_nl, the state-dependent local
projection with its two states in one vector, state 1's 400 entries followed by
state 2's. Their state weight is lpirfs', the two-sided HP cycle at lambda
1600. The other nine are state 1, state 2 and both states under each of the
three other weights of applications/lp_system.h: the causal HP cycle at lambda
1600, the two-sided cycle at the estimated lambda, and the causal cycle at the
estimated lambda, read from the tree of fits carrying that weight's suffix,
out/applications/abm_system_fit_lp_causal_hp/ and so on. The linear projection
does not depend on the weight, so it has no model under the other three. A run
with no usable response, under a both-states model a run missing either state,
is marked by a NaN in its first entry.

float32 for storage and double for the arithmetic, so the only precision lost
is in the stored responses themselves: 2.1 GB for the t-QVARMA, 1.6 GB for each
local projection, 3.2 GB for lp_nl. The cache is a file mapped into memory
rather than read into it, so the kernel can drop pages it needs elsewhere and
read them back, instead of the machine running out of memory with a cache
larger than the memory left.
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
#define RESPONSE_QVARMA_FIT_DIR "out/applications/abm_system_fit_qvarma"
#define RESPONSE_CACHE_DIR_DEFAULT "out/montecarlo/sweep_grid_response_cache"

/* New models go at the end: a cache file's header records the model's index.
   N_LPIRFS_WEIGHT_MODELS counts the first five, the only models the naive-loss
   grid of old/naive_mae_loss/ holds. */
enum {
    MODEL_QVARMA, MODEL_LP_LIN, MODEL_LP_S1, MODEL_LP_S2, MODEL_LP_NL,
    MODEL_LP_S1_CAUSAL_HP, MODEL_LP_S2_CAUSAL_HP, MODEL_LP_NL_CAUSAL_HP,
    MODEL_LP_S1_ESTIMATED_LAMBDA, MODEL_LP_S2_ESTIMATED_LAMBDA, MODEL_LP_NL_ESTIMATED_LAMBDA,
    MODEL_LP_S1_CAUSAL_HP_ESTIMATED_LAMBDA, MODEL_LP_S2_CAUSAL_HP_ESTIMATED_LAMBDA,
    MODEL_LP_NL_CAUSAL_HP_ESTIMATED_LAMBDA,
    N_MODELS
};
#define N_LPIRFS_WEIGHT_MODELS (MODEL_LP_NL + 1)

/* A label is the response's name followed by the suffix of its state weight's
   tree of fits (lp_system_state_weight_suffix). */
static const char *const model_label[N_MODELS] = {
    "qvarma", "lp_lin", "lp_s1", "lp_s2", "lp_nl",
    "lp_s1_causal_hp", "lp_s2_causal_hp", "lp_nl_causal_hp",
    "lp_s1_estimated_lambda", "lp_s2_estimated_lambda", "lp_nl_estimated_lambda",
    "lp_s1_causal_hp_estimated_lambda", "lp_s2_causal_hp_estimated_lambda", "lp_nl_causal_hp_estimated_lambda",
};

/* Longer than any model label, for the buffers that parse one back from a
   csv. */
#define MODEL_LABEL_CAPACITY 64

/* Which response of a local-projection model: the linear one, state 1, state 2
   or both states stacked. */
enum { LP_PART_LIN, LP_PART_S1, LP_PART_S2, LP_PART_NL };

static inline int lp_model_part(int model) {
    assert(model >= MODEL_LP_LIN && model < N_MODELS);
    if (model <= MODEL_LP_NL) return LP_PART_LIN + (model - MODEL_LP_LIN);
    return LP_PART_S1 + (model - MODEL_LP_S1_CAUSAL_HP) % 3;
}

/* Whether a model is state 1 or state 2 of a state-dependent projection on its
   own. Its cache is what the both-states model of its weight is stacked from.
   The programs under montecarlo/ run and report every other model; those of
   old/montecarlo_single_states/ run and report these. */
static inline int model_is_single_state(int model) {
    if (model == MODEL_QVARMA) return 0;
    int part = lp_model_part(model);
    return part == LP_PART_S1 || part == LP_PART_S2;
}

/* The state weight a local-projection model was fitted under. */
static inline LpSystemStateWeight lp_model_state_weight(int model) {
    assert(model >= MODEL_LP_LIN && model < N_MODELS);
    if (model <= MODEL_LP_NL) return lp_system_state_weight(0, 0);
    int weight = (model - MODEL_LP_S1_CAUSAL_HP) / 3;
    return lp_system_state_weight(weight != 1, weight != 0);
}

/* The model with the same state weight and the given part, so a both-states
   model finds its two state models. */
static inline int lp_model_with_part(int model, int part) {
    assert(model >= MODEL_LP_LIN && model < N_MODELS);
    if (model <= MODEL_LP_NL) return MODEL_LP_LIN + (part - LP_PART_LIN);
    assert(part != LP_PART_LIN && "response_cache: the linear projection has no model under another weight");
    return model - (lp_model_part(model) - part);
}

/* Where a local-projection model's fits are: lp_fit_dir, the tree of lpirfs'
   weight, followed by the model's weight suffix. */
static inline void lp_model_fit_dir(int model, const char *lp_fit_dir, char *out, size_t capacity) {
    snprintf(out, capacity, "%s%s", lp_fit_dir, lp_system_state_weight_suffix(lp_model_state_weight(model)));
}

/* A local-projection model's label without its "lp_" prefix, as the
   local-projection programs name their files: lin, s1, s2, nl, s1_causal_hp
   and so on. */
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
    if (model == MODEL_QVARMA) return RESPONSE_QVARMA_DIM;
    return lp_model_part(model) == LP_PART_NL ? 2 * LP_RESPONSE_DIM : LP_RESPONSE_DIM;
}

/* One replicate's response under a local-projection model, response_dim(model)
   entries written to out, from a configuration loaded under the model's state
   weight. Returns 0 when the response, or under a both-states model either
   state's response, is missing or not finite. */
static inline int lp_model_response(const LpConfiguration *c, int replicate, int model, mreal *out) {
    static const int first_loss[] = { LP_LOSS_LIN, LP_LOSS_S1, LP_LOSS_S2, LP_LOSS_S1 };
    int part = lp_model_part(model);
    int first = first_loss[part];
    int last = part == LP_PART_NL ? LP_LOSS_S2 : first;
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
   response comes back non-finite. */
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
                char fit_dir[600];
                lp_model_fit_dir(model, lp_fit_dir, fit_dir, sizeof fit_dir);
                LpConfiguration c = lp_configuration_load_weighted(fit_dir, input_dir, samples[sample].name,
                                                                   n_replicates, lp_model_state_weight(model));
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

/* A both-states model's cache from the two state caches of the same weight,
   each opened or built as any other, so the archives are not read a second
   time for the same responses. A run missing either state is missing. */
static inline ResponseCache response_cache_load_or_build(int model, const LpSample *samples, int n_samples,
                                                         int n_replicates, const char *lp_fit_dir,
                                                         const char *input_dir, const char *dir);

static inline void _response_cache_stack_states(ResponseCache *cache, const LpSample *samples,
                                                const char *lp_fit_dir, const char *input_dir, const char *dir) {
    ResponseCache state[2];
    for (int k = 0; k < 2; k++)
        state[k] = response_cache_load_or_build(lp_model_with_part(cache->model, LP_PART_S1 + k), samples,
                                                cache->n_samples, cache->n_replicates, lp_fit_dir, input_dir, dir);
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
   minutes for a local projection, 10 for the t-QVARMA; a both-states model is
   put together from its two state caches instead). lp_fit_dir is the tree of
   lpirfs' weight; a model under another weight reads the tree with that
   weight's suffix. A build goes to
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
    if (model != MODEL_QVARMA && lp_model_part(model) == LP_PART_NL)
        _response_cache_stack_states(&cache, samples, lp_fit_dir, input_dir, dir);
    else _response_cache_fill(&cache, samples, lp_fit_dir, input_dir);
    memcpy(cache.mapping, expected, sizeof expected);
    int synced = msync(cache.mapping, length, MS_SYNC) == 0 && rename(partial, path) == 0;
    assert(synced && "response_cache: cannot write the cache file");
    (void)synced;
    return cache;
}

#endif /* MONTECARLO_RESPONSE_CACHE_H */
