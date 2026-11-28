/*
Whether the local projections on the differenced series are what they claim
to be: inflation and the interest rate, the pair the t-QVARMA treats as
co-integrated, in first differences, everything else as in the undifferenced
fits.

1. The series (applications/lp_system.h, lp_system_from_stored_differenced), on
   a block of 400 random periods: GDP growth, energy growth and employment
   change equal lp_system_from_stored's exactly; inflation and the interest
   rate at fitted period t equal stored_t - stored_(t-1), computed here as a
   plain loop, exactly; the switching series equals lp_system_from_stored's
   bit for bit; lp_system_fitted_series with the flag off is
   lp_system_from_stored.
2. The model table (montecarlo/response_cache.h): the four differenced models
   are labelled lp_<response>_differenced<weight suffix>, the linear one under
   lpirfs' weight and the three state models under the causal HP cycle at the
   estimated lambda, their fit directories carry _differenced, the state
   models a both-states model stacks are differenced too, and the grids of
   this folder run the linear and the both-states model while no other
   program runs any of them.
3. The stored fits, for runs 0 to 9 of cop_0001, cop_0500 and cop_1000: each
   row of out/applications/abm_system_fit_lp_differenced/ and
   out/applications/abm_system_fit_lp_differenced_causal_hp_estimated_lambda/
   equals, bit for bit, a fit computed here by lp_system_fit_record on the
   differenced series, and does not match the undifferenced series.
4. The response caches the grids read (out/montecarlo/sweep_grid_response_cache/),
   for all 1000 runs of the same three configurations: the differenced linear
   and both-states caches equal the archive's responses cast to float, both
   states being state 1 followed by state 2. Exact equality.

Writes
out/old/montecarlo_differenced_lp/lp_system_differenced_correctness_report.txt.
Exits 1 on any failure.
*/

#include "montecarlo/response_cache.h"
#include <et_al./random/random.h>
#include <cblas.h>
#include <stdio.h>

#define REPORT "out/old/montecarlo_differenced_lp/lp_system_differenced_correctness_report.txt"
#define T_STORED 400
#define N_CHECKED_SAMPLES 3
#define N_REFIT_RUNS 10

static const char *const checked_sample[N_CHECKED_SAMPLES] = { "cop_0001", "cop_0500", "cop_1000" };
static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    fflush(report);
    if (!ok) failures++;
}

static Mat stored_block(Rng *rng, int T) {
    Mat stored = mat_new(STORED_K, T);
    for (int t = 0; t < T; t++) {
        AT(stored, STORED_GDP_GROWTH, t) = (mreal)(0.5 + rng_normal(rng));
        AT(stored, STORED_EN_GROWTH, t) = (mreal)(0.1 + rng_normal(rng));
        AT(stored, STORED_EMPLOYMENT_CHANGE, t) = (mreal)(0.3 * rng_normal(rng));
        AT(stored, STORED_INFLATION, t) = (mreal)(0.5 + 0.2 * rng_normal(rng));
        AT(stored, STORED_INTEREST_RATE, t) = (mreal)(1.0 + 0.3 * rng_normal(rng));
    }
    return stored;
}

static int same_bits(Mat a, Mat b) {
    return a.r == b.r && a.c == b.c && memcmp(a.d, b.d, (size_t)a.r * a.c * sizeof(mreal)) == 0;
}

static void check_series(FILE *report) {
    Rng rng = rng_new(20261126, 0);
    Mat stored = stored_block(&rng, T_STORED);
    Mat y_level, switching_level, y_difference, switching_difference, y_off, switching_off;
    lp_system_from_stored(stored, &y_level, &switching_level);
    lp_system_from_stored_differenced(stored, &y_difference, &switching_difference);
    lp_system_fitted_series(stored, 0, &y_off, &switching_off);

    int first = LP_SYSTEM_SWITCHING_WINDOW - 1;
    int shape = y_difference.r == STORED_K && y_difference.c == T_STORED - first;
    int stationary_rows = shape, differenced_rows = shape;
    for (int t = first; t < T_STORED && shape; t++) {
        for (int k = STORED_GDP_GROWTH; k <= STORED_EMPLOYMENT_CHANGE; k++)
            stationary_rows &= AT(y_difference, k, t - first) == AT(stored, k, t);
        for (int k = STORED_INFLATION; k <= STORED_INTEREST_RATE; k++)
            differenced_rows &= AT(y_difference, k, t - first) == AT(stored, k, t) - AT(stored, k, t - 1);
    }
    check(report, shape, "1. shape: 5 x 397 from 400 stored periods, the undifferenced fit's");
    check(report, stationary_rows, "1. GDP growth, energy growth and employment change as stored, exactly");
    check(report, differenced_rows, "1. inflation and the interest rate are stored_t - stored_(t-1), exactly");
    check(report, same_bits(switching_difference, switching_level),
          "1. the switching series is the undifferenced fit's, bit for bit");
    check(report, same_bits(y_off, y_level) && same_bits(switching_off, switching_level),
          "1. with the flag off the fitted series are lp_system_from_stored's");
    mat_free(stored);
    mat_free(y_level);
    mat_free(switching_level);
    mat_free(y_difference);
    mat_free(switching_difference);
    mat_free(y_off);
    mat_free(switching_off);
}

static void check_model_table(FILE *report) {
    static const char *const part_name[4] = { "lin", "s1", "s2", "nl" };
    int labels_ok = 1, weights_ok = 1, siblings_ok = 1, dirs_ok = 1, selection_ok = 1, n_differenced = 0;
    for (int model = 0; model < N_MODELS; model++) {
        if (!model_is_differenced(model)) {
            selection_ok &= !model_in_differenced_lp_grids(model);
            continue;
        }
        n_differenced++;
        int part = lp_model_part(model);
        LpSystemStateWeight weight = lp_model_state_weight(model);
        char expected[MODEL_LABEL_CAPACITY];
        snprintf(expected, sizeof expected, "lp_%s_differenced%s", part_name[part],
                 lp_system_state_weight_suffix(weight));
        labels_ok &= strcmp(expected, model_label[model]) == 0;
        weights_ok &= part == LP_PART_LIN ? lp_system_state_weight_is_lpirfs(weight)
                                          : weight.causal_hp && weight.estimate_lambda;
        char dir[600], expected_dir[600];
        lp_model_fit_dir(model, LP_FIT_DIR_DEFAULT, dir, sizeof dir);
        snprintf(expected_dir, sizeof expected_dir, "%s_differenced%s", LP_FIT_DIR_DEFAULT,
                 lp_system_state_weight_suffix(weight));
        dirs_ok &= strcmp(dir, expected_dir) == 0;
        if (part != LP_PART_LIN)
            for (int other = LP_PART_S1; other <= LP_PART_NL; other++) {
                int sibling = lp_model_with_part(model, other);
                siblings_ok &= model_is_differenced(sibling) && lp_model_part(sibling) == other;
            }
        selection_ok &= !model_in_main_grid(model) && !model_is_single_state(model)
                     && !model_is_other_state_weight(model) && !model_is_mean_jacobian_qvarma(model)
                     && model_in_differenced_lp_grids(model) == (part == LP_PART_LIN || part == LP_PART_NL);
    }
    check(report, n_differenced == 4, "2. four differenced models");
    check(report, labels_ok, "2. labelled lp_<response>_differenced<weight suffix>");
    check(report, weights_ok, "2. linear under lpirfs' weight, the state models under the causal cycle at the "
                              "estimated lambda");
    check(report, dirs_ok, "2. fit directories carry _differenced before the weight's suffix");
    check(report, siblings_ok, "2. the state models a differenced both-states model stacks are differenced");
    check(report, selection_ok, "2. the grids of old/montecarlo_differenced_lp/ run the differenced linear and "
                                "both-states models, and no other program runs a differenced model");
    int distinct = 1;
    for (int a = 0; a < N_MODELS; a++)
        for (int b = a + 1; b < N_MODELS; b++) distinct &= strcmp(model_label[a], model_label[b]) != 0;
    check(report, distinct, "2. the labels are distinct");
}

static int same_record(const LpSystemRecord *a, const LpSystemRecord *b, int linear) {
    int same = a->replicate == b->replicate && a->has_responses == b->has_responses
            && a->var_ols_status == b->var_ols_status && a->var_rank == b->var_rank
            && a->var_chol_status == b->var_chol_status && a->switching_is_constant == b->switching_is_constant
            && memcmp(a->nl_ols_status, b->nl_ols_status, sizeof a->nl_ols_status) == 0
            && a->data_fingerprint == b->data_fingerprint && a->switching_fingerprint == b->switching_fingerprint
            && memcmp(&a->lambda, &b->lambda, sizeof a->lambda) == 0
            && memcmp(a->s1, b->s1, sizeof a->s1) == 0 && memcmp(a->s2, b->s2, sizeof a->s2) == 0;
    if (linear)
        same = same && memcmp(a->lin_ols_status, b->lin_ols_status, sizeof a->lin_ols_status) == 0
            && memcmp(a->d, b->d, sizeof a->d) == 0 && memcmp(a->lin, b->lin, sizeof a->lin) == 0;
    return same;
}

static void check_archives(FILE *report) {
    LpNlSpec spec = lp_system_spec();
    const int n_weights = 2;
    LpSystemStateWeight weights[2] = { lp_system_state_weight(0, 0), lp_system_state_weight(1, 1) };
    for (int w = 0; w < n_weights; w++) {
        int all_present = 1, all_equal = 1, none_matches_levels = 1, n_compared = 0;
        for (int s = 0; s < N_CHECKED_SAMPLES; s++) {
            char path[640];
            snprintf(path, sizeof path, "%s_differenced%s/%s.npz", LP_FIT_DIR_DEFAULT,
                     lp_system_state_weight_suffix(weights[w]), checked_sample[s]);
            int count = 0;
            LpSystemRecord *stored = lp_system_read_weighted_archive(path, spec, weights[w], &count);
            if (!stored || count < N_REFIT_RUNS) {
                fprintf(report, "  %s: no archive or fewer than %d rows\n", path, N_REFIT_RUNS);
                all_present = 0;
                free(stored);
                continue;
            }
            char dir[560];
            snprintf(dir, sizeof dir, "%s/%s", LP_INPUT_DIR_DEFAULT, checked_sample[s]);
            Mat block[ABM_SYSTEM_BATCH];
            int replicate[ABM_SYSTEM_BATCH];
            int n_block = abm_system_read_batch(dir, 0, block, replicate);
            for (int b = 0; b < n_block; b++) {
                if (replicate[b] < N_REFIT_RUNS && stored[replicate[b]].replicate == replicate[b]) {
                    Mat y, switching, y_level, switching_level;
                    lp_system_from_stored_differenced(block[b], &y, &switching);
                    lp_system_from_stored(block[b], &y_level, &switching_level);
                    LpSystemRecord refit = lp_system_fit_record(y, switching, spec, weights[w], replicate[b]);
                    int equal = same_record(&refit, &stored[replicate[b]], w == 0);
                    if (!equal) {
                        const LpSystemRecord *a = &refit, *c = &stored[replicate[b]];
                        double worst = 0;
                        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) {
                            worst = fmax(worst, fabs((double)a->s1[i] - (double)c->s1[i]));
                            worst = fmax(worst, fabs((double)a->s2[i] - (double)c->s2[i]));
                            if (w == 0) worst = fmax(worst, fabs((double)a->lin[i] - (double)c->lin[i]));
                        }
                        fprintf(report, "  %s run %d differs: fingerprints %d, notes %d, lambda %d, "
                                        "largest response gap %.3e\n",
                                checked_sample[s], replicate[b],
                                a->data_fingerprint == c->data_fingerprint
                                    && a->switching_fingerprint == c->switching_fingerprint,
                                a->var_ols_status == c->var_ols_status && a->var_rank == c->var_rank
                                    && a->var_chol_status == c->var_chol_status,
                                memcmp(&a->lambda, &c->lambda, sizeof a->lambda) == 0, worst);
                    }
                    all_equal &= equal;
                    none_matches_levels &= !lp_system_record_matches(&stored[replicate[b]], y_level, switching_level);
                    n_compared++;
                    mat_free(y);
                    mat_free(switching);
                    mat_free(y_level);
                    mat_free(switching_level);
                }
                mat_free(block[b]);
            }
            free(stored);
        }
        char what[256];
        snprintf(what, sizeof what, "3. %s: archives present for the three configurations",
                 lp_system_state_weight_description(weights[w]));
        check(report, all_present, what);
        snprintf(what, sizeof what, "3. %s: %d stored rows equal a fit on the differenced series bit for bit",
                 lp_system_state_weight_description(weights[w]), n_compared);
        check(report, all_equal && n_compared == N_CHECKED_SAMPLES * N_REFIT_RUNS, what);
        snprintf(what, sizeof what, "3. %s: no stored row matches the undifferenced series",
                 lp_system_state_weight_description(weights[w]));
        check(report, none_matches_levels && n_compared > 0, what);
    }
}

/* The cache file of a model mapped read-only, or data NULL when it is
   missing or its header disagrees. */
static ResponseCache map_cache(int model, int n_samples, int n_replicates) {
    char path[640];
    snprintf(path, sizeof path, "%s/%s.f32", RESPONSE_CACHE_DIR_DEFAULT, model_label[model]);
    int expected[4] = { model, n_samples, n_replicates, response_dim(model) };
    size_t length = sizeof expected + (size_t)n_samples * n_replicates * expected[3] * sizeof(float);
    ResponseCache cache = { model, n_samples, n_replicates, expected[3], NULL, 0, NULL, length };
    struct stat st;
    if (stat(path, &st) != 0 || (size_t)st.st_size != length) return cache;
    cache.mapping = _response_cache_map(path, length, 0);
    if (cache.mapping && memcmp(cache.mapping, expected, sizeof expected) == 0)
        cache.data = (float *)((char *)cache.mapping + sizeof expected);
    return cache;
}

static void check_caches(FILE *report) {
    int n_samples;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    const int checked_model[2] = { MODEL_LP_LIN_DIFFERENCED, MODEL_LP_NL_DIFFERENCED_CAUSAL_HP_ESTIMATED_LAMBDA };
    for (int m = 0; m < 2; m++) {
        int model = checked_model[m];
        ResponseCache cache = map_cache(model, n_samples, n_replicates);
        char what[256];
        snprintf(what, sizeof what, "4. %s: cache present with the expected header", model_label[model]);
        check(report, cache.data != NULL, what);
        if (!cache.data) {
            if (cache.mapping) munmap(cache.mapping, cache.mapping_length);
            continue;
        }
        LpSystemStateWeight weight = lp_model_state_weight(model);
        char fit_dir[600];
        lp_model_fit_dir(model, LP_FIT_DIR_DEFAULT, fit_dir, sizeof fit_dir);
        long n_equal = 0, n_cells = 0;
        for (int s = 0; s < N_CHECKED_SAMPLES; s++) {
            int sample = -1;
            for (int i = 0; i < n_samples; i++)
                if (strcmp(samples[i].name, checked_sample[s]) == 0) sample = i;
            assert(sample >= 0 && "lp_system_differenced_correctness: a checked configuration is not in the design");
            char path[640];
            snprintf(path, sizeof path, "%s/%s.npz", fit_dir, checked_sample[s]);
            int count = 0;
            LpSystemRecord *stored = lp_system_read_weighted_archive(path, lp_system_spec(), weight, &count);
            for (int i = 0; stored && i < count; i++) {
                const LpSystemRecord *r = &stored[i];
                if (r->replicate >= n_replicates || !r->has_responses) continue;
                const float *slot = response_slot(&cache, sample, r->replicate);
                int equal = 1;
                if (m == 0)
                    for (int h = 0; h < LP_RESPONSE_DIM; h++) equal &= slot[h] == (float)r->lin[h];
                else
                    for (int h = 0; h < LP_RESPONSE_DIM; h++)
                        equal &= slot[h] == (float)r->s1[h] && slot[LP_RESPONSE_DIM + h] == (float)r->s2[h];
                n_equal += equal;
                n_cells++;
            }
            free(stored);
        }
        snprintf(what, sizeof what, "4. %s: %ld of %ld runs of the three configurations equal the archive",
                 model_label[model], n_equal, n_cells);
        check(report, n_cells == (long)N_CHECKED_SAMPLES * n_replicates && n_equal == n_cells, what);
        response_cache_free(&cache);
    }
    lp_free_samples(samples, n_samples);
}

int main(void) {
    FILE *report = fopen(REPORT, "w");
    if (!report) {
        fprintf(stderr, "cannot open " REPORT "\n");
        return 1;
    }
    /* The fits were computed one configuration per thread, et_al's own
       parallel loops and OpenBLAS each at one thread. At 16 OpenMP threads the
       same fit differs from the stored one by up to 1e-9 in the responses
       (measured on runs 0 to 9 of cop_0001), so the refit runs at one. */
    openblas_set_num_threads(1);
    omp_set_num_threads(1);
    check_series(report);
    check_model_table(report);
    check_archives(report);
    check_caches(report);
    fprintf(report, "\n%s, %d failure%s\n", failures ? "FAILED" : "passed", failures, failures == 1 ? "" : "s");
    fclose(report);
    return failures ? 1 : 0;
}
