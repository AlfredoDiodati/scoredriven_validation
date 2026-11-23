/*
Whether montecarlo/response_cache.h holds, under each of the four state
weights, exactly the state-dependent responses the fits of that weight stored.

1. The model table: every local-projection model's label is "lp_", its
   response (lin, s1, s2, nl) and its weight's directory suffix, and the state
   models a both-states model stacks are those of its own weight.
2. For three configurations spread over the design (the first, the middle and
   the last by index) and all their runs, the caches of state 1, state 2 and
   both states under every weight are built afresh, and every entry is compared
   with the archive of that weight read directly by
   lp_system_read_weighted_archive: state 1 and state 2 equal the stored
   responses cast to float, both states equal state 1 followed by state 2.
   Exact equality, no tolerance.
3. Under lpirfs' weight the fresh caches equal, byte for byte, the rows of the
   same configurations in the caches the grids read
   (out/montecarlo/sweep_grid_response_cache/).
4. Under each other weight no run's state 1 or state 2 response equals the one
   under lpirfs' weight entry for entry, so the cache did not read lpirfs'
   archives under another weight's name.

The caches are built in a fresh directory named after the start time,
out/tests/response_cache_state_weight_correctness/<YYYYmmdd_HHMMSS>/, about
77 MB, with the report beside them. Exits 1 on any failure.
*/

#include "montecarlo/response_cache.h"
#include <stdio.h>
#include <time.h>

#define OUTPUT_ROOT "out/tests/response_cache_state_weight_correctness"
#define PRODUCTION_CACHE_DIR RESPONSE_CACHE_DIR_DEFAULT
#define N_CHECKED_SAMPLES 3

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    fflush(report);
    if (!ok) failures++;
}

static const char *const part_name[4] = { "lin", "s1", "s2", "nl" };

static void check_model_table(FILE *report) {
    int labels_ok = 1, parts_ok = 1;
    for (int model = MODEL_LP_LIN; model < N_MODELS; model++) {
        char expected[MODEL_LABEL_CAPACITY];
        int part = lp_model_part(model);
        snprintf(expected, sizeof expected, "lp_%s%s", part_name[part],
                 lp_system_state_weight_suffix(lp_model_state_weight(model)));
        if (strcmp(expected, model_label[model]) != 0) {
            fprintf(report, "  model %d is labelled %s, its part and weight say %s\n", model, model_label[model],
                    expected);
            labels_ok = 0;
        }
        if (part == LP_PART_LIN) continue;
        for (int other = LP_PART_S1; other <= LP_PART_NL; other++) {
            int found = lp_model_with_part(model, other);
            LpSystemStateWeight a = lp_model_state_weight(model), b = lp_model_state_weight(found);
            if (lp_model_part(found) != other || a.causal_hp != b.causal_hp || a.estimate_lambda != b.estimate_lambda)
                parts_ok = 0;
        }
    }
    check(report, labels_ok, "every local-projection label is lp_<response><weight suffix>");
    check(report, parts_ok, "a model's sibling of each part has the same weight and that part");
    int distinct = 1;
    for (int a = 0; a < N_MODELS; a++)
        for (int b = a + 1; b < N_MODELS; b++) distinct &= strcmp(model_label[a], model_label[b]) != 0;
    check(report, distinct, "the fourteen labels are distinct");
    int fits = 1;
    for (int model = 0; model < N_MODELS; model++) fits &= strlen(model_label[model]) < MODEL_LABEL_CAPACITY;
    check(report, fits, "every label fits MODEL_LABEL_CAPACITY with its terminator");
}

/* The state models, four weights times state 1, state 2 and both states. */
static const int state_model[4][3] = {
    { MODEL_LP_S1, MODEL_LP_S2, MODEL_LP_NL },
    { MODEL_LP_S1_CAUSAL_HP, MODEL_LP_S2_CAUSAL_HP, MODEL_LP_NL_CAUSAL_HP },
    { MODEL_LP_S1_ESTIMATED_LAMBDA, MODEL_LP_S2_ESTIMATED_LAMBDA, MODEL_LP_NL_ESTIMATED_LAMBDA },
    { MODEL_LP_S1_CAUSAL_HP_ESTIMATED_LAMBDA, MODEL_LP_S2_CAUSAL_HP_ESTIMATED_LAMBDA,
      MODEL_LP_NL_CAUSAL_HP_ESTIMATED_LAMBDA },
};

int main(void) {
    char stamp[32], dir[256], report_path[320];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", localtime(&now));
    mkdir("out", 0755);
    mkdir("out/tests", 0755);
    mkdir(OUTPUT_ROOT, 0755);
    snprintf(dir, sizeof dir, "%s/%s", OUTPUT_ROOT, stamp);
    int created = mkdir(dir, 0755) == 0;
    assert(created && "response_cache_state_weight_correctness: the output directory exists already");
    (void)created;
    snprintf(report_path, sizeof report_path, "%s/response_cache_state_weight_correctness_report.txt", dir);
    FILE *report = fopen(report_path, "w");
    assert(report);

    check_model_table(report);

    int n_all = 0;
    LpSample *all = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_all);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, all[0].name);
    const int pick[N_CHECKED_SAMPLES] = { 0, n_all / 2, n_all - 1 };
    LpSample samples[N_CHECKED_SAMPLES];
    for (int k = 0; k < N_CHECKED_SAMPLES; k++) samples[k] = all[pick[k]];
    fprintf(report, "\nconfigurations %s, %s, %s, %d runs each, caches in %s\n", samples[0].name, samples[1].name,
            samples[2].name, n_replicates, dir);

    ResponseCache cache[4][3];
    for (int weight = 0; weight < 4; weight++)
        for (int part = 0; part < 3; part++)
            cache[weight][part] = response_cache_load_or_build(state_model[weight][part], samples, N_CHECKED_SAMPLES,
                                                               n_replicates, LP_FIT_DIR_DEFAULT,
                                                               LP_INPUT_DIR_DEFAULT, dir);

    for (int weight = 0; weight < 4; weight++) {
        int s1_model = state_model[weight][0];
        LpSystemStateWeight state_weight = lp_model_state_weight(s1_model);
        char fit_dir[600];
        lp_model_fit_dir(s1_model, LP_FIT_DIR_DEFAULT, fit_dir, sizeof fit_dir);
        long rows = 0, missing = 0, mismatched[3] = { 0, 0, 0 };
        for (int k = 0; k < N_CHECKED_SAMPLES; k++) {
            char path[700];
            snprintf(path, sizeof path, "%s/%s.npz", fit_dir, samples[k].name);
            int count = 0;
            LpSystemRecord *records = lp_system_read_weighted_archive(path, lp_system_spec(), state_weight, &count);
            assert(records && "response_cache_state_weight_correctness: an archive did not read");
            for (int replicate = 0; replicate < n_replicates; replicate++) {
                const LpSystemRecord *r = NULL;
                for (int i = 0; i < count && !r; i++)
                    if (records[i].replicate == replicate) r = &records[i];
                if (!r || !r->has_responses) { missing++; continue; }
                rows++;
                const float *s1 = response_slot(&cache[weight][0], k, replicate);
                const float *s2 = response_slot(&cache[weight][1], k, replicate);
                const float *nl = response_slot(&cache[weight][2], k, replicate);
                int bad[3] = { 0, 0, 0 };
                for (int h = 0; h < LP_RESPONSE_DIM; h++) {
                    bad[0] |= s1[h] != (float)r->s1[h];
                    bad[1] |= s2[h] != (float)r->s2[h];
                    bad[2] |= nl[h] != (float)r->s1[h] || nl[LP_RESPONSE_DIM + h] != (float)r->s2[h];
                }
                for (int p = 0; p < 3; p++) mismatched[p] += bad[p];
            }
            free(records);
        }
        char what[1024];
        snprintf(what, sizeof what, "%s: %ld runs read from %s, none missing from the archives (missing %ld)",
                 lp_system_state_weight_description(state_weight), rows, fit_dir, missing);
        check(report, missing == 0 && rows == (long)N_CHECKED_SAMPLES * n_replicates, what);
        for (int p = 0; p < 3; p++) {
            snprintf(what, sizeof what, "  %s cache equals the archive at every entry (runs that differ %ld, "
                                        "cells marked missing %ld)", model_label[state_model[weight][p]],
                     mismatched[p], cache[weight][p].n_missing);
            check(report, mismatched[p] == 0 && cache[weight][p].n_missing == 0, what);
        }
    }

    fprintf(report, "\n");
    for (int part = 0; part < 3; part++) {
        int model = state_model[0][part];
        ResponseCache production = response_cache_load_or_build(model, all, n_all, n_replicates,
                                                                LP_FIT_DIR_DEFAULT, LP_INPUT_DIR_DEFAULT,
                                                                PRODUCTION_CACHE_DIR);
        int equal = 1;
        for (int k = 0; k < N_CHECKED_SAMPLES; k++)
            equal &= memcmp(response_slot(&production, pick[k], 0), response_slot(&cache[0][part], k, 0),
                            (size_t)n_replicates * cache[0][part].dim * sizeof(float)) == 0;
        response_cache_free(&production);
        char what[300];
        snprintf(what, sizeof what, "%s: the fresh cache equals %s/%s.f32 byte for byte on these configurations",
                 model_label[model], PRODUCTION_CACHE_DIR, model_label[model]);
        check(report, equal, what);
    }

    fprintf(report, "\n");
    for (int weight = 1; weight < 4; weight++)
        for (int part = 0; part < 2; part++) {
            long identical = 0;
            for (int k = 0; k < N_CHECKED_SAMPLES; k++)
                for (int replicate = 0; replicate < n_replicates; replicate++)
                    identical += memcmp(response_slot(&cache[weight][part], k, replicate),
                                        response_slot(&cache[0][part], k, replicate),
                                        LP_RESPONSE_DIM * sizeof(float)) == 0;
            char what[300];
            snprintf(what, sizeof what, "%s differs from %s in every run (runs equal entry for entry %ld)",
                     model_label[state_model[weight][part]], model_label[state_model[0][part]], identical);
            check(report, identical == 0, what);
        }

    for (int weight = 0; weight < 4; weight++)
        for (int part = 0; part < 3; part++) response_cache_free(&cache[weight][part]);
    lp_free_samples(all, n_all);
    fprintf(report, "\n%s, %d failures\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    printf("%s, %d failures, report in %s\n", failures ? "FAILED" : "passed", failures, report_path);
    return failures ? 1 : 0;
}
