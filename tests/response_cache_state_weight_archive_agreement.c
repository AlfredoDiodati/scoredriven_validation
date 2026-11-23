/*
Whether the response caches the grids read under the three state weights other
than lpirfs' (out/montecarlo/sweep_grid_response_cache/) hold, at every
configuration and run, exactly the responses the fits of that weight stored.

For each weight, every archive of that weight is read by
lp_system_read_weighted_archive and compared with the mapped caches: state 1
and state 2 equal the stored responses cast to float, both states equal state 1
followed by state 2, a run without responses in the archive is missing (NaN) in
all three caches, and no other entry is NaN. Exact equality, no tolerance.

The caches are only read. The report goes to
out/tests/response_cache_state_weight_archive_agreement/<YYYYmmdd_HHMMSS>/.
Exits 1 on any failure.
*/

#include "montecarlo/response_cache.h"
#include <stdio.h>
#include <time.h>

#define OUTPUT_ROOT "out/tests/response_cache_state_weight_archive_agreement"
#define PRODUCTION_CACHE_DIR RESPONSE_CACHE_DIR_DEFAULT

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    fflush(report);
    if (!ok) failures++;
}

/* State 1, state 2 and both states under the three other weights. */
static const int state_model[3][3] = {
    { MODEL_LP_S1_CAUSAL_HP, MODEL_LP_S2_CAUSAL_HP, MODEL_LP_NL_CAUSAL_HP },
    { MODEL_LP_S1_ESTIMATED_LAMBDA, MODEL_LP_S2_ESTIMATED_LAMBDA, MODEL_LP_NL_ESTIMATED_LAMBDA },
    { MODEL_LP_S1_CAUSAL_HP_ESTIMATED_LAMBDA, MODEL_LP_S2_CAUSAL_HP_ESTIMATED_LAMBDA,
      MODEL_LP_NL_CAUSAL_HP_ESTIMATED_LAMBDA },
};

int main(void) {
    char stamp[32], dir[256], report_path[400];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", localtime(&now));
    mkdir("out", 0755);
    mkdir("out/tests", 0755);
    mkdir(OUTPUT_ROOT, 0755);
    snprintf(dir, sizeof dir, "%s/%s", OUTPUT_ROOT, stamp);
    int created = mkdir(dir, 0755) == 0;
    assert(created && "response_cache_state_weight_archive_agreement: the output directory exists already");
    (void)created;
    snprintf(report_path, sizeof report_path, "%s/response_cache_state_weight_archive_agreement_report.txt", dir);
    FILE *report = fopen(report_path, "w");
    assert(report);

    int n_samples = 0;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    fprintf(report, "%d configurations, %d runs each, caches in %s\n\n", n_samples, n_replicates,
            PRODUCTION_CACHE_DIR);

    for (int weight = 0; weight < 3; weight++) {
        ResponseCache cache[3];
        for (int part = 0; part < 3; part++) {
            int model = state_model[weight][part];
            char path[700];
            snprintf(path, sizeof path, "%s/%s.f32", PRODUCTION_CACHE_DIR, model_label[model]);
            struct stat st;
            int present = stat(path, &st) == 0;
            check(report, present, path);
            if (!present) return 1;
            cache[part] = response_cache_load_or_build(model, samples, n_samples, n_replicates, LP_FIT_DIR_DEFAULT,
                                                       LP_INPUT_DIR_DEFAULT, PRODUCTION_CACHE_DIR);
        }
        int s1_model = state_model[weight][0];
        LpSystemStateWeight state_weight = lp_model_state_weight(s1_model);
        char fit_dir[600];
        lp_model_fit_dir(s1_model, LP_FIT_DIR_DEFAULT, fit_dir, sizeof fit_dir);

        long compared = 0, absent_in_archive = 0, absent_but_cached = 0, cached_missing_but_present = 0;
        long mismatched[3] = { 0, 0, 0 };
        for (int sample = 0; sample < n_samples; sample++) {
            char path[700];
            snprintf(path, sizeof path, "%s/%s.npz", fit_dir, samples[sample].name);
            int count = 0;
            LpSystemRecord *records = lp_system_read_weighted_archive(path, lp_system_spec(), state_weight, &count);
            assert(records && "response_cache_state_weight_archive_agreement: an archive did not read");
            for (int replicate = 0; replicate < n_replicates; replicate++) {
                const LpSystemRecord *r = NULL;
                for (int i = 0; i < count && !r; i++)
                    if (records[i].replicate == replicate) r = &records[i];
                const float *s1 = response_slot(&cache[0], sample, replicate);
                const float *s2 = response_slot(&cache[1], sample, replicate);
                const float *nl = response_slot(&cache[2], sample, replicate);
                if (!r || !r->has_responses) {
                    absent_in_archive++;
                    if (!MISNAN(s1[0]) || !MISNAN(s2[0]) || !MISNAN(nl[0])) absent_but_cached++;
                    continue;
                }
                compared++;
                int bad[3] = { 0, 0, 0 };
                for (int h = 0; h < LP_RESPONSE_DIM; h++) {
                    bad[0] |= s1[h] != (float)r->s1[h];
                    bad[1] |= s2[h] != (float)r->s2[h];
                    bad[2] |= nl[h] != (float)r->s1[h] || nl[LP_RESPONSE_DIM + h] != (float)r->s2[h];
                }
                for (int p = 0; p < 3; p++) mismatched[p] += bad[p];
                if (MISNAN(s1[0]) || MISNAN(s2[0]) || MISNAN(nl[0])) cached_missing_but_present++;
            }
            free(records);
            if (sample % 100 == 0)
                fprintf(stderr, "  %s: configuration %d of %d\n", lp_system_state_weight_description(state_weight),
                        sample, n_samples);
        }

        char what[1200];
        snprintf(what, sizeof what, "%s: %ld runs compared from %s, %ld without responses in the archive",
                 lp_system_state_weight_description(state_weight), compared, fit_dir, absent_in_archive);
        check(report, compared + absent_in_archive == (long)n_samples * n_replicates, what);
        snprintf(what, sizeof what, "  runs without responses are missing in all three caches (runs that are not %ld)",
                 absent_but_cached);
        check(report, absent_but_cached == 0, what);
        snprintf(what, sizeof what, "  runs with responses are not missing in any cache (runs that are %ld)",
                 cached_missing_but_present);
        check(report, cached_missing_but_present == 0, what);
        for (int p = 0; p < 3; p++) {
            snprintf(what, sizeof what, "  %s equals the archives at every entry (runs that differ %ld, cells "
                                        "marked missing %ld)", model_label[state_model[weight][p]], mismatched[p],
                     cache[p].n_missing);
            check(report, mismatched[p] == 0 && cache[p].n_missing == absent_in_archive, what);
        }
        for (int part = 0; part < 3; part++) response_cache_free(&cache[part]);
    }

    lp_free_samples(samples, n_samples);
    fprintf(report, "\n%s, %d failures\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    printf("%s, %d failures, report in %s\n", failures ? "FAILED" : "passed", failures, report_path);
    return failures ? 1 : 0;
}
