/*
Whether qvarma_growth_states (montecarlo/response_cache.h) computes what it
claims, on real runs: the cached t-QVARMA fit of runs 0 to 9 of cop_0001,
cop_0500 and cop_1000, on that run's stored series.

For each run the split and the two conditional means of D_t are recomputed
here from their definitions, with an independent ordering of the conditioning
series:

1. The recession periods are those t = 5..T whose g_(t-1), the mean of GDP
   growth over periods t-4..t-1, is at or below the 20th percentile of those
   g_(t-1), the percentile taken as (1 - w) x_(j) + w x_(j+1) on the sorted
   values with j + w = 0.2 (n - 1), counting from 0. Their number and set equal
   the header's.
2. D_recession and D_expansion equal the plain means of D_t over each set,
   exactly.
3. Pooled back, (n_recession D_recession + n_expansion D_expansion) / (T - 4)
   equals the mean of D_t over t = 5..T, to 1e-12 of its largest entry.
4. The 1050 entries the cache stores for the run are the 525 total responses
   at the recomputed D_expansion followed by the 525 at D_recession, exactly.
   That check reads the production cache and is reported as not run when the
   cache file is absent.

The report goes to out/tests/qvarma_growth_states_correctness/<YYYYmmdd_HHMMSS>/.
Exits 1 on any failure.
*/

#include "montecarlo/response_cache.h"
#include <stdio.h>
#include <time.h>

#define OUTPUT_ROOT "out/tests/qvarma_growth_states_correctness"
#define N_RUNS 10

static const int checked_configuration[] = { 1, 500, 1000 };
static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    fflush(report);
    if (!ok) failures++;
}

static int compare_doubles(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The mapped production cache of qvarma_growth_states, or NULL when it is
   absent or its header disagrees. */
static float *map_cache(int n_samples, int n_replicates, void **mapping, size_t *length) {
    char path[640];
    snprintf(path, sizeof path, "%s/%s.f32", RESPONSE_CACHE_DIR_DEFAULT, model_label[MODEL_QVARMA_GROWTH_STATES]);
    int expected[4] = { MODEL_QVARMA_GROWTH_STATES, n_samples, n_replicates, 2 * RESPONSE_QVARMA_DIM };
    *length = sizeof expected + (size_t)n_samples * n_replicates * expected[3] * sizeof(float);
    struct stat st;
    if (stat(path, &st) != 0 || (size_t)st.st_size != *length) return NULL;
    *mapping = _response_cache_map(path, *length, 0);
    if (!*mapping) return NULL;
    if (memcmp(*mapping, expected, sizeof expected) != 0) {
        munmap(*mapping, *length);
        return NULL;
    }
    return (float *)((char *)*mapping + sizeof expected);
}

int main(void) {
    char stamp[32], dir[256], report_path[400];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", localtime(&now));
    mkdir("out", 0755);
    mkdir("out/tests", 0755);
    mkdir(OUTPUT_ROOT, 0755);
    snprintf(dir, sizeof dir, "%s/%s", OUTPUT_ROOT, stamp);
    int created = mkdir(dir, 0755) == 0;
    assert(created && "qvarma_growth_states_correctness: the output directory exists already");
    (void)created;
    snprintf(report_path, sizeof report_path, "%s/qvarma_growth_states_correctness_report.txt", dir);
    FILE *report = fopen(report_path, "w");
    assert(report);
    openblas_set_num_threads(1);

    int n_samples = 0;
    LpSample *samples = lp_list_samples(LP_FIT_DIR_DEFAULT, &n_samples);
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);
    void *mapping = NULL;
    size_t length = 0;
    float *cache = map_cache(n_samples, n_replicates, &mapping, &length);

    QvarmaParams m = _response_qvarma_shape();
    const int K = m.K, first = RESPONSE_QVARMA_STATE_WINDOW;
    int split_ok = 1, means_ok = 1, pooled_ok = 1, cache_ok = 1, n_checked = 0, n_recession_seen = -1;
    double largest_pooled_gap = 0;
    float *expected = malloc(2 * RESPONSE_QVARMA_DIM * sizeof(float));
    assert(expected);
    int n_configurations = (int)(sizeof checked_configuration / sizeof checked_configuration[0]);
    for (int c = 0; c < n_configurations; c++) {
        int sample = -1;
        for (int s = 0; s < n_samples; s++)
            if (samples[s].index == checked_configuration[c]) sample = s;
        assert(sample >= 0 && "qvarma_growth_states_correctness: a checked configuration is missing");
        char data_dir[560];
        snprintf(data_dir, sizeof data_dir, "%s/%s", LP_INPUT_DIR_DEFAULT, samples[sample].name);
        for (int run = 0; run < N_RUNS; run++) {
            char fit_path[640];
            snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json", RESPONSE_QVARMA_FIT_DIR,
                     samples[sample].name, run, RESPONSE_QVARMA_SPEC_LABEL);
            int loaded = qvarma_load_params(&m, fit_path);
            assert(loaded && m.nu > 2 && "qvarma_growth_states_correctness: a checked run has no usable fit");
            (void)loaded;
            Mat y = abm_system_read_replicate(data_dir, run);
            int T = y.c, n = T - first;

            /* The conditioning series and its percentile, from the definition. */
            double *lagged = malloc((size_t)n * sizeof(double)), *sorted = malloc((size_t)n * sizeof(double));
            assert(lagged && sorted);
            for (int t = first; t < T; t++) {
                double g = 0;
                for (int s = t - first; s < t; s++) g += (double)AT(y, STORED_GDP_GROWTH, s);
                lagged[t - first] = sorted[t - first] = g / first;
            }
            qsort(sorted, (size_t)n, sizeof(double), compare_doubles);
            double position = (n - 1) * RESPONSE_QVARMA_RECESSION_QUANTILE;
            int j = (int)position;
            double w = position - j;
            double threshold = j + 1 < n ? (1 - w) * sorted[j] + w * sorted[j + 1] : sorted[j];

            Mat *D_path = qvarma_score_jacobian(&m, y);
            Mat sum[2] = { mat_new(K, K), mat_new(K, K) }, pooled = mat_new(K, K);
            int count[2] = { 0, 0 };
            for (int t = first; t < T; t++) {
                int recession = lagged[t - first] <= threshold;
                for (int i = 0; i < K * K; i++) {
                    sum[recession].d[i] += D_path[t].d[i];
                    pooled.d[i] += D_path[t].d[i];
                }
                count[recession]++;
            }
            for (int s = 0; s < 2; s++)
                for (int i = 0; i < K * K; i++) sum[s].d[i] /= (mreal)count[s];

            Mat D_expansion, D_recession;
            int n_recession = 0;
            int has_states = response_qvarma_state_jacobians(&m, y, &D_expansion, &D_recession, &n_recession);
            split_ok &= has_states && n_recession == count[1];
            if (n_recession_seen < 0) n_recession_seen = n_recession;
            if (has_states) {
                double largest = 0, gap = 0;
                for (int i = 0; i < K * K; i++)
                    if (fabs((double)pooled.d[i] / n) > largest) largest = fabs((double)pooled.d[i] / n);
                for (int i = 0; i < K * K; i++) {
                    means_ok &= D_expansion.d[i] == sum[0].d[i] && D_recession.d[i] == sum[1].d[i];
                    double back = ((double)count[1] * D_recession.d[i] + (double)count[0] * D_expansion.d[i]) / n;
                    double entry_gap = fabs(back - (double)pooled.d[i] / n) / largest;
                    if (entry_gap > gap) gap = entry_gap;
                }
                if (gap > largest_pooled_gap) largest_pooled_gap = gap;
                pooled_ok &= gap <= 1e-12;
                mat_free(D_expansion);
                mat_free(D_recession);
            }

            if (cache) {
                int ok = _response_qvarma_total(&m, sum[0], expected)
                      && _response_qvarma_total(&m, sum[1], expected + RESPONSE_QVARMA_DIM);
                const float *slot = cache + ((size_t)sample * n_replicates + run) * 2 * RESPONSE_QVARMA_DIM;
                cache_ok &= ok && memcmp(slot, expected, 2 * RESPONSE_QVARMA_DIM * sizeof(float)) == 0;
            }

            qvarma_free_mat_array(D_path, T);
            mat_free(sum[0]);
            mat_free(sum[1]);
            mat_free(pooled);
            free(lagged);
            free(sorted);
            mat_free(y);
            n_checked++;
        }
    }

    fprintf(report, "%d runs: runs 0 to %d of cop_0001, cop_0500 and cop_1000; recession periods in the first %d\n\n",
            n_checked, N_RUNS - 1, n_recession_seen);
    check(report, split_ok, "the recession periods and their number equal the definition's");
    check(report, means_ok, "D_recession and D_expansion equal the plain means over their periods, exactly");
    fprintf(report, "      largest gap of the pooled means, relative to the largest entry of the mean %.3g\n", largest_pooled_gap);
    check(report, pooled_ok, "the two means pooled back give the mean over periods 5..T, to 1e-12 relative");
    if (cache) check(report, cache_ok, "the cache holds the responses at D_expansion then D_recession, exactly");
    else fprintf(report, "not run  the cache comparison, no qvarma_growth_states cache on disk\n");
    fprintf(report, "\n%s\n", failures ? "FAILED" : "all passed");
    fclose(report);

    if (cache) munmap(mapping, length);
    free(expected);
    qvarma_params_free(&m);
    lp_free_samples(samples, n_samples);
    return failures ? 1 : 0;
}
