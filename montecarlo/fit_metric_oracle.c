/*
The oracle version of the fit metric (montecarlo/fit_metric.h,
docs/MONTECARLO_FIT_METRIC.md) for every benchmark of the 1000 by 1000 Monte
Carlo, under each of the five response vectors.

A benchmark is run s of configuration j, the row montecarlo/sweep_grid.c wrote
for it. d is the mean loss of the configuration with the smallest mean loss,
which that row records and which is in the final confidence set for every one
of the 5,000,000 rows; it is recomputed here from the same response cache, over
the same 999 runs, with the same sums. sigma2 is the variance of each response
entry across configuration j's 999 runs other than s, divisor 998, averaged
over the entries. Then s = d / sqrt(sigma2) and v = 1 / (1 + s).

Reads montecarlo/out/sweep_grid.csv.gz and the response caches in
out/sweep_grid_response_cache/, which bin/sweep_grid builds; a missing cache is
built from out/abm_system_fit_qvarma/, out/abm_system_fit_lp/ and
dataset/abm_system/ the same way. Writes montecarlo/out/fit_metric_oracle.csv.gz,
one row per model and benchmark: model, benchmark configuration and run, the
configuration with the smallest mean loss, d, sigma2, s and v. Progress goes to
stderr. About a minute per model on 16 threads.
*/

#include "montecarlo/fit_metric.h"
#include <et_al./frame/gzip.h>
#include <cblas.h>

#define GRID_PATH "montecarlo/out/sweep_grid.csv.gz"
#define TEXT_PATH "montecarlo/out/fit_metric_oracle.csv"
#define RESULT_PATH "montecarlo/out/fit_metric_oracle.csv.gz"

static unsigned char *read_file(const char *path, size_t *length) {
    FILE *f = fopen(path, "rb");
    assert(f && "fit_metric_oracle: cannot open a file");
    fseek(f, 0, SEEK_END);
    *length = (size_t)ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *bytes = malloc(*length + 1);
    assert(bytes && "fit_metric_oracle: out of memory reading a file");
    size_t n_read = fread(bytes, 1, *length, f);
    assert(n_read == *length && "fit_metric_oracle: short read");
    (void)n_read;
    fclose(f);
    bytes[*length] = 0;
    return bytes;
}

/* best[(model * n_samples + sample) * n_replicates + run], the sample index of
   the configuration with the smallest mean loss in that benchmark's row, or -1
   where the grid has no row. */
static int *read_best(const LpSample *samples, int n_samples, int n_replicates) {
    int highest = 0;
    for (int i = 0; i < n_samples; i++) if (samples[i].index > highest) highest = samples[i].index;
    int *sample_of_cop = malloc((size_t)(highest + 1) * sizeof(int));
    for (int c = 0; c <= highest; c++) sample_of_cop[c] = -1;
    for (int i = 0; i < n_samples; i++) sample_of_cop[samples[i].index] = i;

    size_t n_cells = (size_t)N_MODELS * n_samples * n_replicates;
    int *best = malloc(n_cells * sizeof(int));
    assert(best && sample_of_cop && "fit_metric_oracle: out of memory");
    for (size_t i = 0; i < n_cells; i++) best[i] = -1;

    size_t compressed_length, text_length;
    unsigned char *compressed = read_file(GRID_PATH, &compressed_length);
    char *text = (char *)gzip_inflate(compressed, compressed_length, &text_length);
    free(compressed);

    long n_rows = 0;
    char *line = strchr(text, '\n') + 1;
    while (line < text + text_length && *line) {
        char *end = strchr(line, '\n');
        if (end) *end = 0;
        char label[32];
        int cop, run, lowest;
        int matched = sscanf(line, "%31[^,],%d,%d,%*d,%*d,%*f,%*d,%*d,%*f,%*d,%d", label, &cop, &run, &lowest);
        assert(matched == 4 && "fit_metric_oracle: a malformed row in the grid result");
        (void)matched;
        int model = 0;
        while (model < N_MODELS && strcmp(label, model_label[model]) != 0) model++;
        assert(model < N_MODELS && cop <= highest && lowest <= highest && run >= 0 && run < n_replicates
               && sample_of_cop[cop] >= 0 && sample_of_cop[lowest] >= 0 && "fit_metric_oracle: unknown row");
        best[((size_t)model * n_samples + sample_of_cop[cop]) * n_replicates + run] = sample_of_cop[lowest];
        n_rows++;
        if (!end) break;
        line = end + 1;
    }
    fprintf(stderr, "%ld rows read from %s\n", n_rows, GRID_PATH);
    free(text);
    free(sample_of_cop);
    return best;
}

static void compress_result(void) {
    size_t length, compressed_length;
    unsigned char *text = read_file(TEXT_PATH, &length);
    unsigned char *compressed = gzip_deflate_level(text, length, 9, &compressed_length);
    FILE *out = fopen(RESULT_PATH, "wb");
    assert(out && "fit_metric_oracle: cannot open the result for writing");
    size_t n_written = fwrite(compressed, 1, compressed_length, out);
    assert(n_written == compressed_length && "fit_metric_oracle: short write of the result");
    (void)n_written;
    fclose(out);
    remove(TEXT_PATH);
    fprintf(stderr, "wrote %s, %zu bytes from %zu\n", RESULT_PATH, compressed_length, length);
    free(compressed);
    free(text);
}

int main(void) {
    const char *lp_fit_dir = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!lp_fit_dir) lp_fit_dir = LP_FIT_DIR_DEFAULT;
    const char *input_dir = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!input_dir) input_dir = LP_INPUT_DIR_DEFAULT;
    openblas_set_num_threads(1);

    int n_samples;
    LpSample *samples = lp_list_samples(lp_fit_dir, &n_samples);
    assert(n_samples > 0 && "fit_metric_oracle: no configurations in the LP fit cache");
    int n_replicates = lp_count_replicates(input_dir, samples[0].name);
    int *best = read_best(samples, n_samples, n_replicates);

    FILE *text = fopen(TEXT_PATH, "w");
    assert(text && "fit_metric_oracle: cannot open the output");
    fputs("model,benchmark_cop,benchmark_run,lowest_loss_cop,d,sigma2,s,v\n", text);

    size_t n_benchmarks = (size_t)n_samples * n_replicates;
    double *d = malloc(n_benchmarks * sizeof(double));
    double *sigma2 = malloc(n_benchmarks * sizeof(double));
    int *order = malloc(n_benchmarks * sizeof(int));
    int *first_of = malloc((size_t)(n_samples + 1) * sizeof(int));
    assert(d && sigma2 && order && first_of && "fit_metric_oracle: out of memory");

    for (int model = 0; model < N_MODELS; model++) {
        double started = omp_get_wtime();
        ResponseCache cache = response_cache_load_or_build(model, samples, n_samples, n_replicates, lp_fit_dir,
                                                           input_dir, RESPONSE_CACHE_DIR_DEFAULT);
        assert(cache.n_missing == 0 && "fit_metric_oracle: the leave-one-out moments need every response");
        const int dim = cache.dim;
        const int *model_best = best + (size_t)model * n_benchmarks;

        double *mean = malloc((size_t)n_samples * dim * sizeof(double));
        double *squares = malloc((size_t)n_samples * dim * sizeof(double));
        assert(mean && squares && "fit_metric_oracle: out of memory for the moments");
        #pragma omp parallel for schedule(dynamic)
        for (int sample = 0; sample < n_samples; sample++)
            fit_metric_moments(&cache, sample, mean + (size_t)sample * dim, squares + (size_t)sample * dim);

        /* Benchmarks grouped by the configuration their d reads, so the threads
           work through one configuration's runs together while they are in
           cache. */
        for (int c = 0; c <= n_samples; c++) first_of[c] = 0;
        for (size_t i = 0; i < n_benchmarks; i++) if (model_best[i] >= 0) first_of[model_best[i] + 1]++;
        for (int c = 0; c < n_samples; c++) first_of[c + 1] += first_of[c];
        int n_ordered = first_of[n_samples];
        for (size_t i = 0; i < n_benchmarks; i++) if (model_best[i] >= 0) order[first_of[model_best[i]]++] = (int)i;

        #pragma omp parallel
        {
            double *benchmark = malloc((size_t)dim * sizeof(double));
            assert(benchmark && "fit_metric_oracle: out of memory");
            #pragma omp for schedule(dynamic, 16)
            for (int at = 0; at < n_ordered; at++) {
                int i = order[at], sample = i / n_replicates, run = i % n_replicates;
                const float *own = response_slot(&cache, sample, run);
                for (int k = 0; k < dim; k++) benchmark[k] = (double)own[k];
                d[i] = fit_metric_d(&cache, model_best[i], run, benchmark);
                sigma2[i] = fit_metric_oracle_sigma2(mean + (size_t)sample * dim, squares + (size_t)sample * dim,
                                                     benchmark, n_replicates, dim);
            }
            free(benchmark);
        }

        for (size_t i = 0; i < n_benchmarks; i++) {
            if (model_best[i] < 0) continue;
            int sample = (int)(i / n_replicates), run = (int)(i % n_replicates);
            double s = fit_metric_s(d[i], sigma2[i]);
            fprintf(text, "%s,%d,%d,%d,%.10g,%.10g,%.10g,%.10g\n", model_label[model], samples[sample].index, run,
                    samples[model_best[i]].index, d[i], sigma2[i], s, fit_metric_v(s));
        }
        fprintf(stderr, "%s: %d benchmarks in %.0f s\n", model_label[model], n_ordered, omp_get_wtime() - started);
        free(mean);
        free(squares);
        response_cache_free(&cache);
    }
    fclose(text);
    compress_result();

    free(first_of);
    free(order);
    free(sigma2);
    free(d);
    free(best);
    lp_free_samples(samples, n_samples);
    return 0;
}
