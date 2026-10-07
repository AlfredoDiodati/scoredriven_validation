/*
The feasible version of the fit metric (montecarlo/fit_metric.h,
docs/MONTECARLO_FIT_METRIC.md) for one run of every configuration of the 1000
by 1000 Monte Carlo, and the comparison with the oracle version on the same
benchmarks.

A benchmark is run R of configuration j (FIT_METRIC_RUN, default 0). Its own
series is resampled B times by the moving-block bootstrap, every method is
re-estimated on every resample, and sigma2_star is the variance of each
response entry across the resamples, divisor B - 1, averaged over the entries.
d and the oracle sigma2 are read from montecarlo/out/fit_metric_oracle.csv.gz,
so bin/fit_metric_oracle runs first; s_star = d / sqrt(sigma2_star) and
v_star = 1 / (1 + s_star) are formed here.

Settings read from the environment, defaults from
fit_metric_bootstrap_options_default:
  FIT_METRIC_RUN          which run of each configuration is the benchmark, 0
  FIT_METRIC_B            resamples per benchmark, 200
  FIT_METRIC_BLOCK        periods per block, 8
  FIT_METRIC_ITERATIONS   the t-QVARMA solver's budget per resample, 4000
  FIT_METRIC_CONFIGURATIONS  only the first this many configurations, all
The resamples of a benchmark come from et_al's generator at seed 123 and a
stream fixed by the configuration and the run, so they do not depend on which
thread draws them or in what order the benchmarks finish.

One benchmark per thread, every thread busy. Resumable at benchmark
granularity: a benchmark's five rows go to
montecarlo/out/fit_metric_bootstrap_progress.csv together when it finishes,
and a rerun skips the configurations that file holds. When every configuration
is done the rows are sorted into montecarlo/out/fit_metric_bootstrap.csv and
the progress file is removed. Under any setting other than the defaults both
names carry the settings, as in
fit_metric_bootstrap_run0_B50_block25_iterations4000_first100.csv, so a run
that varies them leaves the default result alone. Columns: model, benchmark configuration and run,
B, block length, usable resamples, converged t-QVARMA fits, their mean
iterations, the configuration with the smallest mean loss, d, sigma2,
sigma2_star, s, v, s_star and v_star.

Reads dataset/abm_system/ and out/abm_system_fit_qvarma/ (the benchmark's own
t-QVARMA estimate, where every resample's fit starts) and rebuilds neither.
Progress goes to stderr.
*/

#include "montecarlo/fit_metric.h"
#include <et_al./frame/gzip.h>
#include <cblas.h>

#ifndef ORACLE_PATH
#define ORACLE_PATH "montecarlo/out/fit_metric_oracle.csv.gz"
#endif
#define OUT_STEM "montecarlo/out/fit_metric_bootstrap"
#define SEED 123
#define HEADER "model,benchmark_cop,benchmark_run,n_bootstrap,block_length,n_usable,n_qvarma_converged," \
               "qvarma_mean_iterations,lowest_loss_cop,d,sigma2,sigma2_star,s,v,s_star,v_star\n"

/* Where the rows go: fit_metric_bootstrap.csv under the default settings, and
   a name carrying the settings otherwise, so a run that varies them leaves the
   default result alone. */
static char PROGRESS_PATH[256], RESULT_PATH[256];

static void set_paths(int run, FitMetricBootstrapOptions options, int n_configurations, int all_configurations) {
    FitMetricBootstrapOptions standard = fit_metric_bootstrap_options_default();
    char suffix[160] = "";
    if (run != 0 || options.n_bootstrap != standard.n_bootstrap || options.block_length != standard.block_length
        || options.qvarma_max_iterations != standard.qvarma_max_iterations || !all_configurations)
        snprintf(suffix, sizeof suffix, "_run%d_B%d_block%d_iterations%d_first%d", run, options.n_bootstrap,
                 options.block_length, options.qvarma_max_iterations, n_configurations);
    snprintf(PROGRESS_PATH, sizeof PROGRESS_PATH, "%s%s_progress.csv", OUT_STEM, suffix);
    snprintf(RESULT_PATH, sizeof RESULT_PATH, "%s%s.csv", OUT_STEM, suffix);
}

static int env_int(const char *name, int fallback) {
    const char *text = getenv(name);
    return text ? atoi(text) : fallback;
}

/* What the oracle file holds for the benchmarks of one run: d, sigma2 and the
   configuration with the smallest mean loss, by model and sample. */
typedef struct { double d, sigma2; int lowest_cop; } OracleRow;

static OracleRow *read_oracle(const LpSample *samples, int n_samples, int run) {
    int highest = 0;
    for (int i = 0; i < n_samples; i++) if (samples[i].index > highest) highest = samples[i].index;
    int *sample_of_cop = malloc((size_t)(highest + 1) * sizeof(int));
    for (int c = 0; c <= highest; c++) sample_of_cop[c] = -1;
    for (int i = 0; i < n_samples; i++) sample_of_cop[samples[i].index] = i;
    OracleRow *rows = malloc((size_t)N_MODELS * n_samples * sizeof(OracleRow));
    assert(rows && sample_of_cop && "fit_metric_bootstrap: out of memory");
    for (int i = 0; i < N_MODELS * n_samples; i++) rows[i] = (OracleRow){ NAN, NAN, -1 };

    FILE *f = fopen(ORACLE_PATH, "rb");
    assert(f && "fit_metric_bootstrap: montecarlo/out/fit_metric_oracle.csv.gz is missing, run bin/fit_metric_oracle");
    fseek(f, 0, SEEK_END);
    size_t length = (size_t)ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc(length);
    size_t n_read = fread(compressed, 1, length, f);
    assert(n_read == length && "fit_metric_bootstrap: short read of the oracle file");
    (void)n_read;
    fclose(f);
    size_t text_length;
    char *text = (char *)gzip_inflate(compressed, length, &text_length);
    free(compressed);

    char *line = strchr(text, '\n') + 1;
    while (line < text + text_length && *line) {
        char *end = strchr(line, '\n');
        if (end) *end = 0;
        char label[32];
        int cop, row_run, lowest;
        double d, sigma2;
        int matched = sscanf(line, "%31[^,],%d,%d,%d,%lf,%lf", label, &cop, &row_run, &lowest, &d, &sigma2);
        assert(matched == 6 && "fit_metric_bootstrap: a malformed row in the oracle file");
        (void)matched;
        if (row_run == run && cop <= highest && sample_of_cop[cop] >= 0)
            for (int model = 0; model < N_MODELS; model++)
                if (strcmp(label, model_label[model]) == 0)
                    rows[model * n_samples + sample_of_cop[cop]] = (OracleRow){ d, sigma2, lowest };
        if (!end) break;
        line = end + 1;
    }
    free(text);
    free(sample_of_cop);
    return rows;
}

/* Which configurations the progress file already holds, all five rows of
   them, written under the same settings. */
static unsigned char *read_done(const LpSample *samples, int n_samples, int run,
                                FitMetricBootstrapOptions options, int *n_done) {
    unsigned char *rows_seen = calloc((size_t)n_samples, 1);
    assert(rows_seen && "fit_metric_bootstrap: out of memory");
    *n_done = 0;
    FILE *f = fopen(PROGRESS_PATH, "r");
    if (!f) return rows_seen;
    char line[1024];
    int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        char label[32];
        int cop, row_run, n_bootstrap, block_length;
        if (sscanf(line, "%31[^,],%d,%d,%d,%d,", label, &cop, &row_run, &n_bootstrap, &block_length) != 5) continue;
        assert(row_run == run && n_bootstrap == options.n_bootstrap && block_length == options.block_length
               && "fit_metric_bootstrap: the progress file was written under other settings; move it away first");
        for (int i = 0; i < n_samples; i++)
            if (samples[i].index == cop) rows_seen[i]++;
    }
    fclose(f);
    for (int i = 0; i < n_samples; i++) {
        assert(rows_seen[i] == 0 || rows_seen[i] == N_MODELS);
        rows_seen[i] = rows_seen[i] == N_MODELS;
        *n_done += rows_seen[i];
    }
    return rows_seen;
}

static int compare_lines(const void *a, const void *b) {
    const char *x = *(char *const *)a, *y = *(char *const *)b;
    char label_x[32], label_y[32];
    int cop_x, cop_y, model_x = 0, model_y = 0;
    sscanf(x, "%31[^,],%d", label_x, &cop_x);
    sscanf(y, "%31[^,],%d", label_y, &cop_y);
    while (model_x < N_MODELS && strcmp(label_x, model_label[model_x]) != 0) model_x++;
    while (model_y < N_MODELS && strcmp(label_y, model_label[model_y]) != 0) model_y++;
    if (model_x != model_y) return model_x - model_y;
    return cop_x - cop_y;
}

/* The progress rows sorted by model, then configuration, into the result. */
static void write_result(void) {
    FILE *f = fopen(PROGRESS_PATH, "r");
    assert(f && "fit_metric_bootstrap: the progress file is missing");
    char **lines = NULL, line[1024];
    int n = 0, cap = 0, first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        if (n == cap) {
            cap = cap ? 2 * cap : 1024;
            lines = realloc(lines, (size_t)cap * sizeof(char *));
            assert(lines && "fit_metric_bootstrap: out of memory");
        }
        lines[n] = malloc(strlen(line) + 1);
        strcpy(lines[n++], line);
    }
    fclose(f);
    qsort(lines, (size_t)n, sizeof(char *), compare_lines);
    FILE *out = fopen(RESULT_PATH, "w");
    assert(out && "fit_metric_bootstrap: cannot open the result");
    fputs(HEADER, out);
    for (int i = 0; i < n; i++) {
        fputs(lines[i], out);
        free(lines[i]);
    }
    fclose(out);
    free(lines);
    remove(PROGRESS_PATH);
    fprintf(stderr, "wrote %s, %d rows\n", RESULT_PATH, n);
}

int main(void) {
    const char *input_dir = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!input_dir) input_dir = LP_INPUT_DIR_DEFAULT;
    const char *lp_fit_dir = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!lp_fit_dir) lp_fit_dir = LP_FIT_DIR_DEFAULT;
    FitMetricBootstrapOptions options = fit_metric_bootstrap_options_default();
    options.n_bootstrap = env_int("FIT_METRIC_B", options.n_bootstrap);
    options.block_length = env_int("FIT_METRIC_BLOCK", options.block_length);
    options.qvarma_max_iterations = env_int("FIT_METRIC_ITERATIONS", options.qvarma_max_iterations);
    const int run = env_int("FIT_METRIC_RUN", 0);
    openblas_set_num_threads(1);

    int n_samples;
    LpSample *samples = lp_list_samples(lp_fit_dir, &n_samples);
    assert(n_samples > 0 && "fit_metric_bootstrap: no configurations in the LP fit cache");
    const int n_replicates = lp_count_replicates(input_dir, samples[0].name);
    assert(run >= 0 && run < n_replicates && "fit_metric_bootstrap: FIT_METRIC_RUN is not a run of the design");
    int n_configurations = env_int("FIT_METRIC_CONFIGURATIONS", n_samples);
    assert(n_configurations >= 1 && "fit_metric_bootstrap: FIT_METRIC_CONFIGURATIONS must be at least 1");
    int all_configurations = n_configurations >= n_samples;
    if (!all_configurations) n_samples = n_configurations;
    set_paths(run, options, n_samples, all_configurations);

    struct stat st;
    if (stat(RESULT_PATH, &st) == 0 && stat(PROGRESS_PATH, &st) != 0) {
        fprintf(stderr, "%s already holds every configuration, nothing to do\n", RESULT_PATH);
        lp_free_samples(samples, n_samples);
        return 0;
    }

    OracleRow *oracle = read_oracle(samples, n_samples, run);
    int n_done;
    unsigned char *done = read_done(samples, n_samples, run, options, &n_done);
    int n_to_do = n_samples - n_done;
    fprintf(stderr, "run %d, B %d, blocks of %d, %d iterations: %d of %d configurations already done\n", run,
            options.n_bootstrap, options.block_length, options.qvarma_max_iterations, n_done, n_samples);

    int fresh = stat(PROGRESS_PATH, &st) != 0;
    FILE *progress = fopen(PROGRESS_PATH, "a");
    assert(progress && "fit_metric_bootstrap: cannot open the progress file");
    if (fresh) fputs(HEADER, progress);
    fflush(progress);

    double started = omp_get_wtime();
    int finished = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int sample = 0; sample < n_samples; sample++) {
        if (done[sample]) continue;
        char dir[600], fit_path[700];
        snprintf(dir, sizeof dir, "%s/%s", input_dir, samples[sample].name);
        snprintf(fit_path, sizeof fit_path, "%s/%s/replicate_%03d_%s_fit.json", RESPONSE_QVARMA_FIT_DIR,
                 samples[sample].name, run, RESPONSE_QVARMA_SPEC_LABEL);
        Mat stored = abm_system_read_replicate(dir, run);
        QvarmaParams start = _response_qvarma_shape();
        int loaded = qvarma_load_params(&start, fit_path);
        assert(loaded && "fit_metric_bootstrap: the benchmark's t-QVARMA fit does not load");
        (void)loaded;

        Rng rng = rng_new(SEED, 1 + (uint64_t)samples[sample].index * n_replicates + run);
        FitMetricBootstrap result = fit_metric_bootstrap(stored, &start, options, &rng);

        #pragma omp critical(fit_metric_bootstrap_progress)
        {
            for (int model = 0; model < N_MODELS; model++) {
                const OracleRow *o = &oracle[model * n_samples + sample];
                double s = fit_metric_s(o->d, o->sigma2), s_star = fit_metric_s(o->d, result.sigma2_star[model]);
                fprintf(progress, "%s,%d,%d,%d,%d,%d,%d,%.1f,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                        model_label[model], samples[sample].index, run, options.n_bootstrap, options.block_length,
                        result.n_usable[model], result.n_qvarma_converged, result.qvarma_mean_iterations,
                        o->lowest_cop, o->d, o->sigma2, result.sigma2_star[model], s, fit_metric_v(s), s_star,
                        fit_metric_v(s_star));
            }
            fflush(progress);
            finished++;
            double elapsed = omp_get_wtime() - started;
            fprintf(stderr, "%s done, %d this session, %.1f s each, %.1f h left\n", samples[sample].name, finished,
                    elapsed / finished, elapsed / finished * (n_to_do - finished) / 3600);
        }
        qvarma_params_free(&start);
        mat_free(stored);
    }
    fclose(progress);
    write_result();

    free(done);
    free(oracle);
    lp_free_samples(samples, n_samples);
    return 0;
}
