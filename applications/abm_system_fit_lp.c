/*
Fits the linear and the state-dependent local projections of the collaborator's
R pipeline (_temp/main_code.R, estimate() in _temp/functions.R) to every
replicate under dataset/abm_system/, and stores the fits the way the t-QVARMA
fits are stored: once, beside the dataset, for every later comparison to read
rather than refit, the real data's included.

Each replicate's five stored series, the ones the t-QVARMA is fitted on, go
through applications/lp_system.h's lp_system_from_stored, which drops the first
three periods and adds R's 4-period moving average of log GDP as the series
that decides the state, and then through et_al's lp_lin_and_nl at the settings
lp_system_spec() holds: 4 lags, 15 horizons, unit shocks, the state weight from
the HP cycle of that moving average.

A configuration's fits go into one compressed archive,

    out/abm_system_fit_lp/<sample>.npz

one row per replicate, in the layout lp_system.h describes: the three response
tensors, d, the fit notes, the specification and the fingerprints of the data
each row was fitted on. et_al's own LP cache writes one plain JSON file per fit,
two million files and about 33 GB of text for this design; the archive holds the
same numbers in binary, without the state weights fz, which the data determine.

The estimator is closed form, so a stored fit is final. A rerun reads a
configuration's archive and keeps it when it holds a row for every replicate the
dataset has, under this specification, fitted on exactly that replicate's
series; otherwise the whole configuration is fitted again and its archive
rewritten. A fit without d (the VAR behind the shocks had no Cholesky factor,
or the data held a NaN) is stored as a row with has_responses 0 and counted, so
it is not retried on every run only to fail the same way.

    ./bin/abm_system_fit_lp             fit what is not stored
    ./bin/abm_system_fit_lp --refit     ignore what is stored and fit everything
    ./bin/abm_system_fit_lp --r-levels  the same on R's levels, the collaborator's
                                        transformation, into
                                        out/abm_system_fit_lp_r_levels/

One task per configuration. openblas_set_num_threads(1) for the reason
applications/abm_system_fit_qvarma.c gives: the parallelism is the outer loop.

out/abm_system_fit_lp_manifest.txt records per replicate whether its
configuration was loaded or fitted, the VAR's statuses, how many horizons of
each model went through the rank-deficient fallback, whether the switching
series was constant, and whether it has responses. Nothing printed.
*/

#include "abm_system.h"
#include "lp_system.h"
#include <cblas.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#define INPUT_DIR_DEFAULT "dataset/abm_system"
#define OUTPUT_DIR_DEFAULT "out/abm_system_fit_lp"

static const char *INPUT_DIR;
static const char *OUTPUT_DIR;
static LpLayout LAYOUT;

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_fit_lp: mkdir failed");
}

static int compare_names(const void *a, const void *b) {
    return strcmp(*(char *const *)a, *(char *const *)b);
}

/* The sample directories, sorted so the manifest reads cop_0001 first. */
static char **list_samples(int *count) {
    DIR *handle = opendir(INPUT_DIR);
    assert(handle && "abm_system_fit_lp: cannot open the dataset directory");

    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char path[560];
        snprintf(path, sizeof path, "%s/%s", INPUT_DIR, entry->d_name);
        struct stat st;
        if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            char **grown = realloc(names, (size_t)cap * sizeof(char *));
            assert(grown && "abm_system_fit_lp: out of memory listing samples");
            names = grown;
        }
        size_t length = strlen(entry->d_name);
        names[n] = malloc(length + 1);
        assert(names[n] && "abm_system_fit_lp: out of memory copying a sample name");
        memcpy(names[n], entry->d_name, length + 1);
        n++;
    }
    closedir(handle);
    qsort(names, (size_t)n, sizeof(char *), compare_names);
    *count = n;
    return names;
}

/* One past the highest batch index a sample holds, read off the file names. A
   batch that is absent is skipped by abm_system_read_batch. */
static int count_batches(const char *dir) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_fit_lp: cannot open a sample directory");
    int highest = -1;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int batch;
        if (sscanf(entry->d_name, "batch_%d.npz", &batch) == 1 && batch > highest) highest = batch;
    }
    closedir(handle);
    return highest + 1;
}

/* One configuration's replicates, transformed, in the order the archives hold
   them. */
typedef struct {
    int count;
    int *replicate;
    Mat *y, *switching;
} Series;

static Series read_series(const char *sample) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    int n_batches = count_batches(dir);
    int cap = n_batches * ABM_SYSTEM_BATCH;
    Series s = { 0, malloc((size_t)cap * sizeof(int)), malloc((size_t)cap * sizeof(Mat)),
                 malloc((size_t)cap * sizeof(Mat)) };
    assert(s.replicate && s.y && s.switching && "abm_system_fit_lp: out of memory reading a sample");
    for (int batch = 0; batch < n_batches; batch++) {
        Mat block[ABM_SYSTEM_BATCH];
        int replicate[ABM_SYSTEM_BATCH];
        int count = abm_system_read_batch(dir, batch, block, replicate);
        for (int b = 0; b < count; b++) {
            s.replicate[s.count] = replicate[b];
            lp_system_series(block[b], LAYOUT, &s.y[s.count], &s.switching[s.count]);
            s.count++;
            mat_free(block[b]);
        }
    }
    return s;
}

static void free_series(Series *s) {
    for (int i = 0; i < s->count; i++) {
        mat_free(s->y[i]);
        mat_free(s->switching[i]);
    }
    free(s->replicate);
    free(s->y);
    free(s->switching);
}

/* Whether the stored rows are exactly this configuration's fits: one row per
   replicate the dataset holds, in the same order, each fitted on that
   replicate's own series. */
static int stored_rows_match(const LpSystemRecord *records, int n_records, const Series *s) {
    if (!records || n_records != s->count) return 0;
    for (int i = 0; i < s->count; i++)
        if (records[i].replicate != s->replicate[i]
            || !lp_system_record_matches(&records[i], s->y[i], s->switching[i])) return 0;
    return 1;
}

static int fallback_horizons(const int *ols_status) {
    int n = 0;
    for (int h = 0; h < LP_SYSTEM_HOR; h++) n += ols_status[h] != 0;
    return n;
}

/* What the manifest reports about one replicate. A configuration's full rows
   are about 10 MB and are released as soon as it is done, since holding the
   whole design's would take 10 GB. */
typedef struct {
    int replicate;
    int var_ols_status, var_chol_status;
    int lin_fallback_horizons, nl_fallback_horizons;
    int switching_is_constant, has_responses;
} Summary;

static Summary summarise(const LpSystemRecord *r) {
    return (Summary){ r->replicate, r->var_ols_status, r->var_chol_status,
                      fallback_horizons(r->lin_ols_status), fallback_horizons(r->nl_ols_status),
                      r->switching_is_constant, r->has_responses };
}

int main(int argc, char **argv) {
    int force_refit;
    LAYOUT = lp_parse_flags(argc, argv, &force_refit);

    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = INPUT_DIR_DEFAULT;
    char default_output[256];
    snprintf(default_output, sizeof default_output, "%s%s", OUTPUT_DIR_DEFAULT, lp_layout_suffix(LAYOUT));
    OUTPUT_DIR = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!OUTPUT_DIR) OUTPUT_DIR = default_output;

    openblas_set_num_threads(1);
    make_directory(OUTPUT_DIR);

    int n_samples;
    char **samples = list_samples(&n_samples);
    assert(n_samples > 0 && "abm_system_fit_lp: no sample directories in the dataset");

    LpNlSpec spec = lp_system_spec();
    Summary **records = calloc((size_t)n_samples, sizeof(Summary *));
    int *n_records = calloc((size_t)n_samples, sizeof(int));
    int *loaded = calloc((size_t)n_samples, sizeof(int));
    assert(records && n_records && loaded && "abm_system_fit_lp: out of memory for the results");

    #pragma omp parallel for schedule(dynamic)
    for (int s = 0; s < n_samples; s++) {
        char path[640];
        snprintf(path, sizeof path, "%s/%s.npz", OUTPUT_DIR, samples[s]);
        Series series = read_series(samples[s]);

        int count = 0;
        LpSystemRecord *stored = force_refit ? NULL : lp_system_read_archive(path, spec, &count);
        if (stored_rows_match(stored, count, &series)) {
            loaded[s] = 1;
        } else {
            free(stored);
            count = series.count;
            stored = malloc((size_t)count * sizeof(LpSystemRecord));
            assert(stored && "abm_system_fit_lp: out of memory for a configuration's fits");
            for (int i = 0; i < count; i++) {
                LpLinNlFit fit = lp_lin_and_nl(series.y[i], series.switching[i], spec, (LpBands){0});
                stored[i] = lp_system_record(&fit, series.y[i], series.switching[i], series.replicate[i]);
                lp_lin_nl_fit_free(&fit);
            }
            if (count > 0) lp_system_write_archive(path, stored, count, spec);
        }
        records[s] = malloc((size_t)(count > 0 ? count : 1) * sizeof(Summary));
        assert(records[s] && "abm_system_fit_lp: out of memory for the manifest");
        for (int i = 0; i < count; i++) records[s][i] = summarise(&stored[i]);
        n_records[s] = count;
        free(stored);
        free_series(&series);
    }

    long n_fits = 0, n_loaded = 0, n_without = 0, n_lin_fallback = 0, n_nl_fallback = 0, n_constant = 0;
    for (int s = 0; s < n_samples; s++)
        for (int i = 0; i < n_records[s]; i++) {
            const Summary *r = &records[s][i];
            n_fits++;
            n_loaded += loaded[s];
            n_without += !r->has_responses;
            n_lin_fallback += r->lin_fallback_horizons > 0;
            n_nl_fallback += r->nl_fallback_horizons > 0;
            n_constant += r->switching_is_constant;
        }

    char manifest_path[640];
    snprintf(manifest_path, sizeof manifest_path, "%s_manifest.txt", OUTPUT_DIR);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "abm_system_fit_lp: cannot open the manifest for writing");
    fprintf(manifest, "Local projections of the R pipeline on every replicate of %s\n", INPUT_DIR);
    fprintf(manifest, "data: %s (applications/lp_system.h); state series R's GDP_MA4\n",
            lp_layout_description(LAYOUT));
    fprintf(manifest, "spec: lags_endog_lin %d, lags_endog_nl %d, hor %d, unit shocks, use_hp %d, "
                      "lambda %g, gamma %g, lag_switching %d\n",
            spec.lin.lags_endog_lin, spec.lags_endog_nl, spec.lin.hor, spec.use_hp, spec.lambda,
            spec.gamma, spec.lag_switching);
    fprintf(manifest, "stored in %s/<sample>.npz, one row per replicate\n", OUTPUT_DIR);
    fprintf(manifest, "%s\n\n", force_refit ? "run with --refit: nothing read from the stored fits"
                                            : "stored fits kept where they match, the rest fitted");
    fprintf(manifest, "configurations %d, replicates %ld\n", n_samples, n_fits);
    fprintf(manifest, "  replicates whose configuration was loaded %ld, fitted in this run %ld\n", n_loaded,
            n_fits - n_loaded);
    fprintf(manifest, "  without responses (no d) %ld\n", n_without);
    fprintf(manifest, "  with at least one rank-deficient horizon: linear %ld, state dependent %ld\n",
            n_lin_fallback, n_nl_fallback);
    fprintf(manifest, "  with a constant switching series %ld\n\n", n_constant);

    fprintf(manifest, "origin: loaded, the configuration's archive matched the dataset and was kept; "
                      "fitted, computed in this run.\nvar_ols and var_chol are the VAR's statuses behind d "
                      "(0 is fine; var_ols -1 means the data held a NaN\nor an infinity). lin_fallback and "
                      "nl_fallback count the horizons solved by the minimum-norm\nfallback for a "
                      "rank-deficient design. constant is the switching series' flag. responses says\n"
                      "whether the row has d and the three response tensors.\n\n");
    fprintf(manifest, "%-12s %9s %7s %7s %8s %12s %11s %8s %9s\n", "sample", "replicate", "origin",
            "var_ols", "var_chol", "lin_fallback", "nl_fallback", "constant", "responses");
    for (int s = 0; s < n_samples; s++)
        for (int i = 0; i < n_records[s]; i++) {
            const Summary *r = &records[s][i];
            fprintf(manifest, "%-12s %9d %7s %7d %8d %12d %11d %8d %9s\n", samples[s], r->replicate,
                    loaded[s] ? "loaded" : "fitted", r->var_ols_status, r->var_chol_status,
                    r->lin_fallback_horizons, r->nl_fallback_horizons,
                    r->switching_is_constant, r->has_responses ? "yes" : "no");
        }
    fclose(manifest);

    for (int s = 0; s < n_samples; s++) {
        free(records[s]);
        free(samples[s]);
    }
    free(records);
    free(n_records);
    free(loaded);
    free(samples);
    return 0;
}
