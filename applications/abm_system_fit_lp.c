/*
Fits the linear and the state-dependent local projections, at the settings of
the R code (package lpirfs) they were first run with, to every
replicate under out/dataset/abm_system/, and stores the fits the way the t-QVARMA
fits are stored: once, beside the dataset, for every later comparison to read
rather than refit, the real data's included.

Each replicate's five stored series, the ones the t-QVARMA is fitted on, go
through applications/lp_system.h's lp_system_from_stored, which drops the first
three periods and adds R's 4-period moving average of log GDP as the series
that decides the state, and then through et_al's lp_lin_and_nl at the settings
lp_system_spec() holds: 4 lags, 15 horizons, unit shocks, the state weight from
the HP cycle of that moving average.

A configuration's fits go into one compressed archive,

    out/applications/abm_system_fit_lp/<sample>.npz

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

- ./bin/abm_system_fit_lp fits what is not stored.
- ./bin/abm_system_fit_lp --refit ignores what is stored and fits everything.

The state weight of the state-dependent projections is lpirfs' by default: the
two-sided HP cycle of the switching series at lambda 1600. Two flags change
it, alone or together (applications/lp_system.h, LpSystemStateWeight):

- --causal-hp builds the cycle at t from periods 1..t only, through et_al's
  lp_lin_and_nl_causal.
- --estimate-lambda replaces 1600 with et_al's hp_estimate, the maximum
  likelihood lambda of each replicate's own switching series.

--difference-cointegrated fits both models on the series with inflation and
the interest rate, the pair the t-QVARMA treats as co-integrated, in first
differences (applications/lp_system.h, lp_system_from_stored_differenced); the
switching series is unchanged. Its trees carry _differenced before the
weight's suffix, so out/applications/abm_system_fit_lp_differenced and
out/applications/abm_system_fit_lp_differenced_causal_hp_estimated_lambda.

Each weight writes its own tree, the directory name followed by
_causal_hp, _estimated_lambda or _causal_hp_estimated_lambda, so
out/applications/abm_system_fit_lp_causal_hp for --causal-hp on the full
sample. The linear projection and d do not depend on the weight, so only the
tree of lpirfs' weight stores them; a reader of another weight's tree takes
them from the lpirfs tree of the same sample.

Built with -DESTIMATION_PERIODS=188, each replicate is cut to its first 188
stored periods before lp_system_from_stored, the length of the US sample
applications/us_data.h reads, so the projections see 185 periods as the US fit
does, and everything goes under out/applications/abm_system_fit_lp_t188
instead. The Makefile builds that variant as bin/abm_system_fit_lp_t188.

Adding -DESTIMATION_FIRST_PERIOD=188 moves that window to stored periods 188
to 375 (0-based), the 188 observations right after the first window, and
everything goes under out/applications/abm_system_fit_lp_t188_from188. The
Makefile builds that variant as bin/abm_system_fit_lp_t188_from188.

One task per configuration. openblas_set_num_threads(1) for the reason
applications/abm_system_fit_qvarma.c gives: the parallelism is the outer loop.

out/applications/abm_system_fit_lp_manifest.txt records per replicate whether
its configuration was loaded or fitted, the VAR's statuses, how many horizons of
each model went through the rank-deficient fallback, whether the switching
series was constant, and whether it has responses; under --estimate-lambda
also the estimated lambda and whether its search converged and ended at a
bound of its interval. Each weight's manifest sits beside its tree, named
after it. Nothing printed.
*/

#include "abm_system.h"
#include "lp_system.h"
#include <et_al./stats.h>
#include <cblas.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#define INPUT_DIR_DEFAULT "out/dataset/abm_system"
#define FULL_SAMPLE_OUTPUT_DIR "out/applications/abm_system_fit_lp"
#define QUOTED_VALUE(x) #x
#define QUOTED(x) QUOTED_VALUE(x)
#ifdef ESTIMATION_FIRST_PERIOD
#ifndef ESTIMATION_PERIODS
#error "ESTIMATION_FIRST_PERIOD needs ESTIMATION_PERIODS: it says where the window of that length starts"
#endif
#define OUTPUT_DIR_DEFAULT FULL_SAMPLE_OUTPUT_DIR "_t" QUOTED(ESTIMATION_PERIODS) "_from" QUOTED(ESTIMATION_FIRST_PERIOD)
#elif defined(ESTIMATION_PERIODS)
#define OUTPUT_DIR_DEFAULT FULL_SAMPLE_OUTPUT_DIR "_t" QUOTED(ESTIMATION_PERIODS)
#else
#define OUTPUT_DIR_DEFAULT FULL_SAMPLE_OUTPUT_DIR
#endif
#ifdef ESTIMATION_FIRST_PERIOD
#define ESTIMATION_WINDOW_START ESTIMATION_FIRST_PERIOD
#else
#define ESTIMATION_WINDOW_START 0
#endif

static const char *INPUT_DIR;
static const char *OUTPUT_DIR;
static int DIFFERENCE_COINTEGRATED = 0;

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
#ifdef ESTIMATION_PERIODS
            assert(block[b].c >= ESTIMATION_WINDOW_START + ESTIMATION_PERIODS
                   && "abm_system_fit_lp: replicate shorter than ESTIMATION_FIRST_PERIOD + ESTIMATION_PERIODS");
            Mat head = mat_copy(mat_slice(block[b], 0, block[b].r, ESTIMATION_WINDOW_START,
                                          ESTIMATION_WINDOW_START + ESTIMATION_PERIODS));
            mat_free(block[b]);
            block[b] = head;
#endif
            s.replicate[s.count] = replicate[b];
            lp_system_fitted_series(block[b], DIFFERENCE_COINTEGRATED, &s.y[s.count], &s.switching[s.count]);
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
    double lambda;
    int lambda_is_converged, lambda_at_bound;
} Summary;

static Summary summarise(const LpSystemRecord *r) {
    return (Summary){ r->replicate, r->var_ols_status, r->var_chol_status,
                      fallback_horizons(r->lin_ols_status), fallback_horizons(r->nl_ols_status),
                      r->switching_is_constant, r->has_responses,
                      r->lambda, r->lambda_is_converged, r->lambda_at_bound };
}

int main(int argc, char **argv) {
    /* --difference-cointegrated is taken out here and the other flags left to
       lp_parse_state_weight, which the US fit shares. */
    char **weight_argv = malloc((size_t)(argc + 1) * sizeof(char *));
    assert(weight_argv && "abm_system_fit_lp: out of memory parsing the flags");
    int weight_argc = 0;
    for (int i = 0; i < argc; i++) {
        if (i > 0 && strcmp(argv[i], "--difference-cointegrated") == 0) DIFFERENCE_COINTEGRATED = 1;
        else weight_argv[weight_argc++] = argv[i];
    }
    weight_argv[weight_argc] = NULL;
    int force_refit;
    LpSystemStateWeight weight = lp_parse_state_weight(weight_argc, weight_argv, &force_refit);
    free(weight_argv);

    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = INPUT_DIR_DEFAULT;
    static char weighted_output_dir[512];
    snprintf(weighted_output_dir, sizeof weighted_output_dir, "%s%s%s", OUTPUT_DIR_DEFAULT,
             lp_system_data_suffix(DIFFERENCE_COINTEGRATED), lp_system_state_weight_suffix(weight));
    OUTPUT_DIR = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!OUTPUT_DIR) OUTPUT_DIR = weighted_output_dir;
#ifdef ESTIMATION_PERIODS
    /* A shortened fit does not match the full-sample archives' fingerprints, so
       pointed at them it would refit and rewrite every one. */
    assert(strcmp(OUTPUT_DIR, FULL_SAMPLE_OUTPUT_DIR) != 0
           && "abm_system_fit_lp: a shortened build must not write into the full-sample fits");
#endif
    /* An archive of lpirfs' weight reads as a miss under any other, so pointed
       at those a run would refit and rewrite every one. */
    assert((lp_system_state_weight_is_lpirfs(weight)
            || (strcmp(OUTPUT_DIR, OUTPUT_DIR_DEFAULT) != 0 && strcmp(OUTPUT_DIR, FULL_SAMPLE_OUTPUT_DIR) != 0))
           && "abm_system_fit_lp: another state weight must not write into the fits of lpirfs' weight");
    /* A fit on the differenced series reads as a miss in a tree of the
       undifferenced ones and would rewrite it. */
    if (DIFFERENCE_COINTEGRATED) {
        char undifferenced_output_dir[512];
        snprintf(undifferenced_output_dir, sizeof undifferenced_output_dir, "%s%s", OUTPUT_DIR_DEFAULT,
                 lp_system_state_weight_suffix(weight));
        assert(strcmp(OUTPUT_DIR, undifferenced_output_dir) != 0 && strcmp(OUTPUT_DIR, FULL_SAMPLE_OUTPUT_DIR) != 0
               && "abm_system_fit_lp: a fit on the differenced series must not write into the undifferenced fits");
    }

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
        LpSystemRecord *stored = force_refit ? NULL : lp_system_read_weighted_archive(path, spec, weight, &count);
        if (stored_rows_match(stored, count, &series)) {
            loaded[s] = 1;
        } else {
            free(stored);
            count = series.count;
            stored = malloc((size_t)count * sizeof(LpSystemRecord));
            assert(stored && "abm_system_fit_lp: out of memory for a configuration's fits");
            for (int i = 0; i < count; i++)
                stored[i] = lp_system_fit_record(series.y[i], series.switching[i], spec, weight, series.replicate[i]);
            if (count > 0) lp_system_write_weighted_archive(path, stored, count, spec, weight);
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

    /* The estimated lambdas, as quantiles over every replicate that has one. */
    long n_estimated = 0, n_not_converged = 0, n_at_bound = 0;
    Mat estimated_lambda = mat_new(n_fits > 0 ? (int)n_fits : 1, 1);
    for (int s = 0; s < n_samples && weight.estimate_lambda; s++)
        for (int i = 0; i < n_records[s]; i++) {
            const Summary *r = &records[s][i];
            if (MISNAN(r->lambda)) continue;
            AT(estimated_lambda, (int)n_estimated, 0) = (mreal)r->lambda;
            n_estimated++;
            n_not_converged += !r->lambda_is_converged;
            n_at_bound += r->lambda_at_bound;
        }
    const double lambda_probability[] = { 0, 0.05, 0.25, 0.5, 0.75, 0.95, 1 };
    const int n_lambda_probabilities = (int)(sizeof lambda_probability / sizeof lambda_probability[0]);
    double lambda_quantile[sizeof lambda_probability / sizeof lambda_probability[0]];
    if (n_estimated > 0) {
        Mat finite = mat_slice(estimated_lambda, 0, (int)n_estimated, 0, 1);
        for (int p = 0; p < n_lambda_probabilities; p++)
            lambda_quantile[p] = (double)stats_quantile(finite, (mreal)lambda_probability[p]);
    }
    mat_free(estimated_lambda);

    char manifest_path[640];
    snprintf(manifest_path, sizeof manifest_path, "%s_manifest.txt", OUTPUT_DIR);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "abm_system_fit_lp: cannot open the manifest for writing");
    fprintf(manifest, "Local projections of the R pipeline on every replicate of %s\n", INPUT_DIR);
    fprintf(manifest, "data: %s (applications/lp_system.h); state series R's GDP_MA4\n",
            LP_SYSTEM_SERIES_DESCRIPTION);
    if (DIFFERENCE_COINTEGRATED)
        fprintf(manifest, "Inflation and InterestRate, the series the t-QVARMA treats as co-integrated, in first "
                          "differences; the other three and the state series as stored\n");
#ifdef ESTIMATION_PERIODS
    fprintf(manifest, "each replicate cut to its first %d stored periods, %d after the moving average\n",
            ESTIMATION_PERIODS, lp_system_n_periods(ESTIMATION_PERIODS));
#endif
    char lambda_text[256];
    if (weight.estimate_lambda)
        snprintf(lambda_text, sizeof lambda_text,
                 "lambda estimated per replicate (hp_estimate: log q in [%g, %g], xatol %g, at most %d evaluations)",
                 weight.lambda_options.log_q_lower, weight.lambda_options.log_q_upper,
                 weight.lambda_options.xatol, weight.lambda_options.max_evaluations);
    else
        snprintf(lambda_text, sizeof lambda_text, "lambda %g", spec.lambda);
    fprintf(manifest, "spec: lags_endog_lin %d, lags_endog_nl %d, hor %d, unit shocks, use_hp %d, "
                      "%s, gamma %g, lag_switching %d\n",
            spec.lin.lags_endog_lin, spec.lags_endog_nl, spec.lin.hor, spec.use_hp, lambda_text,
            spec.gamma, spec.lag_switching);
    if (!lp_system_state_weight_is_lpirfs(weight))
        fprintf(manifest, "state weight: %s\n", lp_system_state_weight_description(weight));
    fprintf(manifest, "stored in %s/<sample>.npz, one row per replicate\n", OUTPUT_DIR);
    fprintf(manifest, "%s\n\n", force_refit ? "run with --refit: nothing read from the stored fits"
                                            : "stored fits kept where they match, the rest fitted");
    fprintf(manifest, "configurations %d, replicates %ld\n", n_samples, n_fits);
    fprintf(manifest, "  replicates whose configuration was loaded %ld, fitted in this run %ld\n", n_loaded,
            n_fits - n_loaded);
    fprintf(manifest, "  without responses (no d) %ld\n", n_without);
    /* Only lpirfs' weight stores the linear projection, so only its manifest
       reports on it. */
    int linear = lp_system_state_weight_is_lpirfs(weight);
    if (linear)
        fprintf(manifest, "  with at least one rank-deficient horizon: linear %ld, state dependent %ld\n",
                n_lin_fallback, n_nl_fallback);
    else
        fprintf(manifest, "  with at least one rank-deficient horizon: state dependent %ld (the linear projection "
                          "is stored, and reported, only under lpirfs' weight)\n", n_nl_fallback);
    fprintf(manifest, "  with a constant switching series %ld\n", n_constant);
    if (weight.estimate_lambda) {
        fprintf(manifest, "  with an estimated lambda %ld, of which search not converged %ld, "
                          "ended at a bound of the interval %ld\n",
                n_estimated, n_not_converged, n_at_bound);
        if (n_estimated > 0) {
            fprintf(manifest, "  estimated lambda, quantiles over those replicates:");
            for (int p = 0; p < n_lambda_probabilities; p++)
                fprintf(manifest, " %g: %.6g%s", lambda_probability[p], lambda_quantile[p],
                        p + 1 < n_lambda_probabilities ? "," : "");
            fprintf(manifest, "\n");
        }
    }
    fprintf(manifest, "\n");

    fprintf(manifest, "origin: loaded, the configuration's archive matched the dataset and was kept; "
                      "fitted, computed in this run.\nvar_ols and var_chol are the VAR's statuses behind d "
                      "(0 is fine; var_ols -1 means the data held a NaN\nor an infinity). lin_fallback and "
                      "nl_fallback count the horizons solved by the minimum-norm\nfallback for a "
                      "rank-deficient design, lin_fallback - under a weight that does not store the\n"
                      "linear projection. constant is the switching series' flag. responses says\n"
                      "whether the row has d and its response tensors.\n");
    if (weight.estimate_lambda)
        fprintf(manifest, "lambda is the replicate's estimate, nan where the switching series held a NaN or "
                          "an infinity;\nconverged and at_bound are hp_estimate's flags.\n");
    fprintf(manifest, "\n");
    fprintf(manifest, "%-12s %9s %7s %7s %8s %12s %11s %8s %9s", "sample", "replicate", "origin",
            "var_ols", "var_chol", "lin_fallback", "nl_fallback", "constant", "responses");
    if (weight.estimate_lambda) fprintf(manifest, " %13s %9s %8s", "lambda", "converged", "at_bound");
    fprintf(manifest, "\n");
    for (int s = 0; s < n_samples; s++)
        for (int i = 0; i < n_records[s]; i++) {
            const Summary *r = &records[s][i];
            char lin_fallback[16];
            if (linear) snprintf(lin_fallback, sizeof lin_fallback, "%d", r->lin_fallback_horizons);
            else snprintf(lin_fallback, sizeof lin_fallback, "-");
            fprintf(manifest, "%-12s %9d %7s %7d %8d %12s %11d %8d %9s", samples[s], r->replicate,
                    loaded[s] ? "loaded" : "fitted", r->var_ols_status, r->var_chol_status,
                    lin_fallback, r->nl_fallback_horizons,
                    r->switching_is_constant, r->has_responses ? "yes" : "no");
            if (weight.estimate_lambda)
                fprintf(manifest, " %13.6g %9d %8d", r->lambda, r->lambda_is_converged, r->lambda_at_bound);
            fprintf(manifest, "\n");
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
