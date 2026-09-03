/*
Fits both auxiliary models this project settled on - t-QVARMA(1,1,2) and
t-QVARMA(1,1,4) (the driftless base model, not the drift-carrying variant),
the same K_star 3, K_dagger 2, R 1 partition
applications/us_qvarma_employment_change.c's own grid fits to the real
data - to every replicate dataset/abm_system_extract.c wrote under
dataset/abm_system/, and caches every fitted parameter set the moment it
is estimated, not once at the end: each of the 10,800 replicates times 2
specs writes its own JSON the instant that one fit finishes - the same
caching discipline every other application in this project already uses,
applied per (sample, replicate, spec) rather than once per script, because
here there are 21,600 independent fits rather than one.

A cached fit is reused only when it converged. One that stopped at the
iteration cap is carried on instead: its own stored parameters become the
starting guess for another MAX_ITERATIONS, and the result replaces the
cache. L-BFGS returns the best point it saw and the first point it sees is
the cached one, so a resumed fit cannot land below where it started;
the log-likelihood is compared anyway and the cache is left alone if it
would not improve. That is why this does not call qvarma.h's own
fit_cached, which returns any cache that loads: the load and the refit are
spelled out here so the unconverged case can be told apart from the
converged one.

Resuming is per run, not to convergence. One run gives every unconverged
fit one more budget; running it again gives them another. The consequence
is that the estimates depend on how many times this has been run, which
out/abm_system_fit_qvarma_manifest.txt records per fit in its origin and
gain columns.

Parallelized with OpenMP over the 10,800 (sample, replicate) pairs -
schedule(dynamic) since some fits converge in a handful of iterations and
others need the full budget, the same imbalance every multi-start battery
in this project already schedules this way. Both specs for a given pair
run inside the same task rather than as two separate parallel tasks: they
read the same series, so reading it once and fitting both from the one
in-memory copy avoids opening and reparsing the same file twice - each of
the two fit_cached calls still writes its own JSON independently the
moment that one fit finishes.

openblas_set_num_threads(1) below matters for the same reason: without it,
OpenBLAS spawns its own thread pool inside every OpenMP worker thread,
oversubscribing the machine's own cores - confirmed on this project's own
dev box (15 threads, 333% CPU, on a 4-core machine, for a fourfold-parallel
outer loop that should show close to 400% doing useful work rather than
contending with itself). Each task's own matrices are small (K = 5), so
BLAS-level threading buys nothing here - all the real parallelism belongs
to the outer loop over independent fits.

Naming: this file fits t-QVARMA specifically - applications/abm_system_fit_qvarmad.c
is its exact counterpart for the drift-carrying variant, t-QVARMAd, same
partition and specs, deliberately not sharing any code with this file
(docs/MODEL_TEMPLATE.md entry 16: qvarma.h and qvarma_d.h define the same
names, so a translation unit uses one or the other, and neither script
imports anything from the other). Output mirrors dataset/abm_system/'s own
structure exactly, so which sample and which replicate a cached fit belongs
to is never in question - dataset/abm_system/<sample>/replicate_<NNN>.npz's
own input becomes out/abm_system_fit_qvarma/<sample>/replicate_<NNN>_p1q1r2_fit.json
and ..._p1q1r4_fit.json. Rerunning this script after an interruption resumes
rather than redoes: a cache file with the same data fingerprint is loaded
rather than refitted, and refitted from its own parameters rather than from
build_start when it did not converge.

Starting-value convention identical to
applications/us_qvarma_employment_change.c's own build_start(), r fixed at
each of 2 and 4 rather than grid-searched - the same two specs
abm_system_fit_qvarmad.c uses, chosen for direct comparability against that
file's own output rather than a fresh search on simulated data.

out/abm_system_fit_qvarma_manifest.txt records, per (sample, replicate, spec):
converged, log-likelihood, gradient norm and AIC, plus whether the fit was
reused from the cache, carried on from it, or estimated from build_start,
and what carrying it on was worth. Written once at the end - a summary of
what the individual JSON caches already hold, not a substitute for them.

MAX_ITERATIONS is set to the same 2000 abm_system_fit_qvarmad.c uses, as a
starting point for direct comparability - not because this model's own
convergence behaviour at that cap has been measured yet. That file's own
header comment records what 2000 turned out to mean for the drift-carrying
model (measured 2026-08-15: 61% of its cached fits stopped at the cap,
~12% converged); this driftless model's own numbers are whatever
out/abm_system_fit_qvarma_manifest.txt reports once this has actually run,
and may differ - it has one fewer free parameter and no drift term to
compete with the rest of the co-integrated block for identification.

Not a final implementation: no genuine-maximum curvature check per fit
(21,600 Hessians is a different cost than the handful this project checks
by hand elsewhere), no multi-start battery per replicate - one fixed
starting point, matching the real-data grid's own convention, is what
runs here. Whether any of the 21,600 fits need either is a question for
whoever consumes out/abm_system_fit_qvarma_manifest.txt next.

Not part of `make applications` or EXPERIMENT_STEMS - meant to be run
explicitly, once (or resumed after an interruption), not on every routine
build. Requires dataset/abm_system/ to already exist
(`make bin/abm_system_extract && ./bin/abm_system_extract` first). Nothing
printed.
*/

#include "abm_system.h"
#include <et_al./frame/npz.h>
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define START_NU ((mreal)30)
#define MAX_ITERATIONS 2000
#define P 1
#define Q 1

#define INPUT_DIR "dataset/abm_system"
#define OUTPUT_DIR "out/abm_system_fit_qvarma"

typedef struct { int r; const char *label; } Spec;
static const Spec spec_list[] = { { 2, "p1q1r2" }, { 4, "p1q1r4" } };
#define N_SPECS ((int)(sizeof spec_list / sizeof spec_list[0]))

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_fit_qvarma: mkdir failed");
}

static char **list_subdirs(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_fit_qvarma: cannot open dataset/abm_system/ - run abm_system_extract first");

    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char path[560];
        snprintf(path, sizeof path, "%s/%s", dir, entry->d_name);
        struct stat st;
        if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) continue;
        if (n == cap) { cap = cap ? cap * 2 : 16; names = realloc(names, (size_t)cap * sizeof(char*)); }
        size_t name_len = strlen(entry->d_name);
        names[n] = malloc(name_len + 1);
        memcpy(names[n], entry->d_name, name_len + 1);
        n++;
    }
    closedir(handle);
    *count = n;
    return names;
}

/* The replicate indices one sample actually holds, read from the archives
   themselves rather than counted from file names, since a batch whose runs
   did not all succeed is short. Caller must free. */
static int *list_replicates(const char *sample_dir, int *count) {
    char path[560];
    snprintf(path, sizeof path, "%s/%s", INPUT_DIR, sample_dir);
    return abm_system_list_replicates(path, count);
}

static mreal first_difference_sd(Mat y, int row) {
    int periods = y.c;
    Mat difference = mat_new(1, periods - 1);
    for (int t = 1; t < periods; t++) AT(difference, 0, t - 1) = AT(y, row, t) - AT(y, row, t - 1);
    mreal sd = (mreal)sqrt((double)stats_var(difference));
    mat_free(difference);
    return sd;
}

static QvarmaParams build_start(Mat y, int rlag) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, rlag, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    int K_dag = K - K_STAR;

    for (int k = 0; k < K; k++) {
        if (k < K_STAR) {
            Mat row = mat_slice(y, k, k + 1, 0, y.c);
            AT(m.c, k, 0) = stats_mean(row);
        } else {
            AT(m.c, k, 0) = AT(y, k, 0);
        }
    }
    AT(m.Phi_star, 0, 0) = (mreal)0.3;
    for (int a = 0; a < K_STAR; a++) AT(m.Psi_star[0], a, a) = (mreal)0.05;

    for (int a = 0; a < K; a++)
        for (int b = 0; b <= a; b++)
            AT(m.Omega_inv, a, b) = (mreal)(b == a ? first_difference_sd(y, a) : 0);
    m.nu = START_NU;

    for (int l = 0; l < rlag; l++) {
        mreal decay = (mreal)pow(0.6, (double)l);
        for (int i = 0; i < K_dag; i++) AT(m.alpha[l], i, 0) = decay * (mreal)(0.10 - 0.04 * i);
    }
    for (int j = 0; j < K_dag; j++) AT(m.beta[0], 0, j) = (mreal)(j == 0 ? 1 : 0);
    AT(m.beta[0], 0, K_dag - 1) = (mreal)0.3;

    Vec theta = mat_new(qvarma_n_theta(&m), 1);
    _qvarma_unlink(&m, theta);
    qvarma_params_from_theta(theta, &m);
    mat_free(theta);
    return m;
}

/*
The cache, if it holds a fit of this shape for this data. The shape has to be built before
the load: qvarma_load_params compares every field of it, and
mu_star_stationary_only and phi_star_bound both change the length of theta it
expects, so a cache written by this script only loads into params carrying this
script's own settings. It is spelled out from this file's own constants rather
than copied off the starting guess because K and R are macros here.
*/
static int load_cached(QvarmaFitResult *out, Mat y, int rlag, const char *path) {
    out->params = qvarma_params_new(K, K_STAR, P, Q, rlag, R, SHARED_BETA, WARMUP_LONGEST);
    out->params.phi_star_bound = PHI_STAR_BOUND;
    out->params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    if (qvarma_load_fit(out, y, path)) return 1;
    qvarma_params_free(&out->params);
    return 0;
}

/* What one (sample, replicate, spec) fit did, beyond what the fit result itself
   carries: whether it came from the cache untouched, was carried on from a
   cached point, or was estimated from build_start, and what carrying it on was
   worth in log-likelihood. */
typedef enum { FIT_CACHED, FIT_RESUMED, FIT_FRESH } FitOrigin;

typedef struct { int sample_index, replicate; } Task;

int main(void) {
    /* See this file's own header comment for why. */
    openblas_set_num_threads(1);

    make_directory(OUTPUT_DIR);

    int n_samples;
    char **samples = list_subdirs(INPUT_DIR, &n_samples);
    assert(n_samples > 0 && "abm_system_fit_qvarma: no sample subdirectories under dataset/abm_system/");

    int *n_replicates = malloc((size_t)n_samples * sizeof(int));
    int **replicates = malloc((size_t)n_samples * sizeof(int*));
    int n_pairs = 0;
    for (int s = 0; s < n_samples; s++) {
        replicates[s] = list_replicates(samples[s], &n_replicates[s]);
        char out_dir[560];
        snprintf(out_dir, sizeof out_dir, "%s/%s", OUTPUT_DIR, samples[s]);
        make_directory(out_dir);
        n_pairs += n_replicates[s];
    }
    int total_tasks = n_pairs * N_SPECS;

    /* One task per (sample, replicate) pair, not per (sample, replicate,
       spec) triple: both specs share the same series, so it is read once per
       pair and both fit_cached calls reuse the same y in memory. Each fit_cached call still writes its own JSON
       the moment that one fit finishes, independent of the other spec's
       own progress. */
    Task *tasks = malloc((size_t)n_pairs * sizeof(Task));
    int t_idx = 0;
    for (int s = 0; s < n_samples; s++)
        for (int r = 0; r < n_replicates[s]; r++)
            tasks[t_idx++] = (Task){ s, replicates[s][r] };

    double *log_lik = malloc((size_t)total_tasks * sizeof(double));
    double *gradient_norm = malloc((size_t)total_tasks * sizeof(double));
    double *aic = malloc((size_t)total_tasks * sizeof(double));
    int *converged = malloc((size_t)total_tasks * sizeof(int));
    FitOrigin *origin = malloc((size_t)total_tasks * sizeof(FitOrigin));
    double *gain = malloc((size_t)total_tasks * sizeof(double));

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_pairs; i++) {
        Task task = tasks[i];
        const char *sample = samples[task.sample_index];

        char sample_dir[560];
        snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, sample);
        Mat y = abm_system_read_replicate(sample_dir, task.replicate);

        for (int spec = 0; spec < N_SPECS; spec++) {
            int out_idx = i * N_SPECS + spec;
            char cache_path[560];
            snprintf(cache_path, sizeof cache_path, "%s/%s/replicate_%03d_%s_fit.json",
                     OUTPUT_DIR, sample, task.replicate, spec_list[spec].label);

            QvarmaParams start = build_start(y, spec_list[spec].r);
            QvarmaFitOptions options = qvarma_default_fit_options();
            options.max_iterations = MAX_ITERATIONS;

            QvarmaFitResult result;
            QvarmaFitResult cached;
            origin[out_idx] = FIT_FRESH;
            gain[out_idx] = 0;
            if (load_cached(&cached, y, spec_list[spec].r, cache_path)) {
                if (cached.is_converged) {
                    result = cached;
                    origin[out_idx] = FIT_CACHED;
                } else {
                    /* The solver stopped at the iteration cap rather than at
                       the tolerance, so the cached point is a place it was
                       still descending from. Carrying on from there costs one
                       more budget and cannot land above where it started,
                       since L-BFGS returns the best point it saw and the first
                       one it sees here is the cached point. */
                    result = qvarma_fit(y, &cached.params, options);
                    gain[out_idx] = (double)(result.log_likelihood - cached.log_likelihood);
                    if (result.log_likelihood >= cached.log_likelihood) {
                        qvarma_save_fit(&result, y, cache_path);
                        qvarma_fit_result_free(&cached);
                    } else {
                        qvarma_fit_result_free(&result);
                        result = cached;
                    }
                    origin[out_idx] = FIT_RESUMED;
                }
            } else {
                result = qvarma_fit(y, &start, options);
                qvarma_save_fit(&result, y, cache_path);
            }

            log_lik[out_idx] = (double)result.log_likelihood;
            gradient_norm[out_idx] = (double)result.gradient_norm;
            aic[out_idx] = (double)result.aic;
            converged[out_idx] = result.is_converged;

            qvarma_params_free(&start);
            qvarma_fit_result_free(&result);
        }
        mat_free(y);
    }

    FILE *manifest = fopen("out/abm_system_fit_qvarma_manifest.txt", "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "%d samples, %d total (replicate, spec) fits\n\n", n_samples, total_tasks);
    fprintf(manifest, "origin is what this run did with the cache: cached, a converged fit "
                      "reused untouched;\nresumed, an unconverged one carried on from its own "
                      "parameters for another %d iterations;\nfresh, no usable cache, estimated "
                      "from build_start. gain is what resuming was worth\nin log-likelihood, "
                      "blank where nothing was resumed.\n\n", MAX_ITERATIONS);
    fprintf(manifest, "%-28s %8s %6s %14s %10s %10s %10s %8s %10s\n", "sample", "replicate",
            "spec", "log_lik", "gradient", "aic", "converged", "origin", "gain");
    for (int i = 0; i < n_pairs; i++) {
        Task task = tasks[i];
        for (int spec = 0; spec < N_SPECS; spec++) {
            int out_idx = i * N_SPECS + spec;
            const char *origin_name = origin[out_idx] == FIT_CACHED ? "cached"
                                    : origin[out_idx] == FIT_RESUMED ? "resumed" : "fresh";
            char gain_text[32];
            if (origin[out_idx] == FIT_RESUMED) snprintf(gain_text, sizeof gain_text, "%.4f",
                                                         gain[out_idx]);
            else snprintf(gain_text, sizeof gain_text, "%s", "-");
            fprintf(manifest, "%-28s %8d %6s %14.4f %10.4g %10.4f %10s %8s %10s\n",
                    samples[task.sample_index], task.replicate, spec_list[spec].label,
                    log_lik[out_idx], gradient_norm[out_idx], aic[out_idx],
                    converged[out_idx] ? "yes" : "no", origin_name, gain_text);
        }
    }
    fclose(manifest);

    for (int s = 0; s < n_samples; s++) free(samples[s]);
    free(samples);
    for (int s = 0; s < n_samples; s++) free(replicates[s]);
    free(replicates);
    free(n_replicates);
    free(tasks);
    free(log_lik);
    free(gradient_norm);
    free(aic);
    free(converged);
    free(origin);
    free(gain);
    return 0;
}
