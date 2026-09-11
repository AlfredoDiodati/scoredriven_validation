/*
Fits the auxiliary model t-QVARMA(1,1,2) (the driftless base model, not the
drift-carrying variant), the same K_star 3, K_dagger 2, R 1 partition
applications/us_qvarma_spec_choice.c's own grid fits to the real
data - to every replicate under dataset/abm_system/, and caches every fitted
parameter set the moment it is estimated, not once at the end: each replicate
writes its own JSON the instant its fit finishes - the same caching discipline
every other application in this project already uses, applied per (sample,
replicate, spec) rather than once per script, because here there are a million
independent fits rather than one.

spec_list held t-QVARMA(1,1,4) alongside r = 2 while the dataset was the 10,800
replicates the .Rdata route produced. r = 4 was dropped because it lost the
Model Confidence Set: the joint run over both specs put no r = 4 configuration
in the surviving set at all, which docs/ABM_SYSTEM_MCS_VALIDATION.md records.
Adding it back is one entry in spec_list; applications/abm_system_mse_qvarma.c
asserts on exactly one spec and has to match whatever this list holds.

A cached fit is reused only when it converged. An unconverged one is treated
according to why the solver stopped, which the cache now records.

Stopped at the iteration cap: carried on. Its own stored parameters become the
starting guess for another MAX_ITERATIONS and the result replaces the cache.
L-BFGS returns the best point it saw and the first point it sees is the cached
one, so a resumed fit cannot land below where it started; the log-likelihood is
compared anyway and the cache is left alone if it would not improve.

Stopped because the line search could not lower the objective: held, not
resumed. Over the whole design 504,005 fits reported that and 502,777 of them
gained exactly nothing when they were resumed anyway, so carrying them on is a
full pass of a million fits spent rewriting the point already on disk. They are
left as they are, with the reason in their cache, and the manifest calls them
held. What they need is the cause fixed, not more iterations.

A cache whose reason was never recorded is resumed once, which is what gives it
a reason to be sorted by next time.

That is why this does not call qvarma.h's own fit_cached: the load and the
refit are spelled out here so the three cases can be told apart and counted.

Resuming is per run, not to convergence. One run gives every unconverged
fit one more budget; running it again gives them another. The consequence
is that the estimates depend on how many times this has been run, which
out/abm_system_fit_qvarma_manifest.txt records per fit in its origin and
gain columns.

Parallelized with OpenMP over the million (sample, replicate) pairs -
schedule(dynamic) since some fits converge in a handful of iterations and
others need the full budget, the same imbalance every multi-start battery
in this project already schedules this way. Every spec for a given pair
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
applications/us_qvarma_spec_choice.c's own build_start(), r fixed at
at 2 rather than grid-searched - one of the two specs abm_system_fit_qvarmad.c
uses, chosen for direct comparability against that file's own output rather
than a fresh search on simulated data.

out/abm_system_fit_qvarma_manifest.txt records, per (sample, replicate, spec):
converged, log-likelihood, gradient norm and AIC, plus whether the fit was
reused from the cache, carried on from it, or estimated from build_start,
and what carrying it on was worth. Written once at the end - a summary of
what the individual JSON caches already hold, not a substitute for them.

MAX_ITERATIONS is 1000000. It has been raised four times and every raise is
recorded here, with what it bought, so the next person does not repeat the
experiment.

It started at 2000, matched from the drift-carrying study for comparability
rather than from any measurement of this model. Each run carries on the fits
that ran out of budget and leaves alone the ones that stopped moving, so the
"improvable" column below is the pool the cap can act on at all:

    run   cap      improvable at the start   converged   hit rate
    1     2000                     925,218     194,739        21%
    2     12000                    179,807      21,749        12%
    3     12000                     82,625       8,174        10%
    4     30000                     41,427       5,168        12%
    5     30000                     19,494         173       0.9%
    6     100000                     6,887         131       1.9%
    7     100000                     5,845           7       0.1%

Cumulative convergence went 7.5, 27.0, 31.6, 34.6, 35.1, 35.2, 35.25, 35.25 per
cent. The improvable pool fell from 925,218 to 5,597 while the pool that stops
because no step along the search direction lowers the objective grew to
621,205, and those are held rather than carried on because resuming them is a
full pass spent rewriting the point already on disk.

What the record says plainly is that the cap stopped being the binding
constraint somewhere around run 3. Going from 12000 to 30000 did not raise the
hit rate; going from 30000 to 100000 converted 131 fits and then 7. The 5,597
still running out of budget have had of the order of 286,000 iterations each
and keep moving the objective by more than the function tolerance without the
gradient test ever passing. 1000000 is aimed at them, and the manifest is
where the answer will be.

Raising it is cheap to try because the resume path does not restart a fit: an
unconverged one is carried on from its own stored parameters, and iterations
already spent are kept, which is what the iterations column counts.

Not a final implementation: no genuine-maximum curvature check per fit
(a million Hessians is a different cost than the handful this project checks
by hand elsewhere), no multi-start battery per replicate - one fixed
starting point, matching the real-data grid's own convention, is what
runs here. Whether any of the million fits need either is a question for
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
#define MAX_ITERATIONS 1000000
#define P 1
#define Q 1

#define INPUT_DIR_DEFAULT "dataset/abm_system"
#define OUTPUT_DIR_DEFAULT "out/abm_system_fit_qvarma"

/* Both directories are overridable so the whole procedure can be exercised on a
   copy of a configuration or two before it is turned loose on the real cache,
   which it rewrites in place. */
static const char *INPUT_DIR;
static const char *OUTPUT_DIR;

typedef struct { int r; const char *label; } Spec;
static const Spec spec_list[] = { { 2, "p1q1r2" } };
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
typedef enum { FIT_CACHED, FIT_RESUMED, FIT_FRESH, FIT_HELD } FitOrigin;

typedef struct { int sample_index, replicate; } Task;

/* Where a fit's parameters came from, when they came from a sibling replicate
   rather than from its own solver run. -1 means its own. Kept beside the
   caches rather than inside them because it is this file's own bookkeeping and
   not a property of a QVARMA fit, and read back on the next run because the
   rule below needs to know which donors a fit has already descended from. */
static void read_lineage(const char *sample, int *donor, int n_replicates) {
    for (int r = 0; r < n_replicates; r++) donor[r] = -1;

    char path[600];
    snprintf(path, sizeof path, "%s/%s/lineage.txt", OUTPUT_DIR, sample);
    FILE *f = fopen(path, "r");
    if (!f) return;

    int replicate, from;
    while (fscanf(f, "%d %d", &replicate, &from) == 2)
        if (replicate >= 0 && replicate < n_replicates) donor[replicate] = from;
    fclose(f);
}

static void write_lineage(const char *sample, const int *donor, int n_replicates) {
    int any = 0;
    for (int r = 0; r < n_replicates && !any; r++) any = donor[r] >= 0;
    if (!any) return;

    char path[600];
    snprintf(path, sizeof path, "%s/%s/lineage.txt", OUTPUT_DIR, sample);
    FILE *f = fopen(path, "w");
    assert(f && "abm_system_fit_qvarma: cannot write a lineage file");
    for (int r = 0; r < n_replicates; r++)
        if (donor[r] >= 0) fprintf(f, "%d %d\n", r, donor[r]);
    fclose(f);
}

/* The log-likelihood of one parameter set on one replicate's own data, with
   the gradient at that point, so a cache written from it reports the same two
   numbers a fit would. */
static double likelihood_of(Vec theta, const QvarmaParams *shape, Mat y, double *gradient_norm) {
    QvarmaAnalytic *workspace = qvarma_analytic_new(shape, y.c);
    Vec gradient = mat_new(theta.r, 1);
    double value = (double)qvarma_analytic_log_likelihood(workspace, theta, y, gradient);

    double sum = 0;
    for (int i = 0; i < gradient.r; i++) sum += (double)gradient.d[i] * (double)gradient.d[i];
    *gradient_norm = sqrt(sum);

    mat_free(gradient);
    qvarma_analytic_free(workspace);
    return value;
}

int main(void) {
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = INPUT_DIR_DEFAULT;
    OUTPUT_DIR = getenv("ABM_SYSTEM_FIT_DIR");
    if (!OUTPUT_DIR) OUTPUT_DIR = OUTPUT_DIR_DEFAULT;

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
       spec) triple: every spec shares the same series, so it is read once per
       pair and reused in memory. Each fit writes its own JSON the moment it
       finishes, independent of any other spec's progress. */
    Task *tasks = malloc((size_t)n_pairs * sizeof(Task));
    int t_idx = 0;
    for (int s = 0; s < n_samples; s++)
        for (int r = 0; r < n_replicates[s]; r++)
            tasks[t_idx++] = (Task){ s, replicates[s][r] };

    double *log_lik = malloc((size_t)total_tasks * sizeof(double));
    double *gradient_norm = malloc((size_t)total_tasks * sizeof(double));
    double *aic = malloc((size_t)total_tasks * sizeof(double));
    int *converged = malloc((size_t)total_tasks * sizeof(int));
    int *total_niter = malloc((size_t)total_tasks * sizeof(int));
    /* The solver's own reason, or -1 where the cache never recorded one. */
    int *stop_reason = malloc((size_t)total_tasks * sizeof(int));
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
                } else if (cached.status_is_known && cached.status != LBFGS_MAX_ITERATIONS) {
                    /* The solver did not run out of budget, it ran out of
                       moves: at this point no step along its search direction
                       lowered the objective. Measured over the whole design,
                       504,005 fits reported that and 502,777 of them gained
                       exactly nothing when resumed anyway, so resuming is a
                       full pass of the dataset spent reproducing the point
                       already on disk. It is held instead, and the reason it
                       stopped stays in its cache for whoever fixes the cause. */
                    result = cached;
                    origin[out_idx] = FIT_HELD;
                } else {
                    /* The solver stopped at the iteration cap rather than at
                       the tolerance, so the cached point is a place it was
                       still descending from. Carrying on from there costs one
                       more budget and cannot land above where it started,
                       since L-BFGS returns the best point it saw and the first
                       one it sees here is the cached point. */
                    result = qvarma_fit(y, &cached.params, options);
                    /* qvarma_fit reports its own run only, so the chain the
                       cache carries is extended here, the way
                       qvarma_fit_cached does it, because this file drives the
                       resume itself. total_niter is the real cumulative count
                       either way. A cache written before the chain existed has
                       no earlier reasons, so the chain starts at this run and
                       nruns counts the runs whose reason is recorded, which is
                       what the stored array holds. */
                    result.total_niter = cached.total_niter + result.niter;
                    if (cached.status_is_known) {
                        result.nruns = cached.nruns + 1;
                        free(result.run_status);
                        result.run_status = _qvarma_run_status_extend(cached.run_status,
                                                                     cached.nruns, result.status);
                    }
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
            total_niter[out_idx] = result.total_niter;
            stop_reason[out_idx] = result.status_is_known ? (int)result.status : -1;

            qvarma_params_free(&start);
            qvarma_fit_result_free(&result);
        }
        mat_free(y);
    }

    /* Cross-replicate restart. Replicates of one configuration differ only by
       seed, so a parameter set that fits one of them well is a candidate
       starting point for the rest, and the earlier runs showed the solver's
       outcome depends heavily on where it starts: half the fits stopped
       because the line search could not move, and a fifth of a sample of those
       converged when restarted from a sibling.

       Within each configuration the replicate with the highest log-likelihood
       is the donor. Its parameters are evaluated on every other replicate's
       own data, and adopted wherever they beat that replicate's current fit.
       The comparison is on the recipient's own likelihood, so an adopted set
       is a better optimum of the recipient's own objective and not a borrowed
       estimate.

       An adopted fit is marked unconverged with no recorded reason, which is
       what makes the next run take it up and fit from there. A replicate that
       already descends from this donor is left alone, whether it inherited
       from it just now or inherited and was fitted since: its current
       parameters came from that donor and are at least as good. */
    int *inherited_from = malloc((size_t)total_tasks * sizeof(int));
    int *switched = malloc((size_t)total_tasks * sizeof(int));
    for (int i = 0; i < total_tasks; i++) { inherited_from[i] = -1; switched[i] = 0; }

    int sample_offset = 0;
    int *offset_of = malloc((size_t)n_samples * sizeof(int));
    for (int s = 0; s < n_samples; s++) { offset_of[s] = sample_offset; sample_offset += n_replicates[s]; }

    #pragma omp parallel for schedule(dynamic)
    for (int s = 0; s < n_samples; s++) {
        int first = offset_of[s], count = n_replicates[s];
        int *donor = malloc((size_t)count * sizeof(int));
        read_lineage(samples[s], donor, count);

        char sample_dir[560];
        snprintf(sample_dir, sizeof sample_dir, "%s/%s", INPUT_DIR, samples[s]);

        for (int spec = 0; spec < N_SPECS; spec++) {
            int best = -1;
            double best_value = 0;
            for (int r = 0; r < count; r++) {
                int idx = (first + r) * N_SPECS + spec;
                if (best < 0 || log_lik[idx] > best_value) { best = r; best_value = log_lik[idx]; }
            }
            if (best < 0) continue;

            QvarmaParams shape = qvarma_params_new(K, K_STAR, P, Q, spec_list[spec].r, R,
                                                   SHARED_BETA, WARMUP_LONGEST);
            shape.phi_star_bound = PHI_STAR_BOUND;
            shape.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;

            char donor_path[600];
            snprintf(donor_path, sizeof donor_path, "%s/%s/replicate_%03d_%s_fit.json",
                     OUTPUT_DIR, samples[s], tasks[first + best].replicate, spec_list[spec].label);
            QvarmaFitResult source = qvarma_fit_result_new(&shape);
            Mat donor_y = abm_system_read_replicate(sample_dir, tasks[first + best].replicate);
            int have_donor = qvarma_load_fit(&source, donor_y, donor_path);
            mat_free(donor_y);
            if (!have_donor) {
                qvarma_fit_result_free(&source);
                qvarma_params_free(&shape);
                continue;
            }

            Vec theta = mat_new(qvarma_n_theta(&source.params), 1);
            _qvarma_unlink(&source.params, theta);

            for (int r = 0; r < count; r++) {
                if (r == best) continue;
                if (donor[r] == tasks[first + best].replicate) continue;

                int idx = (first + r) * N_SPECS + spec;
                Mat y = abm_system_read_replicate(sample_dir, tasks[first + r].replicate);

                double gradient_at, value = likelihood_of(theta, &shape, y, &gradient_at);
                if (value > log_lik[idx]) {
                    QvarmaFitResult adopted = qvarma_fit_result_new(&shape);
                    qvarma_params_from_theta(theta, &adopted.params);
                    adopted.log_likelihood = (mreal)value;
                    adopted.gradient_norm = (mreal)gradient_at;
                    mreal k = (mreal)theta.r, periods = (mreal)y.c, mean = (mreal)value / periods;
                    adopted.aic = 2 * k / periods - 2 * mean;
                    adopted.bic = k * (mreal)log((double)periods) / periods - 2 * mean;
                    adopted.hannan_quinn = 2 * k * (mreal)log(log((double)periods)) / periods - 2 * mean;
                    adopted.niter = 0;
                    adopted.total_niter = total_niter[idx];
                    adopted.nruns = 1;
                    adopted.is_converged = 0;
                    adopted.status_is_known = 0;
                    free(adopted.run_status);
                    adopted.run_status = NULL;

                    char path[600];
                    snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json",
                             OUTPUT_DIR, samples[s], tasks[first + r].replicate, spec_list[spec].label);
                    qvarma_save_fit(&adopted, y, path);
                    qvarma_fit_result_free(&adopted);

                    log_lik[idx] = value;
                    gradient_norm[idx] = gradient_at;
                    aic[idx] = 2 * (double)theta.r / y.c - 2 * value / y.c;
                    converged[idx] = 0;
                    stop_reason[idx] = -1;
                    switched[idx] = 1;
                    donor[r] = tasks[first + best].replicate;
                }
                inherited_from[idx] = donor[r];
                mat_free(y);
            }

            mat_free(theta);
            qvarma_fit_result_free(&source);
            qvarma_params_free(&shape);
        }

        write_lineage(samples[s], donor, count);
        free(donor);
    }
    free(offset_of);

    char manifest_path[640];
    snprintf(manifest_path, sizeof manifest_path, "%s_manifest.txt", OUTPUT_DIR);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "%d samples, %d total (replicate, spec) fits\n\n", n_samples, total_tasks);
    fprintf(manifest, "origin is what this run did with the cache: cached, a converged fit "
                      "reused untouched;\nresumed, an unconverged one carried on from its own "
                      "parameters for another %d iterations;\nfresh, no usable cache, estimated "
                      "from build_start;\nheld, an unconverged one left alone because the solver "
                      "reported it could not move rather\nthan that it ran out of budget. gain is "
                      "what resuming was worth in log-likelihood,\nblank where nothing was "
                      "resumed;\ninherited, the parameters of the best-fitting replicate of the same "
                      "configuration, adopted\nbecause they scored higher on this replicate's own "
                      "data than its own fit did.\n\n", MAX_ITERATIONS);
    fprintf(manifest, "iterations is cumulative over every run this fit has had. stopped_because "
                      "is the\nsolver's own reason for the most recent run, and reads unrecorded "
                      "for a fit whose\ncache predates the reason being stored, and inherited from replicate N for one that\njust took another replicate's parameters. An inherited fit is marked unconverged so\nthe next run fits on from where it now sits.\n\n");
    fprintf(manifest, "%-28s %8s %6s %14s %10s %10s %10s %8s %10s %10s %-34s\n", "sample",
            "replicate", "spec", "log_lik", "gradient", "aic", "converged", "origin", "gain",
            "iterations", "stopped_because");
    for (int i = 0; i < n_pairs; i++) {
        Task task = tasks[i];
        for (int spec = 0; spec < N_SPECS; spec++) {
            int out_idx = i * N_SPECS + spec;
            const char *origin_name = switched[out_idx] ? "inherited"
                                    : origin[out_idx] == FIT_CACHED ? "cached"
                                    : origin[out_idx] == FIT_RESUMED ? "resumed"
                                    : origin[out_idx] == FIT_HELD ? "held" : "fresh";
            char gain_text[32];
            if (origin[out_idx] == FIT_RESUMED) snprintf(gain_text, sizeof gain_text, "%.4f",
                                                         gain[out_idx]);
            else snprintf(gain_text, sizeof gain_text, "%s", "-");
            char reason_text[64];
            const char *reason;
            if (switched[out_idx]) {
                snprintf(reason_text, sizeof reason_text, "inherited from replicate %d",
                         inherited_from[out_idx]);
                reason = reason_text;
            } else if (stop_reason[out_idx] < 0) {
                reason = "unrecorded";
            } else {
                reason = lbfgs_status_text((LbfgsStatus)stop_reason[out_idx]);
            }
            fprintf(manifest, "%-28s %8d %6s %14.4f %10.4g %10.4f %10s %8s %10s %10d %-34s\n",
                    samples[task.sample_index], task.replicate, spec_list[spec].label,
                    log_lik[out_idx], gradient_norm[out_idx], aic[out_idx],
                    converged[out_idx] ? "yes" : "no", origin_name, gain_text,
                    total_niter[out_idx], reason);
        }
    }
    fclose(manifest);

    for (int s = 0; s < n_samples; s++) free(samples[s]);
    free(samples);
    for (int s = 0; s < n_samples; s++) free(replicates[s]);
    free(replicates);
    free(n_replicates);
    free(tasks);
    free(inherited_from);
    free(switched);
    free(total_niter);
    free(stop_reason);
    free(log_lik);
    free(gradient_norm);
    free(aic);
    free(converged);
    free(origin);
    free(gain);
    return 0;
}
