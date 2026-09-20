/*
Why resuming an unconverged fit from its own cached parameters leaves half of
them exactly where they were. A study of the points themselves, not a speed
measurement and not a pass-or-fail check.

out/abm_system_fit_qvarma_manifest.txt reports that of 20,497 fits carried on
for a second budget of 2000 iterations, 8,902 came back with the same
log-likelihood and, for 99 percent of those, a bit-identical gradient norm:
L-BFGS took no step at all. A gradient of 1e5 with no step taken means the
line search found nothing downhill, so this walks the steepest-descent ray out
of every cached point and asks whether anything downhill is there to find.

Per cached fit: the objective and its gradient at the stored theta, then
f(theta - s g / ||g||) over a geometric ladder of step lengths, recording the
best decrease found, the step that produced it, and whether the objective went
non-finite further out, which is what a scale matrix at the edge of the
feasible region looks like from inside. The constrained parameters are read
too, since a near-singular Omega_inv or a saturated tanh would explain a
gradient the function does not respond to.

Fits are classified by what the ladder finds, not by what the manifest said,
so the two are independent readings of the same question.

Writes out/qvarma_stuck_fits.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <dirent.h>
#include <string.h>
#include <malloc.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define P 1
#define Q 1

#define INPUT_DIR "dataset/abm_system"
#define CACHE_DIR "out/abm_system_fit_qvarma"

#define N_STEPS 26
#define N_DETAIL 6

typedef struct { int r; const char *label; } Spec;
static const Spec spec_list[] = { { 2, "p1q1r2" }, { 4, "p1q1r4" } };
#define N_SPECS ((int)(sizeof spec_list / sizeof spec_list[0]))

/* What one cached point looks like from the inside. */
typedef struct {
    int loaded, converged;
    double objective, gradient_norm, worst_gradient;
    int worst_index;
    double best_decrease, best_step;
    int went_non_finite;
    double smallest_omega_diagonal, nu, largest_phi;
    LbfgsStatus one_step_status;
    double one_step_gain;
    LbfgsStatus deep_search_status;
    double deep_search_gain;
} Probe;

static char **list_subdirs(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "qvarma_stuck_fits: cannot open dataset/abm_system/");
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
        size_t len = strlen(entry->d_name);
        names[n] = malloc(len + 1);
        memcpy(names[n], entry->d_name, len + 1);
        n++;
    }
    closedir(handle);
    *count = n;
    return names;
}

static int count_replicates(const char *sample) {
    char path[560];
    snprintf(path, sizeof path, "%s/%s", INPUT_DIR, sample);
    int n = 0;
    free(abm_system_list_replicates(path, &n));
    return n;
}

static Mat load_series(const char *sample, int replicate) {
    char path[560];
    snprintf(path, sizeof path, "%s/%s", INPUT_DIR, sample);
    return abm_system_read_replicate(path, replicate);
}

static Probe probe_cached_fit(Mat y, int rlag, const char *cache_path, QvarmaAnalytic *workspace) {
    Probe p;
    memset(&p, 0, sizeof p);

    QvarmaFitResult cached;
    cached.params = qvarma_params_new(K, K_STAR, P, Q, rlag, R, SHARED_BETA, WARMUP_LONGEST);
    cached.params.phi_star_bound = PHI_STAR_BOUND;
    cached.params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    if (!qvarma_load_fit(&cached, y, cache_path)) {
        qvarma_params_free(&cached.params);
        return p;
    }
    p.loaded = 1;
    p.converged = cached.is_converged;

    p.smallest_omega_diagonal = (double)AT(cached.params.Omega_inv, 0, 0);
    for (int k = 1; k < K; k++) {
        double d = (double)AT(cached.params.Omega_inv, k, k);
        if (d < p.smallest_omega_diagonal) p.smallest_omega_diagonal = d;
    }
    p.nu = (double)cached.params.nu;
    p.largest_phi = (double)MABS(AT(cached.params.Phi_star, 0, 0));

    int n = qvarma_n_theta(&cached.params);
    Vec theta = mat_new(n, 1);
    _qvarma_unlink(&cached.params, theta);

    Vec gradient = mat_new(n, 1);
    QvarmaFitContext context = { y, &cached.params, workspace };
    double at_theta = (double)qvarma_negative_log_likelihood(theta, gradient, &context);
    p.objective = at_theta;

    double squared = 0;
    for (int i = 0; i < n; i++) {
        double g = (double)gradient.d[i];
        squared += g * g;
        if (fabs(g) > p.worst_gradient) { p.worst_gradient = fabs(g); p.worst_index = i; }
    }
    p.gradient_norm = sqrt(squared);

    /* The steepest-descent ray, which is the direction L-BFGS falls back to on
       its first iteration after a restart, since the curvature history is
       empty there. */
    Vec trial = mat_new(n, 1);
    Vec no_gradient = { 0, 0, 0, NULL };
    for (int s = 0; s < N_STEPS; s++) {
        double step = pow(10.0, -14.0 + (double)s);
        for (int i = 0; i < n; i++)
            trial.d[i] = theta.d[i] - (mreal)(step * (double)gradient.d[i] / p.gradient_norm);
        double value = (double)qvarma_negative_log_likelihood(trial, no_gradient, &context);
        if (!isfinite(value)) { p.went_non_finite = 1; continue; }
        double decrease = at_theta - value;
        if (decrease > p.best_decrease) { p.best_decrease = decrease; p.best_step = step; }
    }

    /* What L-BFGS itself does from here, one iteration, which is the iteration
       the resume in applications/abm_system_fit_qvarma.c either gets past or
       does not. A fresh run has no curvature history, so its first direction
       is -g and its first trial step moves unit length in parameter space,
       the same ray the ladder above walks. */
    QvarmaFitOptions options = qvarma_default_fit_options();
    options.max_iterations = 1;
    QvarmaFitResult one = qvarma_fit(y, &cached.params, options);
    p.one_step_status = one.status;
    p.one_step_gain = (double)(one.log_likelihood - cached.log_likelihood);
    qvarma_fit_result_free(&one);

    /* The same iteration with a line search allowed to bisect far enough to
       reach the distances the ladder finds the decrease at. max_line_search is
       not a field of QvarmaFitOptions, so this is the one place in this project
       that calls the solver directly rather than through fit: the solver's own
       setting is what is being studied, which is the case et_al.'s model policy
       makes the exception for. */
    QvarmaFitContext deep_context = { y, &cached.params, workspace };
    LbfgsOptions deep = lbfgs_default_options();
    deep.max_iterations = 1;
    deep.max_line_search = 60;
    Vec deep_start = mat_new(n, 1);
    _qvarma_unlink(&cached.params, deep_start);
    LbfgsResult deep_result = lbfgs(qvarma_negative_log_likelihood, &deep_context,
                                    deep_start, deep);
    p.deep_search_status = deep_result.status;
    p.deep_search_gain = (double)(-deep_result.value - cached.log_likelihood);
    mat_free(deep_result.theta);
    mat_free(deep_start);

    mat_free(trial);
    mat_free(gradient);
    mat_free(theta);
    qvarma_fit_result_free(&cached);
    return p;
}

static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

static int compare_double(const void *a, const void *b) {
    double x = *(const double*)a, z = *(const double*)b;
    return x < z ? -1 : x > z ? 1 : 0;
}

static double quantile(double *sorted, int n, double q) {
    if (n == 0) return 0;
    int at = (int)(q * (n - 1));
    return sorted[at];
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);

    int n_samples;
    char **samples = list_subdirs(INPUT_DIR, &n_samples);
    int *n_replicates = malloc((size_t)n_samples * sizeof(int));
    int n_pairs = 0;
    for (int s = 0; s < n_samples; s++) {
        n_replicates[s] = count_replicates(samples[s]);
        n_pairs += n_replicates[s];
    }
    int total = n_pairs * N_SPECS;

    Probe *probe = malloc((size_t)total * sizeof(Probe));
    int *sample_of = malloc((size_t)total * sizeof(int));
    int *replicate_of = malloc((size_t)total * sizeof(int));

    int *pair_sample = malloc((size_t)n_pairs * sizeof(int));
    int *pair_replicate = malloc((size_t)n_pairs * sizeof(int));
    int at = 0;
    for (int s = 0; s < n_samples; s++)
        for (int r = 0; r < n_replicates[s]; r++) { pair_sample[at] = s; pair_replicate[at] = r; at++; }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_pairs; i++) {
        Mat y = load_series(samples[pair_sample[i]], pair_replicate[i]);
        for (int spec = 0; spec < N_SPECS; spec++) {
            char cache_path[560];
            snprintf(cache_path, sizeof cache_path, "%s/%s/replicate_%03d_%s_fit.json",
                     CACHE_DIR, samples[pair_sample[i]], pair_replicate[i], spec_list[spec].label);
            QvarmaParams shape = qvarma_params_new(K, K_STAR, P, Q, spec_list[spec].r, R,
                                                   SHARED_BETA, WARMUP_LONGEST);
            shape.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
            QvarmaAnalytic *workspace = qvarma_analytic_new(&shape, y.c);
            int idx = i * N_SPECS + spec;
            probe[idx] = probe_cached_fit(y, spec_list[spec].r, cache_path, workspace);
            sample_of[idx] = pair_sample[i];
            replicate_of[idx] = pair_replicate[i];
            qvarma_analytic_free(workspace);
            qvarma_params_free(&shape);
        }
        mat_free(y);
    }

    FILE *report = fopen("out/qvarma_stuck_fits.txt", "w");
    assert(report);
    fprintf(report, "Every cached fit under %s, probed at its own stored theta.\n", CACHE_DIR);
    fprintf(report, "The ladder is f(theta - s g/||g||) for s from 1e-14 to 1e11, "
                    "26 steps, one decade apart.\n");
    fprintf(report, "A fit is called stuck when no s on the ladder lowers the objective.\n\n");

    int n_unconverged = 0, n_stuck = 0, n_movable = 0;
    for (int i = 0; i < total; i++) {
        if (!probe[i].loaded || probe[i].converged) continue;
        n_unconverged++;
        if (probe[i].best_decrease > 0) n_movable++; else n_stuck++;
    }
    fprintf(report, "unconverged cached fits          %6d\n", n_unconverged);
    fprintf(report, "  no downhill step on the ladder %6d  (%.1f%%)\n", n_stuck,
            100.0 * n_stuck / n_unconverged);
    fprintf(report, "  a downhill step exists         %6d  (%.1f%%)\n\n", n_movable,
            100.0 * n_movable / n_unconverged);

    static const char *status_name[5] = { "max iterations", "gradient tolerance",
                                          "function tolerance", "no progress", "not finite" };
    int status_tally[5] = { 0, 0, 0, 0, 0 }, moved = 0;
    for (int i = 0; i < total; i++) {
        if (!probe[i].loaded || probe[i].converged) continue;
        status_tally[probe[i].one_step_status]++;
        if (probe[i].one_step_gain > 0) moved++;
    }
    fprintf(report, "why one L-BFGS iteration from the cached point stops\n");
    for (int k = 0; k < 5; k++)
        if (status_tally[k])
            fprintf(report, "  %-22s %6d  (%.1f%%)\n", status_name[k], status_tally[k],
                    100.0 * status_tally[k] / n_unconverged);
    fprintf(report, "  log-likelihood rose    %6d  (%.1f%%)\n\n", moved,
            100.0 * moved / n_unconverged);

    /* How far along the ray the decrease is, against how far the line search
       can look. The first trial step moves unit length in parameter space and
       the bracket is bisected at most max_line_search times, so nothing closer
       than 2^-20, about 9.5e-7, is reachable. */
    {
        double reach = pow(2.0, -20.0);
        int no_progress = 0, downhill = 0, out_of_reach = 0;
        double *step_value = malloc((size_t)n_unconverged * sizeof(double));
        int m = 0;
        for (int i = 0; i < total; i++) {
            if (!probe[i].loaded || probe[i].converged) continue;
            if (probe[i].one_step_status != LBFGS_NO_PROGRESS) continue;
            no_progress++;
            if (probe[i].best_decrease <= 0) continue;
            downhill++;
            step_value[m++] = probe[i].best_step;
            if (probe[i].best_step < reach) out_of_reach++;
        }
        qsort(step_value, (size_t)m, sizeof(double), compare_double);
        fprintf(report, "the %d fits the line search gave up on\n", no_progress);
        fprintf(report, "  a downhill step exists on the ladder for %d of them (%.1f%%)\n",
                downhill, 100.0 * downhill / no_progress);
        fprintf(report, "  its distance in parameter space, quantiles: "
                        "min %.3g  median %.3g  p90 %.3g  max %.3g\n",
                quantile(step_value, m, 0), quantile(step_value, m, 0.5),
                quantile(step_value, m, 0.9), quantile(step_value, m, 1.0));
        fprintf(report, "  closer than the 2^-20 = %.3g the bisection can reach: %d (%.1f%%)\n\n",
                reach, out_of_reach, 100.0 * out_of_reach / downhill);
        free(step_value);
    }

    /* Whether letting the bisection go further is enough on its own. */
    {
        int considered = 0, still_stuck = 0, now_moves = 0;
        double total_gain = 0, best = 0;
        for (int i = 0; i < total; i++) {
            if (!probe[i].loaded || probe[i].converged) continue;
            if (probe[i].one_step_status != LBFGS_NO_PROGRESS) continue;
            considered++;
            if (probe[i].deep_search_status == LBFGS_NO_PROGRESS) still_stuck++;
            if (probe[i].deep_search_gain > 0) {
                now_moves++;
                total_gain += probe[i].deep_search_gain;
                if (probe[i].deep_search_gain > best) best = probe[i].deep_search_gain;
            }
        }
        fprintf(report, "the same one iteration with max_line_search raised from 20 to 60\n");
        fprintf(report, "  fits the default gave up on      %6d\n", considered);
        fprintf(report, "  still no progress                %6d  (%.1f%%)\n", still_stuck,
                100.0 * still_stuck / considered);
        fprintf(report, "  log-likelihood rose              %6d  (%.1f%%)\n", now_moves,
                100.0 * now_moves / considered);
        fprintf(report, "  mean gain over those, one iteration %.6f, largest %.6f\n\n",
                now_moves ? total_gain / now_moves : 0.0, best);
    }

    /* Distributions, stuck against movable, over the quantities that would
       explain a gradient the objective does not respond to. */
    const char *column[5] = { "gradient norm", "worst |g_i|", "min diag Omega_inv",
                              "nu", "|Phi_star|" };
    for (int group = 0; group < 2; group++) {
        int want_stuck = group == 0;
        int count = want_stuck ? n_stuck : n_movable;
        if (count == 0) continue;
        double *value = malloc((size_t)count * sizeof(double));
        fprintf(report, "%s (%d fits), quantiles\n", want_stuck ? "stuck" : "movable", count);
        fprintf(report, "  %-20s %12s %12s %12s %12s\n", "", "min", "median", "p90", "max");
        for (int c = 0; c < 5; c++) {
            int m = 0;
            for (int i = 0; i < total; i++) {
                if (!probe[i].loaded || probe[i].converged) continue;
                int stuck = probe[i].best_decrease <= 0;
                if (stuck != want_stuck) continue;
                value[m++] = c == 0 ? probe[i].gradient_norm
                           : c == 1 ? probe[i].worst_gradient
                           : c == 2 ? probe[i].smallest_omega_diagonal
                           : c == 3 ? probe[i].nu : probe[i].largest_phi;
            }
            qsort(value, (size_t)m, sizeof(double), compare_double);
            fprintf(report, "  %-20s %12.4g %12.4g %12.4g %12.4g\n", column[c],
                    quantile(value, m, 0), quantile(value, m, 0.5),
                    quantile(value, m, 0.9), quantile(value, m, 1.0));
        }
        int non_finite = 0;
        for (int i = 0; i < total; i++) {
            if (!probe[i].loaded || probe[i].converged) continue;
            int stuck = probe[i].best_decrease <= 0;
            if (stuck == want_stuck && probe[i].went_non_finite) non_finite++;
        }
        fprintf(report, "  objective non-finite somewhere on the ladder: %d of %d (%.1f%%)\n\n",
                non_finite, count, 100.0 * non_finite / count);
        free(value);
    }

    /* Which parameter carries the largest gradient, stuck fits only. */
    fprintf(report, "the coordinate carrying the largest gradient, stuck fits\n");
    for (int spec = 0; spec < N_SPECS; spec++) {
        QvarmaParams shape = qvarma_params_new(K, K_STAR, P, Q, spec_list[spec].r, R,
                                               SHARED_BETA, WARMUP_LONGEST);
        shape.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
        int n = qvarma_n_theta(&shape);
        int *tally = calloc((size_t)n, sizeof(int));
        int seen = 0;
        for (int i = 0; i < total; i++) {
            if (!probe[i].loaded || probe[i].converged || probe[i].best_decrease > 0) continue;
            if (i % N_SPECS != spec) continue;
            tally[probe[i].worst_index]++;
            seen++;
        }
        fprintf(report, "  %s, %d stuck fits\n", spec_list[spec].label, seen);
        for (int pass = 0; pass < 4; pass++) {
            int best = 0;
            for (int i = 1; i < n; i++) if (tally[i] > tally[best]) best = i;
            if (tally[best] == 0) break;
            char name[64];
            _qvarma_theta_name(&shape, best, name, (int)sizeof name);
            fprintf(report, "    %-22s %6d  (%.1f%%)\n", name, tally[best],
                    100.0 * tally[best] / seen);
            tally[best] = 0;
        }
        free(tally);
        qvarma_params_free(&shape);
    }

    fprintf(report, "\nthe first %d stuck fits in directory order, in full\n", N_DETAIL);
    fprintf(report, "  %-28s %4s %6s %12s %12s %12s %10s %8s\n", "sample", "rep", "spec",
            "objective", "grad norm", "min diag", "nu", "non-fin");
    int shown = 0;
    for (int i = 0; i < total && shown < N_DETAIL; i++) {
        if (!probe[i].loaded || probe[i].converged || probe[i].best_decrease > 0) continue;
        fprintf(report, "  %-28s %4d %6s %12.4f %12.4g %12.4g %10.4g %8s\n",
                samples[sample_of[i]], replicate_of[i], spec_list[i % N_SPECS].label,
                probe[i].objective, probe[i].gradient_norm, probe[i].smallest_omega_diagonal,
                probe[i].nu, probe[i].went_non_finite ? "yes" : "no");
        shown++;
    }

    fclose(report);
    for (int s = 0; s < n_samples; s++) free(samples[s]);
    free(samples); free(n_replicates); free(probe);
    free(sample_of); free(replicate_of); free(pair_sample); free(pair_replicate);
    return 0;
}
