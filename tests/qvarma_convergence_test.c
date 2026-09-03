/*
Which convergence test tells a fit that has arrived from one that has not,
without the answer depending on the sample size or the units the data is in.
A study on simulated data, where the point the fit is heading for is known
because it can be reached with a budget nothing in production would pay for.

The test in solver/lbfgs.h stops when max_i |g_i| <= gradient_tolerance. That
threshold is absolute on the total log-likelihood, which grows with T and moves
with the scale of the data, so the same statistical situation crosses it at
different points. Its own comment says as much. Three replacements are scored
against it here:

  worst_gradient      max_i |g_i|, what lbfgs uses now
  relative_gradient   max_i |g_i| max(|theta_i|,1) / max(|f|,1), the usual
                      scaled form, dimensionless in both f and theta
  newton_decrement    g' H^-1 g, which near a minimum is twice the objective
                      still to be gained, and is unchanged by any linear
                      reparameterization of theta
  per_observation     max_i |g_i| / T, the cheapest fix, which removes the
                      sample size but not the curvature

Ground truth is not the parameters the data was simulated from: the maximum
likelihood estimate differs from those by sampling error, and a fit that
reaches the estimate has converged whatever the truth was. Ground truth here is
distance from the point the same fit reaches with REFERENCE_BUDGET iterations,
measured in standard errors of the estimate itself, since a movement far below
one standard error cannot change any conclusion drawn from the fit:

  se_movement = max_i |theta_i - theta_star_i| / se_i,  se_i = sqrt((H^-1)_ii)

A fit is called arrived when se_movement <= ARRIVED_AT. Every criterion is then
read at the first budget where that happens, and what matters is not its value
but how much that value moves between cells. A scale-free criterion reads the
same at the same statistical distance whatever T and whatever the units; the
spread of its readings across the cells is the whole result.

Cells vary T and multiply the simulated series by a constant. That constant is
an exact reparameterization of the same model, so every cell is the same
statistical problem seen in different units, and any criterion whose reading
changes across a row is measuring the units.

Writes out/qvarma_convergence_test.txt.
*/
#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./stats.h>
#include <et_al./random/random.h>
#include <et_al./frame/csv.h>
#include <cblas.h>
#include <string.h>
#include <malloc.h>
#include <dirent.h>

#define K 5
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define PHI_STAR_BOUND ((mreal)1)
#define P 1
#define Q 1
#define SPEC_R 2

#define REFERENCE_BUDGET 50000
#define ARRIVED_AT 0.01
#define N_REPLICATES 8

static const int budget[] = { 25, 50, 100, 200, 400, 800, 1600, 3200, 6400 };
#define N_BUDGETS ((int)(sizeof budget / sizeof budget[0]))

static const int period_count[] = { 200, 400, 1600 };
#define N_PERIODS ((int)(sizeof period_count / sizeof period_count[0]))

static const double series_scale[] = { 0.01, 1.0, 100.0 };
#define N_SCALES ((int)(sizeof series_scale / sizeof series_scale[0]))

#define N_CRITERIA 5
static const char *criterion_name[N_CRITERIA] = {
    "worst_gradient", "relative_gradient", "newton_decrement", "per_observation",
    "gradient_in_se"
};

/* The model the data comes from. nu is moderate so the Student-t tail is
   identified and the maximum it is fitted to is a real one, which is what
   makes a distance from that maximum meaningful. */
static QvarmaParams truth(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    m.phi_star_bound = PHI_STAR_BOUND;
    for (int k = 0; k < K; k++) AT(m.c, k, 0) = (mreal)(0.2 * (k + 1));
    AT(m.Phi_star, 0, 0) = (mreal)0.55;
    for (int a = 0; a < K_STAR; a++)
        for (int b = 0; b < K; b++)
            AT(m.Psi_star[0], a, b) = (mreal)(a == b ? 0.18 : 0.03);
    for (int a = 0; a < K; a++)
        for (int b = 0; b <= a; b++)
            AT(m.Omega_inv, a, b) = (mreal)(b == a ? 0.5 + 0.1 * a : 0.05);
    m.nu = (mreal)8;
    for (int l = 0; l < SPEC_R; l++)
        for (int i = 0; i < K - K_STAR; i++)
            AT(m.alpha[l], i, 0) = (mreal)(pow(0.5, l) * (0.12 - 0.05 * i));
    AT(m.beta[0], 0, 0) = 1;
    AT(m.beta[0], 0, 1) = (mreal)0.4;
    Vec theta = mat_new(qvarma_n_theta(&m), 1);
    _qvarma_unlink(&m, theta);
    qvarma_params_from_theta(theta, &m);
    mat_free(theta);
    return m;
}

static mreal first_difference_sd(Mat y, int row) {
    int periods = y.c;
    Mat difference = mat_new(1, periods - 1);
    for (int t = 1; t < periods; t++) AT(difference, 0, t - 1) = AT(y, row, t) - AT(y, row, t - 1);
    mreal sd = (mreal)sqrt((double)stats_var(difference));
    mat_free(difference);
    return sd;
}

/* The same starting rule applications/abm_system_fit_qvarma.c uses, so the
   path measured here is the path that pipeline walks. */
static QvarmaParams build_start(Mat y) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
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
    m.nu = (mreal)30;
    for (int l = 0; l < SPEC_R; l++) {
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
The Hessian with a step proportional to each coordinate rather than the fixed
1e-4 _qvarma_hessian uses. A fixed absolute step is not scale free: multiplying
the series by 100 multiplies theta_c by 100, so the same 1e-4 becomes a
different perturbation of the same statistical quantity, and any criterion read
off the result inherits that. Scoring a criterion for scale-freeness against a
curvature that is not scale free would measure the instrument.
*/
static Mat relative_step_hessian(Vec theta, const QvarmaParams *shape, Mat y) {
    int n = theta.r;
    QvarmaFitContext context = { y, shape, qvarma_analytic_new(shape, y.c) };
    Mat H = mat_new(n, n);
    Vec forward = mat_new(n, 1), backward = mat_new(n, 1), probe = mat_new(n, 1);
    for (int j = 0; j < n; j++) {
        double magnitude = fabs((double)theta.d[j]);
        mreal step = (mreal)(1e-5 * (magnitude > 1 ? magnitude : 1));
        for (int i = 0; i < n; i++) probe.d[i] = theta.d[i];
        probe.d[j] += step;
        qvarma_negative_log_likelihood(probe, forward, &context);
        probe.d[j] -= 2 * step;
        qvarma_negative_log_likelihood(probe, backward, &context);
        for (int i = 0; i < n; i++) AT(H, i, j) = (forward.d[i] - backward.d[i]) / (2 * step);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) {
            mreal mean = (mreal)0.5 * (AT(H, i, j) + AT(H, j, i));
            AT(H, i, j) = AT(H, j, i) = mean;
        }
    qvarma_analytic_free(context.workspace);
    mat_free(forward); mat_free(backward); mat_free(probe);
    return H;
}

/* Cholesky of H if it is positive definite, NULL otherwise. mat_chol asserts
   on an indefinite matrix, and an unconverged point is exactly where one turns
   up, so _potrf's own return code is used instead. */
static int chol_if_positive_definite(Mat H, Mat *out) {
    Mat l = mat_copy(H);
    int info = _potrf(l.d, H.r, l.stride);
    if (info != 0) { mat_free(l); return 0; }
    for (int i = 0; i < l.r; i++)
        for (int j = i + 1; j < l.c; j++) AT(l, i, j) = 0;
    *out = l;
    return 1;
}

static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

static int compare_double(const void *a, const void *b) {
    double x = *(const double*)a, z = *(const double*)b;
    return x < z ? -1 : x > z ? 1 : 0;
}

typedef struct {
    int usable;
    int T;
    double scale;
    int replicate;
    int arrived_budget;              /* first budget with se_movement <= ARRIVED_AT */
    double at_arrival[N_CRITERIA];   /* each criterion read there */
    double before_arrival[N_CRITERIA]; /* and at the budget before it */
    double reference_worst_gradient;
    int reference_converged;
    double se_movement[N_BUDGETS];
    double reading[N_BUDGETS][N_CRITERIA];
    int computable[N_BUDGETS];       /* the Hessian at this point was positive definite */
} Cell;

/* Every budget of every usable cell, flattened, which is what the threshold
   questions below are asked of. */
#define FAR_AWAY 1.0

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);

    int n_cells = N_PERIODS * N_SCALES * N_REPLICATES;
    Cell *cell = malloc((size_t)n_cells * sizeof(Cell));

    #pragma omp parallel for schedule(dynamic)
    for (int index = 0; index < n_cells; index++) {
        int replicate = index % N_REPLICATES;
        int scale_at = (index / N_REPLICATES) % N_SCALES;
        int period_at = index / (N_REPLICATES * N_SCALES);
        int T = period_count[period_at];
        double scale = series_scale[scale_at];

        Cell out;
        memset(&out, 0, sizeof out);
        out.T = T; out.scale = scale; out.replicate = replicate;
        out.arrived_budget = -1;

        /* One draw per (T, replicate), the same series in every scale column,
           so a row differs only in the units. */
        QvarmaParams generating = truth();
        Rng rng = rng_new((uint64_t)(1000 + period_at * 100 + replicate), 7);
        Mat y = qvarma_simulate(&rng, &generating, T);
        for (int i = 0; i < y.r * y.c; i++) y.d[i] = (mreal)(scale * (double)y.d[i]);
        qvarma_params_free(&generating);

        QvarmaParams start = build_start(y);
        int n = qvarma_n_theta(&start);

        QvarmaFitOptions options = qvarma_default_fit_options();
        options.max_iterations = REFERENCE_BUDGET;
        QvarmaFitResult reference = qvarma_fit(y, &start, options);
        out.reference_worst_gradient = (double)reference.gradient_norm;
        out.reference_converged = reference.is_converged;

        Vec theta_star = mat_new(n, 1);
        _qvarma_unlink(&reference.params, theta_star);
        Mat H_star = relative_step_hessian(theta_star, &reference.params, y);

        /* Standard errors from the reference curvature. Without a positive
           definite one there is no yardstick and the cell is dropped rather
           than scored against a bad one. */
        Mat chol_star;
        Vec se = mat_new(n, 1);
        int have_yardstick = chol_if_positive_definite(H_star, &chol_star);
        if (have_yardstick) {
            Vec unit = mat_new(n, 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) unit.d[j] = (mreal)(j == i ? 1 : 0);
                Vec column = vec_chol_solve(chol_star, unit);
                se.d[i] = (mreal)sqrt(fabs((double)column.d[i]));
                mat_free(column);
            }
            mat_free(unit);
            mat_free(chol_star);
            out.usable = 1;
        }

        if (out.usable) {
            for (int b = 0; b < N_BUDGETS; b++) {
                QvarmaParams b_start = build_start(y);
                QvarmaFitOptions b_options = qvarma_default_fit_options();
                b_options.max_iterations = budget[b];
                QvarmaFitResult fit = qvarma_fit(y, &b_start, b_options);

                Vec theta = mat_new(n, 1);
                _qvarma_unlink(&fit.params, theta);

                double se_movement = 0;
                for (int i = 0; i < n; i++) {
                    double move = fabs((double)(theta.d[i] - theta_star.d[i])) / (double)se.d[i];
                    if (move > se_movement) se_movement = move;
                }

                /* The criteria, all read at this same point. */
                QvarmaFitContext context = { y, &fit.params, qvarma_analytic_new(&fit.params, T) };
                Vec gradient = mat_new(n, 1);
                double objective = (double)qvarma_negative_log_likelihood(theta, gradient, &context);
                qvarma_analytic_free(context.workspace);

                double worst = 0, relative = 0;
                for (int i = 0; i < n; i++) {
                    double g = fabs((double)gradient.d[i]);
                    if (g > worst) worst = g;
                    double scaled_theta = fabs((double)theta.d[i]);
                    if (scaled_theta < 1) scaled_theta = 1;
                    double denominator = fabs(objective) < 1 ? 1 : fabs(objective);
                    double value = g * scaled_theta / denominator;
                    if (value > relative) relative = value;
                }

                Mat H = relative_step_hessian(theta, &fit.params, y);
                Mat chol;
                double decrement = INFINITY, in_se = INFINITY;
                if (chol_if_positive_definite(H, &chol)) {
                    Vec step = vec_chol_solve(chol, gradient);
                    double d = 0;
                    for (int i = 0; i < n; i++) d += (double)gradient.d[i] * (double)step.d[i];
                    decrement = fabs(d);
                    mat_free(step);
                    /* Each coordinate's gradient in units of that coordinate's own
                       standard error, which is the diagonal of the same inverse the
                       decrement uses the whole of. */
                    Vec unit = mat_new(n, 1);
                    in_se = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) unit.d[j] = (mreal)(j == i ? 1 : 0);
                        Vec column = vec_chol_solve(chol, unit);
                        double own_se = sqrt(fabs((double)column.d[i]));
                        double value = fabs((double)gradient.d[i]) * own_se;
                        if (value > in_se) in_se = value;
                        mat_free(column);
                    }
                    mat_free(unit);
                    mat_free(chol);
                }
                mat_free(H);

                double reading[N_CRITERIA] = { worst, relative, decrement, worst / T, in_se };
                out.se_movement[b] = se_movement;
                out.computable[b] = isfinite(decrement);
                for (int c = 0; c < N_CRITERIA; c++) out.reading[b][c] = reading[c];
                if (out.arrived_budget < 0 && se_movement <= ARRIVED_AT) {
                    out.arrived_budget = budget[b];
                    for (int c = 0; c < N_CRITERIA; c++) out.at_arrival[c] = reading[c];
                } else if (out.arrived_budget < 0) {
                    for (int c = 0; c < N_CRITERIA; c++) out.before_arrival[c] = reading[c];
                }

                mat_free(gradient);
                mat_free(theta);
                qvarma_fit_result_free(&fit);
                qvarma_params_free(&b_start);
            }
        }

        mat_free(H_star);
        mat_free(se);
        mat_free(theta_star);
        qvarma_fit_result_free(&reference);
        qvarma_params_free(&start);
        mat_free(y);
        cell[index] = out;
    }

    FILE *report = fopen("out/qvarma_convergence_test.txt", "w");
    assert(report);
    fprintf(report, "t-QVARMA(1,1,2), K = %d, K_star = %d, R = %d, simulated from a fixed truth "
                    "with nu = 8.\n", K, K_STAR, R);
    fprintf(report, "Reference point: %d iterations from build_start. Standard errors from the "
                    "Hessian there.\n", REFERENCE_BUDGET);
    fprintf(report, "A fit has arrived when max_i |theta_i - theta_star_i| / se_i <= %.3g.\n",
            ARRIVED_AT);
    fprintf(report, "Budgets: ");
    for (int b = 0; b < N_BUDGETS; b++) fprintf(report, "%d%s", budget[b],
                                                b + 1 < N_BUDGETS ? ", " : "\n");
    fprintf(report, "%d cells: T in {200, 400, 1600} by series scale in {1, 100} by %d "
                    "replicates.\n\n", n_cells, N_REPLICATES);

    int usable = 0, arrived = 0;
    for (int i = 0; i < n_cells; i++) {
        if (cell[i].usable) usable++;
        if (cell[i].usable && cell[i].arrived_budget > 0) arrived++;
    }
    fprintf(report, "cells with a positive definite reference Hessian: %d of %d\n", usable, n_cells);
    fprintf(report, "of those, cells that arrived within the largest budget: %d\n\n", arrived);

    fprintf(report, "per cell, the budget at which it arrived and each criterion read there\n");
    fprintf(report, "  %6s %8s %4s %10s %16s %18s %18s %16s\n", "T", "scale", "rep", "arrived_at",
            "worst_gradient", "relative_gradient", "newton_decrement", "per_observation");
    for (int i = 0; i < n_cells; i++) {
        if (!cell[i].usable || cell[i].arrived_budget < 0) continue;
        fprintf(report, "  %6d %8.0f %4d %10d %16.4g %18.4g %18.4g %16.4g\n", cell[i].T,
                cell[i].scale, cell[i].replicate, cell[i].arrived_budget,
                cell[i].at_arrival[0], cell[i].at_arrival[1], cell[i].at_arrival[2],
                cell[i].at_arrival[3]);
    }

    fprintf(report, "\nspread of each criterion at the moment of arrival, over the %d cells\n",
            arrived);
    fprintf(report, "  a criterion that does not depend on T or on the units reads the same "
                    "number in every cell,\n  so the ratio of its largest reading to its "
                    "smallest is what is being compared here.\n\n");
    fprintf(report, "  %-20s %14s %14s %14s %14s\n", "criterion", "min", "median", "max",
            "max/min");
    double *value = malloc((size_t)n_cells * N_BUDGETS * sizeof(double));
    for (int c = 0; c < N_CRITERIA; c++) {
        int m = 0;
        for (int i = 0; i < n_cells; i++)
            if (cell[i].usable && cell[i].arrived_budget > 0) value[m++] = cell[i].at_arrival[c];
        qsort(value, (size_t)m, sizeof(double), compare_double);
        double lo = value[0], hi = value[m - 1], mid = value[m / 2];
        fprintf(report, "  %-20s %14.4g %14.4g %14.4g %14.4g\n", criterion_name[c], lo, mid, hi,
                lo > 0 ? hi / lo : INFINITY);
    }

    /* A threshold is only useful if one number separates arrived from not
       arrived everywhere. Asking it to separate se_movement 0.009 from 0.011
       is a distinction nothing depends on, so the question is put with a
       margin: never declare convergence while a full standard error of
       movement is still to come, and always declare it once a hundredth of one
       is. Whether a threshold fits between those is what decides a criterion. */
    fprintf(report, "\nis there a single threshold that works in every cell?\n");
    fprintf(report, "  arrived     = se_movement <= %.3g, the threshold must accept all of these\n",
            ARRIVED_AT);
    fprintf(report, "  not arrived = se_movement >  %.3g, the threshold must reject all of these\n",
            FAR_AWAY);
    fprintf(report, "  points between the two are not counted either way.\n\n");
    fprintf(report, "  %-20s %10s %10s %16s %18s %12s\n", "criterion", "arrived",
            "not arrived", "max on arrived", "min on not arrived", "separable");
    for (int c = 0; c < N_CRITERIA; c++) {
        double worst_arrived = 0, best_far = INFINITY;
        int n_arrived = 0, n_far = 0;
        for (int i = 0; i < n_cells; i++) {
            if (!cell[i].usable) continue;
            for (int b = 0; b < N_BUDGETS; b++) {
                if (!cell[i].computable[b] && c != 0 && c != 1 && c != 3) continue;
                double value = cell[i].reading[b][c];
                if (!isfinite(value)) continue;
                if (cell[i].se_movement[b] <= ARRIVED_AT) {
                    n_arrived++;
                    if (value > worst_arrived) worst_arrived = value;
                } else if (cell[i].se_movement[b] > FAR_AWAY) {
                    n_far++;
                    if (value < best_far) best_far = value;
                }
            }
        }
        fprintf(report, "  %-20s %10d %10d %16.4g %18.4g %12s\n", criterion_name[c],
                n_arrived, n_far, worst_arrived, best_far,
                worst_arrived < best_far ? "yes" : "no");
    }

    /* Failing separability, how wrong is the best threshold. Swept over every
       reading any criterion took, keeping the one with the fewest points on
       the wrong side. */
    fprintf(report, "\nthe best single threshold each criterion can manage\n");
    fprintf(report, "  %-20s %12s %14s %16s %18s\n", "criterion", "threshold",
            "errors", "missed arrivals", "false convergences");
    for (int c = 0; c < N_CRITERIA; c++) {
        int m = 0;
        for (int i = 0; i < n_cells; i++) {
            if (!cell[i].usable) continue;
            for (int b = 0; b < N_BUDGETS; b++)
                if (isfinite(cell[i].reading[b][c])) value[m++] = cell[i].reading[b][c];
        }
        qsort(value, (size_t)m, sizeof(double), compare_double);
        double best_threshold = 0;
        int fewest = 1 << 30, best_missed = 0, best_false = 0;
        for (int t = 0; t < m; t++) {
            double threshold = value[t];
            int missed = 0, wrong = 0;
            for (int i = 0; i < n_cells; i++) {
                if (!cell[i].usable) continue;
                for (int b = 0; b < N_BUDGETS; b++) {
                    double reading_here = cell[i].reading[b][c];
                    if (!isfinite(reading_here)) continue;
                    int declares = reading_here <= threshold;
                    if (cell[i].se_movement[b] <= ARRIVED_AT && !declares) missed++;
                    else if (cell[i].se_movement[b] > FAR_AWAY && declares) wrong++;
                }
            }
            if (missed + wrong < fewest) {
                fewest = missed + wrong; best_threshold = threshold;
                best_missed = missed; best_false = wrong;
            }
        }
        fprintf(report, "  %-20s %12.4g %14d %16d %18d\n", criterion_name[c], best_threshold,
                fewest, best_missed, best_false);
    }

    /* What the criteria needing a Hessian cost in availability. */
    {
        int total_rows = 0, computable_rows = 0;
        for (int i = 0; i < n_cells; i++) {
            if (!cell[i].usable) continue;
            for (int b = 0; b < N_BUDGETS; b++) { total_rows++; if (cell[i].computable[b]) computable_rows++; }
        }
        fprintf(report, "\npoints where the Hessian was positive definite, so newton_decrement "
                        "and gradient_in_se\nexist at all: %d of %d (%.1f%%)\n",
                computable_rows, total_rows, 100.0 * computable_rows / total_rows);
    }

    fprintf(report, "\nreference fits, for whether the yardstick itself is sound\n");
    fprintf(report, "  %6s %8s %4s %16s %12s %10s\n", "T", "scale", "rep", "worst_gradient",
            "converged", "usable");
    for (int i = 0; i < n_cells; i++)
        fprintf(report, "  %6d %8.0f %4d %16.4g %12s %10s\n", cell[i].T, cell[i].scale,
                cell[i].replicate, cell[i].reference_worst_gradient,
                cell[i].reference_converged ? "yes" : "no", cell[i].usable ? "yes" : "no");

    /*
    The same question put to the pipeline's own fits, where there is no
    reference point to score against but the test can still be applied and
    counted. What it reports there is the practical consequence of adopting it.
    */
    {
        fprintf(report, "\n\nthe proposed test applied to the pipeline's own cached fits\n");
        fprintf(report, "  out/abm_system_fit_qvarma, spec p1q1r2, every replicate.\n");
        fprintf(report, "  No reference point exists for real data, so these are counts, "
                        "not accuracy.\n\n");

        DIR *handle = opendir("dataset/abm_system");
        if (!handle) {
            fprintf(report, "  dataset/abm_system not present; skipped.\n");
        } else {
            char **name = NULL;
            int n_sample = 0, cap = 0;
            struct dirent *entry;
            while ((entry = readdir(handle)) != NULL) {
                if (entry->d_name[0] == '.') continue;
                if (n_sample == cap) { cap = cap ? cap * 2 : 16;
                                       name = realloc(name, (size_t)cap * sizeof(char*)); }
                size_t len = strlen(entry->d_name);
                name[n_sample] = malloc(len + 1);
                memcpy(name[n_sample], entry->d_name, len + 1);
                n_sample++;
            }
            closedir(handle);

            int loaded = 0, flagged_now = 0, positive_definite = 0, passes_proposed = 0;
            #pragma omp parallel for schedule(dynamic) reduction(+:loaded,flagged_now,positive_definite,passes_proposed)
            for (int s_at = 0; s_at < n_sample; s_at++) {
                for (int rep = 0; rep < 200; rep++) {
                    char csv_path[560], cache_path[560];
                    snprintf(csv_path, sizeof csv_path,
                             "dataset/abm_system/%s/replicate_%03d.csv", name[s_at], rep);
                    FILE *probe = fopen(csv_path, "r");
                    if (!probe) break;
                    fclose(probe);
                    snprintf(cache_path, sizeof cache_path,
                             "out/abm_system_fit_qvarma/%s/replicate_%03d_p1q1r2_fit.json",
                             name[s_at], rep);

                    DataFrame df = df_read_csv(csv_path, csv_read_options_default());
                    Mat y = mat_new(K, df.r);
                    static const char *row_name[K] = { "GDP_growth", "EN_growth",
                        "Employment_change", "Inflation", "InterestRate" };
                    for (int k = 0; k < K; k++) {
                        Mat column = df_col_numeric(&df, row_name[k]);
                        for (int t = 0; t < df.r; t++) AT(y, k, t) = AT(column, t, 0);
                    }
                    df_free(&df);

                    QvarmaFitResult cached;
                    cached.params = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA,
                                                      WARMUP_LONGEST);
                    cached.params.phi_star_bound = PHI_STAR_BOUND;
                    cached.params.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
                    if (!qvarma_load_fit(&cached, y, cache_path)) {
                        qvarma_params_free(&cached.params);
                        mat_free(y);
                        continue;
                    }
                    loaded++;
                    if (cached.is_converged) flagged_now++;

                    int n = qvarma_n_theta(&cached.params);
                    Vec theta = mat_new(n, 1);
                    _qvarma_unlink(&cached.params, theta);
                    Vec gradient = mat_new(n, 1);
                    QvarmaFitContext context = { y, &cached.params,
                                                 qvarma_analytic_new(&cached.params, y.c) };
                    qvarma_negative_log_likelihood(theta, gradient, &context);
                    qvarma_analytic_free(context.workspace);

                    Mat H = relative_step_hessian(theta, &cached.params, y);
                    Mat chol;
                    if (chol_if_positive_definite(H, &chol)) {
                        positive_definite++;
                        Vec step = vec_chol_solve(chol, gradient);
                        double d = 0;
                        for (int i = 0; i < n; i++)
                            d += (double)gradient.d[i] * (double)step.d[i];
                        if (fabs(d) <= 0.0306) passes_proposed++;
                        mat_free(step);
                        mat_free(chol);
                    }
                    mat_free(H);
                    mat_free(gradient);
                    mat_free(theta);
                    qvarma_fit_result_free(&cached);
                    mat_free(y);
                }
            }
            fprintf(report, "  cached fits read                                %6d\n", loaded);
            fprintf(report, "  called converged by the test in use today       %6d  (%.1f%%)\n",
                    flagged_now, 100.0 * flagged_now / loaded);
            fprintf(report, "  Hessian positive definite, so at a maximum      %6d  (%.1f%%)\n",
                    positive_definite, 100.0 * positive_definite / loaded);
            fprintf(report, "  newton_decrement <= 0.0306 as well              %6d  (%.1f%%)\n",
                    passes_proposed, 100.0 * passes_proposed / loaded);
            for (int i = 0; i < n_sample; i++) free(name[i]);
            free(name);
        }
    }

    free(value);
    fclose(report);
    free(cell);
    return 0;
}
