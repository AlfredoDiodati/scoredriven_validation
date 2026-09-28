/*
Which named law fits the marginal of a simulated series, by maximum likelihood.

The sections before this one leave a marginal that is non-Gaussian, has a
finite variance, an excess kurtosis near one, a skewness near zero for the
growth series and no stable behaviour. Several named laws have those
properties. This fits all of them to the same data and reports which one the
likelihood prefers.

The laws are studies/marginal_laws.h, which also says why they are written
there rather than called from et_al: et_al has the Gaussian and the Student t
in their Mat forms and none of the other six, and no Bessel function for the
three that need one. All eight come through the same fitting routine on et_al's
solver/lbfgs.h, so a comparison across them is a comparison of laws and not of
optimisers.

    Gaussian                  2 parameters, the null this is measured against
    Laplace                   2, excess kurtosis 3 whatever the parameters
    Subbotin                  3, exp(-|z|^b), b = 2 Gaussian and b = 1 Laplace
    Student t                 3, tail index nu
    lognormal mixture         3, a normal whose standard deviation is lognormal
    normal inverse Gaussian   4, semi-heavy tails
    variance gamma            4, a normal whose variance is a gamma draw
    generalised hyperbolic    5, which contains the two above it

One fit per series per replicate, on that replicate's own 400 quarters. Nothing
is pooled across replicates: the system is not ergodic, so a sample taken
across them at a fixed date is a mixture over economies in different states.

The comparison is by the Schwarz criterion, p log n - 2 log L, which is what
charges the five-parameter laws for their parameters on a sample of 400. Akaike
is reported beside it because the two disagree about how much to charge, and
the disagreement is the point at which a reader should look at both.

The same eight are fitted to the US series, 187 quarters from out/us_system.csv,
so the simulated answer has something to be read against.

A fit that did not converge is counted and excluded from the comparison rather
than being allowed to win by having stopped somewhere flat.

    ./bin/abm_system_marginal_fit [CONFIGURATIONS [REPLICATES]]

Default 100 and 2, which is 200 series per variable and about an hour. The
sample is a sample because of what the fits cost, not because more is not
available: the three laws with a Bessel function in them evaluate a
sixty-four-node quadrature per observation per likelihood, and a likelihood is
called a few thousand times per fit by a line search that has no analytic
gradient to shorten it. Measured at 4.4 minutes for 96 series on 16 threads, so
the whole million replicates would be years.

Reads dataset/abm_system and out/us_system.csv. Writes, none of it printed:
    out/abm_system_marginal_fit_report.txt
    out/abm_system_marginal_fit_by_law.csv
    out/abm_system_marginal_fit_us.csv
*/

#include "applications/abm_system.h"
#include "applications/us_data.h"
#include "studies/marginal_laws.h"
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define K ABM_SYSTEM_K
#define INPUT_DIR "dataset/abm_system"
#define REPORT_PATH "out/abm_system_marginal_fit_report.txt"
#define BY_LAW_PATH "out/abm_system_marginal_fit_by_law.csv"
#define US_PATH "out/abm_system_marginal_fit_us.csv"

#define PERIODS 400
#define US_PERIODS (ESTIMATION_PERIODS - 1)
#define DEFAULT_CONFIGURATIONS 100
#define DEFAULT_REPLICATES 2

#define N_SERIES (K + 1)
static const char *series_name[N_SERIES] = {
    "GDP growth", "energy growth", "employment change", "inflation", "interest rate",
    "interest rate change"
};
static const char *series_column[N_SERIES] = {
    "gdp_growth", "energy_growth", "employment_change", "inflation", "interest_rate",
    "interest_rate_change"
};

/* Running totals over the fits of one law on one series. */
typedef struct {
    long fits;
    long converged;
    long best_by_bic;
    long best_by_aic;
    double log_likelihood;
    double gain_over_gaussian;
    double shape;            /* the law's own shape parameter, or NAN */
} Tally;

/* The parameter a law's shape lives in, or -1 when it has none. */
static int shape_index(int law) {
    switch (law) {
        case 2: return 2;   /* Subbotin b */
        case 3: return 2;   /* Student t nu */
        case 4: return 2;   /* lognormal mixture spread */
        case 5: return 2;   /* normal inverse Gaussian alpha */
        case 6: return 2;   /* variance gamma kappa */
        case 7: return 4;   /* generalised hyperbolic lambda */
        default: return -1;
    }
}

static void tally_merge(Tally *into, const Tally *from) {
    into->fits += from->fits;
    into->converged += from->converged;
    into->best_by_bic += from->best_by_bic;
    into->best_by_aic += from->best_by_aic;
    into->log_likelihood += from->log_likelihood;
    into->gain_over_gaussian += from->gain_over_gaussian;
    into->shape += from->shape;
}

/* Fits every law to one sample and adds the result to a tally. */
static void fit_all(const double *x, int n, Tally *tally, int stride) {
    MarginalFit fit[MARGINAL_LAWS_COUNT];
    const MarginalLaw *law = marginal_laws_all();

    for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) fit[l] = marginal_laws_fit(&law[l], x, n);

    int best_bic = -1, best_aic = -1;
    for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) {
        if (!fit[l].is_converged) continue;
        if (best_bic < 0 || fit[l].bic < fit[best_bic].bic) best_bic = l;
        if (best_aic < 0 || fit[l].aic < fit[best_aic].aic) best_aic = l;
    }

    for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) {
        Tally *here = &tally[l * stride];
        here->fits++;
        if (!fit[l].is_converged) continue;
        here->converged++;
        here->log_likelihood += fit[l].log_likelihood;
        if (fit[0].is_converged)
            here->gain_over_gaussian += fit[l].log_likelihood - fit[0].log_likelihood;
        int s = shape_index(l);
        if (s >= 0) here->shape += fit[l].theta[s];
        if (l == best_bic) here->best_by_bic++;
        if (l == best_aic) here->best_by_aic++;
    }
}

static int *list_configurations(int *count) {
    DIR *handle = opendir(INPUT_DIR);
    assert(handle && "abm_system_marginal_fit: cannot open dataset/abm_system");

    int *cop = NULL, n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int here;
        if (sscanf(entry->d_name, "cop_%d", &here) != 1) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 1024;
            cop = (int*)realloc(cop, (size_t)cap * sizeof(int));
            assert(cop && "abm_system_marginal_fit: out of memory");
        }
        cop[n++] = here;
    }
    closedir(handle);
    for (int i = 1; i < n; i++) {
        int current = cop[i], j = i;
        while (j > 0 && cop[j - 1] > current) { cop[j] = cop[j - 1]; j--; }
        cop[j] = current;
    }
    *count = n;
    return cop;
}

/* The five US series, built the way applications/us_qvarma_spec_choice.c's own
   block does, which is the route dataset/abm_system's own series mirror. */
static Mat build_us_block(void) {
    Mat original = load_us_system();
    Mat y = mat_new(K, US_PERIODS);
    for (int t = 1; t < ESTIMATION_PERIODS; t++) {
        int c = t - 1;
        AT(y, ROW_GDP_GROWTH, c) = AT(original, LOG_GDP, t) - AT(original, LOG_GDP, t - 1);
        AT(y, ROW_EN_GROWTH, c) = AT(original, LOG_ENERGY_DEMAND, t) - AT(original, LOG_ENERGY_DEMAND, t - 1);
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = AT(original, EMPLOYMENT, t) - AT(original, EMPLOYMENT, t - 1);
        AT(y, ROW_INFLATION, c) = AT(original, LOG_CPI, t) - AT(original, LOG_CPI, t - 1);
        AT(y, ROW_INTEREST_RATE, c) = AT(original, INTEREST_RATE, t);
    }
    mat_free(original);
    return y;
}

int main(int argc, char **argv) {
    int wanted_configurations = argc > 1 ? atoi(argv[1]) : DEFAULT_CONFIGURATIONS;
    int replicates = argc > 2 ? atoi(argv[2]) : DEFAULT_REPLICATES;
    assert(wanted_configurations >= 1 && replicates >= 1 &&
           "abm_system_marginal_fit: nothing to fit");

    int stored;
    int *cop = list_configurations(&stored);
    assert(stored > 0 && "abm_system_marginal_fit: no configurations stored");
    int configurations = wanted_configurations < stored ? wanted_configurations : stored;
    int stride = stored / configurations > 0 ? stored / configurations : 1;

    const MarginalLaw *law = marginal_laws_all();
    Tally *tally = (Tally*)calloc((size_t)N_SERIES * MARGINAL_LAWS_COUNT, sizeof(Tally));
    assert(tally && "abm_system_marginal_fit: out of memory");

    const double started = omp_get_wtime();
    long series_fitted = 0;

    #pragma omp parallel
    {
        Tally *mine = (Tally*)calloc((size_t)N_SERIES * MARGINAL_LAWS_COUNT, sizeof(Tally));
        double *sample = (double*)malloc((size_t)PERIODS * sizeof(double));
        assert(mine && sample && "abm_system_marginal_fit: out of memory");
        long my_series = 0;

        #pragma omp for schedule(dynamic)
        for (int c = 0; c < configurations; c++) {
            char dir[512];
            snprintf(dir, sizeof dir, "%s/cop_%04d", INPUT_DIR, cop[c * stride]);

            int batches = (replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
            int taken = 0;
            for (int batch = 0; batch < batches && taken < replicates; batch++) {
                Mat block[ABM_SYSTEM_BATCH];
                int replicate[ABM_SYSTEM_BATCH];
                int count = abm_system_read_batch(dir, batch, block, replicate);

                for (int b = 0; b < count; b++) {
                    Mat y = block[b];
                    if (taken >= replicates || y.c < PERIODS) { mat_free(y); continue; }
                    taken++;

                    for (int s = 0; s < N_SERIES; s++) {
                        int n;
                        if (s < K) {
                            n = PERIODS;
                            for (int t = 0; t < PERIODS; t++) sample[t] = (double)AT(y, s, t);
                        } else {
                            n = PERIODS - 1;
                            for (int t = 1; t < PERIODS; t++)
                                sample[t - 1] = (double)AT(y, ROW_INTEREST_RATE, t) -
                                                (double)AT(y, ROW_INTEREST_RATE, t - 1);
                        }
                        fit_all(sample, n, &mine[s * MARGINAL_LAWS_COUNT], 1);
                        my_series++;
                    }
                    mat_free(y);
                }
            }
        }

        #pragma omp critical(merge)
        {
            for (int i = 0; i < N_SERIES * MARGINAL_LAWS_COUNT; i++)
                tally_merge(&tally[i], &mine[i]);
            series_fitted += my_series;
        }
        free(mine); free(sample);
    }

    const double elapsed = omp_get_wtime() - started;

    /* The US series, the same eight laws, one fit each. */
    Mat us = build_us_block();
    MarginalFit us_fit[K][MARGINAL_LAWS_COUNT];
    int us_best_bic[K];
    {
        double *sample = (double*)malloc((size_t)US_PERIODS * sizeof(double));
        assert(sample && "abm_system_marginal_fit: out of memory");
        for (int s = 0; s < K; s++) {
            for (int t = 0; t < US_PERIODS; t++) sample[t] = (double)AT(us, s, t);
            int best = -1;
            for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) {
                us_fit[s][l] = marginal_laws_fit(&law[l], sample, US_PERIODS);
                if (!us_fit[s][l].is_converged) continue;
                if (best < 0 || us_fit[s][l].bic < us_fit[s][best].bic) best = l;
            }
            us_best_bic[s] = best;
        }
        free(sample);
    }
    mat_free(us);

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_marginal_fit: cannot open the report path");
    fprintf(report, "Which named law fits the marginal of a simulated series, by maximum\n"
                    "likelihood.\n\n");
    fprintf(report, "%d configurations of the design, %d replicates of each, %ld series fitted,\n"
                    "400 quarters each, %d laws, %.1f minutes on %d threads.\n",
            configurations, replicates, series_fitted, MARGINAL_LAWS_COUNT,
            elapsed / 60.0, omp_get_max_threads());
    fprintf(report, "Nothing is pooled across replicates. 'wins' is the share of series where\n"
                    "that law had the lowest Schwarz criterion, over the series where it\n"
                    "converged. 'gain' is its mean log-likelihood less the Gaussian's on the\n"
                    "same series, so the Gaussian's own gain is zero by construction.\n");

    FILE *table = fopen(BY_LAW_PATH, "w");
    assert(table && "abm_system_marginal_fit: cannot open the table path");
    fprintf(table, "series,law,parameters,fits,converged,share_converged,mean_log_likelihood,"
                   "mean_gain_over_gaussian,share_best_bic,share_best_aic,mean_shape\n");

    for (int s = 0; s < N_SERIES; s++) {
        fprintf(report, "\n%s\n", series_name[s]);
        fprintf(report, "  law                        p     fits   converged      gain   wins BIC"
                        "   wins AIC      shape\n");
        for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) {
            const Tally *t = &tally[s * MARGINAL_LAWS_COUNT + l];
            if (!t->fits) continue;
            double converged = t->converged ? 1.0 / (double)t->converged : NAN;
            int index = shape_index(l);
            fprintf(report, "  %-24s %3d %8ld %11.4f %9.3f %10.4f %10.4f %10.4g\n",
                    law[l].name, law[l].parameters, t->fits,
                    (double)t->converged / (double)t->fits,
                    t->gain_over_gaussian * converged,
                    (double)t->best_by_bic / (double)t->fits,
                    (double)t->best_by_aic / (double)t->fits,
                    index >= 0 ? t->shape * converged : NAN);
            fprintf(table, "%s,%s,%d,%ld,%ld,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    series_column[s], law[l].name, law[l].parameters, t->fits, t->converged,
                    (double)t->converged / (double)t->fits,
                    t->log_likelihood * converged, t->gain_over_gaussian * converged,
                    (double)t->best_by_bic / (double)t->fits,
                    (double)t->best_by_aic / (double)t->fits,
                    index >= 0 ? t->shape * converged : NAN);
        }
    }
    fclose(table);

    FILE *us_table = fopen(US_PATH, "w");
    assert(us_table && "abm_system_marginal_fit: cannot open the US table path");
    fprintf(us_table, "series,law,parameters,converged,log_likelihood,aic,bic,"
                      "p1,p2,p3,p4,p5\n");
    fprintf(report, "\n\nTHE US SERIES, %d quarters, one fit each. The winner by the Schwarz\n"
                    "criterion is marked.\n", US_PERIODS);
    for (int s = 0; s < K; s++) {
        fprintf(report, "\n%s\n", series_name[s]);
        fprintf(report, "  law                        p    log L        BIC     shape\n");
        for (int l = 0; l < MARGINAL_LAWS_COUNT; l++) {
            const MarginalFit *f = &us_fit[s][l];
            int index = shape_index(l);
            fprintf(report, "  %-24s %3d %8.2f %10.2f %9.4g %s\n",
                    law[l].name, law[l].parameters,
                    f->is_converged ? f->log_likelihood : NAN,
                    f->is_converged ? f->bic : NAN,
                    index >= 0 ? f->theta[index] : NAN,
                    l == us_best_bic[s] ? "  <- lowest BIC" :
                    (f->is_converged ? "" : "  did not converge"));
            fprintf(us_table, "%s,%s,%d,%d,%.9g,%.9g,%.9g",
                    series_column[s], law[l].name, law[l].parameters,
                    f->is_converged, f->log_likelihood, f->aic, f->bic);
            for (int p = 0; p < MARGINAL_LAWS_MAX_PARAMETERS; p++)
                fprintf(us_table, ",%.9g", f->theta[p]);
            fprintf(us_table, "\n");
        }
    }
    fclose(us_table);

    fprintf(report, "\nEvery row is in %s and %s.\n", BY_LAW_PATH, US_PATH);
    fclose(report);

    free(tally);
    free(cop);
    return EXIT_SUCCESS;
}
