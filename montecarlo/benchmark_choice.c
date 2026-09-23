/*
Picks the simulated run that stands in for the real data in montecarlo/, and
records why it was picked.

The Monte Carlo experiment asks whether the validation procedure recovers a
configuration it already knows the answer for. One replicate of one
configuration is promoted to the role the US data plays in the main pipeline,
every loss is measured against it, and the question is whether the confidence
set comes back holding the configuration the benchmark was drawn from. That
only tests the procedure if the benchmark is chosen on grounds that have
nothing to do with the answer, so the grounds are fixed here and written down.

Three requirements, in the order they are applied.

The configuration must be one the real-data confidence sets kept. Those are
read from the comparison tables the three real-data procedures wrote, and the
union of their surviving configurations is the candidate pool. A configuration
no procedure kept would make the experiment ask whether a poor benchmark is
recovered, which is a different and less interesting question.

Among those, the configuration whose auxiliary fits converged most often, read
from out/abm_system_fit_qvarma_manifest.txt. Roughly a third of the million
fits in this project converged, so this is the property in shortest supply and
the one most likely to make a benchmark meaningless: a benchmark whose own fit
did not converge is a parameter vector the optimizer stopped at, not an
estimate.

Within that configuration, the converged replicate with the smallest gradient
norm whose observed information matrix is positive definite. The gradient norm
orders the candidates; positive definiteness is a gate, because the weighted
score loss inverts that matrix and a benchmark it cannot be built at would
leave one of the three protocols unrunnable. Candidates are tried in order and
the first that passes is taken. Replicates that failed the gate are listed in
the report, so a reader can see how close the choice was to falling through.

Fits taken from another replicate's parameters are excluded: the manifest marks
them inherited, and a benchmark should be the estimate its own data produced.

Reads out/abm_system_mcs_statistic_comparison*.csv, the fit manifest, the fit
cache and the dataset. Writes montecarlo/out/benchmark_choice.txt, the report,
and montecarlo/out/benchmark.env, two shell assignments the driver reads so the
choice is made once and every downstream step runs against the same run.
Nothing printed.
*/

#include "applications/abm_system.h"
#include <et_al./sd/qvarma.h>
#include <et_al./linalg/decomp.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define K ABM_SYSTEM_K
#define K_STAR 3
#define R 1
#define SHARED_BETA 1
#define WARMUP_LONGEST 0
#define MU_STAR_STATIONARY_ONLY 1
#define P 1
#define Q 1
#define SPEC_R 2
#define SPEC_LABEL "p1q1r2"

#define FIT_DIR "out/abm_system_fit_qvarma"
#define INPUT_DIR "dataset/abm_system"
#define MANIFEST_PATH "out/abm_system_fit_qvarma_manifest.txt"
#define REPORT_PATH "montecarlo/out/benchmark_choice.txt"
#define ENV_PATH "montecarlo/out/benchmark.env"
#define MAX_CANDIDATES 25

static const char *const winner_table[] = {
    "out/abm_system_mcs_statistic_comparison.csv",
    "out/abm_system_mcs_statistic_comparison_score.csv",
    "out/abm_system_mcs_statistic_comparison_score_weighted.csv"
};
#define N_WINNER_TABLES ((int)(sizeof winner_table / sizeof winner_table[0]))

static QvarmaParams spec_shape(void) {
    QvarmaParams m = qvarma_params_new(K, K_STAR, P, Q, SPEC_R, R, SHARED_BETA, WARMUP_LONGEST);
    m.mu_star_stationary_only = MU_STAR_STATIONARY_ONLY;
    return m;
}

/* A model name in a comparison table is "<sample>_qvarma_<spec>"; the sample is
   what everything else here is keyed on. */
static void sample_of_model(const char *model, char *out, size_t n) {
    const char *marker = strstr(model, "_qvarma_");
    assert(marker && "benchmark_choice: a model name does not carry the _qvarma_ marker");
    size_t len = (size_t)(marker - model);
    assert(len < n && "benchmark_choice: a sample name does not fit");
    memcpy(out, model, len);
    out[len] = 0;
}

typedef struct { char name[64]; } Sample;

/* Every configuration any of the three real-data confidence sets kept, without
   repeats. A table that is not there is skipped rather than fatal, so this can
   run before every procedure has been scored. */
static Sample *read_winners(int *count, FILE *report) {
    Sample *winners = NULL;
    int n = 0, cap = 0;
    for (int t = 0; t < N_WINNER_TABLES; t++) {
        FILE *probe = fopen(winner_table[t], "r");
        if (!probe) {
            fprintf(report, "  %s not found, skipped\n", winner_table[t]);
            continue;
        }
        fclose(probe);

        DataFrame table = df_read_csv(winner_table[t], csv_read_options_default());
        char **model = df_col_string(&table, "model");
        Mat in_set = df_col_numeric(&table, "tr_in_set");
        int kept = 0;
        for (int i = 0; i < table.r; i++) {
            if ((int)AT(in_set, i, 0) != 1) continue;
            kept++;
            char sample[64];
            sample_of_model(model[i], sample, sizeof sample);
            int seen = 0;
            for (int j = 0; j < n; j++) if (strcmp(winners[j].name, sample) == 0) seen = 1;
            if (seen) continue;
            if (n == cap) {
                cap = cap ? cap * 2 : 8;
                Sample *grown = realloc(winners, (size_t)cap * sizeof(Sample));
                assert(grown && "benchmark_choice: out of memory collecting winners");
                winners = grown;
            }
            snprintf(winners[n].name, sizeof winners[n].name, "%s", sample);
            n++;
        }
        fprintf(report, "  %s kept %d\n", winner_table[t], kept);
        df_free(&table);
    }
    *count = n;
    return winners;
}

typedef struct {
    int replicate;
    double log_likelihood;
    double gradient;
    int converged;
    int inherited;
} FitRow;

typedef struct {
    char name[64];
    FitRow *row;
    int n_rows;
    int n_converged;
} SampleFits;

/* The manifest's rows for the configurations named, in one pass over a file of
   a million lines. */
static void read_manifest(SampleFits *fits, int n_fits) {
    FILE *f = fopen(MANIFEST_PATH, "r");
    assert(f && "benchmark_choice: out/abm_system_fit_qvarma_manifest.txt is missing - "
                 "make app-abm_system_fit_qvarma writes it");

    char line[512];
    while (fgets(line, sizeof line, f)) {
        if (strncmp(line, "cop_", 4) != 0) continue;
        char sample[64], spec[32], converged[8], origin[32];
        int replicate;
        double log_likelihood, gradient, aic;
        if (sscanf(line, "%63s %d %31s %lf %lf %lf %7s %31s",
                   sample, &replicate, spec, &log_likelihood, &gradient, &aic,
                   converged, origin) != 8) continue;

        for (int i = 0; i < n_fits; i++) {
            if (strcmp(fits[i].name, sample) != 0) continue;
            SampleFits *s = &fits[i];
            s->row = realloc(s->row, (size_t)(s->n_rows + 1) * sizeof(FitRow));
            assert(s->row && "benchmark_choice: out of memory reading the manifest");
            s->row[s->n_rows].replicate = replicate;
            s->row[s->n_rows].log_likelihood = log_likelihood;
            s->row[s->n_rows].gradient = gradient;
            s->row[s->n_rows].converged = strcmp(converged, "yes") == 0;
            s->row[s->n_rows].inherited = strcmp(origin, "inherited") == 0;
            if (s->row[s->n_rows].converged) s->n_converged++;
            s->n_rows++;
            break;
        }
    }
    fclose(f);
}

static int by_gradient(const void *a, const void *b) {
    double ga = ((const FitRow*)a)->gradient, gb = ((const FitRow*)b)->gradient;
    return ga < gb ? -1 : ga > gb ? 1 : 0;
}

/*
Whether the weighted score loss can be built at this fit: the observed
information at the fit's own estimate on its own series, positive definite
against the same floor applications/abm_system_score_loss.c reads a flat
direction against. smallest and condition come back for the report whatever the
verdict is.
*/
static int information_is_usable(const char *sample, int replicate, double *smallest_out,
                                 double *condition_out) {
    QvarmaParams m = spec_shape();
    char path[640];
    snprintf(path, sizeof path, "%s/%s/replicate_%03d_%s_fit.json", FIT_DIR, sample, replicate,
             SPEC_LABEL);
    if (!qvarma_load_params(&m, path)) { qvarma_params_free(&m); return 0; }

    int n = qvarma_n_theta(&m);
    Vec theta = mat_new(n, 1);
    _qvarma_unlink(&m, theta);

    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", INPUT_DIR, sample);
    Mat y = abm_system_read_replicate(dir, replicate);

    Mat information = _qvarma_hessian(theta, &m, y);
    Vec eigenvalues;
    Mat eigenvectors;
    int usable = 0;
    *smallest_out = NAN;
    *condition_out = NAN;
    if (mat_eig_sym_status(information, &eigenvalues, &eigenvectors) == 0) {
        mreal smallest = eigenvalues.d[0], largest = eigenvalues.d[0];
        for (int i = 1; i < n; i++) {
            if (eigenvalues.d[i] < smallest) smallest = eigenvalues.d[i];
            if (eigenvalues.d[i] > largest) largest = eigenvalues.d[i];
        }
        mreal floor_value = (mreal)(n * MEPS) * MABS(largest);
        usable = smallest > floor_value;
        *smallest_out = (double)smallest;
        *condition_out = smallest > 0 ? (double)(largest / smallest) : INFINITY;
        mat_free(eigenvalues);
        mat_free(eigenvectors);
    }

    mat_free(information);
    mat_free(y);
    mat_free(theta);
    qvarma_params_free(&m);
    return usable;
}

int main(void) {
    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "benchmark_choice: cannot open montecarlo/out/benchmark_choice.txt - "
                      "the directory has to exist");

    fprintf(report, "Which simulated run stands in for the real data in montecarlo/\n\n");
    fprintf(report, "step 1, the candidate pool: every configuration a real-data confidence set "
                     "kept, MCS_TR\n");
    int n_winners;
    Sample *winners = read_winners(&n_winners, report);
    assert(n_winners > 0 && "benchmark_choice: no winners found - run the real-data MCS first");

    fprintf(report, "  pool of %d:", n_winners);
    for (int i = 0; i < n_winners; i++) fprintf(report, " %s", winners[i].name);
    fprintf(report, "\n\n");

    SampleFits *fits = calloc((size_t)n_winners, sizeof(SampleFits));
    assert(fits);
    for (int i = 0; i < n_winners; i++) snprintf(fits[i].name, sizeof fits[i].name, "%s",
                                                 winners[i].name);
    read_manifest(fits, n_winners);

    fprintf(report, "step 2, how often each one's fits converged\n");
    fprintf(report, "  %-12s %8s %10s %12s\n", "configuration", "fits", "converged", "rate");
    int best = 0;
    for (int i = 0; i < n_winners; i++) {
        assert(fits[i].n_rows > 0 && "benchmark_choice: a winner has no rows in the fit manifest");
        fprintf(report, "  %-12s %8d %10d %11.1f%%\n", fits[i].name, fits[i].n_rows,
                fits[i].n_converged, 100.0 * fits[i].n_converged / fits[i].n_rows);
        double rate = (double)fits[i].n_converged / fits[i].n_rows;
        double best_rate = (double)fits[best].n_converged / fits[best].n_rows;
        if (rate > best_rate) best = i;
    }
    fprintf(report, "  chosen: %s\n\n", fits[best].name);

    SampleFits *chosen = &fits[best];
    qsort(chosen->row, (size_t)chosen->n_rows, sizeof(FitRow), by_gradient);

    fprintf(report, "step 3, its converged replicates by gradient norm, taking the first whose "
                     "information matrix is positive definite\n");
    fprintf(report, "  %-10s %14s %12s %16s %12s %10s\n", "replicate", "log_lik", "gradient",
            "smallest eigval", "condition", "usable");

    int picked = -1;
    int tried = 0;
    for (int i = 0; i < chosen->n_rows && tried < MAX_CANDIDATES && picked < 0; i++) {
        if (!chosen->row[i].converged || chosen->row[i].inherited) continue;
        tried++;
        double smallest, condition;
        int usable = information_is_usable(chosen->name, chosen->row[i].replicate, &smallest,
                                           &condition);
        fprintf(report, "  %-10d %14.4f %12.3e %16.4e %12.3e %10s\n", chosen->row[i].replicate,
                chosen->row[i].log_likelihood, chosen->row[i].gradient, smallest, condition,
                usable ? "yes" : "no");
        if (usable) picked = i;
    }
    assert(picked >= 0 && "benchmark_choice: no converged replicate of the chosen configuration "
                           "has a usable information matrix");

    fprintf(report, "\nbenchmark: %s replicate %03d\n", chosen->name, chosen->row[picked].replicate);
    fprintf(report, "  log-likelihood %.4f, gradient norm %.6e, converged, its own fit\n",
            chosen->row[picked].log_likelihood, chosen->row[picked].gradient);
    fprintf(report, "  replicate %03d is seed %d of the model, and every configuration's replicate "
                     "%03d is held out of the validation set\n",
            chosen->row[picked].replicate, chosen->row[picked].replicate + 1,
            chosen->row[picked].replicate);

    FILE *env = fopen(ENV_PATH, "w");
    assert(env && "benchmark_choice: cannot open montecarlo/out/benchmark.env for writing");
    fprintf(env, "ABM_SYSTEM_BENCHMARK_SAMPLE=%s\n", chosen->name);
    fprintf(env, "ABM_SYSTEM_BENCHMARK_REPLICATE=%d\n", chosen->row[picked].replicate);
    fclose(env);

    fclose(report);
    for (int i = 0; i < n_winners; i++) free(fits[i].row);
    free(fits);
    free(winners);
    return 0;
}
