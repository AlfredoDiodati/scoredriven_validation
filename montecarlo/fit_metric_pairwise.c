/*
Pairwise tests of whether one auxiliary model's fit metric v is higher than
another's (docs/MONTECARLO_FIT_METRIC_COMPARISON.md, "Pairwise tests"), on the
oracle rows bin/fit_metric_oracle writes: v_ijk for auxiliary model i, CoP j
and benchmark run k.

Across CoPs: v_ij is the mean of v_ijk over the runs of CoP j, and for a pair
of auxiliary models i, i' the null is E[v_ij - v_i'j] <= 0, the expectation
over CoPs, with the CoPs as the observations.

Per CoP: for a pair i, i' and one CoP j the null is E[v_ijk - v_i'jk] <= 0,
the expectation over the runs of that CoP, with its runs as the observations.
A CoP counts as one where i is higher when the two-sided p-value is below
ALPHA / N_COPS, Bonferroni over the CoPs of one pair, and the mean difference
is positive. et_al has no Holm or Benjamini-Hochberg procedure.

Both are et_al's dm_test at its defaults, truncation lag 0, on a table whose
columns are the auxiliary models' v rather than losses: the mean difference
divided by its standard error, against N(0,1). Run k of a CoP is the same
simulated series for every auxiliary model, so both comparisons are paired.

Reads montecarlo/out/fit_metric_oracle.csv.gz. Writes
montecarlo/out/fit_metric_pairwise.txt, the same tables as LaTeX in
montecarlo/out/fit_metric_pairwise.tex, and one row per pair and CoP in
montecarlo/out/fit_metric_pairwise_per_cop.csv. Nothing printed.
*/

#include "montecarlo/fit_metric.h"
#include "montecarlo/fit_metric_table.h"
#include <et_al./frame/gzip.h>

#define OUT_DIR "montecarlo/out"
#define ORACLE_PATH OUT_DIR "/fit_metric_oracle.csv.gz"
#define REPORT_PATH OUT_DIR "/fit_metric_pairwise.txt"
#define LATEX_PATH OUT_DIR "/fit_metric_pairwise.tex"
#define PER_COP_PATH OUT_DIR "/fit_metric_pairwise_per_cop.csv"
#define N_COPS 1000
#define N_RUNS 1000
#define N_PAIRS (N_MODELS * (N_MODELS - 1) / 2)
#define ALPHA 0.05

static int model_of(const char *label) {
    for (int model = 0; model < N_MODELS; model++)
        if (strcmp(label, model_label[model]) == 0) return model;
    assert(0 && "fit_metric_pairwise: unknown model label");
    return -1;
}

static size_t cell_of(int model, int cop, int run) { return ((size_t)model * N_COPS + cop) * N_RUNS + run; }

/* v[cell_of(model, cop, run)], CoPs counted from 0. */
static double *read_oracle_v(void) {
    FILE *f = fopen(ORACLE_PATH, "rb");
    assert(f && "fit_metric_pairwise: montecarlo/out/fit_metric_oracle.csv.gz is missing");
    fseek(f, 0, SEEK_END);
    size_t length = (size_t)ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc(length);
    size_t n_read = fread(compressed, 1, length, f);
    assert(n_read == length && "fit_metric_pairwise: short read of the oracle file");
    (void)n_read;
    fclose(f);
    size_t text_length;
    char *text = (char *)gzip_inflate(compressed, length, &text_length);
    free(compressed);

    size_t n_cells = (size_t)N_MODELS * N_COPS * N_RUNS, n_filled = 0;
    double *v = malloc(n_cells * sizeof *v);
    assert(v && "fit_metric_pairwise: out of memory");
    for (size_t c = 0; c < n_cells; c++) v[c] = NAN;
    char *line = strchr(text, '\n') + 1;
    while (line < text + text_length && *line) {
        char *end = strchr(line, '\n');
        if (end) *end = 0;
        char label[32];
        int cop, run, lowest;
        double d, sigma2, s, value;
        int matched = sscanf(line, "%31[^,],%d,%d,%d,%lf,%lf,%lf,%lf", label, &cop, &run, &lowest, &d, &sigma2, &s,
                             &value);
        assert(matched == 8 && "fit_metric_pairwise: a malformed row in the oracle file");
        (void)matched;
        assert(cop >= 1 && cop <= N_COPS && run >= 0 && run < N_RUNS && "fit_metric_pairwise: CoP or run out of range");
        v[cell_of(model_of(label), cop - 1, run)] = value;
        n_filled++;
        if (!end) break;
        line = end + 1;
    }
    free(text);
    assert(n_filled == n_cells && "fit_metric_pairwise: the oracle file does not have one row per model, CoP and run");
    return v;
}

/* The two one-sided p-values of a dm_test on columns i, i' of v: for the null
   that i is not higher, and for the null that i' is not higher. */
static double p_not_higher(const DieboldMariano *dm) { return special_norm_cdf(-dm->stat); }
static double p_other_not_higher(const DieboldMariano *dm) { return special_norm_cdf(dm->stat); }

static void across_cops_table(Table *t, const double *v) {
    static const char *const names[] = { "model i", "model i'", "mean v_i - v_i'", "std error", "z",
                                         "p, i not higher", "p, i' not higher" };
    set_header(t, 7, names);
    Mat mean_v = mat_new(N_COPS, N_MODELS);
    for (int model = 0; model < N_MODELS; model++)
        for (int cop = 0; cop < N_COPS; cop++) {
            double total = 0;
            for (int run = 0; run < N_RUNS; run++) total += v[cell_of(model, cop, run)];
            AT(mean_v, cop, model) = total / N_RUNS;
        }
    DataFrame frame = df_from_matrix(mean_v, model_label);
    mat_free(mean_v);

    t->n_rows = 0;
    for (int i = 0; i < N_MODELS; i++)
        for (int i_prime = i + 1; i_prime < N_MODELS; i_prime++) {
            DieboldMariano dm = dm_test(&frame, model_label[i], model_label[i_prime], dm_options_default());
            char (*row)[CELL] = t->cell[t->n_rows++];
            snprintf(row[0], CELL, "%s", model_name[i]);
            snprintf(row[1], CELL, "%s", model_name[i_prime]);
            snprintf(row[2], CELL, "%.4f", dm.mean_diff);
            snprintf(row[3], CELL, "%.4f", dm.std_error);
            snprintf(row[4], CELL, "%.2f", dm.stat);
            snprintf(row[5], CELL, "%.3g", p_not_higher(&dm));
            snprintf(row[6], CELL, "%.3g", p_other_not_higher(&dm));
        }
    df_free(&frame);
    snprintf(t->note, sizeof t->note,
             "v_ij is the mean of v over the %d runs of CoP j. Observations: the %d CoPs, treated as independent. "
             "The null of the first p-value is E[v_ij - v_i'j] <= 0, of the second E[v_ij - v_i'j] >= 0.",
             N_RUNS, N_COPS);
}

/* Counts by pair: CoPs where the mean of v_i over the runs is above v_i',
   where i or i' is higher at the Bonferroni level, and the same unadjusted. */
typedef struct { int i_above, i_higher, i_prime_higher, i_higher_unadjusted, i_prime_higher_unadjusted; } PairCounts;

static void per_cop_tests(const double *v, PairCounts *counts) {
    FILE *csv = fopen(PER_COP_PATH, "w");
    assert(csv && "fit_metric_pairwise: cannot open the per-CoP output");
    fprintf(csv, "model_i,model_i_prime,cop,mean_difference,std_error,z,p_two_sided,p_i_not_higher,"
                 "p_i_prime_not_higher\n");
    memset(counts, 0, N_PAIRS * sizeof *counts);
    const double bonferroni = ALPHA / N_COPS;
    Mat runs = mat_new(N_RUNS, N_MODELS);
    for (int cop = 0; cop < N_COPS; cop++) {
        for (int model = 0; model < N_MODELS; model++)
            for (int run = 0; run < N_RUNS; run++) AT(runs, run, model) = v[cell_of(model, cop, run)];
        DataFrame frame = df_from_matrix(runs, model_label);
        int pair = 0;
        for (int i = 0; i < N_MODELS; i++)
            for (int i_prime = i + 1; i_prime < N_MODELS; i_prime++, pair++) {
                DieboldMariano dm = dm_test(&frame, model_label[i], model_label[i_prime], dm_options_default());
                PairCounts *c = &counts[pair];
                c->i_above += dm.mean_diff > 0;
                c->i_higher += dm.pvalue < bonferroni && dm.mean_diff > 0;
                c->i_prime_higher += dm.pvalue < bonferroni && dm.mean_diff < 0;
                c->i_higher_unadjusted += dm.pvalue < ALPHA && dm.mean_diff > 0;
                c->i_prime_higher_unadjusted += dm.pvalue < ALPHA && dm.mean_diff < 0;
                fprintf(csv, "%s,%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", model_label[i], model_label[i_prime],
                        cop + 1, dm.mean_diff, dm.std_error, dm.stat, dm.pvalue, p_not_higher(&dm),
                        p_other_not_higher(&dm));
            }
        df_free(&frame);
    }
    mat_free(runs);
    fclose(csv);
}

static void per_cop_table(Table *t, const PairCounts *counts) {
    static const char *const names[] = { "model i", "model i'", "mean v_i above v_i'", "i higher", "i' higher",
                                         "i higher, unadjusted", "i' higher, unadjusted" };
    set_header(t, 7, names);
    t->n_rows = 0;
    int pair = 0;
    for (int i = 0; i < N_MODELS; i++)
        for (int i_prime = i + 1; i_prime < N_MODELS; i_prime++, pair++) {
            char (*row)[CELL] = t->cell[t->n_rows++];
            snprintf(row[0], CELL, "%s", model_name[i]);
            snprintf(row[1], CELL, "%s", model_name[i_prime]);
            snprintf(row[2], CELL, "%d", counts[pair].i_above);
            snprintf(row[3], CELL, "%d", counts[pair].i_higher);
            snprintf(row[4], CELL, "%d", counts[pair].i_prime_higher);
            snprintf(row[5], CELL, "%d", counts[pair].i_higher_unadjusted);
            snprintf(row[6], CELL, "%d", counts[pair].i_prime_higher_unadjusted);
        }
    snprintf(t->note, sizeof t->note,
             "Number of CoPs out of %d. One test per CoP, observations: its %d runs. Higher: two-sided p-value "
             "below %g / %d = %.0e (Bonferroni over the CoPs of one pair), direction from the sign of the mean "
             "difference. Unadjusted: two-sided p-value below %g.",
             N_COPS, N_RUNS, ALPHA, N_COPS, ALPHA / N_COPS, ALPHA);
}

int main(void) {
    double *v = read_oracle_v();

    static Table table[2];
    table[0] = (Table){ .label = "pairwise_across_cops",
                        .title = "Is v higher under auxiliary model i than under i', across CoPs" };
    across_cops_table(&table[0], v);

    PairCounts counts[N_PAIRS];
    per_cop_tests(v, counts);
    table[1] = (Table){ .label = "pairwise_per_cop",
                        .title = "Is v higher under auxiliary model i than under i', one test per CoP" };
    per_cop_table(&table[1], counts);
    free(v);

    FILE *report = fopen(REPORT_PATH, "w"), *latex = fopen(LATEX_PATH, "w");
    assert(report && latex && "fit_metric_pairwise: cannot open the outputs");
    for (int i = 0; i < 2; i++) {
        write_text(report, &table[i]);
        write_latex(latex, &table[i]);
    }
    fclose(report);
    fclose(latex);
    return 0;
}
