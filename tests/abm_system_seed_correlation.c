/*
Whether reusing the same seeds across configurations made their losses move
together, which is what would make comparing configurations more precise.

Replicate n of every configuration is simulated with seed n. If the same seed
produces similar randomness under different parameters, the losses of two
configurations are positively correlated across replicates and

    Var(L_i - L_j) = Var(L_i) + Var(L_j) - 2 Cov(L_i, L_j)

is smaller than the independent sum. In an agent-based model the parameters
change how many random draws each period consumes, so the streams of two
configurations fall out of step and the correlation can vanish.

Reads out/abm_system_mse_qvarma_joint.csv, one row per replicate and one column
per configuration, and for every pair of configurations computes the Pearson
correlation of their per-replicate losses. The reference for no shared
randomness is the same correlation with the rows of the second configuration
moved down by one replicate, pairing seed n with seed n + 1 (the last row
wraps to the first). Under independence both have mean zero and standard
deviation close to 1/sqrt(1000) = 0.032.

For every pair it also reports the variance ratio

    Var(L_i - L_j) / (Var(L_i) + Var(L_j)),

1 when the seeds bought nothing, below 1 when they reduced the variance of the
comparison. The pair the confidence set separated last is reported on its own,
with the standard error of its mean loss difference computed both ways.

Whether nearby configurations keep more of the correlation is measured by the
rank correlation between a pair's loss correlation and its distance in the
parameter design, each of the nine parameters scaled to [0, 1] by its range in
dataset/abm_system_design.csv.

Reads only. Writes out/abm_system_seed_correlation_report.txt.
*/
#include "applications/abm_system.h"
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LOSS_TABLE_PATH "out/abm_system_mse_qvarma_joint.csv"
#define CONFIDENCE_SET_PATH "out/abm_system_mcs_joint.csv"
#define DESIGN_PATH "dataset/abm_system_design.csv"
#define REPORT_PATH "out/abm_system_seed_correlation_report.txt"

/* Pairs sampled for the rank correlation against design distance. Every pair
   is used for everything else. */
#define DISTANCE_SAMPLE_STEP 7

static double quantile_of(const double *values, long n, double probability) {
    Mat copy = mat_new((int)n, 1);
    for (long i = 0; i < n; i++) copy.d[i] = (mreal)values[i];
    double q = (double)stats_quantile(copy, (mreal)probability);
    mat_free(copy);
    return q;
}

static double mean_of(const double *values, long n) {
    double sum = 0;
    for (long i = 0; i < n; i++) sum += values[i];
    return sum / (double)n;
}

/* Pearson correlation of columns a and b of a replicates x models table, with
   b read `shift` rows further down, wrapping at the end. */
static double column_correlation(const double *table, int rows, int a, int b, int models, int shift,
                                 const double *mean, const double *spread) {
    double sum = 0;
    for (int r = 0; r < rows; r++) {
        int r_b = (r + shift) % rows;
        sum += (table[(long)r * models + a] - mean[a]) * (table[(long)r_b * models + b] - mean[b]);
    }
    return sum / (double)rows / (spread[a] * spread[b]);
}

int main(void) {
    DataFrame losses = df_read_csv(LOSS_TABLE_PATH, csv_read_options_default());
    int rows = losses.r;
    int models = losses.n_cols - 1;
    assert(rows > 2 && models > 1);

    char **name = (char**)malloc((size_t)models * sizeof(char*));
    double *table = (double*)malloc((size_t)rows * models * sizeof(double));
    int at = 0;
    for (int c = 0; c < losses.n_cols; c++) {
        if (strcmp(losses.columns[c].name, "replicate") == 0) continue;
        name[at] = losses.columns[c].name;
        Mat column = df_col_numeric(&losses, name[at]);
        for (int r = 0; r < rows; r++) table[(long)r * models + at] = (double)AT(column, r, 0);
        at++;
    }
    assert(at == models);

    double *mean = (double*)calloc((size_t)models, sizeof(double));
    double *spread = (double*)calloc((size_t)models, sizeof(double));
    for (int m = 0; m < models; m++) {
        for (int r = 0; r < rows; r++) mean[m] += table[(long)r * models + m];
        mean[m] /= rows;
        for (int r = 0; r < rows; r++) {
            double d = table[(long)r * models + m] - mean[m];
            spread[m] += d * d;
        }
        spread[m] = sqrt(spread[m] / rows);
    }

    long n_pairs = (long)models * (models - 1) / 2;
    double *same_seed = (double*)malloc((size_t)n_pairs * sizeof(double));
    double *shifted_seed = (double*)malloc((size_t)n_pairs * sizeof(double));
    double *variance_ratio = (double*)malloc((size_t)n_pairs * sizeof(double));
    long *pair_offset = (long*)malloc((size_t)models * sizeof(long));
    long running = 0;
    for (int a = 0; a < models; a++) { pair_offset[a] = running; running += models - 1 - a; }

    #pragma omp parallel for schedule(dynamic)
    for (int a = 0; a < models; a++)
        for (int b = a + 1; b < models; b++) {
            long k = pair_offset[a] + (b - a - 1);
            double rho = column_correlation(table, rows, a, b, models, 0, mean, spread);
            same_seed[k] = rho;
            shifted_seed[k] = column_correlation(table, rows, a, b, models, 1, mean, spread);
            double var_a = spread[a] * spread[a], var_b = spread[b] * spread[b];
            variance_ratio[k] = (var_a + var_b - 2 * rho * spread[a] * spread[b]) / (var_a + var_b);
        }

    /* The last two configurations of the confidence set: the winner and the
       one eliminated in the final round. */
    DataFrame set = df_read_csv(CONFIDENCE_SET_PATH, csv_read_options_default());
    char **set_model = df_col_string(&set, "model");
    Mat in_set = df_col_numeric(&set, "in_set");
    Mat round = df_col_numeric(&set, "elimination_round");
    Mat set_loss = df_col_numeric(&set, "mean_loss");
    int winner_row = -1, runner_row = -1;
    for (int i = 0; i < set.r; i++) {
        if (AT(in_set, i, 0) > 0 && (winner_row < 0 || AT(set_loss, i, 0) < AT(set_loss, winner_row, 0))) winner_row = i;
        if (AT(in_set, i, 0) == 0 && (runner_row < 0 || AT(round, i, 0) > AT(round, runner_row, 0))) runner_row = i;
    }
    assert(winner_row >= 0 && runner_row >= 0);
    int winner = -1, runner = -1;
    for (int m = 0; m < models; m++) {
        if (strcmp(name[m], set_model[winner_row]) == 0) winner = m;
        if (strcmp(name[m], set_model[runner_row]) == 0) runner = m;
    }
    assert(winner >= 0 && runner >= 0);

    int first = winner < runner ? winner : runner, second = winner < runner ? runner : winner;
    long winner_pair = pair_offset[first] + (second - first - 1);
    double difference_mean = 0, difference_variance = 0;
    for (int r = 0; r < rows; r++) difference_mean += table[(long)r * models + winner] - table[(long)r * models + runner];
    difference_mean /= rows;
    for (int r = 0; r < rows; r++) {
        double d = table[(long)r * models + winner] - table[(long)r * models + runner] - difference_mean;
        difference_variance += d * d;
    }
    difference_variance /= rows;
    double se_independent = sqrt((spread[winner] * spread[winner] + spread[runner] * spread[runner]) / rows);
    double se_paired = sqrt(difference_variance / rows);

    /* Design distance, on every DISTANCE_SAMPLE_STEP-th pair. Configuration
       cop_NNNN is design row NNNN. */
    DataFrame design = df_read_csv(DESIGN_PATH, csv_read_options_default());
    int n_parameters = design.n_cols;
    const char *const *parameter = abm_system_parameter_names();
    assert(n_parameters == ABM_SYSTEM_N_PARAMETERS);
    double low[ABM_SYSTEM_N_PARAMETERS], high[ABM_SYSTEM_N_PARAMETERS];
    Mat parameter_column[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < n_parameters; p++) {
        parameter_column[p] = df_col_numeric(&design, parameter[p]);
        low[p] = INFINITY; high[p] = -INFINITY;
        for (int i = 0; i < design.r; i++) {
            double v = (double)AT(parameter_column[p], i, 0);
            if (v < low[p]) low[p] = v;
            if (v > high[p]) high[p] = v;
        }
    }
    int *design_row = (int*)malloc((size_t)models * sizeof(int));
    for (int m = 0; m < models; m++) {
        int index = 0;
        int parsed = sscanf(name[m], "cop_%d", &index);
        assert(parsed == 1 && index >= 1 && index <= design.r);
        design_row[m] = index - 1;
    }
    long n_sampled = n_pairs / DISTANCE_SAMPLE_STEP;
    Mat sampled_distance = mat_new((int)n_sampled, 1), sampled_correlation = mat_new((int)n_sampled, 1);
    long taken = 0;
    for (int a = 0; a < models && taken < n_sampled; a++)
        for (int b = a + 1; b < models && taken < n_sampled; b++) {
            long k = pair_offset[a] + (b - a - 1);
            if (k % DISTANCE_SAMPLE_STEP != 0) continue;
            double squared = 0;
            for (int p = 0; p < n_parameters; p++) {
                double d = ((double)AT(parameter_column[p], design_row[a], 0)
                            - (double)AT(parameter_column[p], design_row[b], 0)) / (high[p] - low[p]);
                squared += d * d;
            }
            sampled_distance.d[taken] = (mreal)sqrt(squared);
            sampled_correlation.d[taken] = (mreal)same_seed[k];
            taken++;
        }
    Mat distance_view = mat_slice(sampled_distance, 0, (int)taken, 0, 1);
    Mat correlation_view = mat_slice(sampled_correlation, 0, (int)taken, 0, 1);
    double distance_rank_correlation = (double)stats_spearman(distance_view, correlation_view);

    /* Mean correlation by distance quartile. */
    double distance_cut[3];
    for (int q = 0; q < 3; q++) distance_cut[q] = (double)stats_quantile(distance_view, (mreal)(0.25 * (q + 1)));
    double quartile_sum[4] = { 0 };
    long quartile_count[4] = { 0 };
    for (long i = 0; i < taken; i++) {
        double d = (double)sampled_distance.d[i];
        int q = d <= distance_cut[0] ? 0 : (d <= distance_cut[1] ? 1 : (d <= distance_cut[2] ? 2 : 3));
        quartile_sum[q] += (double)sampled_correlation.d[i];
        quartile_count[q]++;
    }

    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "abm_system_seed_correlation: cannot open the report path for writing");
    fprintf(out, "per-replicate losses from %s: %d replicates x %d configurations, %ld pairs\n", LOSS_TABLE_PATH,
            rows, models, n_pairs);
    fprintf(out, "same seed: replicate n of both configurations. shifted: replicate n against replicate n + 1, the\n");
    fprintf(out, "reference for no shared randomness. Under independence the standard deviation is about %.3f.\n\n",
            1.0 / sqrt((double)rows));

    fprintf(out, "%-34s %9s %9s %9s %9s %9s %9s\n", "over all pairs", "mean", "p05", "p25", "p50", "p75", "p95");
    fprintf(out, "%-34s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", "loss correlation, same seed", mean_of(same_seed, n_pairs),
            quantile_of(same_seed, n_pairs, 0.05), quantile_of(same_seed, n_pairs, 0.25), quantile_of(same_seed, n_pairs, 0.5),
            quantile_of(same_seed, n_pairs, 0.75), quantile_of(same_seed, n_pairs, 0.95));
    fprintf(out, "%-34s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", "loss correlation, shifted seed", mean_of(shifted_seed, n_pairs),
            quantile_of(shifted_seed, n_pairs, 0.05), quantile_of(shifted_seed, n_pairs, 0.25),
            quantile_of(shifted_seed, n_pairs, 0.5), quantile_of(shifted_seed, n_pairs, 0.75),
            quantile_of(shifted_seed, n_pairs, 0.95));
    fprintf(out, "%-34s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", "variance ratio, same seed", mean_of(variance_ratio, n_pairs),
            quantile_of(variance_ratio, n_pairs, 0.05), quantile_of(variance_ratio, n_pairs, 0.25),
            quantile_of(variance_ratio, n_pairs, 0.5), quantile_of(variance_ratio, n_pairs, 0.75),
            quantile_of(variance_ratio, n_pairs, 0.95));
    long above_two_sd = 0;
    for (long k = 0; k < n_pairs; k++) above_two_sd += same_seed[k] > 2.0 / sqrt((double)rows);
    long shifted_above_two_sd = 0;
    for (long k = 0; k < n_pairs; k++) shifted_above_two_sd += shifted_seed[k] > 2.0 / sqrt((double)rows);
    fprintf(out, "pairs with correlation above 2/sqrt(%d) = %.3f: same seed %.2f%%, shifted %.2f%%\n\n", rows,
            2.0 / sqrt((double)rows), 100.0 * above_two_sd / n_pairs, 100.0 * shifted_above_two_sd / n_pairs);

    fprintf(out, "last pair of the confidence set: %s (kept) and %s (eliminated last)\n", name[winner], name[runner]);
    fprintf(out, "  loss correlation, same seed %.4f, shifted %.4f\n", same_seed[winner_pair], shifted_seed[winner_pair]);
    fprintf(out, "  mean loss difference %.6f\n", difference_mean);
    fprintf(out, "  standard error of that difference treating replicates as independent %.6f (%.2f standard errors)\n",
            se_independent, fabs(difference_mean) / se_independent);
    fprintf(out, "  standard error from the paired differences, which uses the shared seeds  %.6f (%.2f standard errors)\n\n",
            se_paired, fabs(difference_mean) / se_paired);

    fprintf(out, "correlation against distance in the design, %ld pairs (every %d-th), parameters scaled to [0,1]\n", taken,
            DISTANCE_SAMPLE_STEP);
    fprintf(out, "  Spearman rank correlation, distance against loss correlation: %.4f\n", distance_rank_correlation);
    fprintf(out, "  mean loss correlation by distance quartile (nearest first):");
    for (int q = 0; q < 4; q++) fprintf(out, " %.4f", quartile_sum[q] / (double)quartile_count[q]);
    fprintf(out, "\n");
    fclose(out);

    mat_free(sampled_distance); mat_free(sampled_correlation);
    free(design_row); free(pair_offset); free(same_seed); free(shifted_seed); free(variance_ratio);
    free(mean); free(spread); free(table); free(name);
    df_free(&design); df_free(&set); df_free(&losses);
    return 0;
}
