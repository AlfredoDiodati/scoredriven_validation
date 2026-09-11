/*
Whether the two Model Confidence Set statistics et_al implements pick out the
same ABM configurations on this project's loss matrix.

MCS_TR compares every pair of configurations and rejects when any two look
different. MCS_TMAX compares each configuration with the average of those still
in the running and rejects when the worst of them looks worse than that
average. Both are Hansen, Lunde and Nason's (2011).
applications/abm_system_mcs.c reports MCS_TR; this runs both with every other
setting identical and writes where they agree.

With the bootstrap variance et_al draws the resamples once and reuses them in
every elimination round under both statistics, and the two runs share a seed,
so they are scored on the same resamples and a difference between them comes
from the statistic.

Reads out/abm_system_mse_qvarma_joint.csv. Writes
out/abm_system_mcs_statistic_comparison.txt, a summary followed by the full
report of each run, and out/abm_system_mcs_statistic_comparison.csv, one row
per configuration with each statistic's MCS p-value, membership of the set and
elimination round.
In EXPERIMENT_STEMS. Nothing printed.
*/

#define _POSIX_C_SOURCE 199309L

#include <et_al./inference/mcs.h>
#include <et_al./frame/csv.h>
#include <et_al./stats.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define LOSS_PATH "out/abm_system_mse_qvarma_joint.csv"
#define REPORT_PATH "out/abm_system_mcs_statistic_comparison.txt"
#define TABLE_PATH "out/abm_system_mcs_statistic_comparison.csv"
#define N_LAST_STANDING 10

enum { STATISTIC_TR, STATISTIC_TMAX, N_STATISTICS };
static const MCSStat statistic[N_STATISTICS] = { MCS_TR, MCS_TMAX };
static const char *const statistic_name[N_STATISTICS] = { "MCS_TR", "MCS_TMAX" };

/* Every numeric column except "replicate", which indexes the rows, is a
   configuration. */
static DataFrame read_losses(const char *path) {
    DataFrame raw = df_read_csv(path, csv_read_options_default());
    DataFrame losses = df_new(raw.r);
    for (int j = 0; j < raw.n_cols; j++) {
        if (raw.columns[j].type != COL_NUMERIC) continue;
        if (strcmp(raw.columns[j].name, "replicate") == 0) continue;
        df_add_numeric_col(&losses, raw.columns[j].name, df_col_numeric(&raw, raw.columns[j].name));
    }
    df_free(&raw);
    return losses;
}

/* The round in which each configuration left. The one left at the end has no
   round, and is placed after every round so the order can be ranked. */
static int departure_round(const MCSResult *res, int j) {
    return res->elimination_round[j] > 0 ? res->elimination_round[j] : res->n_rounds + 1;
}

static double seconds_between(struct timespec start, struct timespec end) {
    return (double)(end.tv_sec - start.tv_sec) + 1e-9 * (double)(end.tv_nsec - start.tv_nsec);
}

/* qsort passes no context, so the rounds it orders by are reached through this. */
static int *tr_departure_rounds_for_sort;

static int later_tr_departure_first(const void *a, const void *b) {
    int left = *(const int *)a, right = *(const int *)b;
    return tr_departure_rounds_for_sort[right] - tr_departure_rounds_for_sort[left];
}

int main(void) {
    DataFrame losses = read_losses(LOSS_PATH);
    int n_configurations = mcs_n_models(&losses);

    MCSOptions opt = mcs_options_default();
    opt.bootstrap = 10000;
    opt.block_length = 1;
    opt.variance = MCS_VARIANCE_BOOTSTRAP;

    MCSResult result[N_STATISTICS];
    double wall_seconds[N_STATISTICS];
    for (int s = 0; s < N_STATISTICS; s++) {
        opt.stat = statistic[s];
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        result[s] = mcs(&losses, opt);
        clock_gettime(CLOCK_MONOTONIC, &end);
        wall_seconds[s] = seconds_between(start, end);
    }
    const MCSResult *tr = &result[STATISTIC_TR];
    const MCSResult *tmax = &result[STATISTIC_TMAX];

    int in_both = 0, only_tr = 0, only_tmax = 0;
    Mat tr_round = mat_new(n_configurations, 1), tmax_round = mat_new(n_configurations, 1);
    Mat round_gap = mat_new(n_configurations, 1);
    double largest_pvalue_gap = 0;
    for (int j = 0; j < n_configurations; j++) {
        int in_tr = mcs_in_set(tr, j), in_tmax = mcs_in_set(tmax, j);
        in_both += in_tr && in_tmax;
        only_tr += in_tr && !in_tmax;
        only_tmax += in_tmax && !in_tr;
        AT(tr_round, j, 0) = (mreal)departure_round(tr, j);
        AT(tmax_round, j, 0) = (mreal)departure_round(tmax, j);
        AT(round_gap, j, 0) = (mreal)abs(departure_round(tr, j) - departure_round(tmax, j));
        double pvalue_gap = fabs(tr->pvalue[j] - tmax->pvalue[j]);
        if (pvalue_gap > largest_pvalue_gap) largest_pvalue_gap = pvalue_gap;
    }

    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "abm_system_mcs_statistic_comparison: cannot open the report path for writing");
    fprintf(report, "MCS_TR against MCS_TMAX on %s\n", LOSS_PATH);
    fprintf(report, "%d configurations, %d replicates each\n", n_configurations, losses.r);
    fprintf(report, "shared: alpha = %.3f, %d resamples, block length %d, bootstrap variance, seed %llu, stream %llu\n",
            opt.alpha, opt.bootstrap, opt.block_length, (unsigned long long)opt.seed, (unsigned long long)opt.stream);
    fprintf(report, "\nBoth statistics are scored on the same resamples, drawn once and reused in every round.\n\n");

    fprintf(report, "statistic  set size  decided by an accepted test  deciding round  last p-value  seconds\n");
    for (int s = 0; s < N_STATISTICS; s++)
        fprintf(report, "%-10s %8d %28s %15d %13.4f %8.1f\n", statistic_name[s], result[s].n_surviving,
                result[s].converged ? "yes" : "no", result[s].decided_round, result[s].final_pvalue,
                wall_seconds[s]);

    fprintf(report, "\nsets\n");
    fprintf(report, "  in both sets %d, only in MCS_TR's %d, only in MCS_TMAX's %d\n", in_both, only_tr, only_tmax);
    fprintf(report, "  largest difference between the two MCS p-values of one configuration %.4f\n",
            largest_pvalue_gap);

    fprintf(report, "\nelimination rounds, the configuration left at the end counted as round %d\n",
            n_configurations);
    fprintf(report, "  Spearman rank correlation between the two orders %.4f\n",
            (double)stats_spearman(tr_round, tmax_round));
    fprintf(report, "  absolute difference in round: median %.0f, 90th percentile %.0f, largest %.0f\n",
            (double)stats_median(round_gap), (double)stats_quantile(round_gap, (mreal)0.9),
            (double)stats_quantile(round_gap, (mreal)1));

    int *by_tr_departure = malloc((size_t)n_configurations * sizeof *by_tr_departure);
    int *tr_rounds = malloc((size_t)n_configurations * sizeof *tr_rounds);
    assert(by_tr_departure && tr_rounds);
    for (int j = 0; j < n_configurations; j++) {
        by_tr_departure[j] = j;
        tr_rounds[j] = departure_round(tr, j);
    }
    tr_departure_rounds_for_sort = tr_rounds;
    qsort(by_tr_departure, (size_t)n_configurations, sizeof *by_tr_departure, later_tr_departure_first);

    fprintf(report, "\nthe %d configurations MCS_TR eliminated last, with where MCS_TMAX put them\n",
            N_LAST_STANDING);
    fprintf(report, "  configuration  mean loss  MCS_TR round  MCS_TR p  MCS_TMAX round  MCS_TMAX p\n");
    for (int i = 0; i < N_LAST_STANDING && i < n_configurations; i++) {
        int j = by_tr_departure[i];
        const char *name = mcs_model_name(&losses, j);
        fprintf(report, "  %s  %.5f  %d  %.4f  %d  %.4f\n", name,
                (double)stats_mean(df_col_numeric(&losses, name)), departure_round(tr, j), tr->pvalue[j],
                departure_round(tmax, j), tmax->pvalue[j]);
    }
    free(by_tr_departure);
    free(tr_rounds);

    for (int s = 0; s < N_STATISTICS; s++) {
        char title[160];
        snprintf(title, sizeof title, "MCS over t-QVARMA (driftless) impulse-response MAE, p1q1r2, %s statistic",
                 statistic_name[s]);
        fprintf(report, "\n\n");
        opt.stat = statistic[s];
        mcs_fwrite_options(report, &losses, opt);
        fprintf(report, "\n");
        mcs_fwrite_report(report, title, &losses, &result[s]);
    }
    fclose(report);

    const char **names = malloc((size_t)n_configurations * sizeof *names);
    assert(names);
    Vec mean_loss = vec_new(n_configurations);
    Vec tr_pvalue = vec_new(n_configurations), tr_in_set = vec_new(n_configurations);
    Vec tmax_pvalue = vec_new(n_configurations), tmax_in_set = vec_new(n_configurations);
    Vec tr_elimination_round = vec_new(n_configurations), tmax_elimination_round = vec_new(n_configurations);
    for (int j = 0; j < n_configurations; j++) {
        names[j] = mcs_model_name(&losses, j);
        AT(mean_loss, j, 0) = stats_mean(df_col_numeric(&losses, names[j]));
        AT(tr_pvalue, j, 0) = (mreal)tr->pvalue[j];
        AT(tr_in_set, j, 0) = (mreal)mcs_in_set(tr, j);
        AT(tmax_pvalue, j, 0) = (mreal)tmax->pvalue[j];
        AT(tmax_in_set, j, 0) = (mreal)mcs_in_set(tmax, j);
        AT(tr_elimination_round, j, 0) = (mreal)tr->elimination_round[j];
        AT(tmax_elimination_round, j, 0) = (mreal)tmax->elimination_round[j];
    }
    DataFrame table = df_new(n_configurations);
    df_add_string_col(&table, "model", names);
    df_add_numeric_col(&table, "mean_loss", mean_loss);
    df_add_numeric_col(&table, "tr_pvalue", tr_pvalue);
    df_add_numeric_col(&table, "tr_in_set", tr_in_set);
    df_add_numeric_col(&table, "tr_elimination_round", tr_elimination_round);
    df_add_numeric_col(&table, "tmax_pvalue", tmax_pvalue);
    df_add_numeric_col(&table, "tmax_in_set", tmax_in_set);
    df_add_numeric_col(&table, "tmax_elimination_round", tmax_elimination_round);
    df_write_csv(&table, TABLE_PATH, csv_write_options_default());

    df_free(&table);
    free(names);
    mat_free(mean_loss);
    mat_free(tr_pvalue); mat_free(tr_in_set); mat_free(tr_round);
    mat_free(tmax_pvalue); mat_free(tmax_in_set); mat_free(tmax_round);
    mat_free(round_gap);
    mat_free(tr_elimination_round); mat_free(tmax_elimination_round);
    for (int s = 0; s < N_STATISTICS; s++) mcs_free(&result[s]);
    df_free(&losses);
    return 0;
}
