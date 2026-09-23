/*
The headline confidence set of the Monte Carlo experiment: MCS_TR over
montecarlo/out/irf_loss.csv, the impulse-response distance between one chosen
simulated run and every configuration's replicates.

This is applications/abm_system_mcs.c's procedure over a different loss matrix.
Every numeric column except "replicate" is a model, one per configuration, and
every row is an observation, one per replicate kept - 999 of them, the
benchmark's own seed being held out of all of them.

The question here has a known answer. The benchmark is a replicate of a
particular configuration, so that configuration should come back in the set.
docs/MONTECARLO_VALIDATION.md reports whether it does.

Settings identical to applications/abm_system_mcs.c so the two runs are
comparable: MCS_TR rather than the library default MCS_TMAX, 10000 resamples,
block length 1 for the iid bootstrap the independent replicates call for, and
the bootstrap variance of the resampled mean. That file's own header comment
explains each choice; none of them is revisited here.

Requires montecarlo/out/irf_loss.csv. Writes montecarlo/out/mcs.txt and
montecarlo/out/mcs.csv. Nothing printed.
*/

#include <et_al./inference/mcs.h>
#include <et_al./frame/csv.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define LOSS_PATH "montecarlo/out/irf_loss.csv"
#define REPORT_PATH "montecarlo/out/mcs.txt"
#define TABLE_PATH "montecarlo/out/mcs.csv"

int main(void) {
    const char *loss_path = LOSS_PATH;
    const char *report_path = REPORT_PATH;
    const char *table_path = TABLE_PATH;
    /* Every numeric column except "replicate" (the row index, not a
       model) becomes a model - no name lookup needed, unlike the
       per-spec version this replaced, since the joint csv's own column
       names are already unique and self-describing (spec-suffixed). */
    DataFrame raw = df_read_csv(loss_path, csv_read_options_default());
    DataFrame losses = df_new(raw.r);
    for (int j = 0; j < raw.n_cols; j++) {
        if (raw.columns[j].type != COL_NUMERIC) continue;
        if (strcmp(raw.columns[j].name, "replicate") == 0) continue;
        Mat col = df_col_numeric(&raw, raw.columns[j].name);
        df_add_numeric_col(&losses, raw.columns[j].name, col);
    }
    df_free(&raw);

    MCSOptions opt = mcs_options_default();
    /* Five times the library default of 2000, so a round p-value is resolved to
       1/10000 rather than 1/2000. */
    opt.bootstrap = 10000;
    opt.block_length = 1;                    /* iid bootstrap - see this file's own header comment */
    opt.variance = MCS_VARIANCE_BOOTSTRAP;   /* also the library default, set explicitly - see
                                                 this file's own header comment */
    /* MCS_TR ("range") instead of the default MCS_TMAX: every pairwise loss
       differential d_ij(t) = L(t,i) - L(t,j), max_{i!=j} |t_ij| - rejects
       when ANY two models look different from each other, not only when
       one model looks worse than the field average. O(M^2) t-statistics
       per round instead of TMAX's O(M), which is the whole reason TMAX is
       the default - requested directly, to see whether the two statistics
       actually disagree on this data rather than assuming they would. */
    opt.stat = MCS_TR;
    MCSResult res = mcs(&losses, opt);

    FILE *report = fopen(report_path, "w");
    assert(report && "abm_system_mcs: cannot open the report path for writing");
    mcs_fwrite_options(report, &losses, opt);
    fprintf(report, "\n");
    mcs_fwrite_report(report, "MCS over t-QVARMA (driftless) impulse-response MAE, "
                              "p1q1r2, one model per configuration, "
                              "iid bootstrap, bootstrap variance, MCS_TR statistic",
                      &losses, &res);
    fclose(report);

    DataFrame pvalues = mcs_pvalue_frame(&losses, &res);
    df_write_csv(&pvalues, table_path, csv_write_options_default());
    df_free(&pvalues);

    mcs_free(&res);
    df_free(&losses);
    return 0;
}
