/*
Model Confidence Set (Hansen, Lunde and Nason 2011), over
out/abm_system_irf_loss.csv: the mean absolute error between the
driftless t-QVARMA real-data fit's own impulse response function (p1q1r2)
and every simulated replicate's own IRF. Every numeric column except
"replicate" is a "model", one per configuration of the design, and every row
is an "observation", one per replicate. Over the design experiment that is
1000 models and 1000 replicates in one MCS run. Which configurations cannot be
statistically
distinguished from the one with the smallest average loss - using MCS_TR
(every pairwise differential, rejects when any two models look different
from each other), not et_al.'s own default MCS_TMAX (each model against the
field average, rejects only when the single worst model looks worse than
everyone else) - see the opt.stat assignment below for why this file asks
for TR specifically.

The comparison object is the impulse response function, not the fitted
model's own parameters - see abm_system_irf_loss.c's own header comment
for why (in short: two fits with different-looking coefficients can imply
nearly identical dynamics and vice versa, so a parameter distance does not
answer the question this test needs answered, "do these two models behave
alike"). This changed from an earlier version that fed this file a
constrained-parameter distance instead; the switch is in
abm_system_irf_loss.c, not here - this file only ever consumes whatever
loss the joint csv holds, whatever it is computed from.

The drift-carrying t-QVARMAd is deliberately excluded. Both model families
were run through one joint MCS once, all 400 columns together, back when
both used constrained-parameter loss, and every one of them survived - a
result driven by measurement noise (see the absolute-vs-squared error point
below) rather than any real inability to tell the models apart. Its loss was
never switched from parameter distance to IRF distance and its code is not
part of this repository, so the two are not comparable side by side even if
it were wanted. docs/ABM_SYSTEM_MCS_VALIDATION.md records that run.

Absolute error, not squared: abm_system_irf_loss.c computes loss as mean
absolute error (et_al.'s stats_mae), not mean squared error (stats_mse,
used originally, for both the earlier parameter-distance version and the
first version of the IRF distance). A squared difference turns one
badly-fit replicate's IRF into a term that can outweigh every other
replicate combined by many orders of magnitude, which then inflates not
just that model's mean loss but the bootstrap variance this file's own
mcs() call estimates from - exactly what made an earlier squared-error run
over both model families unable to reject anything. Absolute error still
counts that same replicate, linearly rather than squared, so it can no
longer single-handedly swamp every other observation.

The replicates are independent Monte Carlo draws, not a time series -
there is no serial dependence for a block bootstrap to protect against, so
it is switched off:

  - block_length = 1. et_al.'s own mcs.h documents this as the literal
    iid-bootstrap degenerate case, not an approximation to one: a block of
    length 1 draws a single independent row at a time, which is what an iid
    bootstrap is. block_length = 0 is not a meaningful request (a block
    bootstrap on zero-length blocks draws nothing) and mcs()'s own assert
    (block_length >= 1) would reject it; block_length = 1 is how this
    library expresses exactly the request made - no blocking at all.

variance = MCS_VARIANCE_BOOTSTRAP (mcs.h's own default, set explicitly here
anyway since the choice matters and should not depend on silently inheriting
whatever the library defaults to): the variance of each t-statistic's own
denominator comes from the spread of the resampled means themselves, across
the bootstrap draws - Hansen, Lunde and Nason's own estimator, not a HAC
long-run variance computed once from the data. hac_lag is not set here
because MCS_VARIANCE_BOOTSTRAP ignores it entirely (mcs.h's own comment on
MCSOptions.hac_lag). An earlier version of this file used
MCS_VARIANCE_HAC with hac_lag = 0, the library's only option before this
project's own et_al update added MCS_VARIANCE_BOOTSTRAP - that forced "no
HAC adjustment" to be expressed as a zero-lag HAC estimate, functionally
the plain sample variance, rather than what independent replicates actually
call for.

One joint run over both specs, not one per spec: out/abm_system_irf_loss.csv
already carries one column per configuration, each labeled with the spec it
was fitted under (applications/abm_system_irf_loss.c builds it that way), so
every column of it except "replicate" is a model this file's own mcs() call
treats uniformly. It reads however many columns the file holds, so restoring a
second spec upstream needs no change here: the columns simply double and the
MCS runs jointly across both.

Requires out/abm_system_irf_loss.csv to exist
(applications/abm_system_irf_loss.c). Output:
out/abm_system_mcs.txt
(mcs_fwrite_report's own table: every model's mean loss, MCS p-value, and
whether it survived) and out/abm_system_mcs.csv (the same as data,
via mcs_pvalue_frame). In EXPERIMENT_STEMS. Nothing printed.
*/

#include <et_al./inference/mcs.h>
#include <et_al./frame/csv.h>
#include <string.h>
#include <assert.h>

#define LOSS_PATH "out/abm_system_irf_loss.csv"

int main(void) {
    /* Every numeric column except "replicate" (the row index, not a
       model) becomes a model - no name lookup needed, unlike the
       per-spec version this replaced, since the joint csv's own column
       names are already unique and self-describing (spec-suffixed). */
    DataFrame raw = df_read_csv(LOSS_PATH, csv_read_options_default());
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

    FILE *report = fopen("out/abm_system_mcs.txt", "w");
    assert(report && "abm_system_mcs: cannot open the report path for writing");
    mcs_fwrite_options(report, &losses, opt);
    fprintf(report, "\n");
    mcs_fwrite_report(report, "MCS over t-QVARMA (driftless) impulse-response MAE, "
                              "p1q1r2, one model per configuration, "
                              "iid bootstrap, bootstrap variance, MCS_TR statistic",
                      &losses, &res);
    fclose(report);

    DataFrame pvalues = mcs_pvalue_frame(&losses, &res);
    df_write_csv(&pvalues, "out/abm_system_mcs.csv", csv_write_options_default());
    df_free(&pvalues);

    mcs_free(&res);
    df_free(&losses);
    return 0;
}
