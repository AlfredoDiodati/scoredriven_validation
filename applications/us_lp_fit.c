/*
Fits the linear and the state-dependent local projections to the US data, the
real-data benchmark the simulated fits of applications/abm_system_fit_lp.c are
compared against.

The series are the ones the t-QVARMA is fitted on: applications/us_data.h's
us_system_growth_block turns out/us_system.csv into GDP growth, energy growth,
employment change, inflation and the interest rate, 187 quarters from 1973Q2
to 2019Q4. applications/lp_system.h's lp_system_from_stored then drops the
first three, where R's 4-period moving average of log GDP, the series that
decides the state, does not exist: 184 quarters from 1974Q1. The settings are
lp_system_spec()'s, the same as for every simulated run: 4 lags, 15 horizons,
unit shocks in the order above, the state weight from the HP cycle (lambda
1600) of that moving average through a logistic of slope 2, lagged one period.
The simulated fits went through exactly this route, so the two sides differ
only in their data.

The fit is stored as a one-row archive in lp_system.h's layout,
out/us_lp_fit.npz, the layout every configuration's fits are stored in. A rerun
keeps it when it was fitted under this specification on exactly these series,
and fits again otherwise.

    ./bin/us_lp_fit              fit unless stored
    ./bin/us_lp_fit --refit      fit and overwrite
    ./bin/us_lp_fit --r-levels   the same on R's levels, the collaborator's
                                 transformation, applied to the US growth
                                 rates exactly as to the simulated ones

With --r-levels every file below carries _r_levels after us_lp_fit. The
collaborator's own preparation of the US data (rw_data_prep.R) was not
available; the thesis states the same transformation code runs on the real and
the simulated data, so R's transform_data is applied to the US growth rates as
it is to every simulated run.

Writes out/us_lp_fit.npz, out/us_lp_fit.txt (the data, the settings and the fit
notes) and out/us_lp_fit_irf.csv (every response: model, response, horizon,
shock, value). Nothing printed.
*/

#include "us_data.h"
#include "lp_system.h"
#include <cblas.h>
#include <string.h>
#include <assert.h>

static const char *const growth_name[LP_SYSTEM_K] = {
    "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
};
static const char *const r_level_name[LP_SYSTEM_K] = { "GDP", "Employment", "Inflation", "InterestRate", "Energy" };

int main(int argc, char **argv) {
    int force_refit;
    LpLayout layout = lp_parse_flags(argc, argv, &force_refit);
    openblas_set_num_threads(1);
    const char *const *series_name = layout == LP_LAYOUT_R_LEVELS ? r_level_name : growth_name;
    char archive_path[256], report_path[256], irf_path[256];
    snprintf(archive_path, sizeof archive_path, "out/us_lp_fit%s.npz", lp_layout_suffix(layout));
    snprintf(report_path, sizeof report_path, "out/us_lp_fit%s.txt", lp_layout_suffix(layout));
    snprintf(irf_path, sizeof irf_path, "out/us_lp_fit%s_irf.csv", lp_layout_suffix(layout));

    Mat original = load_us_system();
    Mat stored = us_system_growth_block(original);
    Mat y, switching;
    lp_system_series(stored, layout, &y, &switching);
    LpNlSpec spec = lp_system_spec();

    int count = 0;
    LpSystemRecord *record = force_refit ? NULL : lp_system_read_archive(archive_path, spec, &count);
    int loaded = record && count == 1 && lp_system_record_matches(&record[0], y, switching);
    if (!loaded) {
        free(record);
        record = malloc(sizeof(LpSystemRecord));
        assert(record && "us_lp_fit: out of memory");
        LpLinNlFit fit = lp_lin_and_nl(y, switching, spec, (LpBands){0});
        record[0] = lp_system_record(&fit, y, switching, 0);
        lp_lin_nl_fit_free(&fit);
        lp_system_write_archive(archive_path, record, 1, spec);
    }
    const LpSystemRecord *r = &record[0];

    FILE *report = fopen(report_path, "w");
    assert(report && "us_lp_fit: cannot open the report");
    fprintf(report, "Local projections on the US data\n\n");
    fprintf(report, "data: out/us_system.csv through us_system_growth_block, then %s; %d quarters from "
                    "1974Q1 to 2019Q4, in the order of the recursive identification\n",
            lp_layout_description(layout), y.c);
    fprintf(report, "state series: R's 4-period moving average of 100 log(100 + cumsum(GDP_growth))\n");
    fprintf(report, "spec: lags %d, horizons 0 to %d, unit shocks, HP lambda %g, logistic gamma %g, "
                    "switching lagged %d period\n", spec.lin.lags_endog_lin, spec.lin.hor, spec.lambda,
            spec.gamma, spec.lag_switching);
    fprintf(report, "%s %s\n\n", loaded ? "loaded, fitted on these data under this spec, from"
                                        : "fitted in this run and written to", archive_path);
    fprintf(report, "responses: %s\n", r->has_responses ? "yes" : "NO, the VAR behind the shocks failed");
    fprintf(report, "VAR behind the shocks: ols status %d, rank %d of %d, Cholesky status %d\n",
            r->var_ols_status, r->var_rank, 1 + LP_SYSTEM_K * LP_SYSTEM_LAGS, r->var_chol_status);
    int lin_fallback = 0, nl_fallback = 0;
    for (int h = 0; h < LP_SYSTEM_HOR; h++) {
        lin_fallback += r->lin_ols_status[h] != 0;
        nl_fallback += r->nl_ols_status[h] != 0;
    }
    fprintf(report, "horizons solved by the minimum-norm fallback for a rank-deficient design: linear %d, "
                    "state dependent %d, of %d each\n", lin_fallback, nl_fallback, LP_SYSTEM_HOR);
    fprintf(report, "switching series constant: %s\n", r->switching_is_constant ? "yes" : "no");
    fclose(report);

    FILE *irf = fopen(irf_path, "w");
    assert(irf && "us_lp_fit: cannot open the response table");
    fprintf(irf, "model,response,horizon,shock,value\n");
    const mreal *part[3] = { r->lin, r->s1, r->s2 };
    static const char *const part_name[3] = { "lin", "s1", "s2" };
    for (int p = 0; p < 3; p++)
        for (int k = 0; k < LP_SYSTEM_K; k++)
            for (int h = 0; h <= LP_SYSTEM_HOR; h++)
                for (int j = 0; j < LP_SYSTEM_K; j++)
                    fprintf(irf, "%s,%s,%d,%s,%.17g\n", part_name[p], series_name[k], h, series_name[j],
                            (double)part[p][(k * (LP_SYSTEM_HOR + 1) + h) * LP_SYSTEM_K + j]);
    fclose(irf);

    free(record);
    mat_free(y);
    mat_free(switching);
    mat_free(stored);
    mat_free(original);
    return 0;
}
