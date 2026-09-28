/*
The local-projection loss tables against the US data: how far each
configuration's simulated responses sit from the US responses, one table per
model,

    lin   the linear local projection, irf_lin_mean
    s1    the state-dependent one in state 1 (lags times 1 - F, mostly periods
          with GDP above trend), irf_s1_mean
    s2    the same in state 2 (lags times F, mostly below trend), irf_s2_mean

the local-projection counterpart of applications/abm_system_irf_loss.c.

A cell's loss is the mean absolute difference between the US response vector
(out/us_lp_fit.npz, applications/us_lp_fit.c) and the run's, over all
K x (hor + 1) x K = 400 values of one model. The collaborator's R pipeline uses
the mean squared difference; the absolute one is kept for the reason
applications/abm_system_irf_loss.c gives and for comparability with the
t-QVARMA results. The Monte Carlo experiment with these models and this loss
recovered a known configuration in 913 (linear), 858 (state 1) and 781 (state
2) of 1000 repetitions (docs/MONTECARLO_LP_VALIDATION.md).

Nothing is estimated. The US fit is checked against the US series and every
simulated fit against its own series before it is used (applications/
abm_system_lp.h), so a stale file counts as missing rather than as data. A run
with a missing or non-finite response is a hole, and et_al's mcs refuses a
table with one, so a replicate with a hole in any column is dropped from that
model's table and counted in the manifest.

Requires out/us_lp_fit.npz, out/abm_system_fit_lp/ and dataset/abm_system/.
Writes out/abm_system_lp_irf_loss_<model>.csv and
out/abm_system_lp_irf_loss_manifest.txt. Nothing printed.

With --r-levels, the same on R's levels, the collaborator's transformation:
reads out/us_lp_fit_r_levels.npz and out/abm_system_fit_lp_r_levels/, writes
out/abm_system_lp_r_levels_irf_loss_<model>.csv and its manifest.
*/

#include "us_data.h"
#include "abm_system_lp.h"
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <math.h>


static Mat as_row(mreal *values) { return (Mat){ 1, LP_RESPONSE_DIM, LP_RESPONSE_DIM, values }; }

int main(int argc, char **argv) {
    LpLayout layout = lp_parse_flags(argc, argv, NULL);
    openblas_set_num_threads(1);
    char us_fit_path[256], output_stem[256], fit_dir[256];
    snprintf(us_fit_path, sizeof us_fit_path, "out/us_lp_fit%s.npz", lp_layout_suffix(layout));
    snprintf(output_stem, sizeof output_stem, "out/abm_system_lp%s_irf_loss", lp_layout_suffix(layout));
    snprintf(fit_dir, sizeof fit_dir, "%s%s", LP_FIT_DIR_DEFAULT, lp_layout_suffix(layout));

    /* The US responses, refused unless fitted on exactly the US series. */
    static mreal reference[LP_N_LOSSES][LP_RESPONSE_DIM];
    {
        Mat original = load_us_system();
        Mat stored = us_system_growth_block(original);
        Mat y, switching;
        lp_system_series(stored, layout, &y, &switching);
        int count = 0;
        LpSystemRecord *us = lp_system_read_archive(us_fit_path, lp_system_spec(), &count);
        assert(us && count == 1 && lp_system_record_matches(&us[0], y, switching) && us[0].has_responses
               && "abm_system_lp_irf_loss: out/us_lp_fit.npz is missing or stale - run bin/us_lp_fit");
        for (int loss = 0; loss < LP_N_LOSSES; loss++)
            memcpy(reference[loss], lp_response(&us[0], loss), LP_RESPONSE_DIM * sizeof(mreal));
        free(us);
        mat_free(y); mat_free(switching); mat_free(stored); mat_free(original);
    }

    int n_samples;
    LpSample *samples = lp_list_samples(fit_dir, &n_samples);
    assert(n_samples > 0 && "abm_system_lp_irf_loss: no configurations in out/abm_system_fit_lp/");
    int n_replicates = lp_count_replicates(LP_INPUT_DIR_DEFAULT, samples[0].name);

    Mat values[LP_N_LOSSES];
    for (int loss = 0; loss < LP_N_LOSSES; loss++) values[loss] = mat_new(n_replicates, n_samples);

    #pragma omp parallel for schedule(dynamic)
    for (int col = 0; col < n_samples; col++) {
        LpConfiguration c = lp_configuration_load(fit_dir, LP_INPUT_DIR_DEFAULT, samples[col].name,
                                                  n_replicates, layout);
        for (int row = 0; row < n_replicates; row++)
            for (int loss = 0; loss < LP_N_LOSSES; loss++)
                AT(values[loss], row, col) = lp_configuration_ok(&c, row, loss)
                    ? stats_mae(as_row(reference[loss]), as_row((mreal *)lp_response(c.row[row], loss)))
                    : (mreal)NAN;
        lp_configuration_free(&c);
    }

    char manifest_path[320];
    snprintf(manifest_path, sizeof manifest_path, "%s_manifest.txt", output_stem);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "abm_system_lp_irf_loss: cannot open the manifest");
    fprintf(manifest, "%d configurations, %d replicates each, loss = mean absolute difference between the US "
                      "local-projection responses and each run's (%d values, horizons 0 to %d)\n",
            n_samples, n_replicates, LP_RESPONSE_DIM, LP_SYSTEM_HOR);
    fprintf(manifest, "series: %s\n\n", lp_layout_description(layout));

    for (int loss = 0; loss < LP_N_LOSSES; loss++) {
        int *usable = malloc((size_t)n_replicates * sizeof(int));
        long missing = 0;
        int keep = 0;
        for (int row = 0; row < n_replicates; row++) {
            usable[row] = 1;
            for (int col = 0; col < n_samples; col++)
                if (MISNAN(AT(values[loss], row, col))) {
                    missing++;
                    usable[row] = 0;
                    fprintf(manifest, "%s missing: %s replicate %03d\n", lp_loss_name(loss), samples[col].name, row);
                }
            keep += usable[row];
        }
        assert(keep > 0 && "abm_system_lp_irf_loss: every replicate has a missing cell");

        Mat table = mat_new(keep, n_samples + 1);
        for (int row = 0, at = 0; row < n_replicates; row++) {
            if (!usable[row]) continue;
            AT(table, at, 0) = (mreal)row;
            for (int col = 0; col < n_samples; col++) AT(table, at, col + 1) = AT(values[loss], row, col);
            at++;
        }
        char **names = malloc((size_t)(n_samples + 1) * sizeof(char *));
        names[0] = frame_strdup("replicate");
        for (int col = 0; col < n_samples; col++) {
            char buffer[128];
            snprintf(buffer, sizeof buffer, "%s_lp_%s", samples[col].name, lp_loss_name(loss));
            names[col + 1] = frame_strdup(buffer);
        }
        DataFrame frame = df_from_matrix(table, (const char *const *)names);
        char path[320];
        snprintf(path, sizeof path, "%s_%s.csv", output_stem, lp_loss_name(loss));
        df_write_csv(&frame, path, csv_write_options_default());
        fprintf(manifest, "%s: %ld of %ld cells missing, %d of %d replicates kept\n\n", lp_loss_name(loss),
                missing, (long)n_replicates * n_samples, keep, n_replicates);

        df_free(&frame);
        for (int col = 0; col <= n_samples; col++) free(names[col]);
        free(names);
        mat_free(table);
        free(usable);
        mat_free(values[loss]);
    }
    fclose(manifest);
    lp_free_samples(samples, n_samples);
    return 0;
}
