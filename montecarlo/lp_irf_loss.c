/*
The impulse-response loss matrices of the Monte Carlo experiment under the
local projections of the collaborator's R pipeline, one per model:

    lin   the linear local projection, irf_lin_mean
    s1    the state-dependent one in state 1 (lags times 1 - F), irf_s1_mean
    s2    the same in state 2 (lags times F), irf_s2_mean

The procedure is montecarlo/irf_loss.c's with the auxiliary model changed.
The benchmark is the simulated run montecarlo/benchmark_choice.c chose, read
back through montecarlo/benchmark.h, so both auxiliary models are asked to
recover the same known answer from the same series. Its replicate index is
held out of every configuration's column, for the reason benchmark.h gives:
replicate n of every configuration is seed n+1 of the model.

A cell's loss is the mean absolute error between the benchmark's response
vector and the cell's, all K x (hor + 1) x K = 400 entries of one model, as in
the t-QVARMA protocol. The R pipeline uses the mean squared error instead; the
absolute error is kept for comparability with the t-QVARMA results and for the
reason montecarlo/irf_loss.c gives.

Nothing is estimated here. Every fit is read from the per-configuration
archives in out/abm_system_fit_lp/, which applications/abm_system_fit_lp.c
writes, and checked against the cell's own data on the way in
(applications/abm_system_lp.h).

A cell whose fit is missing, stale or non-finite is a hole, and the
confidence set refuses a matrix with one, so a replicate with a hole in any
column is dropped from that model's matrix and counted in the manifest.
Dropping the row rather than the configuration keeps every configuration in
the comparison.

Requires out/abm_system_fit_lp/, dataset/abm_system/ and
montecarlo/out/benchmark.env. Writes montecarlo/out/lp_irf_loss_<model>.csv,
one per model, and montecarlo/out/lp_irf_loss_manifest.txt. Nothing printed.
*/

#include "applications/abm_system_lp.h"
#include "montecarlo/benchmark.h"
#include <et_al./stats.h>
#include <et_al./frame/csv.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <math.h>

#define OUTPUT_STEM "montecarlo/out/lp_irf_loss"

static const char *FIT_DIR;
static const char *INPUT_DIR;

/* A vector of LP_RESPONSE_DIM entries as a 1 x n view, for stats_mae. */
static Mat as_row(mreal *values) { return (Mat){ 1, LP_RESPONSE_DIM, LP_RESPONSE_DIM, values }; }

int main(void) {
    FIT_DIR = getenv("ABM_SYSTEM_LP_FIT_DIR");
    if (!FIT_DIR) FIT_DIR = LP_FIT_DIR_DEFAULT;
    INPUT_DIR = getenv("ABM_SYSTEM_INPUT_DIR");
    if (!INPUT_DIR) INPUT_DIR = LP_INPUT_DIR_DEFAULT;
    openblas_set_num_threads(1);

    Benchmark benchmark = benchmark_read();

    int n_samples;
    LpSample *samples = lp_list_samples(FIT_DIR, &n_samples);
    assert(n_samples > 0 && "lp_irf_loss: no configurations in the LP fit cache");
    int n_replicates = lp_count_replicates(INPUT_DIR, samples[0].name);
    assert(benchmark.replicate < n_replicates && "lp_irf_loss: the benchmark replicate is out of range");

    static mreal reference[LP_N_LOSSES][LP_RESPONSE_DIM];
    {
        LpConfiguration own = lp_configuration_load(FIT_DIR, INPUT_DIR, benchmark.sample, n_replicates, LP_LAYOUT_GROWTH);
        for (int loss = 0; loss < LP_N_LOSSES; loss++) {
            assert(lp_configuration_ok(&own, benchmark.replicate, loss)
                   && "lp_irf_loss: the benchmark's own LP fit is missing or not finite");
            memcpy(reference[loss], lp_response(own.row[benchmark.replicate], loss),
                   LP_RESPONSE_DIM * sizeof(mreal));
        }
        lp_configuration_free(&own);
    }

    /* values[loss] is replicate x configuration, NaN where the cell has no loss. */
    Mat values[LP_N_LOSSES];
    for (int loss = 0; loss < LP_N_LOSSES; loss++) values[loss] = mat_new(n_replicates, n_samples);
    long n_missing[LP_N_LOSSES] = { 0, 0, 0 };

    /* One task per configuration. Each task writes only its own column, so no
       cell is shared between threads. */
    #pragma omp parallel for schedule(dynamic)
    for (int col = 0; col < n_samples; col++) {
        LpConfiguration c = lp_configuration_load(FIT_DIR, INPUT_DIR, samples[col].name, n_replicates, LP_LAYOUT_GROWTH);
        for (int row = 0; row < n_replicates; row++)
            for (int loss = 0; loss < LP_N_LOSSES; loss++) {
                mreal value = (mreal)NAN;
                if (row != benchmark.replicate && lp_configuration_ok(&c, row, loss))
                    value = stats_mae(as_row(reference[loss]), as_row((mreal *)lp_response(c.row[row], loss)));
                AT(values[loss], row, col) = value;
            }
        lp_configuration_free(&c);
    }

    char manifest_path[256];
    snprintf(manifest_path, sizeof manifest_path, "%s_manifest.txt", OUTPUT_STEM);
    FILE *manifest = fopen(manifest_path, "w");
    assert(manifest && "lp_irf_loss: cannot open the manifest for writing");
    fprintf(manifest, "%d configurations, %d replicates each, loss = MAE between stacked "
                      "local-projection impulse responses (%d entries, horizons 0 to %d)\n",
            n_samples, n_replicates, LP_RESPONSE_DIM, LP_SYSTEM_HOR);
    fprintf(manifest, "benchmark: %s replicate %03d, its own cached LP fits and its own series; "
                      "replicate %03d is held out of every column\n\n",
            benchmark.sample, benchmark.replicate, benchmark.replicate);

    for (int loss = 0; loss < LP_N_LOSSES; loss++) {
        int *usable = malloc((size_t)n_replicates * sizeof(int));
        int keep = 0;
        for (int row = 0; row < n_replicates; row++) {
            usable[row] = row != benchmark.replicate;
            if (!usable[row]) continue;
            for (int col = 0; col < n_samples; col++)
                if (MISNAN(AT(values[loss], row, col))) {
                    n_missing[loss]++;
                    fprintf(manifest, "%s missing: %s replicate %03d\n", lp_loss_name(loss),
                            samples[col].name, row);
                    usable[row] = 0;
                }
            keep += usable[row];
        }
        assert(keep > 0 && "lp_irf_loss: every replicate has a missing cell");

        Mat table = mat_new(keep, n_samples + 1);
        int at = 0;
        for (int row = 0; row < n_replicates; row++) {
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
        char path[256];
        snprintf(path, sizeof path, "%s_%s.csv", OUTPUT_STEM, lp_loss_name(loss));
        df_write_csv(&frame, path, csv_write_options_default());

        fprintf(manifest, "%s: %ld of %ld cells missing, %d of %d replicates kept "
                          "(the benchmark's own held out, %d dropped for a missing cell)\n\n",
                lp_loss_name(loss), n_missing[loss], (long)(n_replicates - 1) * n_samples, keep,
                n_replicates, n_replicates - 1 - keep);

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
