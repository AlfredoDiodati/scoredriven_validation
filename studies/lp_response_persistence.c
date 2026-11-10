/*
Whether the impulse responses of the local projections die out within the 15
quarters they are estimated over, or are still large at the end.

For one response, meaning one model (lin, s1 or s2), one responding series k
and one shocked series j, the persistence ratio is

    mean over h = 12..15 of |irf(k, h, j)| / max over h = 0..15 of |irf(k, h, j)|

so 0 means the response is gone by the last year and 1 means it is as large in
the last year as at its peak. A response counts as persistent when the ratio
is above PERSISTENT_ABOVE. A response whose peak is exactly 0 has no ratio and
is left out and counted.

The responses are the ones applications/abm_system_fit_lp.c and
applications/abm_system_fit_lp_t188.c stored for every replicate of
out/dataset/abm_system, the full-length simulations and the same simulations
cut to their first 188 periods, and the one applications/us_lp_fit.c stored for
the US data. Nothing is refitted; a replicate without responses or with a
non-finite one is left out and counted.

The series are growth rates (GDP, energy, employment), inflation and the
interest rate in levels, so a persistent response of a growth rate is a
permanent drift in the corresponding level.

Writes out/studies/lp_response_persistence.txt.
*/

#include "applications/abm_system_lp.h"
#include <et_al./stats.h>
#include <omp.h>

#define PERSISTENT_ABOVE 0.5
#define TAIL_FIRST_HORIZON 12
#define N_CONFIGURATIONS 1000
#define N_REPLICATES 1000
#define N_PAIRS (LP_SYSTEM_K * LP_SYSTEM_K)

static const char *const series_name[LP_SYSTEM_K] = {
    "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
};

/* The ratio of one response, or NAN when its peak is 0. */
static double persistence_ratio(const mreal *irf, int k, int j) {
    double peak = 0, tail = 0;
    for (int h = 0; h <= LP_SYSTEM_HOR; h++) {
        double value = fabs((double)irf[(k * (LP_SYSTEM_HOR + 1) + h) * LP_SYSTEM_K + j]);
        if (value > peak) peak = value;
        if (h >= TAIL_FIRST_HORIZON) tail += value;
    }
    if (peak == 0) return NAN;
    return tail / (LP_SYSTEM_HOR + 1 - TAIL_FIRST_HORIZON) / peak;
}

/* Share above the threshold and quartiles of the finite entries of values,
   which is reordered. */
static void summarize(FILE *out, const char *label, mreal *values, size_t n) {
    size_t kept = 0, persistent = 0;
    for (size_t i = 0; i < n; i++)
        if (!MISNAN(values[i])) {
            if (values[i] > PERSISTENT_ABOVE) persistent++;
            values[kept++] = values[i];
        }
    if (kept == 0) {
        fprintf(out, "%-40s no responses\n", label);
        return;
    }
    double q25 = (double)stats_quantile_inplace(values, (int)kept, (mreal)0.25);
    double q50 = (double)stats_quantile_inplace(values, (int)kept, (mreal)0.50);
    double q75 = (double)stats_quantile_inplace(values, (int)kept, (mreal)0.75);
    fprintf(out, "%-40s persistent %5.1f%%  ratio quartiles %.2f %.2f %.2f  (n %zu, no ratio %zu)\n",
            label, 100.0 * (double)persistent / (double)kept, q25, q50, q75, kept, n - kept);
}

static void study_archive(FILE *out, const char *fit_dir, const char *title) {
    size_t n_slots = (size_t)N_CONFIGURATIONS * N_REPLICATES;
    /* ratio[((model * N_PAIRS + pair) * n_slots) + slot], NAN where missing */
    mreal *ratio = malloc((size_t)LP_N_LOSSES * N_PAIRS * n_slots * sizeof(mreal));
    assert(ratio && "lp_response_persistence: out of memory");
    for (size_t i = 0; i < (size_t)LP_N_LOSSES * N_PAIRS * n_slots; i++) ratio[i] = (mreal)NAN;

    size_t missing = 0, non_finite = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:missing, non_finite)
    for (int c = 0; c < N_CONFIGURATIONS; c++) {
        char path[640];
        snprintf(path, sizeof path, "%s/cop_%04d.npz", fit_dir, c + 1);
        int n_records = 0;
        LpSystemRecord *records = lp_system_read_archive(path, lp_system_spec(), &n_records);
        unsigned char seen[N_REPLICATES] = { 0 };
        for (int i = 0; i < n_records; i++) {
            const LpSystemRecord *r = &records[i];
            if (r->replicate < 0 || r->replicate >= N_REPLICATES || !r->has_responses) continue;
            seen[r->replicate] = 1;
            size_t slot = (size_t)c * N_REPLICATES + r->replicate;
            for (int model = 0; model < LP_N_LOSSES; model++) {
                const mreal *irf = lp_response(r, model);
                if (!_lp_all_finite(irf)) {
                    non_finite++;
                    continue;
                }
                for (int k = 0; k < LP_SYSTEM_K; k++)
                    for (int j = 0; j < LP_SYSTEM_K; j++)
                        ratio[((size_t)model * N_PAIRS + k * LP_SYSTEM_K + j) * n_slots + slot] =
                            (mreal)persistence_ratio(irf, k, j);
            }
        }
        for (int i = 0; i < N_REPLICATES; i++) missing += !seen[i];
        free(records);
    }

    fprintf(out, "%s\n", title);
    fprintf(out, "fits read from %s, %d configurations x %d replicates; without responses %zu, "
                 "model fits with a non-finite response %zu\n\n", fit_dir, N_CONFIGURATIONS, N_REPLICATES,
            missing, non_finite);

    size_t pooled_n = (size_t)N_PAIRS * n_slots;
    mreal *pooled = malloc(pooled_n * sizeof(mreal));
    assert(pooled && "lp_response_persistence: out of memory");
    char label[128];
    for (int model = 0; model < LP_N_LOSSES; model++) {
        mreal *block = ratio + (size_t)model * N_PAIRS * n_slots;
        memcpy(pooled, block, pooled_n * sizeof(mreal));
        snprintf(label, sizeof label, "%s, all 25 responses pooled", lp_loss_name(model));
        summarize(out, label, pooled, pooled_n);

        size_t own_n = 0;
        for (int k = 0; k < LP_SYSTEM_K; k++) {
            memcpy(pooled + own_n, block + (size_t)(k * LP_SYSTEM_K + k) * n_slots, n_slots * sizeof(mreal));
            own_n += n_slots;
        }
        snprintf(label, sizeof label, "%s, own responses pooled", lp_loss_name(model));
        summarize(out, label, pooled, own_n);

        for (int k = 0; k < LP_SYSTEM_K; k++) {
            memcpy(pooled, block + (size_t)(k * LP_SYSTEM_K + k) * n_slots, n_slots * sizeof(mreal));
            snprintf(label, sizeof label, "%s, %s to own shock", lp_loss_name(model), series_name[k]);
            summarize(out, label, pooled, n_slots);
        }
        fprintf(out, "\n");
    }
    free(pooled);
    free(ratio);
}

static void study_us(FILE *out) {
    int count = 0;
    LpSystemRecord *record = lp_system_read_archive("out/applications/us_lp_fit.npz", lp_system_spec(), &count);
    assert(record && count == 1 && record->has_responses && "lp_response_persistence: run app-us_lp_fit first");
    fprintf(out, "US data, out/applications/us_lp_fit.npz, one fit\n\n");
    for (int model = 0; model < LP_N_LOSSES; model++) {
        const mreal *irf = lp_response(record, model);
        int persistent = 0, own_persistent = 0;
        for (int k = 0; k < LP_SYSTEM_K; k++)
            for (int j = 0; j < LP_SYSTEM_K; j++) {
                double value = persistence_ratio(irf, k, j);
                if (value > PERSISTENT_ABOVE) {
                    persistent++;
                    if (k == j) own_persistent++;
                }
            }
        fprintf(out, "%s: persistent %d of 25 responses, %d of 5 own responses; own ratios", lp_loss_name(model),
                persistent, own_persistent);
        for (int k = 0; k < LP_SYSTEM_K; k++) fprintf(out, " %s %.2f", series_name[k], persistence_ratio(irf, k, k));
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
    free(record);
}

int main(void) {
    frame_mkdir_p("out/studies");
    const char *path = "out/studies/lp_response_persistence.txt";
    FILE *out = fopen(path, "w");
    assert(out && "lp_response_persistence: cannot open the report");

    fprintf(out, "Persistence of the local projections' impulse responses\n\n");
    fprintf(out, "ratio of one response: mean |irf| over horizons %d to %d divided by max |irf| over horizons 0 to %d.\n"
                 "persistent: ratio above %.2f. Quartiles are over replicates (and responses, where pooled).\n"
                 "Own response: series k to a unit shock in series k, which is 1 at horizon 0.\n\n",
            TAIL_FIRST_HORIZON, LP_SYSTEM_HOR, LP_SYSTEM_HOR, PERSISTENT_ABOVE);

    study_us(out);
    study_archive(out, "out/applications/abm_system_fit_lp", "Simulations, full length");
    fflush(out);
    study_archive(out, "out/applications/abm_system_fit_lp_t188",
                  "Simulations cut to their first 188 periods, the US sample length");
    fclose(out);
    printf("wrote %s\n", path);
    return 0;
}
