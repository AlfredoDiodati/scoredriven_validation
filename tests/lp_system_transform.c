/*
Whether applications/lp_system.h builds the series it claims to: R's
transform_data() in lp_system_r_transform, and the layout the fits use in
lp_system_from_stored.

Every local-projection fit, on the simulated data now and on the real data
later, sees its series through that header, so a wrong row order, a
moving average aligned one period off, or a cumulative sum started in the
wrong place changes every response without any downstream output noticing.

lp_system_r_transform was compared once against transform_data() itself, run
unmodified in R on four stored replicates, and agreed to 3.7e-15 relative.
This test is what keeps that agreement without R: it compares the header
against a reference written here as plain loops, the R code read one line at a
time.

1. The reference on a series with random growth rates, every element of every
   row and of the switching series, including the three periods dropped. Then
   the fitted layout: the stored rows from period 3 on, exactly, and R's
   switching series, exactly.
2. A strided input, a column window of a wider block, against its copy, since
   the header hands the stored block to mat_cumsum as it is.
3. A cumulative GDP growth below -init: R's log gives NaN there, and so must
   the header, so that lp.h reports the replicate rather than fitting it.
4. The specification: the values R's estimate() passes to lpirfs.
5. The archive the fits are stored in: rows written and read back bit for bit,
   a fit without d stored as a row without responses, a different
   specification and a missing file both read as a miss, and a row refused for
   series other than the ones it was fitted on.

Writes out/lp_system_transform_report.txt, and a scratch archive under out/
that it removes on the way out. Exits 1 on any failure.
*/

#include "applications/lp_system.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define REPORT "out/lp_system_transform_report.txt"
#define SCRATCH_ARCHIVE "out/lp_system_transform_scratch.npz"
#define T_STORED 400

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) failures++;
}

/* transform_data() written out one period at a time, in double. */
static void reference(Mat stored, Mat *y, Mat *switching) {
    int T = stored.c, first = LP_SYSTEM_SWITCHING_WINDOW - 1, n = T - first;
    double *gdp = malloc((size_t)T * sizeof(double));
    double gdp_sum = 0, energy_sum = 0, employment_sum = 0;
    *y = mat_new(LP_SYSTEM_K, n);
    *switching = mat_new(n, 1);
    for (int t = 0; t < T; t++) {
        gdp_sum += (double)AT(stored, STORED_GDP_GROWTH, t);
        energy_sum += (double)AT(stored, STORED_EN_GROWTH, t);
        employment_sum += (double)AT(stored, STORED_EMPLOYMENT_CHANGE, t);
        gdp[t] = log(100.0 + gdp_sum) * 100.0;
        if (t < first) continue;
        int c = t - first;
        AT(*y, 0, c) = (mreal)gdp[t];
        AT(*y, 1, c) = (mreal)(0.90 + employment_sum);
        AT(*y, 2, c) = AT(stored, STORED_INFLATION, t);
        AT(*y, 3, c) = AT(stored, STORED_INTEREST_RATE, t);
        AT(*y, 4, c) = (mreal)(log(100.0 + energy_sum) * 100.0);
        AT(*switching, c, 0) = (mreal)((gdp[t] + gdp[t - 1] + gdp[t - 2] + gdp[t - 3]) / 4.0);
    }
    free(gdp);
}

static double largest_relative_gap(Mat a, Mat b) {
    double worst = 0;
    for (int i = 0; i < a.r; i++)
        for (int j = 0; j < a.c; j++) {
            double gap = fabs((double)AT(a, i, j) - (double)AT(b, i, j))
                       / fmax(1.0, fabs((double)AT(b, i, j)));
            if (gap > worst) worst = gap;
        }
    return worst;
}

/* Growth rates on the stored scales: GDP and energy near 0.5 per period with
   unit noise, employment changes near 0, inflation near 0.5, the rate near 1. */
static Mat stored_block(Rng *rng, int T) {
    Mat stored = mat_new(STORED_K, T);
    for (int t = 0; t < T; t++) {
        AT(stored, STORED_GDP_GROWTH, t) = (mreal)(0.5 + rng_normal(rng));
        AT(stored, STORED_EN_GROWTH, t) = (mreal)(0.1 + rng_normal(rng));
        AT(stored, STORED_EMPLOYMENT_CHANGE, t) = (mreal)(0.3 * rng_normal(rng));
        AT(stored, STORED_INFLATION, t) = (mreal)(0.5 + 0.2 * rng_normal(rng));
        AT(stored, STORED_INTEREST_RATE, t) = (mreal)(1.0 + 0.3 * rng_normal(rng));
    }
    return stored;
}

int main(void) {
    FILE *report = fopen(REPORT, "w");
    if (!report) { fprintf(stderr, "cannot open " REPORT "\n"); return 1; }
    Rng rng = rng_new(2026, 0);
    const double tolerance = 64 * MEPS;

    Mat stored = stored_block(&rng, T_STORED);
    Mat y, switching, y_ref, switching_ref;
    lp_system_r_transform(stored, &y, &switching);
    reference(stored, &y_ref, &switching_ref);
    check(report, y.r == 5 && y.c == T_STORED - 3 && switching.r == T_STORED - 3 && switching.c == 1,
          "1. shape: 5 x 397 and 397 x 1 from 400 stored periods");
    double gap_y = largest_relative_gap(y, y_ref), gap_s = largest_relative_gap(switching, switching_ref);
    fprintf(report, "      largest relative gap: series %.2e, switching %.2e, tolerance %.2e\n", gap_y,
            gap_s, tolerance);
    check(report, gap_y <= tolerance && gap_s <= tolerance, "1. every element against the reference");

    Mat fitted_layout, fitted_switching_series;
    lp_system_from_stored(stored, &fitted_layout, &fitted_switching_series);
    int layout_exact = fitted_layout.r == STORED_K && fitted_layout.c == T_STORED - 3;
    for (int k = 0; layout_exact && k < STORED_K; k++)
        for (int c = 0; c < fitted_layout.c; c++)
            if (AT(fitted_layout, k, c) != AT(stored, k, c + 3)) { layout_exact = 0; break; }
    check(report, layout_exact && largest_relative_gap(fitted_switching_series, switching) == 0,
          "1. fitted layout: the stored rows from period 3 on and R's switching series, exactly");
    mat_free(fitted_layout);
    mat_free(fitted_switching_series);

    Mat wide = mat_new(STORED_K, T_STORED + 7);
    for (int k = 0; k < STORED_K; k++)
        for (int t = 0; t < T_STORED; t++) AT(wide, k, t + 5) = AT(stored, k, t);
    Mat view = mat_slice(wide, 0, STORED_K, 5, 5 + T_STORED);
    Mat y_view, switching_view;
    lp_system_r_transform(view, &y_view, &switching_view);
    check(report, largest_relative_gap(y_view, y) == 0 && largest_relative_gap(switching_view, switching) == 0,
          "2. a strided input gives exactly what its copy gives");

    Mat falling = mat_copy(stored);
    for (int t = 0; t < T_STORED; t++) AT(falling, STORED_GDP_GROWTH, t) = (mreal)-1;
    Mat y_fall, switching_fall;
    lp_system_r_transform(falling, &y_fall, &switching_fall);
    /* Cumulated growth reaches -100 at period 99 (log of 0) and passes it at 100. */
    int nan_after = 1, finite_before = 1;
    for (int c = 0; c < y_fall.c; c++) {
        int t = c + 3;
        mreal v = AT(y_fall, LP_R_ROW_GDP, c);
        if (t >= 100 && !MISNAN(v)) nan_after = 0;
        if (t < 99 && (MISNAN(v) || MISINF(v))) finite_before = 0;
    }
    check(report, nan_after && finite_before && !mat_all_finite(switching_fall),
          "3. cumulated GDP growth below -100 gives NaN, as R's log does, and reaches the switching series");

    LpNlSpec spec = lp_system_spec();
    check(report, spec.lin.K == 5 && spec.lin.lags_endog_lin == 4 && spec.lags_endog_nl == 4
                  && spec.lin.hor == 15 && spec.lin.shock_type == LP_SHOCK_UNIT && spec.use_logistic == 1
                  && spec.use_hp == 1 && spec.lambda == 1600 && spec.gamma == 2 && spec.lag_switching == 1,
          "4. spec: 4 lags in both parts, hor 15, unit shock, HP at 1600, gamma 2, lagged switching");

    /* Two fitted replicates and one whose data hold a NaN, so it has no d. */
    LpSystemRecord written[3];
    Mat fitted_y[3], fitted_switching[3];
    for (int i = 0; i < 3; i++) {
        Mat block = stored_block(&rng, T_STORED);
        if (i == 2) AT(block, STORED_INFLATION, 50) = (mreal)NAN;
        lp_system_from_stored(block, &fitted_y[i], &fitted_switching[i]);
        LpLinNlFit fit = lp_lin_and_nl(fitted_y[i], fitted_switching[i], spec, (LpBands){0});
        written[i] = lp_system_record(&fit, fitted_y[i], fitted_switching[i], 10 + i);
        lp_lin_nl_fit_free(&fit);
        mat_free(block);
    }
    check(report, written[0].has_responses && written[1].has_responses && !written[2].has_responses,
          "5. a fit on data holding a NaN has no responses, the others do");
    lp_system_write_archive(SCRATCH_ARCHIVE, written, 3, spec);
    int n_read = 0;
    LpSystemRecord *read_back = lp_system_read_archive(SCRATCH_ARCHIVE, spec, &n_read);
    int identical = read_back && n_read == 3;
    for (int i = 0; identical && i < 3; i++) {
        const LpSystemRecord *a = &written[i], *b = &read_back[i];
        identical = a->replicate == b->replicate && a->has_responses == b->has_responses
                 && a->var_ols_status == b->var_ols_status && a->var_rank == b->var_rank
                 && a->var_chol_status == b->var_chol_status
                 && a->switching_is_constant == b->switching_is_constant
                 && memcmp(a->lin_ols_status, b->lin_ols_status, sizeof a->lin_ols_status) == 0
                 && memcmp(a->nl_ols_status, b->nl_ols_status, sizeof a->nl_ols_status) == 0
                 && a->data_fingerprint == b->data_fingerprint
                 && a->switching_fingerprint == b->switching_fingerprint;
        /* A NaN never equals itself, so the responses are compared by bit
           pattern, which also holds the finite values to exact equality. */
        identical = identical && memcmp(a->d, b->d, sizeof a->d) == 0
                 && memcmp(a->lin, b->lin, sizeof a->lin) == 0
                 && memcmp(a->s1, b->s1, sizeof a->s1) == 0
                 && memcmp(a->s2, b->s2, sizeof a->s2) == 0;
    }
    check(report, identical, "5. every field of every row read back bit for bit");
    check(report, read_back && lp_system_record_matches(&read_back[0], fitted_y[0], fitted_switching[0])
                  && !lp_system_record_matches(&read_back[0], fitted_y[1], fitted_switching[1]),
          "5. a row matches the series it was fitted on and refuses another replicate's");
    free(read_back);

    LpNlSpec other = spec;
    other.gamma = 3;
    int n_other = -1;
    LpSystemRecord *other_read = lp_system_read_archive(SCRATCH_ARCHIVE, other, &n_other);
    check(report, other_read == NULL && n_other == 0, "5. an archive written for another gamma reads as a miss");
    free(other_read);
    remove(SCRATCH_ARCHIVE);
    int n_missing = -1;
    LpSystemRecord *missing = lp_system_read_archive(SCRATCH_ARCHIVE, spec, &n_missing);
    check(report, missing == NULL && n_missing == 0, "5. a missing archive reads as a miss");
    for (int i = 0; i < 3; i++) { mat_free(fitted_y[i]); mat_free(fitted_switching[i]); }

    fprintf(report, "\n%s\n", failures ? "FAILED" : "all passed");
    fclose(report);

    mat_free(stored); mat_free(y); mat_free(switching); mat_free(y_ref); mat_free(switching_ref);
    mat_free(wide); mat_free(y_view); mat_free(switching_view);
    mat_free(falling); mat_free(y_fall); mat_free(switching_fall);
    if (failures) fprintf(stderr, "lp_system_transform: %d check(s) failed, see " REPORT "\n", failures);
    return failures ? 1 : 0;
}
