/*
Discarded from both pipelines with the local projections on R's levels:
docs/OLD_EXPERIMENTS.md. Whether old/lp_r_levels/lp_system_r_levels.h builds
R's transform_data(), and whether the state series the live fits use,
applications/lp_system.h's lp_system_switching, is exactly the GDP_MA4 of that
transformation.

lp_system_r_transform was compared once against transform_data() itself, run
unmodified in R on four stored replicates, and agreed to 3.7e-15 relative.
This check is what keeps that agreement without R: it compares the header
against a reference written here as plain loops, the R code read one line at a
time.

1. The reference on a series with random growth rates, every element of every
   row and of the switching series, including the three periods dropped. Then
   the fitted series: the stored rows from period 3 on, exactly, and R's
   switching series, exactly.
2. A strided input, a column window of a wider block, against its copy, since
   the header hands the stored block to mat_cumsum as it is.
3. A cumulative GDP growth below -init: R's log gives NaN there, and so must
   the header.

Writes out/old/lp_system_r_levels_transform_report.txt. Exits 1 on any failure.
*/

#include "old/lp_r_levels/lp_system_r_levels.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define REPORT "out/old/lp_system_r_levels_transform_report.txt"
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

    fprintf(report, "\n%s\n", failures ? "FAILED" : "all passed");
    fclose(report);

    mat_free(stored); mat_free(y); mat_free(switching); mat_free(y_ref); mat_free(switching_ref);
    mat_free(wide); mat_free(y_view); mat_free(switching_view);
    mat_free(falling); mat_free(y_fall); mat_free(switching_fall);
    if (failures) fprintf(stderr, "lp_system_r_levels_transform: %d check(s) failed, see " REPORT "\n", failures);
    return failures ? 1 : 0;
}
