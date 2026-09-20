/*
Whether applications/abm_system.h stores and returns what it says it does.

That header is the one place the layout of the simulated dataset is defined.
Two writers go through it - applications/abm_system_simulate.c from the model
and applications/abm_system_convert_rdata.c from the .Rdata files - and four readers
come back out of it: the fitting stage, the loss table, the winner's impulse
responses and every robustness study. A disagreement there about units, about
which row is which series, or about which period anchors the first difference
changes every number downstream, and no other test in this repository would see
it: the pipeline's own outputs are all computed through the same header, so they
agree with each other whatever it does.

Five questions, one section each. Nothing here runs the simulator or fits
anything; it takes a few hundredths of a second.

1. Does abm_system_transform compute the transformation its comment claims?
   Built on a series whose answer is known in closed form: a level growing at a
   constant factor g per period has 100 ln(g) as its growth row at every period,
   exactly, and a linear employment path has a constant difference.

2. Does it agree with the real data's own convention? The US series reach the
   auxiliary model through a different route - applications/us_prepare_data.c
   writes 100 ln(level) and applications/us_qvarma_spec_choice.c's build_block
   differences that - while the simulated series are transformed here in one
   step. The two are the same function composed in the opposite order, and the
   whole comparison the project makes rests on them agreeing. Checked directly,
   on the same raw levels, rather than argued.

3. Does the burn-in keep the periods it claims? ABM_SYSTEM_BURN_IN is 200 and a
   600-period run must come back with 400 periods, not 399: period 200 is spent
   anchoring period 201's difference rather than being reported.

4. Does an archive give back what was put into it? Write, read, compare every
   element. Including a short batch, which is what a configuration whose runs
   did not all succeed actually writes, since that is the case a reader is most
   likely to get wrong.

5. Does abm_system_list_replicates find every replicate, in ascending order,
   across several archives and with a whole batch missing?

Writes out/abm_system_layout_report.txt and scratches in
out/abm_system_layout_scratch/, which it removes on the way out.
*/

#include "applications/abm_system.h"

#include <et_al./linalg/mat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#define REPORT "out/abm_system_layout_report.txt"
#define SCRATCH "out/abm_system_layout_scratch"

static int failures = 0;
static FILE *report;

#define CHECK(cond, ...) do { \
    if (!(cond)) { fprintf(report, "  FAIL %s:%d: ", __FILE__, __LINE__); \
                   fprintf(report, __VA_ARGS__); fprintf(report, "\n"); \
                   printf("  FAIL %s:%d: ", __FILE__, __LINE__); \
                   printf(__VA_ARGS__); printf("\n"); failures++; } \
} while (0)

#define CHECK_NEAR(got, want, tol, label) do { \
    double got_value = (double)(got), want_value = (double)(want); \
    if (!(fabs(got_value - want_value) <= (double)(tol))) { \
        CHECK(0, "%s: got %.12g, want %.12g, tolerance %g", label, got_value, want_value, \
              (double)(tol)); } \
} while (0)

/* float32 carries about seven digits, so a difference of two logs of numbers
   near 1e5 keeps far fewer. The identities below are exact in real arithmetic
   and the tolerance is only about how much of that survives the build. */
#ifdef MAT_DOUBLE
#define TOLERANCE 1e-9
#else
#define TOLERANCE 2e-3
#endif

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_layout: mkdir failed");
}

/* A levels matrix in the units the model writes: GDP and energy growing at a
   constant factor, employment a 0-1 proportion on a straight line, the price
   level compounding, the interest rate a decimal proportion. Every row's
   transform is then known in closed form. */
#define GDP_GROWTH_FACTOR 1.004
#define ENERGY_GROWTH_FACTOR 0.9985
#define PRICE_GROWTH_FACTOR 1.0125
#define EMPLOYMENT_STEP 0.0002
#define EMPLOYMENT_FIRST 0.9
#define INTEREST_LEVEL 0.0175

static Mat known_levels(int periods) {
    Mat levels = mat_new(ABM_SYSTEM_K, periods);
    for (int t = 0; t < periods; t++) {
        AT(levels, LEVEL_GDP, t) = (mreal)(220000.0 * pow(GDP_GROWTH_FACTOR, t));
        AT(levels, LEVEL_ENERGY, t) = (mreal)(180000.0 * pow(ENERGY_GROWTH_FACTOR, t));
        AT(levels, LEVEL_EMPLOYMENT, t) = (mreal)(EMPLOYMENT_FIRST + EMPLOYMENT_STEP * t);
        AT(levels, LEVEL_PRICE, t) = (mreal)(1.3 * pow(PRICE_GROWTH_FACTOR, t));
        AT(levels, LEVEL_INTEREST, t) = (mreal)INTEREST_LEVEL;
    }
    return levels;
}

static void test_transform_is_what_it_claims(void) {
    const int periods = 40, start = 5;
    fprintf(report, "1. abm_system_transform against the closed form\n");
    printf("the transformation against its closed form\n");

    Mat levels = known_levels(periods);
    Mat y = abm_system_transform(levels, start);

    CHECK(y.r == ABM_SYSTEM_K, "transform must return %d rows, got %d", ABM_SYSTEM_K, y.r);
    CHECK(y.c == periods - start, "transform must return %d periods, got %d", periods - start, y.c);

    /* A constant growth factor g gives 100 ln(g) every period, and the level's
       own size drops out. Employment is a difference of a 0-1 proportion
       rescaled by 100, so a step of 0.0002 is 0.02 points. The interest rate is
       not differenced at all: it is 100 times the level. */
    double gdp = 100.0 * log(GDP_GROWTH_FACTOR);
    double energy = 100.0 * log(ENERGY_GROWTH_FACTOR);
    double price = 100.0 * log(PRICE_GROWTH_FACTOR);
    double employment = 100.0 * EMPLOYMENT_STEP;
    double interest = 100.0 * INTEREST_LEVEL;

    for (int t = 0; t < y.c; t++) {
        CHECK_NEAR(AT(y, ROW_GDP_GROWTH, t), gdp, TOLERANCE, "GDP growth");
        CHECK_NEAR(AT(y, ROW_EN_GROWTH, t), energy, TOLERANCE, "energy growth");
        CHECK_NEAR(AT(y, ROW_EMPLOYMENT_CHANGE, t), employment, TOLERANCE, "employment change");
        CHECK_NEAR(AT(y, ROW_INFLATION, t), price, TOLERANCE, "inflation");
        CHECK_NEAR(AT(y, ROW_INTEREST_RATE, t), interest, TOLERANCE, "interest rate");
    }

    fprintf(report, "   %d periods, constant-growth levels: GDP %.6f, energy %.6f, "
                    "inflation %.6f,\n   employment %.6f, interest %.6f, every period\n",
            y.c, gdp, energy, price, employment, interest);

    /* The first reported period differences against the one before it, so the
       period at index start - 1 is the anchor and is not itself reported. */
    double anchored = 100.0 * (log((double)AT(levels, LEVEL_GDP, start))
                               - log((double)AT(levels, LEVEL_GDP, start - 1)));
    CHECK_NEAR(AT(y, ROW_GDP_GROWTH, 0), anchored, TOLERANCE,
               "the first reported period anchors on the period before it");

    mat_free(y);
    mat_free(levels);
}

/* The route the US data takes: 100 ln(level) written by us_prepare_data.c, then
   differenced by us_qvarma_spec_choice.c's build_block. Employment there is
   already on a 0-100 scale and the interest rate already in percentage points,
   which is what the 100s below stand in for. */
static Mat real_data_route(Mat levels, int start) {
    int periods = levels.c;
    Mat prepared = mat_new(ABM_SYSTEM_K, periods);
    for (int t = 0; t < periods; t++) {
        AT(prepared, LEVEL_GDP, t) = (mreal)(100.0 * log((double)AT(levels, LEVEL_GDP, t)));
        AT(prepared, LEVEL_ENERGY, t) = (mreal)(100.0 * log((double)AT(levels, LEVEL_ENERGY, t)));
        AT(prepared, LEVEL_EMPLOYMENT, t) = (mreal)(100.0 * (double)AT(levels, LEVEL_EMPLOYMENT, t));
        AT(prepared, LEVEL_PRICE, t) = (mreal)(100.0 * log((double)AT(levels, LEVEL_PRICE, t)));
        AT(prepared, LEVEL_INTEREST, t) = (mreal)(100.0 * (double)AT(levels, LEVEL_INTEREST, t));
    }

    Mat y = mat_new(ABM_SYSTEM_K, periods - start);
    for (int t = start; t < periods; t++) {
        int c = t - start;
        AT(y, ROW_GDP_GROWTH, c) = AT(prepared, LEVEL_GDP, t) - AT(prepared, LEVEL_GDP, t - 1);
        AT(y, ROW_EN_GROWTH, c) = AT(prepared, LEVEL_ENERGY, t) - AT(prepared, LEVEL_ENERGY, t - 1);
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = AT(prepared, LEVEL_EMPLOYMENT, t)
                                        - AT(prepared, LEVEL_EMPLOYMENT, t - 1);
        AT(y, ROW_INFLATION, c) = AT(prepared, LEVEL_PRICE, t) - AT(prepared, LEVEL_PRICE, t - 1);
        AT(y, ROW_INTEREST_RATE, c) = AT(prepared, LEVEL_INTEREST, t);
    }
    mat_free(prepared);
    return y;
}

static void test_agrees_with_the_real_data_route(void) {
    const int periods = 60, start = 9;
    fprintf(report, "\n2. the simulated route against the real data's route\n");
    printf("the simulated route against the real data's route\n");

    Mat levels = known_levels(periods);
    Mat simulated = abm_system_transform(levels, start);
    Mat real = real_data_route(levels, start);

    CHECK(simulated.r == real.r && simulated.c == real.c,
          "the two routes must produce the same shape, got %dx%d and %dx%d",
          simulated.r, simulated.c, real.r, real.c);

    double worst = 0;
    for (int k = 0; k < ABM_SYSTEM_K; k++)
        for (int t = 0; t < simulated.c; t++) {
            double difference = fabs((double)AT(simulated, k, t) - (double)AT(real, k, t));
            if (difference > worst) worst = difference;
        }
    CHECK(worst <= TOLERANCE,
          "the two routes must agree: worst |difference| %.3g over %d elements",
          worst, ABM_SYSTEM_K * simulated.c);

    fprintf(report, "   %d elements, worst |difference| %.3g\n", ABM_SYSTEM_K * simulated.c, worst);

    mat_free(real);
    mat_free(simulated);
    mat_free(levels);
}

static void test_burn_in_keeps_the_periods_it_claims(void) {
    const int model_periods = 600;
    fprintf(report, "\n3. the burn-in\n");
    printf("the burn-in\n");

    Mat levels = known_levels(model_periods);
    Mat y = abm_system_transform(levels, ABM_SYSTEM_BURN_IN);

    /* 600 periods with the first 200 discarded leaves 400, not 399: period 200
       is the anchor for period 201's difference rather than a reported period. */
    CHECK(y.c == model_periods - ABM_SYSTEM_BURN_IN,
          "a %d-period run must give %d periods, got %d", model_periods,
          model_periods - ABM_SYSTEM_BURN_IN, y.c);

    fprintf(report, "   %d model periods, burn-in %d, reported %d\n", model_periods,
            ABM_SYSTEM_BURN_IN, y.c);

    mat_free(y);
    mat_free(levels);
}

/* One block per replicate, each distinguishable from the others so that a
   reader returning the wrong one is visible rather than plausible. */
static Mat marked_block(int replicate, int periods) {
    Mat y = mat_new(ABM_SYSTEM_K, periods);
    for (int k = 0; k < ABM_SYSTEM_K; k++)
        for (int t = 0; t < periods; t++)
            AT(y, k, t) = (mreal)(1000.0 * replicate + 10.0 * k + 0.25 * t);
    return y;
}

static int blocks_differ(Mat a, Mat b, double *worst) {
    *worst = 0;
    if (a.r != b.r || a.c != b.c) return 1;
    for (int k = 0; k < a.r; k++)
        for (int t = 0; t < a.c; t++) {
            double difference = fabs((double)AT(a, k, t) - (double)AT(b, k, t));
            if (difference > *worst) *worst = difference;
        }
    return *worst > 0;
}

/* count blocks starting at first_replicate, written as one archive. A count
   below ABM_SYSTEM_BATCH is a short batch, which is what a configuration whose
   runs did not all succeed leaves on disk. */
static void write_one_batch(const char *dir, int first_replicate, int count, int periods) {
    Mat block[ABM_SYSTEM_BATCH];
    int replicate[ABM_SYSTEM_BATCH];
    for (int b = 0; b < count; b++) {
        replicate[b] = first_replicate + b;
        block[b] = marked_block(replicate[b], periods);
    }
    abm_system_write_batch(dir, block, replicate, count);
    for (int b = 0; b < count; b++) mat_free(block[b]);
}

static void test_archive_round_trip(void) {
    const int periods = 12;
    fprintf(report, "\n4. write an archive, read it back\n");
    printf("the archive round trip\n");

    char dir[512];
    snprintf(dir, sizeof dir, "%s/round_trip", SCRATCH);
    make_directory(dir);

    /* A full batch, then a short one: the second is the case a reader is most
       likely to get wrong, since its archive holds fewer rows than its name
       suggests. */
    write_one_batch(dir, 0, ABM_SYSTEM_BATCH, periods);
    const int short_count = 3;
    write_one_batch(dir, ABM_SYSTEM_BATCH, short_count, periods);

    int checked = 0;
    double worst_overall = 0;
    for (int replicate = 0; replicate < ABM_SYSTEM_BATCH + short_count; replicate++) {
        Mat expected = marked_block(replicate, periods);
        Mat got = abm_system_read_replicate(dir, replicate);
        double worst;
        int differ = blocks_differ(expected, got, &worst);
        CHECK(!differ, "replicate %d must come back unchanged, worst |difference| %.3g "
                       "(shape %dx%d against %dx%d)",
              replicate, worst, got.r, got.c, expected.r, expected.c);
        if (worst > worst_overall) worst_overall = worst;
        checked++;
        mat_free(got);
        mat_free(expected);
    }

    fprintf(report, "   %d replicates over a full batch of %d and a short batch of %d, "
                    "%d periods each\n   worst |difference| %.3g\n",
            checked, ABM_SYSTEM_BATCH, short_count, periods, worst_overall);

    /* The archive a replicate lands in is its index divided by the batch size,
       so the short batch above is batch_001 and nothing else. */
    char path[512];
    abm_system_batch_path(path, sizeof path, dir, ABM_SYSTEM_BATCH + 1);
    char expected_path[560];
    snprintf(expected_path, sizeof expected_path, "%s/batch_001.npz", dir);
    CHECK(strcmp(path, expected_path) == 0, "batch path must be %s, got %s", expected_path, path);
}

static void test_listing_replicates(void) {
    const int periods = 6;
    fprintf(report, "\n5. listing the replicates an archive directory holds\n");
    printf("listing the replicates\n");

    char dir[512];
    snprintf(dir, sizeof dir, "%s/listing", SCRATCH);
    make_directory(dir);

    /* Batches 0 and 2 written, batch 1 left out: a configuration whose middle
       batch failed entirely, which the header says is simply absent rather than
       a hole to step over. */
    write_one_batch(dir, 0, ABM_SYSTEM_BATCH, periods);
    write_one_batch(dir, 2 * ABM_SYSTEM_BATCH, 4, periods);

    int count = 0;
    int *replicate = abm_system_list_replicates(dir, &count);

    CHECK(count == ABM_SYSTEM_BATCH + 4, "must list %d replicates, got %d",
          ABM_SYSTEM_BATCH + 4, count);

    int ascending = 1;
    for (int i = 1; i < count; i++) if (replicate[i] <= replicate[i - 1]) ascending = 0;
    CHECK(ascending, "the listing must be ascending");

    for (int i = 0; i < ABM_SYSTEM_BATCH && i < count; i++)
        CHECK(replicate[i] == i, "replicate %d of the first batch must be %d, got %d",
              i, i, replicate[i]);
    for (int i = 0; i < 4 && ABM_SYSTEM_BATCH + i < count; i++)
        CHECK(replicate[ABM_SYSTEM_BATCH + i] == 2 * ABM_SYSTEM_BATCH + i,
              "replicate %d of the third batch must be %d, got %d", i,
              2 * ABM_SYSTEM_BATCH + i, replicate[ABM_SYSTEM_BATCH + i]);

    fprintf(report, "   batches 0 and 2 written, batch 1 absent: %d replicates listed, "
                    "ascending %s\n", count, ascending ? "yes" : "NO");
    fprintf(report, "   first %d, last %d\n", count ? replicate[0] : -1,
            count ? replicate[count - 1] : -1);

    free(replicate);
}

static void remove_tree(const char *path) {
    char command[1024];
    snprintf(command, sizeof command, "rm -rf \"%s\"", path);
    if (system(command) != 0) fprintf(stderr, "abm_system_layout: could not remove %s\n", path);
}

int main(void) {
    make_directory(SCRATCH);

    report = fopen(REPORT, "w");
    assert(report && "abm_system_layout: cannot open the report path");
    fprintf(report, "applications/abm_system.h: the stored layout of the simulated dataset\n");
    fprintf(report, "%s build\n\n", sizeof(mreal) == sizeof(double) ? "float64" : "float32");

    printf("abm_system.h layout, %s build\n\n",
           sizeof(mreal) == sizeof(double) ? "float64" : "float32");

    test_transform_is_what_it_claims();
    test_agrees_with_the_real_data_route();
    test_burn_in_keeps_the_periods_it_claims();
    test_archive_round_trip();
    test_listing_replicates();

    fprintf(report, "\n%s, %d failure%s\n", failures ? "FAILED" : "PASSED", failures,
            failures == 1 ? "" : "s");
    fclose(report);

    remove_tree(SCRATCH);

    printf("\n%s, %d failure%s\n", failures ? "FAILED" : "PASSED", failures,
           failures == 1 ? "" : "s");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
