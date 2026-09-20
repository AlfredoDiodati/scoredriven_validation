/*
Whether counting the consumption good in a bigger unit leaves the model alone.

model/dsk_sfc/dsk_sfc_good_unit.h divides every quantity of the good by a power
of two once real GDP grows large, and multiplies the money price of a unit by
the same power, so that output and capacity, which grow without limit, stay
inside the range of a double. Nothing in the model depends on whether the good
is counted in tonnes or in hundreds of tonnes, so a run that changes the unit
should be the run that does not, with quantities and prices carrying different
exponents and the same digits.

That is what this checks, and it is also what finds a mistake in the list of
what is measured in the good. A quantity left out keeps the old unit while
everything around it changes, and something wrongly included changes when it
should not. Either way firms compare quantities against capacities that no
longer match, and the run diverges.

The comparison needs no list of its own. Two runs, identical but for the
ceiling, are compared column by column, and every value must be either
unchanged, divided or multiplied by exactly the factor applied by that period:
quantities divide, prices multiply, money stocks and rates stay unchanged, and
a column that is none of the three is the failure. The factor is read from the
log the model writes when it changes the unit.

One column is expected to miss, and it is the model's own doing rather than the
rescaling's. Upstream's mean productivity adds the consumption firms' output per
worker to the capital firms' machines per worker (`module_macro_sfc.cpp`, the
loops that build Am), which are amounts of different things. Changing what a
unit of the good is therefore moves that average by a fraction of a percent, the
wage rule reads its growth, and R&D staffing follows the wage. Measured at the
first change: every economy-wide column scales exactly and R&D labour moves by
2e-4. The test allows one column to miss by up to 1e-3 and requires the rest to
be exact.

After that first period the two runs are different paths, because the model
amplifies any perturbation - its own dsk_ulp_sensitivity test shows a change of
one part in 1e15 reaching the output within a few periods - so from there the
comparison is of the process: mean growth of real GDP within four standard
errors of the reference's.

Exactness: multiplying by a power of two changes only a double's exponent, so
"scaled" means the values match to the last digit printed, not approximately.
The results file is written at ten decimal places, so the comparison allows
half a unit in the last printed place, scaled up with the factor.

Setup: the model's own parameter file, 600 periods, seeds 1 and 2. The
reference run leaves the ceiling at its default, which no 600-period run
reaches, so it never changes the unit. The other sets the ceiling at 2^18 and
the step at 2^1, small enough that real GDP grows back through the ceiling
within the run, so the comparison covers several changes rather than one.

Writes out/dsk_good_unit_invariance.txt.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>
#include "tests/dsk_scratch.h"

#define MODEL "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_good_unit_invariance.txt"
#define RUN_NAME "goodunit"
#define MODEL_COLUMNS 83
#define PERIODS 2000
#define DEFAULT_SEEDS 2
#define CEILING_EXPONENT 18
#define STEP_EXPONENT 1
#define PRINTED_DECIMALS 10
/* The model's own mean productivity adds output per worker to machines per
   worker, so changing the good's unit moves it, the wage rule reads its growth,
   and the wage moves with it. These four columns are the ones that reach: the
   two mean productivities, the real wage, and R&D staffing, which is paid for
   out of wages. Nothing else may move. */
static const int mixed_unit_column[] = {11, 13, 14, 36};
static const char *const mixed_unit_name[] = {"R&D labour", "mean energy efficiency",
                                              "mean productivity", "real wage"};
#define MIXED_UNIT_COLUMNS (int)(sizeof mixed_unit_column / sizeof mixed_unit_column[0])
#define WORST_ALLOWED 2e-3
#define STANDARD_ERRORS 4.0

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_good_unit_invariance: mkdir failed");
}

/* The model's parameter file with the horizon fixed and, for the redenominating
   run, a ceiling low enough to fire. Written by text substitution rather than a
   JSON library: the file is machine-written, one key per line. */
static void write_inputs(const char *source, const char *target, int redenominating) {
    FILE *in = fopen(source, "rb");
    assert(in && "dsk_good_unit_invariance: cannot read the model's parameter file");
    FILE *out = fopen(target, "wb");
    assert(out && "dsk_good_unit_invariance: cannot write a parameter file");

    char line[4096];
    while (fgets(line, sizeof line, in)) {
        if (strstr(line, "\"T\"")) {
            fprintf(out, "      \"T\": %d,\n", PERIODS);
            if (redenominating)
                fprintf(out, "      \"good_unit_ceiling_exponent\": %d,\n"
                             "      \"good_unit_step_exponent\": %d,\n",
                        CEILING_EXPONENT, STEP_EXPONENT);
        } else {
            fputs(line, out);
        }
    }
    fclose(out);
    fclose(in);
}

static int run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    return system(command);
}

static double *read_results(const char *dir, int seed, int *periods) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, seed);
    FILE *f = fopen(path, "r");
    if (!f) { *periods = -1; return NULL; }

    int capacity = PERIODS + 16, rows = 0;
    double *value = malloc((size_t)capacity * MODEL_COLUMNS * sizeof(double));
    assert(value && "dsk_good_unit_invariance: out of memory");

    while (rows < capacity) {
        int read = 0;
        for (int c = 0; c < MODEL_COLUMNS; c++)
            read += fscanf(f, "%lf", &value[rows * MODEL_COLUMNS + c]) == 1;
        if (read < MODEL_COLUMNS) break;
        rows++;
    }
    fclose(f);
    *periods = rows;
    return value;
}

/* The power of two the good's unit has grown by, period by period, from the log
   the model writes. Returns 0 if the run never changed the unit. */
static int read_exponents(const char *dir, int seed, int *exponent, int periods, int *first_change) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d_good_units.txt", dir, RUN_NAME, seed);
    for (int t = 0; t < periods; t++) exponent[t] = 0;

    FILE *f = fopen(path, "r");
    if (!f) return 0;

    int period, step, cumulative, events = 0;
    *first_change = 0;
    while (fscanf(f, "%d %d %d", &period, &step, &cumulative) == 3) {
        if (!events) *first_change = period;
        /* The period's row is written before the redenomination happens, so the
           new unit first shows up in the row after it. */
        for (int t = period; t < periods; t++) exponent[t] = cumulative;
        events++;
    }
    fclose(f);
    return events;
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;

    char model[PATH_MAX], inputs[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "dsk_good_unit_invariance: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, inputs) && "dsk_good_unit_invariance: the parameter file is missing");

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], plain_dir[640], scaled_dir[640];
    snprintf(root, sizeof root, "%s/dsk_good_unit_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(plain_dir, sizeof plain_dir, "%s/plain", root);
    snprintf(scaled_dir, sizeof scaled_dir, "%s/scaled", root);
    make_directory(plain_dir);
    make_directory(scaled_dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", plain_dir);
    assert(symlink(model, link) == 0 && "dsk_good_unit_invariance: cannot link the model");
    snprintf(link, sizeof link, "%s/dsk_SFC", scaled_dir);
    assert(symlink(model, link) == 0 && "dsk_good_unit_invariance: cannot link the model");

    char plain_json[1024], scaled_json[1024];
    snprintf(plain_json, sizeof plain_json, "%s/plain.json", root);
    snprintf(scaled_json, sizeof scaled_json, "%s/scaled.json", root);
    write_inputs(inputs, plain_json, 0);
    write_inputs(inputs, scaled_json, 1);

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_good_unit_invariance: cannot open the report path");
    fprintf(report, "Counting the good in a bigger unit, against a run that does not.\n\n");
    fprintf(report, "%d periods, seeds 1 to %d, ceiling 2^%d, step 2^%d, %d columns compared per period.\n",
            PERIODS, n_seeds, CEILING_EXPONENT, STEP_EXPONENT, MODEL_COLUMNS);
    fprintf(report, "Each value must be unchanged, divided or multiplied by exactly the factor applied.\n\n");

    int failures = 0;
    for (int seed = 1; seed <= n_seeds; seed++) {
        int plain_status = run(plain_dir, plain_json, seed);
        int scaled_status = run(scaled_dir, scaled_json, seed);

        int plain_periods, scaled_periods;
        double *plain = read_results(plain_dir, seed, &plain_periods);
        double *scaled = read_results(scaled_dir, seed, &scaled_periods);

        if (plain_status != 0 || scaled_status != 0 || !plain || !scaled ||
            plain_periods != PERIODS || scaled_periods != PERIODS) {
            fprintf(report, "seed %d: a run did not complete (exit %d and %d, %d and %d periods)  FAILED\n",
                    seed, plain_status, scaled_status, plain_periods, scaled_periods);
            failures++;
            free(plain); free(scaled);
            continue;
        }

        int *exponent = malloc((size_t)PERIODS * sizeof(int));
        assert(exponent && "dsk_good_unit_invariance: out of memory");
        int first_change = 0;
        int events = read_exponents(scaled_dir, seed, exponent, PERIODS, &first_change);

        if (events == 0) {
            fprintf(report, "seed %d: the run never changed the unit, so nothing was tested  FAILED\n", seed);
            failures++;
            free(exponent); free(plain); free(scaled);
            continue;
        }

        /* Before the first change the two runs are the same run. */
        int differing = 0;
        for (int row = 0; row < first_change; row++)
            for (int c = 0; c < MODEL_COLUMNS; c++)
                differing += plain[row * MODEL_COLUMNS + c] != scaled[row * MODEL_COLUMNS + c];

        /* The first period written in the new unit. */
        const int row = first_change;
        const double factor = ldexp(1.0, -exponent[row]);
        const double printed = pow(10.0, -PRINTED_DECIMALS);
        int exact = 0, missed = 0, moved = 0, worst_column = 0;
        double worst_gap = 0, worst_moved = 0, worst_plain = 0, worst_scaled = 0;

        for (int c = 0; c < MODEL_COLUMNS; c++) {
            double a = plain[row * MODEL_COLUMNS + c];
            double b = scaled[row * MODEL_COLUMNS + c];
            double as_is = fabs(b - a);
            double divided = fabs(b - a * factor);
            double multiplied = fabs(b - a / factor);
            double closest = fmin(as_is, fmin(divided, multiplied));
            double allowed = 0.5 * printed * (1.0 + 1.0 / factor) + fabs(a) / factor * 1e-12;

            if (closest <= allowed) {
                exact++;
                continue;
            }

            int expected_to_move = 0;
            for (int k = 0; k < MIXED_UNIT_COLUMNS; k++) expected_to_move |= mixed_unit_column[k] == c + 1;

            double scale = closest == divided ? fabs(a * factor) : (closest == multiplied ? fabs(a / factor) : fabs(a));
            double gap = scale != 0 ? closest / scale : closest;
            if (expected_to_move) {
                moved++;
                if (gap > worst_moved) worst_moved = gap;
            } else {
                missed++;
                if (gap > worst_gap) {
                    worst_gap = gap; worst_column = c + 1; worst_plain = a; worst_scaled = b;
                }
            }
        }

        int seed_failures = 0;
        seed_failures += differing != 0;
        seed_failures += missed != 0;
        seed_failures += worst_moved > WORST_ALLOWED;
        failures += seed_failures;

        fprintf(report, "seed %d: %d unit changes, first at period %d\n", seed, events, first_change);
        fprintf(report, "  periods 1 to %d, before it: %d values differ  %s\n",
                first_change, differing, differing ? "FAILED" : "");
        fprintf(report, "  period %d, the first written in the new unit: %d of %d columns exactly unchanged,\n"
                        "    divided or multiplied\n", row + 1, exact, MODEL_COLUMNS);
        fprintf(report, "    %d of the %d columns the mixed-unit mean productivity reaches moved, worst by %.2e  %s\n",
                moved, MIXED_UNIT_COLUMNS, worst_moved, worst_moved > WORST_ALLOWED ? "FAILED" : "");
        if (missed)
            fprintf(report, "    column %d moved and should not have: %.10g against %.10g, %.2e away from a clean factor  FAILED\n",
                    worst_column, worst_plain, worst_scaled, worst_gap);

        /* Afterwards the paths differ, so the process is what is compared. */
        {
            const int from = first_change + 30;
            double mean[2] = {0, 0}, spread[2] = {0, 0};
            int n = 0;
            for (int pass = 0; pass < 2; pass++) {
                const double *series = pass ? scaled : plain;
                double sum = 0, square = 0;
                n = 0;
                for (int r = from; r < PERIODS; r++) {
                    /* Put the rescaled run back in the unit it started in, which
                       is exact, so that a change of unit is not read as a fall
                       in output. */
                    double back = pass ? ldexp(1.0, exponent[r]) : 1.0;
                    double back_previous = pass ? ldexp(1.0, exponent[r - 1]) : 1.0;
                    double previous = series[(r - 1) * MODEL_COLUMNS + 1] * back_previous;
                    double now = series[r * MODEL_COLUMNS + 1] * back;
                    if (previous <= 0 || now <= 0) continue;
                    double growth = log(now / previous);
                    sum += growth; square += growth * growth; n++;
                }
                mean[pass] = sum / n;
                spread[pass] = sqrt(square / n - mean[pass] * mean[pass]);
            }
            double standard_error = sqrt((spread[0] * spread[0] + spread[1] * spread[1]) / n);
            double distance = fabs(mean[1] - mean[0]) / standard_error;
            int ok = distance <= STANDARD_ERRORS;
            failures += !ok;
            fprintf(report, "  real GDP growth over periods %d to %d: mean %.4f%% against %.4f%%, "
                            "%.1f standard errors apart  %s\n\n",
                    from + 1, PERIODS, 100 * mean[0], 100 * mean[1], distance, ok ? "" : "FAILED");
        }

        free(exponent);
        free(plain);
        free(scaled);
    }

    fprintf(report, "\n%s\n", failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    printf("good unit invariance: %d failures\n", failures);
    printf("%s\n", failures == 0 ? "PASSED, 0 failures" : "FAILED");
    dsk_scratch_finish(root, failures == 0);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
