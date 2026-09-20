/*
Whether dividing the money side of the model by a power of two leaves the model
alone.

model/dsk_sfc/dsk_sfc_redenomination.h divides every money quantity by a power
of two once the wage grows large, so that prices and money stocks, which grow
without limit, stay inside the range of a double. Nothing in the model depends
on the unit money is counted in, so a run that redenominates should be the run
that does not, with every money number carrying a different exponent and the
same digits.

That is what this checks, and it is also what finds a mistake in the list of
what counts as money. A money variable left out of the list keeps the old unit
while everything around it changes, and a variable wrongly called money changes
when it should not. Either way the model's own behaviour changes: firms compare
prices against costs that no longer match, and the run diverges.

The comparison needs no list of its own. Two runs, identical but for the
ceiling, are compared column by column, and every value must be either
unchanged or scaled by exactly the factor the run applied by that period:
real quantities, employment and rates come back unchanged, money columns come
back scaled, and a column that is neither is the failure. The factor is read
from the log the model writes when it redenominates, so the test does not need
to know which columns are money.

Exactness: multiplying by a power of two changes only a double's exponent, so
"scaled" means the values match to the last digit printed, not approximately.
The results file is written at ten decimal places, so the comparison allows
half a unit in the last printed place, scaled up with the factor.

Setup: the model's own parameter file, 600 periods, seeds 1 and 2. The
reference run leaves the ceiling at its default, which no 600-period run
reaches, so it never redenominates. The other sets the ceiling at 2^4 and the
step at 2^2, small enough that the wage grows back through the ceiling within
the run, so the comparison covers several redenominations rather than one.

Writes out/dsk_redenomination_invariance.txt.
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
#define REPORT "out/dsk_redenomination_invariance.txt"
#define RUN_NAME "reden"
#define MODEL_COLUMNS 83
#define PERIODS 600
#define DEFAULT_SEEDS 2
#define CEILING_EXPONENT 4
#define STEP_EXPONENT 2
#define PRINTED_DECIMALS 10

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_redenomination_invariance: mkdir failed");
}

/* The model's parameter file with the horizon fixed and, for the redenominating
   run, a ceiling low enough to fire. Written by text substitution rather than a
   JSON library: the file is machine-written, one key per line. */
static void write_inputs(const char *source, const char *target, int redenominating) {
    FILE *in = fopen(source, "rb");
    assert(in && "dsk_redenomination_invariance: cannot read the model's parameter file");
    FILE *out = fopen(target, "wb");
    assert(out && "dsk_redenomination_invariance: cannot write a parameter file");

    char line[4096];
    while (fgets(line, sizeof line, in)) {
        if (strstr(line, "\"T\"")) {
            fprintf(out, "      \"T\": %d,\n", PERIODS);
            if (redenominating)
                fprintf(out, "      \"redenomination_ceiling_exponent\": %d,\n"
                             "      \"redenomination_step_exponent\": %d,\n",
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
    assert(value && "dsk_redenomination_invariance: out of memory");

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

/* The factor the run's money side has been divided by, period by period, from
   the log the model writes. Returns 0 if the run never redenominated. */
static int read_exponents(const char *dir, int seed, int *exponent, int periods) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d_redenominations.txt", dir, RUN_NAME, seed);
    for (int t = 0; t < periods; t++) exponent[t] = 0;

    FILE *f = fopen(path, "r");
    if (!f) return 0;

    int period, step, cumulative, events = 0;
    while (fscanf(f, "%d %d %d", &period, &step, &cumulative) == 3) {
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
           "dsk_redenomination_invariance: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, inputs) && "dsk_redenomination_invariance: the parameter file is missing");

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], plain_dir[640], scaled_dir[640];
    snprintf(root, sizeof root, "%s/dsk_redenomination_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(plain_dir, sizeof plain_dir, "%s/plain", root);
    snprintf(scaled_dir, sizeof scaled_dir, "%s/scaled", root);
    make_directory(plain_dir);
    make_directory(scaled_dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", plain_dir);
    assert(symlink(model, link) == 0 && "dsk_redenomination_invariance: cannot link the model");
    snprintf(link, sizeof link, "%s/dsk_SFC", scaled_dir);
    assert(symlink(model, link) == 0 && "dsk_redenomination_invariance: cannot link the model");

    char plain_json[1024], scaled_json[1024];
    snprintf(plain_json, sizeof plain_json, "%s/plain.json", root);
    snprintf(scaled_json, sizeof scaled_json, "%s/scaled.json", root);
    write_inputs(inputs, plain_json, 0);
    write_inputs(inputs, scaled_json, 1);

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_redenomination_invariance: cannot open the report path");
    fprintf(report, "Dividing the money side of the model by powers of two, against a run that does not.\n\n");
    fprintf(report, "%d periods, seeds 1 to %d, ceiling 2^%d, step 2^%d, %d columns compared per period.\n",
            PERIODS, n_seeds, CEILING_EXPONENT, STEP_EXPONENT, MODEL_COLUMNS);
    fprintf(report, "Each value must be unchanged or scaled by exactly the factor applied by that period.\n\n");

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
        assert(exponent && "dsk_redenomination_invariance: out of memory");
        int events = read_exponents(scaled_dir, seed, exponent, PERIODS);

        int unchanged = 0, rescaled = 0, wrong = 0;
        int worst_period = 0, worst_column = 0;
        double worst_plain = 0, worst_scaled = 0;

        for (int t = 0; t < PERIODS; t++) {
            double factor = ldexp(1.0, exponent[t]);
            for (int c = 0; c < MODEL_COLUMNS; c++) {
                double a = plain[t * MODEL_COLUMNS + c];
                double b = scaled[t * MODEL_COLUMNS + c];
                double printed = pow(10.0, -PRINTED_DECIMALS);
                double same_tolerance = 0.5 * printed;
                double scaled_tolerance = 0.5 * printed * (1.0 + factor) + fabs(a) * factor * 1e-12;

                if (fabs(a - b) <= same_tolerance) unchanged++;
                else if (fabs(b - a * factor) <= scaled_tolerance) rescaled++;
                else {
                    if (!wrong) { worst_period = t + 1; worst_column = c + 1; worst_plain = a; worst_scaled = b; }
                    wrong++;
                }
            }
        }

        fprintf(report, "seed %d: %d redenominations, %d values unchanged, %d scaled, %d neither  %s\n",
                seed, events, unchanged, rescaled, wrong, wrong ? "FAILED" : "passed");
        if (events == 0) {
            fprintf(report, "  the run never redenominated, so nothing was tested  FAILED\n");
            failures++;
        }
        if (wrong) {
            fprintf(report, "  first mismatch at period %d, column %d: %.10g against %.10g, ratio %.12g\n",
                    worst_period, worst_column, worst_plain, worst_scaled,
                    worst_plain != 0 ? worst_scaled / worst_plain : 0.0);
            failures++;
        }

        free(exponent);
        free(plain);
        free(scaled);
    }

    fprintf(report, "\n%s\n", failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    printf("redenomination invariance: %d failures\n", failures);
    printf("%s\n", failures == 0 ? "PASSED, 0 failures" : "FAILED");
    dsk_scratch_finish(root, failures == 0);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
