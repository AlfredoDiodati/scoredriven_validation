/*
Whether counting machines in bigger lots leaves the economy where it was, and
whether it gets a run past the period where machine counts break.

model/dsk_sfc/dsk_sfc_machine_lots.h doubles the output one machine makes and
halves every count once a firm's holding passes a ceiling, because a double
counts whole numbers exactly only up to 2^53 and machine counts grow with real
output. Halving an odd count rounds, so this cannot be exact the way
redenominating money is. What it can be is small: at the ceiling a holding is
around 1e12 machines, so rounding moves it by about one part in 1e12, and the
economy that comes out is the same economy.

Three runs, seed 1, the model's own parameter file:

  reference   lots switched off, 11,000 periods
  lots        lots at their defaults, 11,000 periods
  long        lots at their defaults, 14,000 periods

What is required:

  - the lots run completes with an empty error log, which is also the model's own
    accounting check: it writes there when net worth and assets disagree by more
    than its tolerance, so a rebase that left the stocks inconsistent shows up;
  - it rebases at least once, and each rebase doubles the output a machine makes,
    read from the log the model writes;
  - every period before the first rebase is identical to the reference;
  - in the period the first rebase happens, real GDP, consumption, investment,
    employment and the capital stock are within 1e-9 of the reference: the
    rebase itself is a rounding;
  - afterwards the two runs are different paths, so what has to agree is the
    process rather than the path. The mean growth of real GDP and of the price
    level over the remaining periods must be within four standard errors of the
    reference's, and their spreads within a fifth of each other;
  - the long run completes all 14,000 periods, past period 12,908, where the
    model without lots stops advancing and never exits.

The paths separate because the model amplifies any perturbation: its own
dsk_ulp_sensitivity test shows a change of one part in 1e15 in a parameter
reaching the printed output within a few periods. A rebase rounds holdings by
about one part in 1e12, so a different path afterwards is what a correct
implementation produces, and a level comparison is only meaningful in the
rebase period itself.

The reference run stops at 11,000 periods because further on it hangs.

Writes out/dsk_machine_lot_rebase.txt.
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
#define REPORT "out/dsk_machine_lot_rebase.txt"
#define RUN_NAME "lots"
#define MODEL_COLUMNS 83
#define SEED 1
#define REFERENCE_PERIODS 11000
#define LONG_PERIODS 14000
#define HANGS_WITHOUT_LOTS 12908
#define AGREEMENT 1e-9
#define STANDARD_ERRORS 4.0
#define SPREAD_TOLERANCE 0.2
#define MACHINE_SIZE 40.0

/* Columns of the results file: the economy-wide quantities a rebase must leave
   where they were. */
static const int compared_column[] = {2, 3, 4, 5, 10};
static const char *const compared_name[] = {"real GDP", "consumption", "investment",
                                            "employment", "capital stock"};
#define COMPARED (int)(sizeof compared_column / sizeof compared_column[0])

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_machine_lot_rebase: mkdir failed");
}

static void write_inputs(const char *source, const char *target, int periods, int lots) {
    FILE *in = fopen(source, "rb");
    assert(in && "dsk_machine_lot_rebase: cannot read the model's parameter file");
    FILE *out = fopen(target, "wb");
    assert(out && "dsk_machine_lot_rebase: cannot write a parameter file");

    char line[4096];
    while (fgets(line, sizeof line, in)) {
        if (strstr(line, "\"T\"")) {
            fprintf(out, "      \"T\": %d,\n", periods);
            if (!lots) fprintf(out, "      \"machine_lot_ceiling_exponent\": 0,\n");
        } else {
            fputs(line, out);
        }
    }
    fclose(out);
    fclose(in);
}

static int run(const char *dir, const char *inputs, const char *name) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, name, SEED);
    return system(command);
}

static double *read_results(const char *dir, const char *name, int expected, int *periods) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, name, SEED);
    FILE *f = fopen(path, "r");
    if (!f) { *periods = -1; return NULL; }

    int capacity = expected + 16, rows = 0;
    double *value = malloc((size_t)capacity * MODEL_COLUMNS * sizeof(double));
    assert(value && "dsk_machine_lot_rebase: out of memory");

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

static long error_log_size(const char *dir, const char *name) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/errors/Errors_%s_%d.txt", dir, name, SEED);
    struct stat info;
    if (stat(path, &info) != 0) return -1;
    return (long)info.st_size;
}

/* The periods at which the run rebased, and the output one machine makes after
   each one, from the log the model writes. */
static int read_rebases(const char *dir, const char *name, int *period, double *machine_size, int room) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d_machine_lots.txt", dir, name, SEED);
    FILE *f = fopen(path, "r");
    if (!f) return 0;

    int events = 0, at, step, cumulative;
    double size;
    while (events < room && fscanf(f, "%d %d %d %lf", &at, &step, &cumulative, &size) == 4) {
        period[events] = at;
        machine_size[events] = size;
        events++;
    }
    fclose(f);
    return events;
}

int main(void) {
    char model[PATH_MAX], inputs[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "dsk_machine_lot_rebase: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, inputs) && "dsk_machine_lot_rebase: the parameter file is missing");

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], dir[640];
    snprintf(root, sizeof root, "%s/dsk_machine_lots_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(dir, sizeof dir, "%s/run", root);
    make_directory(dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    assert(symlink(model, link) == 0 && "dsk_machine_lot_rebase: cannot link the model");

    char reference_json[1024], lots_json[1024], long_json[1024];
    snprintf(reference_json, sizeof reference_json, "%s/reference.json", root);
    snprintf(lots_json, sizeof lots_json, "%s/lots.json", root);
    snprintf(long_json, sizeof long_json, "%s/long.json", root);
    write_inputs(inputs, reference_json, REFERENCE_PERIODS, 0);
    write_inputs(inputs, lots_json, REFERENCE_PERIODS, 1);
    write_inputs(inputs, long_json, LONG_PERIODS, 1);

    int reference_status = run(dir, reference_json, "ref");
    int lots_status = run(dir, lots_json, RUN_NAME);
    int long_status = run(dir, long_json, "long");

    int reference_periods, lots_periods, long_periods;
    double *reference = read_results(dir, "ref", REFERENCE_PERIODS, &reference_periods);
    double *lots = read_results(dir, RUN_NAME, REFERENCE_PERIODS, &lots_periods);
    double *longer = read_results(dir, "long", LONG_PERIODS, &long_periods);

    int rebase_period[64];
    double machine_size[64];
    int rebases = read_rebases(dir, RUN_NAME, rebase_period, machine_size, 64);
    int long_rebases = read_rebases(dir, "long", rebase_period, machine_size, 0);
    (void)long_rebases;

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_machine_lot_rebase: cannot open the report path");
    fprintf(report, "Counting machines in lots, against a run that does not, seed %d.\n\n", SEED);

    int failures = 0;

    fprintf(report, "reference  %d periods, exit %d, %ld bytes of errors\n",
            reference_periods, reference_status, error_log_size(dir, "ref"));
    fprintf(report, "lots       %d periods, exit %d, %ld bytes of errors, %d rebases\n",
            lots_periods, lots_status, error_log_size(dir, RUN_NAME), rebases);
    fprintf(report, "long       %d periods, exit %d, %ld bytes of errors\n\n",
            long_periods, long_status, error_log_size(dir, "long"));

    if (lots_status != 0 || lots_periods != REFERENCE_PERIODS || error_log_size(dir, RUN_NAME) != 0) {
        fprintf(report, "the run with lots did not complete cleanly  FAILED\n");
        failures++;
    }
    if (rebases < 1) {
        fprintf(report, "the run never rebased, so nothing was tested  FAILED\n");
        failures++;
    }

    /* Each rebase doubles the output one machine makes. */
    for (int k = 0; k < rebases; k++) {
        double expected = MACHINE_SIZE * ldexp(1.0, k + 1);
        if (fabs(machine_size[k] - expected) > 0) {
            fprintf(report, "rebase %d at period %d: a machine makes %.10g, expected %.10g  FAILED\n",
                    k + 1, rebase_period[k], machine_size[k], expected);
            failures++;
        }
    }
    if (rebases > 0)
        fprintf(report, "first rebase at period %d, a machine makes %.10g afterwards\n\n",
                rebase_period[0], machine_size[0]);

    /* Everything before the first rebase must be the run that never rebased. */
    if (rebases > 0 && reference && lots) {
        int differing = 0;
        for (int row = 0; row < rebase_period[0] - 1 && row < reference_periods && row < lots_periods; row++)
            for (int c = 0; c < MODEL_COLUMNS; c++)
                differing += reference[row * MODEL_COLUMNS + c] != lots[row * MODEL_COLUMNS + c];
        fprintf(report, "periods 1 to %d, before the first rebase: %d values differ  %s\n\n",
                rebase_period[0] - 1, differing, differing ? "FAILED" : "");
        failures += differing != 0;
    }

    /* The rebase period itself, against the run that never rebased. */
    if (rebases > 0 && reference && lots && reference_periods >= rebase_period[0] &&
        lots_periods >= rebase_period[0]) {
        int row = rebase_period[0] - 1;
        fprintf(report, "period %d, in which the first rebase happens:\n", rebase_period[0]);
        for (int k = 0; k < COMPARED; k++) {
            int c = compared_column[k] - 1;
            double a = reference[row * MODEL_COLUMNS + c];
            double b = lots[row * MODEL_COLUMNS + c];
            double gap = a != 0 ? fabs(b - a) / fabs(a) : fabs(b);
            int ok = gap <= AGREEMENT;
            failures += !ok;
            fprintf(report, "  %-14s %.10g against %.10g, relative gap %.3e  %s\n",
                    compared_name[k], a, b, gap, ok ? "" : "FAILED");
        }
        fprintf(report, "\n");
    }

    /* Afterwards the paths differ, so what is compared is the process: how fast
       the economy and prices grow, and how much that varies. */
    if (rebases > 0 && reference && lots) {
        const int from = rebase_period[0] + 30;
        const int columns[2] = {2, 34};
        const char *const names[2] = {"real GDP", "price level"};
        for (int k = 0; k < 2; k++) {
            int c = columns[k] - 1;
            double mean[2] = {0, 0}, spread[2] = {0, 0};
            int n = 0;
            for (int pass = 0; pass < 2; pass++) {
                const double *series = pass ? lots : reference;
                int periods = pass ? lots_periods : reference_periods;
                double sum = 0, square = 0;
                n = 0;
                for (int row = from; row < periods; row++) {
                    double previous = series[(row - 1) * MODEL_COLUMNS + c];
                    double now = series[row * MODEL_COLUMNS + c];
                    if (previous <= 0 || now <= 0) continue;
                    double growth = log(now / previous);
                    sum += growth;
                    square += growth * growth;
                    n++;
                }
                mean[pass] = sum / n;
                spread[pass] = sqrt(square / n - mean[pass] * mean[pass]);
            }
            double standard_error = sqrt((spread[0] * spread[0] + spread[1] * spread[1]) / n);
            double distance = fabs(mean[1] - mean[0]) / standard_error;
            double spread_gap = fabs(spread[1] - spread[0]) / spread[0];
            int ok = distance <= STANDARD_ERRORS && spread_gap <= SPREAD_TOLERANCE;
            failures += !ok;
            fprintf(report,
                    "%-12s growth over periods %d to %d: mean %.4f%% against %.4f%%, %.1f standard errors apart, "
                    "spread %.4f%% against %.4f%%  %s\n",
                    names[k], from + 1, reference_periods, 100 * mean[0], 100 * mean[1], distance,
                    100 * spread[0], 100 * spread[1], ok ? "" : "FAILED");
        }
        fprintf(report, "\n");
    }

    if (long_status != 0 || long_periods != LONG_PERIODS || error_log_size(dir, "long") != 0) {
        fprintf(report, "the %d-period run did not complete, so lots did not get past period %d  FAILED\n",
                LONG_PERIODS, HANGS_WITHOUT_LOTS);
        failures++;
    } else {
        fprintf(report, "the %d-period run completed, past period %d where the model without lots\n"
                        "stops advancing\n", LONG_PERIODS, HANGS_WITHOUT_LOTS);
    }

    fprintf(report, "\n%s\n", failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    printf("machine lot rebase: %d failures\n", failures);
    printf("%s\n", failures == 0 ? "PASSED, 0 failures" : "FAILED");

    free(reference);
    free(lots);
    free(longer);
    dsk_scratch_finish(root, failures == 0);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
