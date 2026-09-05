/*
Whether this project's build of the DSK simulator computes what the upstream
build computes.

model/dsk_sfc is upstream's code, compiled optimised where upstream compiles it
unoptimised, and with the rewrites docs/DSK_MODEL_CHANGES.md records: hot loops
reading the machine vintage arrays along the direction they are stored in,
quantities carried forward instead of rescanned, arrays that nothing reads
deleted, and the fixed 64-byte output filename buffers widened
(tests/dsk_long_path.c is that one's own test). All of it is there to make a
million runs affordable, and none of it is allowed to change a single number:
the compiler is told not to contract multiplies and adds, and every rewritten
loop accumulates the same terms in the same order.

So the check is equality, not similarity. Two runs that agree byte for byte are
indistinguishable under any test that could be applied to them, and no test
statistic is needed while that holds. If a future change breaks it, the worst
absolute and relative gap is reported with the period and column it sits at,
which is what separates a slightly different arithmetic (a reassociated sum)
from a different model (a different sequence of random draws).

Three things are required of each seed. The upstream build has to produce the
same results file twice, since every claim here rests on it being a fixed
reference. The two builds have to agree on the results file. And they have to
agree on the error log, so a run that ends early ends early in both for the
same stated reason.

tests/dsk_full_output_equivalence.c compares every file a run writes, not only
these two, and tests/dsk_ulp_sensitivity.c measures how small a difference in
the arithmetic these comparisons can still see.

Both builds come from the same vendored tree: `make model` builds the one this
project runs, and `make model-upstream` builds the reference from
model/dsk_sfc/upstream's copies of the changed files at upstream's own flags.
Nothing here reaches the network.

Seeds 1 to N, N = 3 unless given as the first argument. Writes
out/dsk_build_equivalence.txt.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "tests/dsk_upstream_scratch.h"

#define MODIFIED "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_build_equivalence.txt"
#define RUN_NAME "equiv"
#define MODEL_COLUMNS 83
#define DEFAULT_SEEDS 3

static void run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    int status = system(command);
    assert(status == 0 && "dsk_build_equivalence: a build did not complete");
}

static char *read_output(const char *dir, const char *name, long *size) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/%.200s", dir, name);

    FILE *f = fopen(path, "rb");
    if (!f) { *size = -1; return NULL; }
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size + 1);
    assert(text && "dsk_build_equivalence: cannot hold an output file");
    if (*size > 0) assert(fread(text, 1, (size_t)*size, f) == (size_t)*size);
    text[*size] = '\0';
    fclose(f);
    unlink(path);
    return text;
}

static char *read_results(const char *dir, int seed, long *size) {
    char name[256];
    snprintf(name, sizeof name, "results_%s_%d.txt", RUN_NAME, seed);
    char *text = read_output(dir, name, size);
    assert(text && "dsk_build_equivalence: a run produced no results file");
    return text;
}

/* The model writes what it complains about here, and a run that ends early
   says why here rather than in the results file. */
static char *read_errors(const char *dir, int seed, long *size) {
    char name[256];
    snprintf(name, sizeof name, "errors/Errors_%s_%d.txt", RUN_NAME, seed);
    return read_output(dir, name, size);
}

/* The worst gap between two outputs that are not identical, and where it is.
   Periods and columns are 1-indexed, as the model writes them. */
static void report_gap(FILE *report, const char *upstream_text, const char *modified_text) {
    const char *a = upstream_text, *b = modified_text;
    double worst_absolute = 0.0, worst_relative = 0.0;
    long index = 0, worst_index = -1;

    for (;;) {
        char *a_end, *b_end;
        double x = strtod(a, &a_end);
        double y = strtod(b, &b_end);
        if (a_end == a || b_end == b) break;

        double absolute = fabs(x - y);
        double relative = fabs(x) > 0 ? absolute / fabs(x) : absolute;
        if (absolute > worst_absolute) { worst_absolute = absolute; worst_index = index; }
        if (relative > worst_relative) worst_relative = relative;

        a = a_end;
        b = b_end;
        index++;
    }

    fprintf(report, "      worst absolute gap %.3g, worst relative gap %.3g\n",
            worst_absolute, worst_relative);
    if (worst_index >= 0)
        fprintf(report, "      worst gap at period %ld, column %ld, of %ld values\n",
                worst_index / MODEL_COLUMNS + 1, worst_index % MODEL_COLUMNS + 1, index);
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;
    assert(n_seeds >= 1 && "dsk_build_equivalence: at least one seed");

    char modified[PATH_MAX], upstream[PATH_MAX], inputs[PATH_MAX];
    assert(realpath(MODIFIED, modified) &&
           "dsk_build_equivalence: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(UPSTREAM, upstream) &&
           "dsk_build_equivalence: bin/dsk_SFC_upstream is not built - run make model-upstream");
    assert(realpath(INPUTS, inputs) && "dsk_build_equivalence: the model's own inputs JSON is missing");

    DskScratch scratch = dsk_scratch_open("eq", RUN_NAME, inputs, upstream, modified, n_seeds);
    const char *upstream_dir = scratch.upstream, *modified_dir = scratch.modified;

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_build_equivalence: cannot open the report path");
    fprintf(report, "upstream build: %s\n", upstream);
    fprintf(report, "this project's: %s\n", modified);
    fprintf(report, "inputs:         %s\n", inputs);
    fprintf(report, "scratch:        %s, %d characters of padding so upstream writes\n",
            scratch.root, scratch.padding);
    fprintf(report, "seeds 1 to %d; each seed checks that upstream repeats itself, then compares\nthe results file and the error log byte for byte\n\n", n_seeds);

    int matched = 0, repeatable = 0;
    for (int seed = 1; seed <= n_seeds; seed++) {
        run(upstream_dir, inputs, seed);

        /* The reference is only a reference if it says the same thing twice.
           Nothing else here would notice an upstream build that varied between
           runs, and every equality this project claims rests on it. */
        long first_size, second_size, first_error_size;
        char *first = read_results(upstream_dir, seed, &first_size);
        /* The model opens its error log with ios::app, so the first run's copy has
           to go before the second run adds to it. */
        free(read_errors(upstream_dir, seed, &first_error_size));

        run(upstream_dir, inputs, seed);
        char *second = read_results(upstream_dir, seed, &second_size);
        long upstream_error_size;
        char *upstream_errors = read_errors(upstream_dir, seed, &upstream_error_size);

        if (first_size == second_size && memcmp(first, second, (size_t)first_size) == 0) {
            repeatable++;
        } else {
            fprintf(report, "  seed %-4d the upstream build did not repeat itself\n", seed);
            report_gap(report, first, second);
        }
        free(first);

        run(modified_dir, inputs, seed);
        long modified_size, modified_error_size;
        char *modified_text = read_results(modified_dir, seed, &modified_size);
        char *modified_errors = read_errors(modified_dir, seed, &modified_error_size);

        int identical = (second_size == modified_size &&
                         memcmp(second, modified_text, (size_t)second_size) == 0);
        int errors_identical = (upstream_error_size == modified_error_size &&
                                (upstream_error_size < 0 ||
                                 memcmp(upstream_errors, modified_errors, (size_t)upstream_error_size) == 0));

        if (identical && errors_identical) {
            matched++;
            if (upstream_error_size < 0)
                fprintf(report, "  seed %-4d identical, %ld bytes, neither build wrote an error log\n",
                        seed, second_size);
            else
                fprintf(report, "  seed %-4d identical, %ld bytes and an error log of %ld\n",
                        seed, second_size, upstream_error_size);
        } else if (!identical) {
            fprintf(report, "  seed %-4d results DIFFER\n", seed);
            report_gap(report, second, modified_text);
        } else {
            fprintf(report, "  seed %-4d results identical but the error log DIFFERS,"
                            " %ld bytes upstream and %ld here\n",
                    seed, upstream_error_size, modified_error_size);
        }
        fflush(report);

        free(second);
        free(modified_text);
        free(upstream_errors);
        free(modified_errors);
    }

    fprintf(report, "\n%d of %d seeds identical; the upstream build repeated itself on %d of %d\n",
            matched, n_seeds, repeatable, n_seeds);
    fprintf(report, "%s\n", matched == n_seeds && repeatable == n_seeds ? "PASSED" : "FAILED");
    fclose(report);

    int passed = matched == n_seeds && repeatable == n_seeds;
    printf("%d of %d seeds identical between the upstream build and this project's\n", matched, n_seeds);
    printf("the upstream build repeated itself on %d of %d seeds\n", repeatable, n_seeds);
    printf("%s\n", passed ? "PASSED, 0 failures" : "FAILED");

    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
