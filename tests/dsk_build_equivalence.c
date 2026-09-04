/*
Whether this project's build of the DSK simulator computes what the upstream
build computes.

model/dsk_sfc is upstream's code with three changes: it is compiled optimised
instead of unoptimised, two loops in MACH() were rewritten to read the machine
vintage arrays along the direction they are stored in, and the fixed 64-byte
output filename buffers were widened (tests/dsk_long_path.c is that one's own
test). The first two are there to make a million runs affordable, and neither
is allowed to change a single number: the compiler is told not to contract
multiplies and adds, and the rewritten loops accumulate the same terms in the
same order, only reading them in a different sequence.

So the check is equality, not similarity. Two runs that agree byte for byte
are indistinguishable under any test that could be applied to them, and no
test statistic is needed while that holds. If a future change breaks it, the
worst absolute and relative gap is reported with the period and column it sits
at, which is what separates a slightly different arithmetic (a reassociated
sum) from a different model (a different sequence of random draws).

Both builds come from the same vendored tree: `make model` builds the one this
project runs, and `make model-upstream` builds the reference from
model/dsk_sfc/upstream's copies of the three changed files at upstream's own
flags. Nothing here reaches the network.

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

#define MODIFIED "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_build_equivalence.txt"
#define RUN_NAME "equiv"
#define MODEL_COLUMNS 83
#define DEFAULT_SEEDS 3

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_build_equivalence: mkdir failed");
}

/* Each build runs through a symlink in its own directory: the model writes its
   output beside the path it was invoked as, so two builds sharing one
   directory would overwrite each other. */
static void prepare(const char *dir, const char *executable) {
    make_directory(dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(executable, link) == 0 && "dsk_build_equivalence: cannot link a build into place");
}

static void run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    int status = system(command);
    assert(status == 0 && "dsk_build_equivalence: a build did not complete");
}

static char *read_results(const char *dir, int seed, long *size) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, seed);

    FILE *f = fopen(path, "rb");
    assert(f && "dsk_build_equivalence: a run produced no results file");
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size + 1);
    assert(text && fread(text, 1, (size_t)*size, f) == (size_t)*size &&
           "dsk_build_equivalence: cannot read a results file");
    text[*size] = '\0';
    fclose(f);
    unlink(path);
    return text;
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

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], upstream_dir[640], modified_dir[640];
    snprintf(root, sizeof root, "%s/dsk_equivalence_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(upstream_dir, sizeof upstream_dir, "%s/upstream", root);
    snprintf(modified_dir, sizeof modified_dir, "%s/modified", root);

    prepare(upstream_dir, upstream);
    prepare(modified_dir, modified);

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_build_equivalence: cannot open the report path");
    fprintf(report, "upstream build: %s\n", upstream);
    fprintf(report, "this project's: %s\n", modified);
    fprintf(report, "inputs:         %s\n", inputs);
    fprintf(report, "seeds 1 to %d, comparing every byte of each run's results file\n\n", n_seeds);

    int matched = 0;
    for (int seed = 1; seed <= n_seeds; seed++) {
        run(upstream_dir, inputs, seed);
        run(modified_dir, inputs, seed);

        long upstream_size, modified_size;
        char *upstream_text = read_results(upstream_dir, seed, &upstream_size);
        char *modified_text = read_results(modified_dir, seed, &modified_size);

        int identical = (upstream_size == modified_size &&
                         memcmp(upstream_text, modified_text, (size_t)upstream_size) == 0);
        if (identical) {
            matched++;
            fprintf(report, "  seed %-4d identical, %ld bytes\n", seed, upstream_size);
        } else {
            fprintf(report, "  seed %-4d DIFFERS\n", seed);
            report_gap(report, upstream_text, modified_text);
        }
        fflush(report);

        free(upstream_text);
        free(modified_text);
    }

    fprintf(report, "\n%d of %d seeds identical\n", matched, n_seeds);
    fprintf(report, "%s\n", matched == n_seeds ? "PASSED" : "FAILED");
    fclose(report);

    printf("%d of %d seeds identical between the upstream build and this project's\n", matched, n_seeds);
    printf("%s\n", matched == n_seeds ? "PASSED, 0 failures" : "FAILED");

    return matched == n_seeds ? EXIT_SUCCESS : EXIT_FAILURE;
}
