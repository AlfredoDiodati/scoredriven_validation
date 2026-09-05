/*
Whether this project's build of the DSK simulator agrees with the upstream
build away from the baseline calibration.

tests/dsk_build_equivalence.c runs the model at one parameter vector, the one
model/dsk_sfc/dsk_sfc_inputs.json ships with. The experiment runs it at a
thousand others: dataset/abm_system_design.csv is a Latin hypercube over nine
parameters, and every one of the million runs uses a point of that design. A
change can be exact at the baseline and wrong elsewhere, because the parameters
decide how many firms enter and exit, how many machines a firm holds, and
therefore which branches the model takes and how sparse the machine arrays are.
Several changes in docs/DSK_MODEL_CHANGES.md turn on exactly that sparsity.

Four points are used, all inside the design's own box:

- the design's column-wise minima, and its column-wise maxima, which are the
  corners of the box the design actually reached and the most demanding points
  in it;
- two rows of the design itself, the first and the middle one, which are points
  the experiment will really simulate.

For each point both builds are run and required to agree on three things: the
exit status, the results file byte for byte, and the error log byte for byte. A
parameter vector the model refuses to simulate is a valid outcome and still a
comparison, as long as both builds refuse it the same way at the same period.

Writes out/dsk_design_equivalence.txt. One seed per point unless a seed count
is given as the first argument.
*/

#define _XOPEN_SOURCE 700

#include "applications/abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./json.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#include "tests/dsk_upstream_scratch.h"

#define MODIFIED "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define DESIGN "dataset/abm_system_design.csv"
#define REPORT "out/dsk_design_equivalence.txt"
#define RUN_NAME "design"
#define N_POINTS 4
#define DEFAULT_SEEDS 1

/* The baseline parameter file with the nine design parameters overwritten,
   which is what applications/abm_system_simulate.c writes for a run. */
static void write_inputs(const char *base_json, const char *path, const double *point) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 && "dsk_design_equivalence: no params block in the base JSON");

    JsonValue *block = json_array_get(params, 0);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
        assert(json_object_get(block, abm_system_parameter_names()[p]) &&
               "dsk_design_equivalence: the base JSON has no such parameter");
        json_object_set(block, abm_system_parameter_names()[p], json_number(point[p]));
    }

    json_write_file(inputs, path);
    json_free(inputs);
}

static int run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    return system(command);
}

/* A file the run may or may not have written: a point the model refuses leaves
   no results, and an untroubled run leaves an empty error log. */
static char *read_optional(const char *path, long *size) {
    FILE *f = fopen(path, "rb");
    if (!f) { *size = -1; return NULL; }

    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size + 1);
    assert(text && "dsk_design_equivalence: cannot hold a file");
    if (*size > 0) assert(fread(text, 1, (size_t)*size, f) == (size_t)*size);
    text[*size] = '\0';
    fclose(f);
    unlink(path);
    return text;
}

static int compare(FILE *report, const char *label, const char *what,
                   const char *upstream_dir, const char *modified_dir,
                   const char *relative) {
    char upstream_path[1400], modified_path[1400];
    snprintf(upstream_path, sizeof upstream_path, "%s/%s", upstream_dir, relative);
    snprintf(modified_path, sizeof modified_path, "%s/%s", modified_dir, relative);

    long upstream_size, modified_size;
    char *upstream_text = read_optional(upstream_path, &upstream_size);
    char *modified_text = read_optional(modified_path, &modified_size);

    int bad = 0;
    if (upstream_size != modified_size) {
        fprintf(report, "    %s %s: %ld bytes upstream, %ld here\n", label, what, upstream_size, modified_size);
        bad = 1;
    } else if (upstream_size > 0 &&
               memcmp(upstream_text, modified_text, (size_t)upstream_size) != 0) {
        fprintf(report, "    %s %s DIFFERS over %ld bytes\n", label, what, upstream_size);
        bad = 1;
    }

    free(upstream_text);
    free(modified_text);
    return bad;
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;
    assert(n_seeds >= 1 && "dsk_design_equivalence: at least one seed");

    char modified[PATH_MAX], upstream[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(MODIFIED, modified) &&
           "dsk_design_equivalence: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(UPSTREAM, upstream) &&
           "dsk_design_equivalence: bin/dsk_SFC_upstream is not built - run make model-upstream");
    assert(realpath(INPUTS, base_json) && "dsk_design_equivalence: the model's own inputs JSON is missing");

    DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
    assert(design.r >= 1 && "dsk_design_equivalence: the design is empty - run make app-abm_system_design");

    Mat column[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        column[p] = df_col_numeric(&design, abm_system_parameter_names()[p]);

    double point[N_POINTS][ABM_SYSTEM_N_PARAMETERS];
    const char *point_name[N_POINTS] = {
        "lowest of each parameter", "highest of each parameter",
        "design row 1", "design row middle"
    };
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
        double low = (double)AT(column[p], 0, 0), high = low;
        for (int r = 1; r < design.r; r++) {
            double v = (double)AT(column[p], r, 0);
            if (v < low) low = v;
            if (v > high) high = v;
        }
        point[0][p] = low;
        point[1][p] = high;
        point[2][p] = (double)AT(column[p], 0, 0);
        point[3][p] = (double)AT(column[p], design.r / 2, 0);
    }

    DskScratch scratch = dsk_scratch_open("de", RUN_NAME, base_json, upstream, modified, n_seeds);
    const char *upstream_dir = scratch.upstream, *modified_dir = scratch.modified;
    const char *root = scratch.root;
    char point_json[720];

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_design_equivalence: cannot open the report path");
    fprintf(report, "upstream build: %s\n", upstream);
    fprintf(report, "this project's: %s\n", modified);
    fprintf(report, "design:         %s, %d configurations\n", DESIGN, design.r);
    fprintf(report, "scratch:        %s, %d characters of padding so upstream writes\n",
            scratch.root, scratch.padding);
    fprintf(report, "%d parameter points, seeds 1 to %d, comparing exit status, results and error log\n\n",
            N_POINTS, n_seeds);

    int failures = 0, comparisons = 0;
    for (int k = 0; k < N_POINTS; k++) {
        snprintf(point_json, sizeof point_json, "%s/point_%d.json", root, k);
        write_inputs(base_json, point_json, point[k]);

        fprintf(report, "  %s\n", point_name[k]);
        for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
            fprintf(report, "    %-10s %.6f\n", abm_system_parameter_names()[p], point[k][p]);

        for (int seed = 1; seed <= n_seeds; seed++) {
            char label[64];
            snprintf(label, sizeof label, "seed %d", seed);

            int upstream_status = run(upstream_dir, point_json, seed);
            int modified_status = run(modified_dir, point_json, seed);

            int bad = 0;
            if (upstream_status != modified_status) {
                fprintf(report, "    %s exit status differs: upstream %d, this build %d\n",
                        label, upstream_status, modified_status);
                bad = 1;
            }

            char results[128], errors[128];
            snprintf(results, sizeof results, "output/results_%s_%d.txt", RUN_NAME, seed);
            snprintf(errors, sizeof errors, "output/errors/Errors_%s_%d.txt", RUN_NAME, seed);
            bad |= compare(report, label, "results", upstream_dir, modified_dir, results);
            bad |= compare(report, label, "error log", upstream_dir, modified_dir, errors);

            if (!bad)
                fprintf(report, "    %s identical, exit status %d\n", label, upstream_status);
            failures += bad;
            comparisons++;
            fflush(report);
        }
        fprintf(report, "\n");
        unlink(point_json);
    }

    fprintf(report, "%d point-seed comparisons\n%s\n", comparisons, failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    /* The columns are views into the frame, so freeing the frame is enough. */
    df_free(&design);

    printf("%d parameter points times %d seeds compared\n", N_POINTS, n_seeds);
    printf("%s, %d failures\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
