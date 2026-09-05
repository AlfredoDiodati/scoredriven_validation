/*
Whether this project's build of the DSK simulator reads and writes only memory
it owns.

The other three DSK tests compare output. That cannot find the mistake this one
looks for. The speed work in docs/DSK_MODEL_CHANGES.md replaced array subscripts
with pointer arithmetic in several hot loops - `A_s[(vintage-1)*N1+(supplier-1)]`
in COSTPROD(), `bank_ratio_s[(j-1)*NB+(i-1)]` in LOANRATES(), a walk down g_c by
a fixed stride - and it carries two sentinel indices, `chosen` and `dearest`,
which are -1 until a candidate is found. An index that runs off the end of an
array reads whatever happens to sit next to it. On the seeds the equivalence
tests cover, that could read a plausible number and produce output which
compares equal to upstream's, and the defect would surface only somewhere in the
million runs the experiment actually performs.

So this asks the compiler instead of the output. bin/dsk_SFC_sanitized is the
same source built under AddressSanitizer and UndefinedBehaviorSanitizer, which
check every access against the bounds of the object it belongs to, and also
report signed overflow, bad shifts and invalid conversions. Any report is a
failure here even when the run's numbers are right.

The run's numbers are checked too: the sanitized build has to produce the same
results file as the ordinary one. That catches a build whose sanitizers are
silent because they are not actually enabled.

Seeds 1 to N, N = 3 unless given as the first argument, at the baseline
calibration and at the two corners of the parameter design. Writes
out/dsk_memory_safety.txt.
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

#define SANITIZED "bin/dsk_SFC_sanitized"
#define ORDINARY "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define DESIGN "dataset/abm_system_design.csv"
#define REPORT "out/dsk_memory_safety.txt"
#define RUN_NAME "safety"
#define DEFAULT_SEEDS 3
#define N_POINTS 3

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_memory_safety: mkdir failed");
}

static void prepare(const char *dir, const char *executable) {
    make_directory(dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(executable, link) == 0 && "dsk_memory_safety: cannot link a build into place");
}

/* The parameter file with the nine design parameters overwritten, as
   applications/abm_system_simulate.c writes it for a run. */
static void write_inputs(const char *base_json, const char *path, const double *point) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 && "dsk_memory_safety: no params block in the base JSON");

    JsonValue *block = json_array_get(params, 0);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        json_object_set(block, abm_system_parameter_names()[p], json_number(point[p]));

    json_write_file(inputs, path);
    json_free(inputs);
}

/* Everything the sanitized run wrote to stderr, which is where both sanitizers
   report. An empty file is the passing case. */
static int sanitizer_reports(const char *dir, const char *inputs, int seed, char *first, size_t first_size) {
    char log[1024], command[8192];
    snprintf(log, sizeof log, "%s/sanitizer.log", dir);
    snprintf(command, sizeof command,
             "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>\"%s\"",
             dir, inputs, RUN_NAME, seed, log);
    int status = system(command);
    (void)status;

    FILE *f = fopen(log, "r");
    if (!f) return 0;

    int reports = 0;
    char line[4096];
    first[0] = '\0';
    while (fgets(line, sizeof line, f)) {
        if (strstr(line, "ERROR: AddressSanitizer") || strstr(line, "ERROR: LeakSanitizer") ||
            strstr(line, "runtime error:")) {
            reports++;
            if (first[0] == '\0') snprintf(first, first_size, "%s", line);
        }
    }
    fclose(f);
    unlink(log);
    return reports;
}

static char *read_results(const char *dir, int seed, long *size) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, seed);

    FILE *f = fopen(path, "rb");
    if (!f) { *size = -1; return NULL; }
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size + 1);
    assert(text && "dsk_memory_safety: cannot hold a results file");
    if (*size > 0) assert(fread(text, 1, (size_t)*size, f) == (size_t)*size);
    text[*size] = '\0';
    fclose(f);
    unlink(path);
    return text;
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;
    assert(n_seeds >= 1 && "dsk_memory_safety: at least one seed");

    char sanitized[PATH_MAX], ordinary[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(SANITIZED, sanitized) &&
           "dsk_memory_safety: bin/dsk_SFC_sanitized is not built - run make model-sanitized");
    assert(realpath(ORDINARY, ordinary) &&
           "dsk_memory_safety: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, base_json) && "dsk_memory_safety: the model's own inputs JSON is missing");

    /* The baseline, and the two corners of the box the design actually reached:
       the parameters decide how many firms enter and exit and how sparse the
       machine arrays are, which is what the pointer walks index through. */
    double point[N_POINTS][ABM_SYSTEM_N_PARAMETERS];
    const char *point_name[N_POINTS] = {
        "baseline calibration", "lowest of each parameter", "highest of each parameter"
    };
    int have_design = 0;

    DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
    if (design.r >= 1) {
        have_design = 1;
        for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
            Mat column = df_col_numeric(&design, abm_system_parameter_names()[p]);
            double low = (double)AT(column, 0, 0), high = low;
            for (int r = 1; r < design.r; r++) {
                double v = (double)AT(column, r, 0);
                if (v < low) low = v;
                if (v > high) high = v;
            }
            point[1][p] = low;
            point[2][p] = high;
        }
    }

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], sanitized_dir[640], ordinary_dir[640], point_json[720];
    snprintf(root, sizeof root, "%s/dsk_safety_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(sanitized_dir, sizeof sanitized_dir, "%s/sanitized", root);
    snprintf(ordinary_dir, sizeof ordinary_dir, "%s/ordinary", root);

    prepare(sanitized_dir, sanitized);
    prepare(ordinary_dir, ordinary);

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_memory_safety: cannot open the report path");
    fprintf(report, "sanitized build: %s\n", sanitized);
    fprintf(report, "ordinary build:  %s\n", ordinary);
    fprintf(report, "AddressSanitizer and UndefinedBehaviorSanitizer, seeds 1 to %d\n", n_seeds);
    fprintf(report, "at %d parameter points; any report is a failure\n\n",
            have_design ? N_POINTS : 1);

    int failures = 0, runs = 0;
    int n_points = have_design ? N_POINTS : 1;

    for (int k = 0; k < n_points; k++) {
        const char *inputs = base_json;
        if (k > 0) {
            snprintf(point_json, sizeof point_json, "%s/point_%d.json", root, k);
            write_inputs(base_json, point_json, point[k]);
            inputs = point_json;
        }
        fprintf(report, "  %s\n", point_name[k]);

        for (int seed = 1; seed <= n_seeds; seed++) {
            char first[4096];
            int reports = sanitizer_reports(sanitized_dir, inputs, seed, first, sizeof first);
            runs++;

            if (reports > 0) {
                fprintf(report, "    seed %-4d %d sanitizer report(s)\n", seed, reports);
                fprintf(report, "      %s", first);
                failures++;
                continue;
            }

            /* Silent sanitizers prove nothing if the build did not really run
               the model, so its numbers are compared with the ordinary build's. */
            char command[8192];
            snprintf(command, sizeof command,
                     "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
                     ordinary_dir, inputs, RUN_NAME, seed);
            int status = system(command);
            assert(status == 0 && "dsk_memory_safety: the ordinary build did not complete");

            long sanitized_size, ordinary_size;
            char *sanitized_text = read_results(sanitized_dir, seed, &sanitized_size);
            char *ordinary_text = read_results(ordinary_dir, seed, &ordinary_size);

            int identical = (sanitized_size == ordinary_size && sanitized_size > 0 &&
                             memcmp(sanitized_text, ordinary_text, (size_t)sanitized_size) == 0);
            if (identical) {
                fprintf(report, "    seed %-4d clean, %ld bytes matching the ordinary build\n",
                        seed, sanitized_size);
            } else {
                fprintf(report, "    seed %-4d clean but the sanitized build's output differs\n", seed);
                failures++;
            }

            free(sanitized_text);
            free(ordinary_text);
            fflush(report);
        }
        if (k > 0) unlink(point_json);
        fprintf(report, "\n");
    }

    if (!have_design)
        fprintf(report, "dataset/abm_system_design.csv is missing, so only the baseline was used\n");
    fprintf(report, "%d runs under the sanitizers\n%s\n", runs, failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    df_free(&design);

    printf("%d sanitized runs over %d parameter points\n", runs, n_points);
    printf("%s, %d failures\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
