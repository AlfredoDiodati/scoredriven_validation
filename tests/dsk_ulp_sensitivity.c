/*
How small a difference in the DSK model's arithmetic the equivalence tests can
still see.

The other DSK tests compare the files the model writes. Those files are rounded:
the aggregate results file prints ten decimal places, the per-firm files four.
"Byte for byte identical" therefore means identical at the precision printed,
not identical in every bit of every double. Somebody may reasonably ask whether
the speed work changed the arithmetic by an amount too small to print, which
would still be a change to the model.

This measures how small that amount would have to be. One input parameter at a
time is multiplied by 1+eps, and eps is walked up a ladder from 1e-16, the size
of a rounding error in a double, to 1e-6. The smallest eps whose run no longer
matches the unperturbed one, and the period at which the difference first
appears, are what the test reports. A parameter whose threshold sits at 1e-15
means the printed output reacts to a change of that relative size somewhere in
that parameter's path, so a difference that hides behind the printed digits has
to be smaller than that.

The change reaches the model through the parameter JSON, which is written with
seventeen significant digits and so carries a change of this size intact.

What the bar is on is the results file, not any one parameter. Gamma reaches
the model in exactly one place, multiplied by a client count and rounded to an
integer, so a change below half a client produces the identical integer and
the identical run: no output format could resolve it, and its coarse threshold
says nothing about the file. The test therefore requires that two thirds of the
parameter and seed pairs resolve 1e-13 or finer, which a file printed at fewer
digits would fail outright while one quantised parameter leaves intact.

Seeds 1 to N, N = 2 unless given as the first argument, over all nine design
parameters. Writes out/dsk_ulp_sensitivity.txt.
*/

#define _XOPEN_SOURCE 700

#include "applications/abm_system.h"

#include <et_al./json.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#define MODEL "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_ulp_sensitivity.txt"
#define RUN_NAME "ulp"
#define DEFAULT_SEEDS 2
#define N_RUNGS 15
#define FINE_THRESHOLD 1e-13
#define REQUIRED_SHARE (2.0 / 3.0)

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_ulp_sensitivity: mkdir failed");
}

/* The baseline parameter file with one parameter replaced. A NULL name writes
   the baseline unchanged. */
static void write_perturbed(const char *base_json, const char *path, const char *name, double value) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 && "dsk_ulp_sensitivity: no params block in the base JSON");

    JsonValue *block = json_array_get(params, 0);
    if (name) {
        assert(json_object_get(block, name) && "dsk_ulp_sensitivity: the base JSON has no such parameter");
        json_object_set(block, name, json_number(value));
    }

    json_write_file(inputs, path);
    json_free(inputs);
}

static double baseline_value(const char *base_json, const char *name) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *block = json_array_get(json_object_get(inputs, "params"), 0);
    JsonValue *field = json_object_get(block, name);
    assert(field && "dsk_ulp_sensitivity: the base JSON has no such parameter");
    double value = json_as_number(field);
    json_free(inputs);
    return value;
}

static void run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command,
             "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    assert(system(command) == 0 && "dsk_ulp_sensitivity: a run did not complete");
}

static char *take_results(const char *dir, int seed) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, seed);

    FILE *f = fopen(path, "rb");
    assert(f && "dsk_ulp_sensitivity: a run wrote no results file");
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)size + 1);
    assert(text && "dsk_ulp_sensitivity: cannot hold a results file");
    if (size > 0) assert(fread(text, 1, (size_t)size, f) == (size_t)size);
    text[size] = 0;
    fclose(f);
    unlink(path);
    return text;
}

/* The line number, counting from one, of the first place the two texts differ,
   or 0 if they are the same. Line number is period number here: the results
   file carries one line per simulated period. */
static int first_differing_period(const char *a, const char *b) {
    int period = 1;
    for (long i = 0;; i++) {
        if (a[i] != b[i]) return period;
        if (a[i] == '\0') return 0;
        if (a[i] == '\n') period++;
    }
}

static int count_periods(const char *text) {
    int lines = 0;
    for (long i = 0; text[i]; i++)
        if (text[i] == '\n') lines++;
    return lines;
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;
    assert(n_seeds >= 1 && "dsk_ulp_sensitivity: at least one seed");

    char model[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(MODEL, model) && "dsk_ulp_sensitivity: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, base_json) && "dsk_ulp_sensitivity: the model's own inputs JSON is missing");

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], dir[640], json_base[720], json_test[720], link[800];
    snprintf(root, sizeof root, "%s/dsk_ulp_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(dir, sizeof dir, "%s/run", root);
    make_directory(dir);
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    assert(symlink(model, link) == 0 && "dsk_ulp_sensitivity: cannot link the model into place");

    snprintf(json_base, sizeof json_base, "%s/unperturbed.json", root);
    snprintf(json_test, sizeof json_test, "%s/perturbed.json", root);
    write_perturbed(base_json, json_base, NULL, 0);

    double rung[N_RUNGS];
    for (int r = 0; r < N_RUNGS; r++) rung[r] = pow(10.0, -16.0 + r);

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_ulp_sensitivity: cannot open the report path");
    fprintf(report, "build: %s\n", model);
    fprintf(report, "one parameter at a time multiplied by 1+eps, eps walked from %g to %g\n",
            rung[0], rung[N_RUNGS - 1]);
    fprintf(report, "the threshold is the smallest eps whose results file stops matching the\n");
    fprintf(report, "unperturbed run; the results file prints ten decimal places per column\n");
    fprintf(report, "at least %.0f%% of the parameter and seed pairs must resolve %g or finer\n\n",
            100 * REQUIRED_SHARE, FINE_THRESHOLD);

    int fine = 0, pairs = 0, inert = 0, runs = 0;
    double worst = 0;

    for (int seed = 1; seed <= n_seeds; seed++) {
        run(dir, json_base, seed);
        runs++;
        char *unperturbed = take_results(dir, seed);
        int periods = count_periods(unperturbed);

        fprintf(report, "  seed %d, %d periods\n", seed, periods);

        for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
            const char *name = abm_system_parameter_names()[p];
            double base = baseline_value(base_json, name);

            int found = 0;
            for (int r = 0; r < N_RUNGS && !found; r++) {
                double value = base * (1.0 + rung[r]);
                if (value == base) continue;

                write_perturbed(base_json, json_test, name, value);
                run(dir, json_test, seed);
                runs++;

                char *perturbed = take_results(dir, seed);
                int period = first_differing_period(unperturbed, perturbed);
                free(perturbed);

                if (period == 0) continue;

                found = 1;
                fprintf(report, "    %-8s threshold %-8g first differs at period %d\n",
                        name, rung[r], period);
                if (rung[r] > worst) worst = rung[r];
                if (rung[r] <= FINE_THRESHOLD) fine++;
            }

            if (!found) {
                fprintf(report, "    %-8s unchanged output up to eps %g: inert under these flags\n",
                        name, rung[N_RUNGS - 1]);
                inert++;
            }
            pairs++;
            fflush(report);
        }

        free(unperturbed);
        fprintf(report, "\n");
    }

    double share = pairs > 0 ? (double)fine / pairs : 0;
    int passed = share >= REQUIRED_SHARE;

    fprintf(report, "%d model runs; %d of %d pairs resolve %g or finer, a share of %.2f\n",
            runs, fine, pairs, FINE_THRESHOLD, share);
    fprintf(report, "coarsest threshold %g; %d pairs unmoved by the whole ladder\n", worst, inert);
    fprintf(report, "%s\n", passed ? "PASSED" : "FAILED");
    fclose(report);

    unlink(json_base);
    unlink(json_test);

    printf("%d runs, %d of %d parameter and seed pairs resolve %g or finer\n",
           runs, fine, pairs, FINE_THRESHOLD);
    printf("%s, share %.2f against a bar of %.2f\n",
           passed ? "PASSED" : "FAILED", share, REQUIRED_SHARE);
    return passed ? 0 : 1;
}
