/*
Runs the DSK stock-flow-consistent model over the parameter design in
dataset/abm_system_design.csv and writes, for every (configuration,
replication) pair, the five series the auxiliary model is fitted on:

    dataset/abm_system/cop_0001/batch_000.npz

That is the layout and the five series applications/abm_system_extract.c
produces from the older .Rdata dataset, so applications/abm_system_fit_qvarma.c
reads either without knowing which wrote it and no extraction pass is needed
after this one. Both writers go through abm_system_write_batch, so the layout
is defined in one place.

Runs are written ABM_SYSTEM_BATCH at a time rather than one at a time, which
is what lets one deflate stream see ten replications of a series instead of
one and is worth about a third of the stored size. What it costs is the unit
of loss: an interrupted run loses the batch being filled, at most
ABM_SYSTEM_BATCH runs, rather than the single run in flight. The model's own
83-column output is still deleted as soon as it has been read, and nothing is
held beyond the batch being filled.

Rerunning skips every batch already on disk, which is what makes resuming and
re-submitting the same thing. A batch is skipped whole, failures included: the
model is deterministic in its seed, so a run that failed once fails the same
way again, and the failure log is where it is recorded.

The model reads one JSON parameter file and one seed per execution and writes
its output into the directory holding the executable it was invoked as. Each
process therefore reaches it through a symlink in its own scratch directory,
which is what lets several of these run at once against one build.

The design is read, not generated here: applications/abm_system_design.c
draws it, and separating the two is what lets a thousand simulation processes
share one design without any of them being able to redraw it.

Usage:

    ./bin/abm_system_simulate DSK_EXECUTABLE BASE_JSON [FIRST_COP LAST_COP [N_MC]]

FIRST_COP and LAST_COP are 1-indexed design rows and default to the whole
design; N_MC is the replications per configuration and defaults to 1000. One
process covers the whole design on one machine; one process per configuration
is what a cluster array wants.

Writes dataset/, not out/, apart from the failure logs and the manifest.
*/

#define _XOPEN_SOURCE 700

#include "abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./frame/npz.h>
#include <et_al./frame/txt.h>
#include <et_al./json.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>

#define DESIGN_PATH "dataset/abm_system_design.csv"
#define OUTPUT_DIR "dataset/abm_system"
#define REPORT_DIR "out/abm_system_simulate"
#define MANIFEST_PATH "out/abm_system_simulate_manifest.txt"

#define N_STEPS 600
#define N_MODEL_COLUMNS 83
#define DEFAULT_N_MC 1000

/* 1-indexed columns of the model's own results file, which has no header.
   Column 34 is the price level cpi(1), not column 6, which is the
   four-period gross inflation factor cpi(1)/cpi(5); column 50 is the
   central bank's policy rate, not the average commercial loan rate, which
   is column 82. */
static const int model_column[ABM_SYSTEM_K] = {
    [LEVEL_GDP] = 2,
    [LEVEL_ENERGY] = 8,
    [LEVEL_EMPLOYMENT] = 5,
    [LEVEL_PRICE] = 34,
    [LEVEL_INTEREST] = 50
};

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_simulate: mkdir failed");
}

static int file_exists(const char *path) { return access(path, F_OK) == 0; }

static long file_size(const char *path) {
    struct stat info;
    if (stat(path, &info) != 0) return -1;
    return (long)info.st_size;
}

static void copy_file(const char *from, const char *to) {
    FILE *in = fopen(from, "rb");
    if (!in) return;
    FILE *out = fopen(to, "wb");
    if (!out) { fclose(in); return; }

    char buffer[4096];
    size_t got;
    while ((got = fread(buffer, 1, sizeof buffer, in)) > 0) fwrite(buffer, 1, got, out);

    fclose(out);
    fclose(in);
}

/* Appended and closed on every call, so what is on disk is what has actually
   happened up to now rather than what a buffer would eventually hold. */
static void append_line(const char *path, const char *line) {
    FILE *f = fopen(path, "a");
    assert(f && "abm_system_simulate: cannot append to a log");
    fprintf(f, "%s\n", line);
    fclose(f);
}

static const char *scratch_root(void) {
    const char *root = getenv("SLURM_TMPDIR");
    if (!root) root = getenv("TMPDIR");
    if (!root) root = "/tmp";
    return root;
}

/* The model's parameter file for one design row: the baseline with the nine
   sampled parameters and the horizon overwritten. */
static void write_inputs(const char *base_json, const char *path,
                         const Mat *parameter, int design_row) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 && "abm_system_simulate: no params block in the base JSON");

    JsonValue *block = json_array_get(params, 0);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
        assert(json_object_get(block, abm_system_parameter_names()[p]) &&
               "abm_system_simulate: the base JSON has no such parameter");
        json_object_set(block, abm_system_parameter_names()[p], json_number((double)AT(parameter[p], design_row, 0)));
    }
    json_object_set(block, "T", json_number(N_STEPS));

    json_write_file(inputs, path);
    json_free(inputs);
}

/* The model's 83 columns for one run, reduced to the five level series
   abm_system_transform consumes. Returns 0 and leaves levels untouched if
   the run did not produce what it should have. */
static int read_levels(const char *raw_path, Mat *levels) {
    DataFrame raw = df_read_txt(raw_path, txt_read_options_default());

    int usable = (raw.r == N_STEPS && raw.n_cols == N_MODEL_COLUMNS && raw.n_string == 0);
    if (usable) {
        *levels = mat_new(ABM_SYSTEM_K, N_STEPS);
        for (int series = 0; series < ABM_SYSTEM_K; series++)
            for (int t = 0; t < N_STEPS; t++)
                AT(*levels, series, t) = AT(raw.numeric, t, model_column[series] - 1);
    }

    df_free(&raw);
    return usable;
}

static int all_finite(const Mat m) {
    for (int i = 0; i < m.r; i++)
        for (int j = 0; j < m.c; j++)
            if (MISNAN(AT(m, i, j)) || MISINF(AT(m, i, j))) return 0;
    return 1;
}

int main(int argc, char **argv) {
    assert(argc >= 3 &&
           "usage: abm_system_simulate DSK_EXECUTABLE BASE_JSON [FIRST_COP LAST_COP [N_MC]]");

    char executable[PATH_MAX];
    assert(realpath(argv[1], executable) && "abm_system_simulate: cannot resolve the executable path");
    char base_json[PATH_MAX];
    assert(realpath(argv[2], base_json) && "abm_system_simulate: cannot resolve the base JSON path");

    DataFrame design = df_read_csv(DESIGN_PATH, csv_read_options_default());
    int n_cop = design.r;
    int first_cop = argc > 3 ? atoi(argv[3]) : 1;
    int last_cop = argc > 4 ? atoi(argv[4]) : n_cop;
    int n_mc = argc > 5 ? atoi(argv[5]) : DEFAULT_N_MC;

    assert(first_cop >= 1 && last_cop <= n_cop && first_cop <= last_cop &&
           "abm_system_simulate: the requested range is not inside the design");
    assert(n_mc >= 1 && "abm_system_simulate: at least one replication per configuration");

    Mat parameter[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) parameter[p] = df_col_numeric(&design, abm_system_parameter_names()[p]);

    char scratch[256];
    snprintf(scratch, sizeof scratch, "%s/dsk_%d", scratch_root(), (int)getpid());
    make_directory(scratch);

    char local_exe[512];
    snprintf(local_exe, sizeof local_exe, "%s/dsk_SFC", scratch);
    unlink(local_exe);
    assert(symlink(executable, local_exe) == 0 && "abm_system_simulate: cannot link the executable into scratch");

    make_directory(OUTPUT_DIR);
    make_directory(REPORT_DIR);

    for (int cop = first_cop; cop <= last_cop; cop++) {
        char cop_dir[256];
        snprintf(cop_dir, sizeof cop_dir, "%s/cop_%04d", OUTPUT_DIR, cop);
        make_directory(cop_dir);

        char json_path[512];
        snprintf(json_path, sizeof json_path, "%s/inputs_cop%04d.json", scratch, cop);
        write_inputs(base_json, json_path, parameter, cop - 1);

        char failure_path[512];
        snprintf(failure_path, sizeof failure_path, "%s/cop_%04d_failures.txt", REPORT_DIR, cop);

        char run_name[32];
        snprintf(run_name, sizeof run_name, "cop%04d", cop);

        Mat block[ABM_SYSTEM_BATCH];
        int replicate[ABM_SYSTEM_BATCH];
        int held = 0;

        int completed = 0, failed = 0, skipped = 0;
        for (int mc = 0; mc < n_mc; mc++) {
            if (mc % ABM_SYSTEM_BATCH == 0) {
                char batch_path[512];
                abm_system_batch_path(batch_path, sizeof batch_path, cop_dir, mc);

                if (file_exists(batch_path)) {
                    int in_batch = n_mc - mc < ABM_SYSTEM_BATCH ? n_mc - mc : ABM_SYSTEM_BATCH;
                    skipped += in_batch;
                    mc += in_batch - 1;
                    continue;
                }
            }

            /* The same seeds for every configuration: comparisons across
               configurations are what the design exists for, and shared
               randomness reduces the variance of a difference. */
            int seed = mc + 1;

            char command[2048];
            snprintf(command, sizeof command, "\"%s\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
                     local_exe, json_path, run_name, seed);
            int status = system(command);

            char raw_path[640], error_path[640];
            snprintf(raw_path, sizeof raw_path, "%s/output/results_%s_%d.txt", scratch, run_name, seed);
            snprintf(error_path, sizeof error_path, "%s/output/errors/Errors_%s_%d.txt", scratch, run_name, seed);

            const char *reason = NULL;
            Mat levels = {0};
            int have_levels = 0;

            if (status != 0) reason = "nonzero exit status";
            else if (file_size(error_path) > 0) reason = "non-empty model error log";
            else if (!file_exists(raw_path)) reason = "no results file";
            else if (!read_levels(raw_path, &levels)) reason = "results file has the wrong shape";
            else {
                have_levels = 1;
                if (!all_finite(levels)) reason = "non-finite model output";
            }

            if (!reason) {
                Mat y = abm_system_transform(levels, ABM_SYSTEM_BURN_IN);
                if (all_finite(y)) {
                    block[held] = y;
                    replicate[held] = mc;
                    held++;
                    completed++;
                } else {
                    reason = "non-finite transformed series";
                    mat_free(y);
                }
            }

            if (have_levels) mat_free(levels);

            if (reason) {
                failed++;
                char line[256];
                snprintf(line, sizeof line, "seed %d: %s", seed, reason);
                append_line(failure_path, line);

                if (file_size(error_path) > 0) {
                    char kept[512];
                    snprintf(kept, sizeof kept, "%s/cop%04d_seed%04d_error.txt", REPORT_DIR, cop, seed);
                    copy_file(error_path, kept);
                }
            }

            /* About 3 MB of 83-column output per run, none of it needed once
               the five series are in hand. */
            unlink(raw_path);
            unlink(error_path);

            if ((mc % ABM_SYSTEM_BATCH == ABM_SYSTEM_BATCH - 1 || mc == n_mc - 1) && held > 0) {
                abm_system_write_batch(cop_dir, block, replicate, held);
                for (int b = 0; b < held; b++) mat_free(block[b]);
                held = 0;
            }
        }

        char line[256];
        snprintf(line, sizeof line, "cop_%04d completed %d failed %d skipped %d of %d",
                 cop, completed, failed, skipped, n_mc);
        append_line(MANIFEST_PATH, line);

        unlink(json_path);
    }

    unlink(local_exe);
    df_free(&design);
    return 0;
}
