/*
The per-firm output quantities behind GDP, for a few of the runs the experiment
already holds.

dataset/abm_system stores five macro series per replicate and nothing below
them. GDP is Q1.Sum() * dim_mach + Q2.Sum() over the K-firms and the C-firms,
so the cross-section those sums run over is what a central limit theorem would
act on, and it is not in the archives. The simulator writes it when given -m 1,
a switch this project added; docs/DSK_MODEL_CHANGES.md records it and
tests/dsk_build_equivalence.c is what says the switch changed no number.

Configurations are spread evenly over dataset/abm_system_design.csv, first and
last included. Seeds are 1 to SEEDS, which are replicates 0 to SEEDS-1 of those
same configurations in dataset/abm_system, so a run here is a run the archives
already hold and its macro series can be checked against the stored one.

Each run writes three files into dataset/abm_system_micro/cop_XXXX/:

    Q1all_<seed>.txt    600 rows, one per period, N1 columns of K-firm output
    Q2all_<seed>.txt    600 rows, one per period, N2 columns of C-firm output
    results_<seed>.txt  the 83-column aggregate file, kept so the per-firm rows
                        can be summed back against the GDP column they make up

Periods 1 to 600 are all written. Which of them to read is the reader's
choice; studies/abm_system_micro_clt.c takes 201 to 600, the window the
archives keep.

    ./bin/abm_system_micro_simulate

Needs `make model`. Existing files are left alone, so a second run is a no-op
and deleting a configuration's directory is how to redo it.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#include "applications/abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./json.h>

#define MODEL "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define DESIGN "dataset/abm_system_design.csv"
#define OUTPUT_DIR "dataset/abm_system_micro"
#define RUN_NAME "m"
#define N_STEPS 600
#define CONFIGURATIONS 5
#define SEEDS 5

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0)
        assert(errno == EEXIST && "abm_system_micro_simulate: mkdir failed");
}

static int file_exists(const char *path) { return access(path, F_OK) == 0; }

static void write_inputs(const char *base_json, const char *path, const Mat *parameter,
                         int design_row) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 &&
           "abm_system_micro_simulate: no params block");

    JsonValue *block = json_array_get(params, 0);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        json_object_set(block, abm_system_parameter_names()[p],
                        json_number((double)AT(parameter[p], design_row, 0)));
    json_object_set(block, "T", json_number(N_STEPS));

    json_write_file(inputs, path);
    json_free(inputs);
}

/* The scratch directory is usually tmpfs and the dataset is not, so a rename
   across the two fails and the bytes have to be copied. */
static void take(const char *from, const char *to) {
    if (rename(from, to) == 0) return;
    assert(errno == EXDEV && "abm_system_micro_simulate: cannot move a run's output");

    FILE *in = fopen(from, "rb");
    FILE *out = fopen(to, "wb");
    assert(in && out && "abm_system_micro_simulate: cannot copy a run's output");

    char block[1 << 16];
    size_t read;
    while ((read = fread(block, 1, sizeof block, in)) > 0)
        assert(fwrite(block, 1, read, out) == read &&
               "abm_system_micro_simulate: a copy did not complete");

    fclose(in);
    assert(fclose(out) == 0 && "abm_system_micro_simulate: a copy did not complete");
    unlink(from);
}

int main(void) {
    char model[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "abm_system_micro_simulate: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, base_json) &&
           "abm_system_micro_simulate: the parameter file is missing");

    DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
    Mat parameter[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        parameter[p] = df_col_numeric(&design, abm_system_parameter_names()[p]);
    assert(design.r >= CONFIGURATIONS && "abm_system_micro_simulate: the design is too short");

    make_directory(OUTPUT_DIR);

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    /* Short directory name: upstream builds every output filename in a 64-byte
       buffer, and this build runs in the same places as the tests that check
       against it. */
    char root[128];
    snprintf(root, sizeof root, "%s/mi%d", scratch, (int)getpid());
    make_directory(root);

    char link[192];
    snprintf(link, sizeof link, "%s/dsk_SFC", root);
    unlink(link);
    assert(symlink(model, link) == 0 && "abm_system_micro_simulate: cannot link the build into place");

    int written = 0, skipped = 0;
    for (int k = 0; k < CONFIGURATIONS; k++) {
        int cop = 1 + (int)((long)k * (design.r - 1) / (CONFIGURATIONS - 1));

        char cop_dir[512];
        snprintf(cop_dir, sizeof cop_dir, "%s/cop_%04d", OUTPUT_DIR, cop);
        make_directory(cop_dir);

        char json_path[256];
        snprintf(json_path, sizeof json_path, "%s/c%d.json", root, cop);
        write_inputs(base_json, json_path, parameter, cop - 1);

        for (int seed = 1; seed <= SEEDS; seed++) {
            char kept_q1[640], kept_q2[640], kept_results[640];
            snprintf(kept_q1, sizeof kept_q1, "%s/Q1all_%d.txt", cop_dir, seed);
            snprintf(kept_q2, sizeof kept_q2, "%s/Q2all_%d.txt", cop_dir, seed);
            snprintf(kept_results, sizeof kept_results, "%s/results_%d.txt", cop_dir, seed);
            if (file_exists(kept_q1) && file_exists(kept_q2) && file_exists(kept_results)) {
                skipped++;
                continue;
            }

            char command[8192];
            snprintf(command, sizeof command,
                     "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -m 1 -c 0 -v 0 >/dev/null 2>&1",
                     root, json_path, RUN_NAME, seed);
            assert(system(command) == 0 && "abm_system_micro_simulate: a run did not complete");

            char produced[640];
            snprintf(produced, sizeof produced, "%s/output/Q1all_%s_%d.txt", root, RUN_NAME, seed);
            take(produced, kept_q1);
            snprintf(produced, sizeof produced, "%s/output/Q2all_%s_%d.txt", root, RUN_NAME, seed);
            take(produced, kept_q2);
            snprintf(produced, sizeof produced, "%s/output/results_%s_%d.txt", root, RUN_NAME, seed);
            take(produced, kept_results);

            snprintf(produced, sizeof produced, "%s/output/errors/Errors_%s_%d.txt", root, RUN_NAME, seed);
            unlink(produced);
            written++;
        }
    }

    printf("micro simulate: %d configurations x %d seeds, %d runs written, %d already there\n",
           CONFIGURATIONS, SEEDS, written, skipped);

    df_free(&design);
    return EXIT_SUCCESS;
}
