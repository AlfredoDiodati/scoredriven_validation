/*
Whether this project's build of the DSK simulator agrees with the upstream
build on everything the model can report, not only on the aggregate file.

tests/dsk_build_equivalence.c compares the 83-column results file, which holds
one row per period of economy-wide totals. Run with -f 1 the model writes
twelve more files instead, and those hold one column per firm: each K-firm's
and C-firm's productivity, energy efficiency, environmental friendliness and
net worth, each bank's net worth, and each C-firm's debt. A change that moved
one firm's machines to another firm while leaving the totals alone would pass
the aggregate comparison and fail this one.

That is the case worth guarding. Several of the changes in
docs/DSK_MODEL_CHANGES.md rest on an argument about which values are ever read
- ages of machines nobody holds, working copies that only a disabled branch
consumes - and the per-firm files are where a mistake in such an argument would
first become visible.

Every regular file both runs write is compared byte for byte, and the two runs
are required to have written the same set of names. The error logs are compared
too, so a run that fails has to fail the same way in both builds.

Seeds 1 to N, N = 2 unless given as the first argument. Writes
out/dsk_full_output_equivalence.txt.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <assert.h>

#include "tests/dsk_upstream_scratch.h"

#define MODIFIED "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_full_output_equivalence.txt"
#define RUN_NAME "full"
#define DEFAULT_SEEDS 2
#define MAX_FILES 64

static int run(const char *dir, const char *inputs, int seed) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 1 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, seed);
    return system(command);
}

/* The files one run wrote under its output directory, sorted, so two runs'
   listings can be compared name by name. One level of subdirectory is walked
   because the error log sits in output/errors/, and a run that fails has to
   fail the same way in both builds. */
static int list_output(const char *dir, char names[MAX_FILES][256]) {
    char root[1024];
    snprintf(root, sizeof root, "%s/output", dir);

    DIR *handle = opendir(root);
    assert(handle && "dsk_full_output_equivalence: a run produced no output directory");

    int count = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        if (entry->d_name[0] == '.') continue;

        char full[1400];
        snprintf(full, sizeof full, "%s/%.255s", root, entry->d_name);
        struct stat info;
        if (stat(full, &info) != 0) continue;

        if (S_ISREG(info.st_mode)) {
            assert(count < MAX_FILES && "dsk_full_output_equivalence: more output files than expected");
            snprintf(names[count], sizeof names[0], "%s", entry->d_name);
            count++;
            continue;
        }
        if (!S_ISDIR(info.st_mode)) continue;

        DIR *inner = opendir(full);
        if (!inner) continue;
        struct dirent *below;
        while ((below = readdir(inner)) != NULL) {
            if (below->d_name[0] == '.') continue;

            char deep[1700];
            snprintf(deep, sizeof deep, "%s/%.255s", full, below->d_name);
            struct stat below_info;
            if (stat(deep, &below_info) != 0 || !S_ISREG(below_info.st_mode)) continue;

            assert(count < MAX_FILES && "dsk_full_output_equivalence: more output files than expected");
            snprintf(names[count], sizeof names[0], "%.120s/%.120s", entry->d_name, below->d_name);
            count++;
        }
        closedir(inner);
    }
    closedir(handle);

    for (int a = 0; a < count; a++)
        for (int b = a + 1; b < count; b++)
            if (strcmp(names[b], names[a]) < 0) {
                char swap[256];
                snprintf(swap, sizeof swap, "%s", names[a]);
                snprintf(names[a], sizeof names[0], "%s", names[b]);
                snprintf(names[b], sizeof names[0], "%s", swap);
            }
    return count;
}

static char *read_file(const char *dir, const char *name, long *size) {
    char path[1400];
    snprintf(path, sizeof path, "%s/output/%.255s", dir, name);

    FILE *f = fopen(path, "rb");
    assert(f && "dsk_full_output_equivalence: cannot open an output file");
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size + 1);
    assert(text && "dsk_full_output_equivalence: cannot hold an output file");
    assert(fread(text, 1, (size_t)*size, f) == (size_t)*size ||
           *size == 0);
    text[*size] = '\0';
    fclose(f);
    return text;
}

static void clear_output(const char *dir, char names[MAX_FILES][256], int count) {
    for (int k = 0; k < count; k++) {
        char path[1400];
        snprintf(path, sizeof path, "%s/output/%.255s", dir, names[k]);
        unlink(path);
    }
}

int main(int argc, char **argv) {
    int n_seeds = argc > 1 ? atoi(argv[1]) : DEFAULT_SEEDS;
    assert(n_seeds >= 1 && "dsk_full_output_equivalence: at least one seed");

    char modified[PATH_MAX], upstream[PATH_MAX], inputs[PATH_MAX];
    assert(realpath(MODIFIED, modified) &&
           "dsk_full_output_equivalence: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(UPSTREAM, upstream) &&
           "dsk_full_output_equivalence: bin/dsk_SFC_upstream is not built - run make model-upstream");
    assert(realpath(INPUTS, inputs) && "dsk_full_output_equivalence: the model's own inputs JSON is missing");

    DskScratch scratch = dsk_scratch_open("fo", RUN_NAME, inputs, upstream, modified, n_seeds);
    const char *upstream_dir = scratch.upstream, *modified_dir = scratch.modified;

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_full_output_equivalence: cannot open the report path");
    fprintf(report, "upstream build: %s\n", upstream);
    fprintf(report, "this project's: %s\n", modified);
    fprintf(report, "inputs:         %s\n", inputs);
    fprintf(report, "scratch:        %s, %d characters of padding so upstream writes\n",
            scratch.root, scratch.padding);
    fprintf(report, "seeds 1 to %d, run with -f 1, comparing every file both runs write\n\n", n_seeds);

    int failures = 0, files_compared = 0;
    for (int seed = 1; seed <= n_seeds; seed++) {
        int upstream_status = run(upstream_dir, inputs, seed);
        int modified_status = run(modified_dir, inputs, seed);

        if (upstream_status != modified_status) {
            fprintf(report, "  seed %-4d exit status differs: upstream %d, this build %d\n",
                    seed, upstream_status, modified_status);
            failures++;
            continue;
        }

        char upstream_names[MAX_FILES][256], modified_names[MAX_FILES][256];
        int upstream_count = list_output(upstream_dir, upstream_names);
        int modified_count = list_output(modified_dir, modified_names);

        if (upstream_count != modified_count) {
            fprintf(report, "  seed %-4d wrote %d files upstream and %d here\n",
                    seed, upstream_count, modified_count);
            failures++;
            clear_output(upstream_dir, upstream_names, upstream_count);
            clear_output(modified_dir, modified_names, modified_count);
            continue;
        }

        int seed_failures = 0;
        long seed_bytes = 0;
        for (int k = 0; k < upstream_count; k++) {
            if (strcmp(upstream_names[k], modified_names[k]) != 0) {
                fprintf(report, "  seed %-4d file names differ: %s against %s\n",
                        seed, upstream_names[k], modified_names[k]);
                seed_failures++;
                continue;
            }

            long upstream_size, modified_size;
            char *upstream_text = read_file(upstream_dir, upstream_names[k], &upstream_size);
            char *modified_text = read_file(modified_dir, modified_names[k], &modified_size);

            int identical = (upstream_size == modified_size &&
                             memcmp(upstream_text, modified_text, (size_t)upstream_size) == 0);
            if (!identical) {
                fprintf(report, "  seed %-4d %s DIFFERS, %ld bytes upstream, %ld here\n",
                        seed, upstream_names[k], upstream_size, modified_size);
                seed_failures++;
            }
            seed_bytes += upstream_size;
            files_compared++;

            free(upstream_text);
            free(modified_text);
        }

        if (seed_failures == 0)
            fprintf(report, "  seed %-4d %d files identical, %ld bytes\n", seed, upstream_count, seed_bytes);
        failures += seed_failures;

        clear_output(upstream_dir, upstream_names, upstream_count);
        clear_output(modified_dir, modified_names, modified_count);
        fflush(report);
    }

    fprintf(report, "\n%d files compared over %d seeds\n%s\n",
            files_compared, n_seeds, failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    printf("%d files compared over %d seeds, %d differing\n", files_compared, n_seeds, failures);
    printf("%s, %d failures\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
