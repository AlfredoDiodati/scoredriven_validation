/*
Converts every replicate of every file under dataset/simulated/ into the
five-series layout abm_system.h's own header comment describes, and writes
each one into a compressed .npz archive - a data-preparation step, not a
fit: nothing here fits anything. The fitting stage reads these back rather
than re-parsing the .Rdata files itself, so the two stages can be run,
checked and rerun independently, and so which converted series belongs to
which original sample and which replicate is never in question.

This is the older of the two ways to fill a dataset directory. The design
experiment the project now reports is filled by
applications/abm_system_simulate.c instead, which runs the simulator over
dataset/abm_system_design.csv and writes cop_NNNN directories. The two
cannot share an output directory: applications/abm_system_fit_qvarma.c fits
every subdirectory it finds whatever wrote it, so two datasets sitting
there at once would be fitted together with nothing in the results to say
so. This script refuses to start when the directory it is about to write
into already holds cop_* directories, which is the same refusal
applications/abm_system_simulate_all.sh makes in the other direction.

ABM_SYSTEM_EXTRACT_DIR overrides where the archives go. It defaults to
dataset/abm_system_rdata, which is where this route's output was moved when
the design runs took dataset/abm_system over, so a plain run now rebuilds
that copy in place rather than landing on top of the design experiment.

Naming, so there is no doubt which simulation a given archive came from:
dataset/simulated/EstimationSeriesSample1_<N>.Rdata's own <N> - the
original file's own base name, unmodified - becomes the output
subdirectory name, and each of that file's own replicates (0-indexed,
matching the array position abm_system_slice reads them at) goes into
batch_<B>.npz inside it, ABM_SYSTEM_BATCH replicates to an archive, B
zero-padded to 3 digits:

    dataset/abm_system_rdata/EstimationSeriesSample1_1/batch_000.npz
    dataset/abm_system_rdata/EstimationSeriesSample1_1/batch_001.npz
    ...
    dataset/abm_system_rdata/EstimationSeriesSample1_100/batch_010.npz

Each archive carries the five named series (GDP_growth, EN_growth,
Employment_change, Inflation, InterestRate), one row per period, plus the
replicate index each row belongs to, so a reader takes one replicate out of
it through abm_system_read_replicate without needing to know how the
archives were filled. tests/abm_system_layout.c checks that round trip.

out/abm_system_extract_manifest.txt records, per source file: its own
sample number, how many replicates it held, and how many periods each
replicate has after the burn-in and differencing - the actual counts, not
an assumption, since dataset/simulated/'s own 100 files were not all
checked to hold exactly 108 replicates before this ran.

Not a final implementation: no parallelism, no resumability (a rerun
overwrites everything it already wrote), no handling for a file whose
"estimation" object does not match the expected shape beyond abm_system.h's
own asserts. The output directory is regenerated from dataset/simulated/
every time this runs, never edited by hand.

In EXPERIMENT_STEMS, so `make app-abm_system_extract` runs it. Not part of
`make applications`: it writes to dataset/, not out/, and is meant to be
run explicitly, once, before the fitting stage, not on every routine build.
Nothing printed.
*/

#include "abm_system.h"

#include <et_al./frame/npz.h>
#include <et_al./frame/csv.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#define SIMULATED_DIR "dataset/simulated"
#define OUTPUT_DIR_DEFAULT "dataset/abm_system_rdata"

static const char *OUTPUT_DIR;

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_extract: mkdir failed");
}

/* Refuses a directory the design experiment already wrote into. The fitting
   stage reads every subdirectory it finds, so a cop_* tree and this script's
   EstimationSeriesSample1_* tree sharing one directory would be fitted together
   and nothing downstream would say so. Returns 0 and leaves the reason on
   stderr, since there is no output file to write it to yet. */
static int output_dir_is_free(const char *dir) {
    DIR *handle = opendir(dir);
    if (!handle) return 1;

    int occupied = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL)
        if (strncmp(entry->d_name, "cop_", 4) == 0) { occupied = 1; break; }
    closedir(handle);

    if (occupied) {
        fprintf(stderr,
                "%s already holds cop_* directories, which the parameter design wrote.\n"
                "applications/abm_system_fit_qvarma.c fits every subdirectory it finds and\n"
                "cannot tell two datasets apart, so this script will not add to that one.\n\n"
                "Write somewhere else instead:\n\n"
                "  ABM_SYSTEM_EXTRACT_DIR=dataset/abm_system_rdata ./bin/abm_system_extract\n",
                dir);
    }
    return !occupied;
}

/* Every ".Rdata" file's own base name (without the extension), in
   whatever order readdir returns them - the manifest records what was
   actually found, so the order here does not matter. Caller must free
   each entry and the array itself. */
static char **list_rdata_basenames(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_extract: cannot open dataset/simulated/");

    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        size_t len = strlen(entry->d_name);
        if (len < 6 || strcmp(entry->d_name + len - 6, ".Rdata") != 0) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            char **grown = realloc(names, (size_t)cap * sizeof(char*));
            assert(grown && "abm_system_extract: out of memory listing source files");
            names = grown;
        }
        names[n] = malloc(len - 6 + 1);
        assert(names[n] && "abm_system_extract: out of memory copying a source name");
        memcpy(names[n], entry->d_name, len - 6);
        names[n][len - 6] = '\0';
        n++;
    }
    closedir(handle);
    *count = n;
    return names;
}

int main(void) {
    OUTPUT_DIR = getenv("ABM_SYSTEM_EXTRACT_DIR");
    if (!OUTPUT_DIR) OUTPUT_DIR = OUTPUT_DIR_DEFAULT;

    if (!output_dir_is_free(OUTPUT_DIR)) return 1;
    make_directory(OUTPUT_DIR);

    int n_files;
    char **basenames = list_rdata_basenames(SIMULATED_DIR, &n_files);
    assert(n_files > 0 && "abm_system_extract: no .Rdata files found under dataset/simulated/");

    FILE *manifest = fopen("out/abm_system_extract_manifest.txt", "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "source file -> %s/<source>/batch_<000..>.npz\n", OUTPUT_DIR);
    fprintf(manifest, "%d source files found under %s\n\n", n_files, SIMULATED_DIR);
    fprintf(manifest, "%-32s %12s %12s\n", "sample", "replicates", "periods");

    int total_replicates = 0;
    for (int f = 0; f < n_files; f++) {
        char rdata_path[512], out_dir[512];
        snprintf(rdata_path, sizeof rdata_path, "%s/%s.Rdata", SIMULATED_DIR, basenames[f]);
        snprintf(out_dir, sizeof out_dir, "%s/%s", OUTPUT_DIR, basenames[f]);
        make_directory(out_dir);

        RData d = abm_system_read(rdata_path);
        int n_rep = abm_system_n_replicates(&d);
        int n_periods = abm_system_n_periods(&d);

        Mat block[ABM_SYSTEM_BATCH];
        int replicate[ABM_SYSTEM_BATCH];
        int held = 0;

        for (int r = 0; r < n_rep; r++) {
            block[held] = abm_system_slice(&d, r);
            replicate[held] = r;
            held++;

            if (held == ABM_SYSTEM_BATCH || r == n_rep - 1) {
                abm_system_write_batch(out_dir, block, replicate, held);
                for (int b = 0; b < held; b++) mat_free(block[b]);
                held = 0;
            }
        }
        abm_system_free(&d);

        fprintf(manifest, "%-32s %12d %12d\n", basenames[f], n_rep, n_periods);
        total_replicates += n_rep;
        free(basenames[f]);
    }
    free(basenames);

    fprintf(manifest, "\ntotal: %d source files, %d replicates written\n", n_files, total_replicates);
    fclose(manifest);
    return 0;
}
