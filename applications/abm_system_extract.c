/*
Converts every replicate of every file under dataset/simulated/ into the
t-QVARMAd(1,1,2)/(1,1,4) layout abm_system.h's own header comment
describes, and writes each one to its own CSV - a data-preparation step,
not a fit: nothing here calls qvarma_d.h's own fit(). The next pipeline
stage (indirect-inference refitting across every simulated dataset) reads
these back rather than re-parsing the .Rdata files itself, so the two
stages can be run, checked and rerun independently, and so which converted
series belongs to which original sample and which replicate is never in
question.

Naming, so there is no doubt which simulation a given CSV came from:
dataset/simulated/EstimationSeriesSample1_<N>.Rdata's own <N> - the
original file's own base name, unmodified - becomes the output
subdirectory name under dataset/abm_system/, and each of that file's own
replicates (0-indexed, matching the array position abm_system_slice reads
them at) becomes replicate_<R>.csv inside it, R zero-padded to 3 digits
since 108 replicates need three:

    dataset/abm_system/EstimationSeriesSample1_1/replicate_000.csv
    dataset/abm_system/EstimationSeriesSample1_1/replicate_001.csv
    ...
    dataset/abm_system/EstimationSeriesSample1_100/replicate_107.csv

Each CSV has the same five named columns
applications/us_qvarmad_employment_change.c's own residual CSVs use
(GDP_growth, EN_growth, Employment_change, Inflation, InterestRate), one
row per period, so the next stage can read it with the same df_read_csv
call it would use for real data.

out/abm_system_extract_manifest.txt records, per source file: its own
sample number, how many replicates it held, and how many periods each
replicate has after the burn-in and differencing - the actual counts, not
an assumption, since dataset/simulated/'s own 100 files were not all
checked to hold exactly 108 replicates before this ran.

Not a final implementation: no parallelism, no resumability (a rerun
overwrites everything it already wrote), no handling for a file whose
"estimation" object does not match the expected shape beyond abm_system.h's
own asserts. dataset/abm_system/ is regenerated from dataset/simulated/
every time this runs, never edited by hand.

Not part of `make applications` or EXPERIMENT_STEMS - it writes to
dataset/, not out/, and is meant to be run explicitly, once, before the
fitting stage, not on every routine build. Output aside from the manifest:
dataset/abm_system/ itself. Nothing printed.
*/

#include "abm_system.h"
#include <frame/csv.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define SIMULATED_DIR "dataset/simulated"
#define OUTPUT_DIR "dataset/abm_system"

static const char *row_name[ABM_SYSTEM_K] = {
    "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
};

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_extract: mkdir failed");
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
        if (n == cap) { cap = cap ? cap * 2 : 16; names = realloc(names, (size_t)cap * sizeof(char*)); }
        names[n] = malloc(len - 6 + 1);
        memcpy(names[n], entry->d_name, len - 6);
        names[n][len - 6] = '\0';
        n++;
    }
    closedir(handle);
    *count = n;
    return names;
}

int main(void) {
    make_directory(OUTPUT_DIR);

    int n_files;
    char **basenames = list_rdata_basenames(SIMULATED_DIR, &n_files);
    assert(n_files > 0 && "abm_system_extract: no .Rdata files found under dataset/simulated/");

    FILE *manifest = fopen("out/abm_system_extract_manifest.txt", "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "source file -> dataset/abm_system/<source>/replicate_<000..>.csv\n");
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
        int n_periods = abm_system_n_periods(&d) - ABM_SYSTEM_BURN_IN;

        for (int r = 0; r < n_rep; r++) {
            Mat y = abm_system_slice(&d, r);

            Mat transposed = mat_new(y.c, y.r);
            for (int t = 0; t < y.c; t++)
                for (int a = 0; a < y.r; a++) AT(transposed, t, a) = AT(y, a, t);
            DataFrame df = df_from_matrix(transposed, row_name);

            char csv_path[560];
            snprintf(csv_path, sizeof csv_path, "%s/replicate_%03d.csv", out_dir, r);
            df_write_csv(&df, csv_path, csv_write_options_default());

            df_free(&df);
            mat_free(transposed);
            mat_free(y);
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
