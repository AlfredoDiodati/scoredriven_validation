/*
Builds the 500 x 1000 dataset the throughput test of
applications/abm_system_scale_fit_qvarma.c runs on: 500 folders of 1000
five-variable series each, 500,000 in total, in the same
t-QVARMA(1,1,2) layout applications/abm_system_extract.c writes and read
back by the same df_read_csv call.

The conversion itself is abm_system.h's own abm_system_slice, unchanged -
same burn-in, same transforms, same five columns. What this file adds is
the count. dataset/simulated/ holds 100 files of 108 replicates, 10,800
series, and the test needs 500,000, so the originals are padded out with
perturbed copies of themselves.

Layout, so which source a given series came from is never in question.
Each of the 100 source files gets five folders under
dataset/abm_system_scale/:

    EstimationSeriesSample1_7/          series_000..107 original
                                        series_108..999 perturbed
    EstimationSeriesSample1_7_noise1/   series_000..999 perturbed
    EstimationSeriesSample1_7_noise2/   ...
    EstimationSeriesSample1_7_noise3/
    EstimationSeriesSample1_7_noise4/

100 samples x 5 folders = 500 folders, 1000 series each. Every original
replicate is written exactly once, in its own sample's own first folder,
before any perturbed series is written. A perturbed series in position i
is built from original replicate i mod 108 of that same sample, so each
original is the parent of roughly nine perturbed series per folder.

The perturbation exists to make the 500,000 fits land in 500,000
different places, not to model anything. Each perturbed series draws, per
variable k, one level shift m_k ~ N(0, MEAN_SHIFT_SD^2) applied to all
400 periods, plus one N(0, NOISE_SD^2) draw per period:

    y_noisy[k][t] = y[k][t] + m_k + e[k][t]

Both standard deviations are absolute and the same for every variable,
which is defensible here only because the five variables sit within one
order of magnitude of each other: measured over replicates 0 and 1 of
samples 1, 50 and 100, the per-variable standard deviations run from
0.063 (InterestRate, sample 100) to 1.89 (EN_growth, sample 100), with
the three growth variables around 1.2 to 1.9 and Inflation and
InterestRate around 0.06 to 0.32. At 0.01 the perturbation is under one
per cent of a growth variable's own spread and at most sixteen per cent
of the tightest InterestRate series, so no variable is either left
untouched or swamped.

Seeding is by (sample index, folder copy, series index), with the sample
index taken from the sorted file listing rather than readdir order, so a
rerun reproduces the identical 500,000 files regardless of how the
directory happens to enumerate or how OpenMP happens to schedule.

Parallelized over source files. A worker slices its file's own 108
replicates once and releases the RData immediately, holding 1.7 MB of
converted series rather than the 4.1 MB parsed array while it writes the
5000 CSVs those series generate.

Cost, before running it: 500,000 files of about 42 KB (five %.17g fields
per row, 400 rows) is roughly 21 GB on disk and 500,000 inodes. It is a
full rewrite every time, not a resume.

out/abm_system_scale_extract_manifest.txt records, per source file, how
many replicates it held, how many periods each has after the burn-in and
differencing, and how many folders and series were written from it.

Not part of `make applications` or of the abm_system chain: it writes to
dataset/, not out/, it takes tens of minutes, and abm_system_scale_fit_qvarma
deliberately does not depend on it in the Makefile so that rerunning the
fits does not rewrite 21 GB. Run it explicitly, once. Nothing printed.
*/

#include "abm_system.h"
#include <et_al./random.h>
#include <frame/csv.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define SIMULATED_DIR "dataset/simulated"
#define OUTPUT_DIR "dataset/abm_system_scale"

#define FOLDERS_PER_SAMPLE 5
#define SERIES_PER_FOLDER 1000
#define NOISE_SD ((mreal)0.01)
#define MEAN_SHIFT_SD ((mreal)0.01)
#define NOISE_SEED ((uint64_t)20260830)

static const char *row_name[ABM_SYSTEM_K] = {
    "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate"
};

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "abm_system_scale_extract: mkdir failed");
}

static int compare_basenames(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

/* Every ".Rdata" file's own base name, sorted, so that a file's own
   position in this array - which the noise seeds are derived from - does
   not depend on readdir order. Caller must free each entry and the array. */
static char **list_rdata_basenames(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_scale_extract: cannot open dataset/simulated/");

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
    qsort(names, (size_t)n, sizeof(char*), compare_basenames);
    *count = n;
    return names;
}

/* The K x T series as the T x K matrix df_from_matrix wants, with the
   level shift and the per-period noise added when add_noise is set. The
   two are drawn in a fixed order, all K shifts first, so the stream a
   given series consumes depends only on its own seed. */
static Mat transposed_with_noise(Mat y, int add_noise, Rng *rng) {
    mreal shift[ABM_SYSTEM_K];
    for (int k = 0; k < ABM_SYSTEM_K; k++)
        shift[k] = add_noise ? MEAN_SHIFT_SD * (mreal)rng_normal(rng) : (mreal)0;

    Mat out = mat_new(y.c, y.r);
    for (int t = 0; t < y.c; t++)
        for (int k = 0; k < y.r; k++) {
            mreal noise = add_noise ? NOISE_SD * (mreal)rng_normal(rng) : (mreal)0;
            AT(out, t, k) = AT(y, k, t) + shift[k] + noise;
        }
    return out;
}

int main(void) {
    make_directory(OUTPUT_DIR);

    int n_files;
    char **basenames = list_rdata_basenames(SIMULATED_DIR, &n_files);
    assert(n_files > 0 && "abm_system_scale_extract: no .Rdata files found under dataset/simulated/");

    int *n_replicates = malloc((size_t)n_files * sizeof(int));
    int *n_periods = malloc((size_t)n_files * sizeof(int));

    #pragma omp parallel for schedule(dynamic, 1)
    for (int f = 0; f < n_files; f++) {
        char rdata_path[512];
        snprintf(rdata_path, sizeof rdata_path, "%s/%s.Rdata", SIMULATED_DIR, basenames[f]);

        RData d = abm_system_read(rdata_path);
        int n_rep = abm_system_n_replicates(&d);
        assert(n_rep > 0 && "abm_system_scale_extract: a source file holds no replicates");
        n_replicates[f] = n_rep;
        n_periods[f] = abm_system_n_periods(&d) - ABM_SYSTEM_BURN_IN;

        /* Sliced once and the parsed array released before any writing, so a
           worker carries the converted series rather than both. */
        Mat *original = malloc((size_t)n_rep * sizeof(Mat));
        for (int r = 0; r < n_rep; r++) original[r] = abm_system_slice(&d, r);
        abm_system_free(&d);

        for (int copy = 0; copy < FOLDERS_PER_SAMPLE; copy++) {
            char out_dir[560];
            if (copy == 0) snprintf(out_dir, sizeof out_dir, "%s/%s", OUTPUT_DIR, basenames[f]);
            else snprintf(out_dir, sizeof out_dir, "%s/%s_noise%d", OUTPUT_DIR, basenames[f], copy);
            make_directory(out_dir);

            for (int i = 0; i < SERIES_PER_FOLDER; i++) {
                int is_original = (copy == 0 && i < n_rep);
                Mat source = original[is_original ? i : i % n_rep];

                uint64_t stream = ((uint64_t)f * FOLDERS_PER_SAMPLE + (uint64_t)copy)
                                  * SERIES_PER_FOLDER + (uint64_t)i;
                Rng rng = rng_new(NOISE_SEED, stream);

                Mat transposed = transposed_with_noise(source, !is_original, &rng);
                DataFrame df = df_from_matrix(transposed, row_name);
                mat_free(transposed);

                char csv_path[600];
                snprintf(csv_path, sizeof csv_path, "%s/series_%03d.csv", out_dir, i);
                df_write_csv(&df, csv_path, csv_write_options_default());
                df_free(&df);
            }
        }

        for (int r = 0; r < n_rep; r++) mat_free(original[r]);
        free(original);
    }

    FILE *manifest = fopen("out/abm_system_scale_extract_manifest.txt", "w");
    assert(manifest && "cannot open the manifest path for writing");
    fprintf(manifest, "source file -> dataset/abm_system_scale/<source>[_noise1..%d]/series_<000..%d>.csv\n",
            FOLDERS_PER_SAMPLE - 1, SERIES_PER_FOLDER - 1);
    fprintf(manifest, "%d source files found under %s\n\n", n_files, SIMULATED_DIR);
    fprintf(manifest, "series_i is an original replicate when it is in the source's own first\n"
                      "folder and i is below that source's replicate count; every other series is\n"
                      "original replicate i mod (replicates) with a N(0, %g^2) level shift per\n"
                      "variable and N(0, %g^2) added per variable per period, seed %llu, stream\n"
                      "(sample index, folder, series).\n\n",
            (double)MEAN_SHIFT_SD, (double)NOISE_SD, (unsigned long long)NOISE_SEED);
    fprintf(manifest, "%-32s %12s %10s %10s %10s %10s\n", "sample", "replicates", "periods",
            "folders", "original", "perturbed");

    long total_original = 0, total_perturbed = 0;
    for (int f = 0; f < n_files; f++) {
        long original_here = n_replicates[f] < SERIES_PER_FOLDER ? n_replicates[f] : SERIES_PER_FOLDER;
        long perturbed_here = (long)FOLDERS_PER_SAMPLE * SERIES_PER_FOLDER - original_here;
        fprintf(manifest, "%-32s %12d %10d %10d %10ld %10ld\n", basenames[f], n_replicates[f],
                n_periods[f], FOLDERS_PER_SAMPLE, original_here, perturbed_here);
        total_original += original_here;
        total_perturbed += perturbed_here;
        free(basenames[f]);
    }
    free(basenames);
    free(n_replicates);
    free(n_periods);

    fprintf(manifest, "\ntotal: %d folders, %ld series (%ld original, %ld perturbed)\n",
            n_files * FOLDERS_PER_SAMPLE, total_original + total_perturbed,
            total_original, total_perturbed);
    fclose(manifest);
    return 0;
}
