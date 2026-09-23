/*
Whether every value stored in dataset/abm_system is a finite number, and
whether the archives hold the runs the experiment claims to have stored.

The 1000 x 1000 experiment writes 5 x 400 numbers per (configuration,
replication) pair, grouped ABM_SYSTEM_BATCH to an archive. The only finiteness
check it makes is at write time: applications/abm_system_simulate.c's own
all_finite refuses a replication whose model output or transformed series
carries a NaN or an infinity, logs the seed in out/abm_system_simulate/ and
stores nothing for it. Everything downstream reads the archives back and
carries its own guards against a non-finite loss, which is a different
quantity. Nothing reads the stored series themselves and asks whether they are
numbers, so a value damaged after the check - by the transform, by the npz
round trip, by a truncated file from an interrupted run - would first be seen
as a fit that behaves strangely. This is that pass.

Four questions, all answered from the archives alone. Nothing here runs the
simulator or fits anything.

1. Is every stored value finite? NaN and infinity are counted separately, and
   the first offending (replication, series, period) of a configuration is
   named, because a single damaged period and a whole damaged series are
   different accidents.

2. Does every block have the shape the layout promises, ABM_SYSTEM_K rows by
   the period count the first archive read establishes?

3. Does every configuration hold every replication once? A replication index
   appearing twice, appearing outside 0 .. N_MC-1, or not appearing at all is
   reported per configuration.

4. What range does each series actually cover? Reported as context rather than
   as a verdict: the check has no view on how large a growth rate may be, and a
   bound invented here would fail runs the experiment meant to keep.

    ./bin/abm_system_dataset_finiteness [FIRST_COP [LAST_COP [N_MC]]]

Defaults are the whole design and 1000 replications per configuration.
Configurations are independent, so the sweep is divided over the threads
OpenMP gives it; the work is reading and decompressing 15 GB of archives.

Writes out/abm_system_dataset_finiteness_report.txt. Needs dataset/abm_system.
*/

#include "applications/abm_system.h"

#include <et_al./linalg/mat.h>
#include <et_al./frame/csv.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <assert.h>
#include <omp.h>

#define DATASET "dataset/abm_system"
#define DESIGN "dataset/abm_system_design.csv"
#define REPORT "out/abm_system_dataset_finiteness_report.txt"
#define DEFAULT_N_MC 1000

typedef struct {
    int present;
    int archives_missing;
    int blocks;
    int duplicated;
    int out_of_range;
    int wrong_shape;
    long nan_values;
    long infinite_values;
    int first_bad_replicate, first_bad_series, first_bad_period;
    double smallest[ABM_SYSTEM_K], largest[ABM_SYSTEM_K];
} Configuration;

static int directory_exists(const char *path) {
    struct stat info;
    return stat(path, &info) == 0 && S_ISDIR(info.st_mode);
}

static void configuration_path(char *out, size_t n, int cop) {
    int written = snprintf(out, n, "%s/cop_%04d", DATASET, cop);
    assert(written > 0 && (size_t)written < n && "abm_system_dataset_finiteness: path does not fit");
}

/* The period count every block is required to have, taken from the first
   archive that can be read rather than recomputed from the burn-in, so the
   check compares the archives against themselves and a design run at another
   horizon still passes on its own terms. Returns 0 if nothing could be read. */
static int established_periods(int first_cop, int last_cop) {
    Mat block[ABM_SYSTEM_BATCH];
    int replicate[ABM_SYSTEM_BATCH];
    for (int cop = first_cop; cop <= last_cop; cop++) {
        char dir[512];
        configuration_path(dir, sizeof dir, cop);
        if (!directory_exists(dir)) continue;
        int count = abm_system_read_batch(dir, 0, block, replicate);
        if (count > 0) {
            int periods = block[0].c;
            for (int b = 0; b < count; b++) mat_free(block[b]);
            return periods;
        }
    }
    return 0;
}

static void scan_configuration(int cop, int n_mc, int periods, Configuration *out, char *seen) {
    char dir[512];
    configuration_path(dir, sizeof dir, cop);

    out->first_bad_replicate = -1;
    for (int k = 0; k < ABM_SYSTEM_K; k++) {
        out->smallest[k] = INFINITY;
        out->largest[k] = -INFINITY;
    }
    if (!directory_exists(dir)) return;
    out->present = 1;

    memset(seen, 0, (size_t)n_mc);
    int n_batches = (n_mc + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;

    Mat block[ABM_SYSTEM_BATCH];
    int replicate[ABM_SYSTEM_BATCH];
    for (int batch = 0; batch < n_batches; batch++) {
        int count = abm_system_read_batch(dir, batch, block, replicate);
        if (count == 0) { out->archives_missing++; continue; }

        for (int b = 0; b < count; b++) {
            out->blocks++;
            if (replicate[b] < 0 || replicate[b] >= n_mc) out->out_of_range++;
            else if (seen[replicate[b]]) out->duplicated++;
            else seen[replicate[b]] = 1;

            if (block[b].r != ABM_SYSTEM_K || block[b].c != periods) out->wrong_shape++;

            for (int series = 0; series < block[b].r; series++) {
                for (int t = 0; t < block[b].c; t++) {
                    mreal value = AT(block[b], series, t);
                    if (MISNAN(value) || MISINF(value)) {
                        if (MISNAN(value)) out->nan_values++; else out->infinite_values++;
                        if (out->first_bad_replicate < 0) {
                            out->first_bad_replicate = replicate[b];
                            out->first_bad_series = series;
                            out->first_bad_period = t;
                        }
                        continue;
                    }
                    if ((double)value < out->smallest[series]) out->smallest[series] = (double)value;
                    if ((double)value > out->largest[series]) out->largest[series] = (double)value;
                }
            }
            mat_free(block[b]);
        }
    }
}

int main(int argc, char **argv) {
    int first_cop = argc > 1 ? atoi(argv[1]) : 1;
    int last_cop = argc > 2 ? atoi(argv[2]) : 0;
    int n_mc = argc > 3 ? atoi(argv[3]) : DEFAULT_N_MC;

    if (last_cop == 0) {
        DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
        last_cop = design.r;
        df_free(&design);
    }
    assert(first_cop >= 1 && first_cop <= last_cop && "abm_system_dataset_finiteness: empty range");
    assert(n_mc >= 1 && "abm_system_dataset_finiteness: at least one replication per configuration");

    if (!directory_exists(DATASET)) {
        fprintf(stderr, "abm_system_dataset_finiteness: %s is not there - the experiment has not been run "
                        "on this machine\n", DATASET);
        return EXIT_FAILURE;
    }

    int periods = established_periods(first_cop, last_cop);
    if (periods == 0) {
        fprintf(stderr, "abm_system_dataset_finiteness: no archive in cop_%04d .. cop_%04d could be read\n",
                first_cop, last_cop);
        return EXIT_FAILURE;
    }

    int n_cop = last_cop - first_cop + 1;
    Configuration *result = calloc((size_t)n_cop, sizeof(Configuration));
    assert(result && "abm_system_dataset_finiteness: out of memory");

    double started = omp_get_wtime();
    #pragma omp parallel
    {
        char *seen = malloc((size_t)n_mc);
        assert(seen && "abm_system_dataset_finiteness: out of memory");
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n_cop; i++) scan_configuration(first_cop + i, n_mc, periods, &result[i], seen);
        free(seen);
    }
    double elapsed = omp_get_wtime() - started;

    FILE *report = fopen(REPORT, "w");
    assert(report && "abm_system_dataset_finiteness: cannot open the report");
    fprintf(report, "Finiteness and completeness of %s\n\n", DATASET);
    fprintf(report, "configurations cop_%04d .. cop_%04d, %d replications each, %d series x %d periods\n\n",
            first_cop, last_cop, n_mc, ABM_SYSTEM_K, periods);

    long nan_values = 0, infinite_values = 0, blocks = 0;
    int missing_configurations = 0, archives_missing = 0, duplicated = 0, out_of_range = 0;
    int wrong_shape = 0, incomplete = 0, damaged = 0;
    double smallest[ABM_SYSTEM_K], largest[ABM_SYSTEM_K];
    for (int k = 0; k < ABM_SYSTEM_K; k++) { smallest[k] = INFINITY; largest[k] = -INFINITY; }

    for (int i = 0; i < n_cop; i++) {
        const Configuration *c = &result[i];
        int cop = first_cop + i;
        int stored = c->blocks - c->duplicated - c->out_of_range;

        if (!c->present) {
            missing_configurations++;
            fprintf(report, "cop_%04d: no directory  FAILED\n", cop);
            continue;
        }
        blocks += c->blocks;
        archives_missing += c->archives_missing;
        duplicated += c->duplicated;
        out_of_range += c->out_of_range;
        wrong_shape += c->wrong_shape;
        nan_values += c->nan_values;
        infinite_values += c->infinite_values;
        for (int k = 0; k < ABM_SYSTEM_K; k++) {
            if (c->smallest[k] < smallest[k]) smallest[k] = c->smallest[k];
            if (c->largest[k] > largest[k]) largest[k] = c->largest[k];
        }

        if (c->nan_values || c->infinite_values) {
            damaged++;
            fprintf(report, "cop_%04d: %ld NaN and %ld infinite values, first at replication %d, %s, period %d"
                            "  FAILED\n",
                    cop, c->nan_values, c->infinite_values, c->first_bad_replicate,
                    abm_system_column_names()[c->first_bad_series], c->first_bad_period);
        }
        if (c->wrong_shape)
            fprintf(report, "cop_%04d: %d blocks are not %d x %d  FAILED\n",
                    cop, c->wrong_shape, ABM_SYSTEM_K, periods);
        if (c->duplicated || c->out_of_range)
            fprintf(report, "cop_%04d: %d replication indices repeated, %d outside 0 .. %d  FAILED\n",
                    cop, c->duplicated, c->out_of_range, n_mc - 1);
        if (stored != n_mc) {
            incomplete++;
            fprintf(report, "cop_%04d: %d replications stored of %d, %d archives absent\n",
                    cop, stored, n_mc, c->archives_missing);
        }
    }

    fprintf(report, "\nreplications read               %ld of %d\n", blocks, n_cop * n_mc);
    fprintf(report, "values read                     %ld\n", blocks * ABM_SYSTEM_K * periods);
    fprintf(report, "NaN values                      %ld\n", nan_values);
    fprintf(report, "infinite values                 %ld\n", infinite_values);
    fprintf(report, "configurations with either      %d\n", damaged);
    fprintf(report, "configurations absent           %d\n", missing_configurations);
    fprintf(report, "configurations short            %d\n", incomplete);
    fprintf(report, "archives absent                 %d\n", archives_missing);
    fprintf(report, "blocks of the wrong shape       %d\n", wrong_shape);
    fprintf(report, "indices repeated                %d\n", duplicated);
    fprintf(report, "indices out of range            %d\n", out_of_range);

    fprintf(report, "\nrange of each series over everything read, context rather than a verdict\n");
    for (int k = 0; k < ABM_SYSTEM_K; k++)
        fprintf(report, "  %-18s %14.6g  %14.6g\n", abm_system_column_names()[k], smallest[k], largest[k]);

    int failed = nan_values || infinite_values || missing_configurations || wrong_shape ||
                 duplicated || out_of_range || blocks == 0;
    fprintf(report, "\n%s\n", failed ? "FAILED" : "PASSED");
    fclose(report);

    printf("dataset finiteness: %ld replications, %ld values, %ld NaN, %ld infinite, "
           "%d configurations short, %.1f minutes\n",
           blocks, blocks * ABM_SYSTEM_K * periods, nan_values, infinite_values, incomplete, elapsed / 60.0);
    printf("%s\n", failed ? "FAILED" : "PASSED, 0 failures");

    free(result);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
