/*
The local projections fitted with inflation and the interest rate in first
differences against the same projections on the stored series, under the two
Seri grids: unadjusted (squared error, every run a benchmark) and adjusted
(squared error minus the own deviation from the mean of runs 900 to 999,
benchmarks runs 0 to 899), seri_mcs at alpha 0.05 in both.

Pairs: the linear projection (lp_lin against lp_lin_differenced) and both
states under the causal HP cycle at the estimated lambda (lp_nl_... against
lp_nl_differenced_...). Per grid and model, over every benchmark the grid
holds: the share whose set holds the benchmark's configuration, the share
where it has the smallest mean loss, the median and mean set size. Per pair,
over the benchmarks both models hold: in both sets, only in the stored one's,
only in the differenced one's, and how often the two give the same
configuration the smallest mean loss.

The stored projections' rows are read from the grid's result in
out/montecarlo/, the differenced ones' from its copy's in
out/old/montecarlo_differenced_lp/. Each with its progress file, a benchmark
counted once, so the report can be written while a grid is still running.
Writes
out/old/montecarlo_differenced_lp/differenced_lp_against_stored_seri_report.txt.
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#define OUTPUT_DIR "out/old/montecarlo_differenced_lp"
#define REPORT_PATH OUTPUT_DIR "/differenced_lp_against_stored_seri_report.txt"
#define MAX_COP 1024
#define MAX_RUN 1000

enum { SERI_UNADJUSTED, SERI_ADJUSTED, N_GRIDS };
static const char *const grid_label[N_GRIDS] = {
    "Seri et al., unadjusted (squared error, runs 0 to 999 as benchmarks, 999 loss rows)",
    "Seri et al., adjusted (own deviation from the mean of runs 900 to 999 subtracted, runs 0 to 899 as benchmarks)",
};
static const char *const grid_name[N_GRIDS] = {
    "seri_uncorrected_true_configuration_recovery",
    "seri_held_out_mean_true_configuration_recovery_100_mean_runs",
};
/* Where the stored projections' rows are and where the differenced ones'. */
#define N_SOURCES 2
static const char *const source_dir[N_SOURCES] = { "out/montecarlo", OUTPUT_DIR };

#define N_PAIRS 2
static const int stored_model[N_PAIRS] = { MODEL_LP_LIN, MODEL_LP_NL_CAUSAL_HP_ESTIMATED_LAMBDA };
static const int differenced_model[N_PAIRS] = { MODEL_LP_LIN_DIFFERENCED,
                                                MODEL_LP_NL_DIFFERENCED_CAUSAL_HP_ESTIMATED_LAMBDA };
static const char *const pair_label[N_PAIRS] = { "LP (linear)",
                                                 "NLP (both states, causal HP cycle, estimated lambda)" };
#define N_SLOTS (2 * N_PAIRS)

/* One benchmark's outcome: set_size 0 means not seen. */
typedef struct {
    unsigned char in_set, lowest_is_true;
    short lowest;
    int set_size;
} Cell;

static Cell *cells[N_GRIDS][N_SLOTS];
static long unparsed[N_GRIDS];

static int slot_of(const char *label) {
    for (int p = 0; p < N_PAIRS; p++) {
        if (strcmp(label, model_label[stored_model[p]]) == 0) return 2 * p;
        if (strcmp(label, model_label[differenced_model[p]]) == 0) return 2 * p + 1;
    }
    return -1;
}

static void add_row(int grid, const char *line) {
    char label[MODEL_LABEL_CAPACITY];
    int cop, run, in_set, rank, set_size, converged, dropped, lowest;
    double pvalue, final_pvalue;
    if (sscanf(line, "%63[^,],%d,%d,%d,%d,%lf,%d,%d,%lf,%d,%d,", label, &cop, &run, &in_set, &rank, &pvalue,
               &set_size, &converged, &final_pvalue, &dropped, &lowest) != 11) {
        unparsed[grid]++;
        return;
    }
    int slot = slot_of(label);
    if (slot < 0) return;
    if (cop < 0 || cop >= MAX_COP || run < 0 || run >= MAX_RUN || set_size < 1) {
        unparsed[grid]++;
        return;
    }
    Cell *c = &cells[grid][slot][(size_t)cop * MAX_RUN + run];
    c->in_set = (unsigned char)(in_set == 1);
    c->lowest_is_true = (unsigned char)(rank == 1);
    c->lowest = (short)lowest;
    c->set_size = set_size;
}

static void read_stream(int grid, FILE *rows) {
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, rows) > 0) {
        if (header) { header = 0; continue; }
        /* A progress file cut by an unclean shutdown can hold NUL bytes. */
        if (line[0] == '\0') continue;
        add_row(grid, line);
    }
    free(line);
}

static void read_source(int grid, const char *dir, FILE *report) {
    char result_path[512], progress_path[512];
    snprintf(result_path, sizeof result_path, "%s/%s.csv.gz", dir, grid_name[grid]);
    snprintf(progress_path, sizeof progress_path, "%s/%s_progress.csv", dir, grid_name[grid]);
    struct stat st;
    if (stat(result_path, &st) == 0) {
        FILE *f = fopen(result_path, "rb");
        assert(f && "differenced_lp_against_stored_seri_report: cannot open a result");
        unsigned char *compressed = malloc((size_t)st.st_size);
        assert(compressed);
        size_t n_read = fread(compressed, 1, (size_t)st.st_size, f);
        assert(n_read == (size_t)st.st_size);
        (void)n_read;
        fclose(f);
        size_t text_length;
        unsigned char *text = gzip_inflate(compressed, (size_t)st.st_size, &text_length);
        assert(text && "differenced_lp_against_stored_seri_report: a result did not inflate");
        free(compressed);
        FILE *rows = fmemopen(text, text_length, "r");
        assert(rows);
        read_stream(grid, rows);
        fclose(rows);
        free(text);
        fprintf(report, "  read %s\n", result_path);
    }
    if (stat(progress_path, &st) == 0) {
        FILE *rows = fopen(progress_path, "r");
        assert(rows);
        read_stream(grid, rows);
        fclose(rows);
        fprintf(report, "  read %s\n", progress_path);
    }
}

static void read_grid(int grid, FILE *report) {
    for (int source = 0; source < N_SOURCES; source++) read_source(grid, source_dir[source], report);
    fprintf(report, "  rows that did not parse: %ld\n", unparsed[grid]);
}

static void summarise_model(FILE *report, int grid, int slot, int model) {
    const Cell *c = cells[grid][slot];
    long n = 0, in_set = 0, lowest = 0;
    double size_sum = 0;
    for (size_t i = 0; i < (size_t)MAX_COP * MAX_RUN; i++) n += c[i].set_size > 0;
    if (!n) {
        fprintf(report, "    %-46s no benchmarks yet\n", model_label[model]);
        return;
    }
    Mat sizes = mat_new((int)n, 1);
    long at = 0;
    for (size_t i = 0; i < (size_t)MAX_COP * MAX_RUN; i++) {
        if (c[i].set_size <= 0) continue;
        in_set += c[i].in_set;
        lowest += c[i].lowest_is_true;
        size_sum += c[i].set_size;
        AT(sizes, at, 0) = (mreal)c[i].set_size;
        at++;
    }
    fprintf(report, "    %-46s benchmarks %7ld, in set %6.2f%%, smallest mean loss %6.2f%%, set size median %g, "
                    "mean %.2f\n",
            model_label[model], n, 100.0 * in_set / n, 100.0 * lowest / n, (double)stats_median(sizes),
            size_sum / n);
    mat_free(sizes);
}

static void compare_pair(FILE *report, int grid, int pair) {
    const Cell *a = cells[grid][2 * pair], *b = cells[grid][2 * pair + 1];
    long n = 0, both = 0, only_stored = 0, only_differenced = 0, same_lowest = 0;
    for (size_t i = 0; i < (size_t)MAX_COP * MAX_RUN; i++) {
        if (a[i].set_size <= 0 || b[i].set_size <= 0) continue;
        n++;
        both += a[i].in_set && b[i].in_set;
        only_stored += a[i].in_set && !b[i].in_set;
        only_differenced += !a[i].in_set && b[i].in_set;
        same_lowest += a[i].lowest == b[i].lowest;
    }
    if (!n) {
        fprintf(report, "    paired: no benchmark held by both models yet\n");
        return;
    }
    fprintf(report, "    paired over %ld benchmarks both hold: in both sets %.2f%%, only stored %.2f%%, only "
                    "differenced %.2f%%, in neither %.2f%%; same smallest-loss configuration %.2f%%\n",
            n, 100.0 * both / n, 100.0 * only_stored / n, 100.0 * only_differenced / n,
            100.0 * (n - both - only_stored - only_differenced) / n, 100.0 * same_lowest / n);
}

int main(void) {
    mkdir("out/old", 0755);
    mkdir(OUTPUT_DIR, 0755);
    FILE *report = fopen(REPORT_PATH, "w");
    assert(report && "differenced_lp_against_stored_seri_report: cannot open the report");
    fprintf(report, "Local projections with inflation and the interest rate in first differences against the "
                    "stored series, seri_mcs at alpha 0.05.\n"
                    "in set: the benchmark's own configuration is in the final set; smallest mean loss: it has "
                    "rank 1 by mean loss.\n\n");
    for (int grid = 0; grid < N_GRIDS; grid++) {
        for (int slot = 0; slot < N_SLOTS; slot++) {
            cells[grid][slot] = calloc((size_t)MAX_COP * MAX_RUN, sizeof(Cell));
            assert(cells[grid][slot]);
        }
        fprintf(report, "%s\n", grid_label[grid]);
        read_grid(grid, report);
        for (int pair = 0; pair < N_PAIRS; pair++) {
            fprintf(report, "  %s\n", pair_label[pair]);
            summarise_model(report, grid, 2 * pair, stored_model[pair]);
            summarise_model(report, grid, 2 * pair + 1, differenced_model[pair]);
            compare_pair(report, grid, pair);
        }
        fprintf(report, "\n");
        for (int slot = 0; slot < N_SLOTS; slot++) free(cells[grid][slot]);
    }
    fclose(report);
    return 0;
}
