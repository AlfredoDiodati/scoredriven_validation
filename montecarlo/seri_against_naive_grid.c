/*
The grid of montecarlo/seri_true_configuration_recovery.c against the
discarded naive-loss grid of old/naive_mae_loss/montecarlo/sweep_grid.c,
benchmark by benchmark. Both hold one row per model and benchmark with the
same columns, and both hold out the benchmark's run index from every
configuration.

- Seri: squared error minus each run's own deviation from its configuration's
  mean response, scaled by n/(n-1), then et_al's seri_mcs at alpha 0.05.
- Naive: mean absolute error between response vectors, no correction, then
  et_al's mcs, MCS_TR, alpha 0.05, block length 1, bootstrap variance, 2000
  resamples, seed 123 stream 0.

Writes out/montecarlo/seri_against_naive_grid_report.txt:

1. Per model and grid: the benchmark's configuration in the set, alone in it,
   first by mean loss, its mean rank by mean loss, and the mean and largest
   set size.
2. Per model: whether the benchmark's configuration is in each grid's set, as
   a two by two table over the benchmarks both grids hold, and how often the
   two losses pick the same configuration as the one with the smallest mean
   loss.
3. Per model and configuration: the share of its own benchmarks whose set
   holds it, summarised per grid as how many configurations reach at least
   half, at least 90% and none, and the correlation across configurations of
   the two grids' shares.
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#define SERI_PATH "out/montecarlo/seri_true_configuration_recovery.csv.gz"
#define NAIVE_PATH "out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz"
#define REPORT_PATH "out/montecarlo/seri_against_naive_grid_report.txt"
#define MAX_COP 1024
#define MAX_RUN 1000
#define N_GRIDS 2

enum { GRID_SERI, GRID_NAIVE };
static const char *const grid_label[N_GRIDS] = { "Seri", "naive" };

/* One benchmark's row, absent when in_set is -1. */
typedef struct {
    signed char in_set;
    short rank, lowest;
    int set_size;
} Cell;

static Cell *cells[N_GRIDS][N_MODELS];
static long n_unparsed[N_GRIDS], n_duplicated[N_GRIDS];

static size_t cell_index(int cop, int run) {
    return (size_t)cop * MAX_RUN + run;
}

static void add_row(int grid, const char *line) {
    char label[32];
    int cop, run, in_set, rank, set_size, decided, dropped, lowest;
    double pvalue, final_pvalue;
    if (sscanf(line, "%31[^,],%d,%d,%d,%d,%lf,%d,%d,%lf,%d,%d", label, &cop, &run, &in_set, &rank, &pvalue,
               &set_size, &decided, &final_pvalue, &dropped, &lowest) != 11
        || cop < 0 || cop >= MAX_COP || run < 0 || run >= MAX_RUN || lowest < 0 || lowest >= MAX_COP) {
        n_unparsed[grid]++;
        return;
    }
    int model = -1;
    for (int m = 0; m < N_MODELS; m++)
        if (strcmp(label, model_label[m]) == 0) model = m;
    if (model < 0) {
        n_unparsed[grid]++;
        return;
    }
    Cell *cell = &cells[grid][model][cell_index(cop, run)];
    if (cell->in_set >= 0) n_duplicated[grid]++;
    cell->in_set = (signed char)in_set;
    cell->rank = (short)rank;
    cell->lowest = (short)lowest;
    cell->set_size = set_size;
}

static void read_grid(int grid, const char *path) {
    FILE *f = fopen(path, "rb");
    assert(f && "seri_against_naive_grid: cannot open a result file");
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc((size_t)length);
    assert(compressed);
    size_t n_read = fread(compressed, 1, (size_t)length, f);
    assert(n_read == (size_t)length);
    (void)n_read;
    fclose(f);
    size_t text_length;
    unsigned char *text = gzip_inflate(compressed, (size_t)length, &text_length);
    free(compressed);

    FILE *rows = fmemopen(text, text_length, "r");
    assert(rows);
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, rows) > 0) {
        if (header) { header = 0; continue; }
        add_row(grid, line);
    }
    free(line);
    fclose(rows);
    free(text);
}

static void write_per_grid(FILE *out) {
    fprintf(out, "1. Per model and grid, over every benchmark the grid holds. Percentages are of those "
                 "benchmarks.\n\n");
    fprintf(out, "model, grid, benchmarks, own configuration in the set, alone in the set, first by mean loss, "
                 "mean rank, mean set size, median set size, largest set\n");
    for (int model = 0; model < N_MODELS; model++)
        for (int grid = 0; grid < N_GRIDS; grid++) {
            long n = 0, in_set = 0, alone = 0, first = 0, largest = 0;
            double rank_sum = 0, size_sum = 0;
            Mat sizes = mat_new(MAX_COP * MAX_RUN, 1);
            for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) {
                const Cell *c = &cells[grid][model][at];
                if (c->in_set < 0) continue;
                n++;
                in_set += c->in_set;
                alone += c->in_set && c->set_size == 1;
                first += c->rank == 1;
                rank_sum += c->rank;
                size_sum += c->set_size;
                AT(sizes, n - 1, 0) = c->set_size;
                if (c->set_size > largest) largest = c->set_size;
            }
            if (!n) {
                mat_free(sizes);
                continue;
            }
            double median_size = stats_median(mat_slice(sizes, 0, (int)n, 0, 1));
            mat_free(sizes);
            fprintf(out, "%s, %s, %ld, %.1f%%, %.1f%%, %.1f%%, %.1f, %.2f, %g, %ld\n", model_label[model],
                    grid_label[grid], n, 100.0 * in_set / n, 100.0 * alone / n, 100.0 * first / n, rank_sum / n,
                    size_sum / n, median_size, largest);
        }
    fprintf(out, "\n");
}

static void write_paired(FILE *out) {
    fprintf(out, "2. Per model, over the benchmarks both grids hold: whether the benchmark's configuration is "
                 "in each set, and whether both losses give the same configuration the smallest mean loss.\n\n");
    fprintf(out, "model, benchmarks in both, in both sets, Seri set only, naive set only, neither, "
                 "same smallest-loss configuration\n");
    for (int model = 0; model < N_MODELS; model++) {
        long n = 0, both = 0, seri_only = 0, naive_only = 0, neither = 0, same_lowest = 0;
        for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) {
            const Cell *s = &cells[GRID_SERI][model][at], *v = &cells[GRID_NAIVE][model][at];
            if (s->in_set < 0 || v->in_set < 0) continue;
            n++;
            both += s->in_set && v->in_set;
            seri_only += s->in_set && !v->in_set;
            naive_only += !s->in_set && v->in_set;
            neither += !s->in_set && !v->in_set;
            same_lowest += s->lowest == v->lowest;
        }
        if (!n) continue;
        fprintf(out, "%s, %ld, %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%\n", model_label[model], n, 100.0 * both / n,
                100.0 * seri_only / n, 100.0 * naive_only / n, 100.0 * neither / n, 100.0 * same_lowest / n);
    }
    fprintf(out, "\n");
}

static void write_per_configuration(FILE *out) {
    fprintf(out, "3. Per model and configuration: the share of its own benchmarks whose set holds it, over the "
                 "configurations both grids hold, with the Pearson correlation of the two grids' shares across "
                 "configurations.\n\n");
    fprintf(out, "model, configurations, Seri at least half, Seri at least 90%%, Seri none, naive at least half, "
                 "naive at least 90%%, naive none, correlation\n");
    for (int model = 0; model < N_MODELS; model++) {
        int n_configurations = 0, half[N_GRIDS] = { 0, 0 }, ninety[N_GRIDS] = { 0, 0 }, none[N_GRIDS] = { 0, 0 };
        double sum[N_GRIDS] = { 0, 0 }, square[N_GRIDS] = { 0, 0 }, cross = 0;
        for (int cop = 0; cop < MAX_COP; cop++) {
            long n[N_GRIDS] = { 0, 0 }, held[N_GRIDS] = { 0, 0 };
            for (int run = 0; run < MAX_RUN; run++)
                for (int grid = 0; grid < N_GRIDS; grid++) {
                    const Cell *c = &cells[grid][model][cell_index(cop, run)];
                    if (c->in_set < 0) continue;
                    n[grid]++;
                    held[grid] += c->in_set;
                }
            if (!n[GRID_SERI] || !n[GRID_NAIVE]) continue;
            n_configurations++;
            double share[N_GRIDS];
            for (int grid = 0; grid < N_GRIDS; grid++) {
                share[grid] = (double)held[grid] / n[grid];
                half[grid] += share[grid] >= 0.5;
                ninety[grid] += share[grid] >= 0.9;
                none[grid] += held[grid] == 0;
                sum[grid] += share[grid];
                square[grid] += share[grid] * share[grid];
            }
            cross += share[GRID_SERI] * share[GRID_NAIVE];
        }
        if (!n_configurations) continue;
        double k = n_configurations;
        double covariance = cross / k - sum[GRID_SERI] / k * sum[GRID_NAIVE] / k;
        double variance_seri = square[GRID_SERI] / k - sum[GRID_SERI] / k * sum[GRID_SERI] / k;
        double variance_naive = square[GRID_NAIVE] / k - sum[GRID_NAIVE] / k * sum[GRID_NAIVE] / k;
        fprintf(out, "%s, %d, %d, %d, %d, %d, %d, %d, %.3f\n", model_label[model], n_configurations,
                half[GRID_SERI], ninety[GRID_SERI], none[GRID_SERI], half[GRID_NAIVE], ninety[GRID_NAIVE],
                none[GRID_NAIVE], covariance / sqrt(variance_seri * variance_naive));
    }
}

int main(void) {
    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) {
            cells[grid][model] = malloc((size_t)MAX_COP * MAX_RUN * sizeof(Cell));
            assert(cells[grid][model]);
            for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) cells[grid][model][at].in_set = -1;
        }
    read_grid(GRID_SERI, SERI_PATH);
    read_grid(GRID_NAIVE, NAIVE_PATH);

    mkdir("out", 0755);
    mkdir("out/montecarlo", 0755);
    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "seri_against_naive_grid: cannot open the report");
    fprintf(out, "%s against %s.\n", SERI_PATH, NAIVE_PATH);
    for (int grid = 0; grid < N_GRIDS; grid++)
        fprintf(out, "%s: %ld rows that did not parse, %ld benchmarks found twice.\n", grid_label[grid],
                n_unparsed[grid], n_duplicated[grid]);
    fprintf(out, "\n");
    write_per_grid(out);
    write_paired(out);
    write_per_configuration(out);
    fclose(out);

    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) free(cells[grid][model]);
    return 0;
}
