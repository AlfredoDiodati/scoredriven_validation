/*
The grid of montecarlo/seri_true_configuration_recovery.c against the grid of
montecarlo/seri_uncorrected_true_configuration_recovery.c, benchmark by
benchmark. Both run et_al's seri_mcs at alpha 0.05 on the same benchmarks with
the same run held out; they differ only in the loss.

- Corrected: squared error between response vectors averaged over the entries,
  minus each run's own squared deviation from its configuration's mean
  response averaged over the entries, scaled by n/(n-1).
- Uncorrected: the same squared error, nothing subtracted.

Every count below is over the benchmarks both grids hold, so the report can be
run while the uncorrected grid is still going. The uncorrected rows are read
from its result, or from its progress file when the result does not exist yet.

Writes out/seri_corrected_against_uncorrected_grid_report.txt:

1. Per model and grid: the benchmark's configuration in the set, alone in it,
   first by mean loss, its mean rank by mean loss, and the mean, median and
   largest set size.
2. Per model: whether the benchmark's configuration is in each grid's set, as
   a two by two table, and how often the two losses give the same
   configuration the smallest mean loss.
3. Per model and configuration: the share of its own benchmarks whose set
   holds it, summarised per grid as how many configurations reach at least
   half, at least 90% and none, and the Pearson correlation across
   configurations of the two grids' shares.
4. Per model and grid, against each configuration's noise, its mean response
   variance over all its runs averaged over the entries, read from
   out/response_variance_by_configuration.csv: the Spearman correlation across
   configurations between noise and the own rate (share of its own sets
   holding it) and the off-target rate (share of the other configurations'
   sets holding it), and the mean noise rank, 1 the least noisy, of the
   configuration with the smallest mean loss and of the configurations in the
   set. A loss that does not depend on noise puts both mean ranks near the
   middle, (n + 1) / 2 for n configurations.
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#ifndef CORRECTED_PATH
#define CORRECTED_PATH "out/montecarlo/seri_true_configuration_recovery.csv.gz"
#endif
#ifndef UNCORRECTED_RESULT_PATH
#define UNCORRECTED_RESULT_PATH "out/montecarlo/seri_uncorrected_true_configuration_recovery.csv.gz"
#endif
#ifndef UNCORRECTED_PROGRESS_PATH
#define UNCORRECTED_PROGRESS_PATH "out/montecarlo/seri_uncorrected_true_configuration_recovery_progress.csv"
#endif
#ifndef REPORT_PATH
#define REPORT_PATH "out/seri_corrected_against_uncorrected_grid_report.txt"
#endif
#define NOISE_PATH "out/response_variance_by_configuration.csv"
#define MAX_COP 1024
#define MAX_RUN 1000
#define N_GRIDS 2

enum { GRID_CORRECTED, GRID_UNCORRECTED };
static const char *const grid_label[N_GRIDS] = { "corrected", "uncorrected" };

/* One benchmark's row, absent when in_set is -1. */
typedef struct {
    signed char in_set;
    short rank, lowest;
    int set_size;
} Cell;

static Cell *cells[N_GRIDS][N_MODELS];
static long n_unparsed[N_GRIDS], n_duplicated[N_GRIDS];

/* Each configuration's noise and its noise rank among the configurations of
   the same model, 1 the least noisy; NaN and 0 when absent. */
static double noise[N_MODELS][MAX_COP];
static int noise_rank[N_MODELS][MAX_COP];

/* Filled by the second pass over each grid, paired benchmarks only: how many
   sets of other configurations hold each configuration, and the sums of the
   noise ranks of the smallest-loss configuration and of the set members. */
static long off_target_held[N_GRIDS][N_MODELS][MAX_COP];
static double lowest_noise_rank_sum[N_GRIDS][N_MODELS], member_noise_rank_sum[N_GRIDS][N_MODELS];
static long member_count[N_GRIDS][N_MODELS];

static size_t cell_index(int cop, int run) {
    return (size_t)cop * MAX_RUN + run;
}

static int model_of(const char *label) {
    for (int model = 0; model < N_MODELS; model++)
        if (strcmp(label, model_label[model]) == 0) return model;
    return -1;
}

static int paired(int model, size_t at) {
    return cells[GRID_CORRECTED][model][at].in_set >= 0 && cells[GRID_UNCORRECTED][model][at].in_set >= 0;
}

static void store_row(int grid, const char *line) {
    char label[32];
    int cop, run, in_set, rank, set_size, decided, dropped, lowest;
    double pvalue, final_pvalue;
    int model = -1;
    if (sscanf(line, "%31[^,],%d,%d,%d,%d,%lf,%d,%d,%lf,%d,%d", label, &cop, &run, &in_set, &rank, &pvalue,
               &set_size, &decided, &final_pvalue, &dropped, &lowest) == 11)
        model = model_of(label);
    if (model < 0 || cop < 0 || cop >= MAX_COP || run < 0 || run >= MAX_RUN || lowest < 0 || lowest >= MAX_COP) {
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

/* The second pass: the set members, which the first pass does not keep. */
static void count_members(int grid, const char *line) {
    char label[32];
    int cop, run, model = -1;
    if (sscanf(line, "%31[^,],%d,%d,", label, &cop, &run) == 3) model = model_of(label);
    if (model < 0 || cop < 0 || cop >= MAX_COP || run < 0 || run >= MAX_RUN) return;
    size_t at = cell_index(cop, run);
    if (!paired(model, at)) return;

    lowest_noise_rank_sum[grid][model] += noise_rank[model][cells[grid][model][at].lowest];
    const char *members = line;
    for (int comma = 0; comma < 11 && members; comma++) {
        members = strchr(members, ',');
        if (members) members++;
    }
    if (!members) return;
    char *end;
    for (long member = strtol(members, &end, 10); end != members; member = strtol(members, &end, 10)) {
        members = end;
        if (member < 0 || member >= MAX_COP) continue;
        member_noise_rank_sum[grid][model] += noise_rank[model][member];
        member_count[grid][model]++;
        if (member != cop) off_target_held[grid][model][member]++;
    }
}

static void read_lines(FILE *rows, int grid, void (*add)(int, const char *)) {
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, rows) > 0) {
        if (header) { header = 0; continue; }
        add(grid, line);
    }
    free(line);
}

/* Every row of a grid, from a gzip-compressed result when the path ends in
   .gz and from a plain progress file otherwise. */
static void read_grid(int grid, const char *path, void (*add)(int, const char *)) {
    size_t path_length = strlen(path);
    if (path_length < 3 || strcmp(path + path_length - 3, ".gz") != 0) {
        FILE *rows = fopen(path, "r");
        assert(rows && "seri_corrected_against_uncorrected_grid: cannot open a progress file");
        read_lines(rows, grid, add);
        fclose(rows);
        return;
    }
    FILE *f = fopen(path, "rb");
    assert(f && "seri_corrected_against_uncorrected_grid: cannot open a result file");
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
    read_lines(rows, grid, add);
    fclose(rows);
    free(text);
}

static void read_noise(void) {
    for (int model = 0; model < N_MODELS; model++)
        for (int cop = 0; cop < MAX_COP; cop++) noise[model][cop] = NAN;
    FILE *f = fopen(NOISE_PATH, "r");
    assert(f && "seri_corrected_against_uncorrected_grid: run bin/response_variance_by_configuration first");
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, f) > 0) {
        if (header) { header = 0; continue; }
        char label[32];
        int cop;
        double value;
        if (sscanf(line, "%31[^,],%d,%lf", label, &cop, &value) != 3) continue;
        int model = model_of(label);
        if (model < 0 || cop < 0 || cop >= MAX_COP) continue;
        noise[model][cop] = value;
    }
    free(line);
    fclose(f);

    for (int model = 0; model < N_MODELS; model++)
        for (int cop = 0; cop < MAX_COP; cop++) {
            noise_rank[model][cop] = 0;
            if (isnan(noise[model][cop])) continue;
            int rank = 1;
            for (int other = 0; other < MAX_COP; other++) rank += noise[model][other] < noise[model][cop];
            noise_rank[model][cop] = rank;
        }
}

static void write_per_grid(FILE *out) {
    fprintf(out, "1. Per model and grid. Percentages are of the benchmarks both grids hold.\n\n");
    fprintf(out, "model, grid, benchmarks, own configuration in the set, alone in the set, first by mean loss, "
                 "mean rank, mean set size, median set size, largest set\n");
    for (int model = 0; model < N_MODELS; model++)
        for (int grid = 0; grid < N_GRIDS; grid++) {
            long n = 0, in_set = 0, alone = 0, first = 0, largest = 0;
            double rank_sum = 0, size_sum = 0;
            Mat sizes = mat_new(MAX_COP * MAX_RUN, 1);
            for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) {
                if (!paired(model, at)) continue;
                const Cell *c = &cells[grid][model][at];
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
    fprintf(out, "2. Per model: whether the benchmark's configuration is in each set, and whether both losses "
                 "give the same configuration the smallest mean loss.\n\n");
    fprintf(out, "model, benchmarks in both, in both sets, corrected set only, uncorrected set only, neither, "
                 "same smallest-loss configuration\n");
    for (int model = 0; model < N_MODELS; model++) {
        long n = 0, both = 0, corrected_only = 0, uncorrected_only = 0, neither = 0, same_lowest = 0;
        for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) {
            if (!paired(model, at)) continue;
            const Cell *c = &cells[GRID_CORRECTED][model][at], *u = &cells[GRID_UNCORRECTED][model][at];
            n++;
            both += c->in_set && u->in_set;
            corrected_only += c->in_set && !u->in_set;
            uncorrected_only += !c->in_set && u->in_set;
            neither += !c->in_set && !u->in_set;
            same_lowest += c->lowest == u->lowest;
        }
        if (!n) continue;
        fprintf(out, "%s, %ld, %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%\n", model_label[model], n, 100.0 * both / n,
                100.0 * corrected_only / n, 100.0 * uncorrected_only / n, 100.0 * neither / n,
                100.0 * same_lowest / n);
    }
    fprintf(out, "\n");
}

/* Per configuration, the own rate of each grid and the off-target rate, over
   the paired benchmarks. Returns 0 when the configuration has none. */
static int configuration_rates(int model, int cop, long n_paired, double own[N_GRIDS], double off[N_GRIDS]) {
    long n = 0, held[N_GRIDS] = { 0, 0 };
    for (int run = 0; run < MAX_RUN; run++) {
        size_t at = cell_index(cop, run);
        if (!paired(model, at)) continue;
        n++;
        for (int grid = 0; grid < N_GRIDS; grid++) held[grid] += cells[grid][model][at].in_set;
    }
    if (!n) return 0;
    for (int grid = 0; grid < N_GRIDS; grid++) {
        own[grid] = (double)held[grid] / n;
        off[grid] = n_paired > n ? (double)off_target_held[grid][model][cop] / (n_paired - n) : NAN;
    }
    return 1;
}

static long count_paired(int model) {
    long n = 0;
    for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) n += paired(model, at);
    return n;
}

static void write_per_configuration(FILE *out) {
    fprintf(out, "3. Per model and configuration: the share of its own benchmarks whose set holds it, with the "
                 "Pearson correlation of the two grids' shares across configurations.\n\n");
    fprintf(out, "model, configurations, corrected at least half, corrected at least 90%%, corrected none, "
                 "uncorrected at least half, uncorrected at least 90%%, uncorrected none, correlation\n");
    for (int model = 0; model < N_MODELS; model++) {
        long n_paired = count_paired(model);
        int n_configurations = 0, half[N_GRIDS] = { 0, 0 }, ninety[N_GRIDS] = { 0, 0 }, none[N_GRIDS] = { 0, 0 };
        Mat shares = mat_new(MAX_COP, N_GRIDS);
        for (int cop = 0; cop < MAX_COP; cop++) {
            double own[N_GRIDS], off[N_GRIDS];
            if (!configuration_rates(model, cop, n_paired, own, off)) continue;
            for (int grid = 0; grid < N_GRIDS; grid++) {
                half[grid] += own[grid] >= 0.5;
                ninety[grid] += own[grid] >= 0.9;
                none[grid] += own[grid] == 0;
                AT(shares, n_configurations, grid) = own[grid];
            }
            n_configurations++;
        }
        if (n_configurations) {
            Mat corrected = mat_slice(shares, 0, n_configurations, GRID_CORRECTED, GRID_CORRECTED + 1);
            Mat uncorrected = mat_slice(shares, 0, n_configurations, GRID_UNCORRECTED, GRID_UNCORRECTED + 1);
            fprintf(out, "%s, %d, %d, %d, %d, %d, %d, %d, %.3f\n", model_label[model], n_configurations,
                    half[GRID_CORRECTED], ninety[GRID_CORRECTED], none[GRID_CORRECTED], half[GRID_UNCORRECTED],
                    ninety[GRID_UNCORRECTED], none[GRID_UNCORRECTED], (double)stats_corr(corrected, uncorrected));
        }
        mat_free(shares);
    }
    fprintf(out, "\n");
}

static void write_noise(FILE *out) {
    fprintf(out, "4. Per model and grid, against each configuration's noise (%s). Correlations are Spearman, "
                 "across the configurations with paired benchmarks and a noise value. Noise ranks run from 1, "
                 "the least noisy.\n\n", NOISE_PATH);
    fprintf(out, "model, grid, configurations, noise against own rate, noise against off-target rate, "
                 "mean noise rank of the smallest-loss configuration, mean noise rank of the set members, "
                 "middle rank\n");
    for (int model = 0; model < N_MODELS; model++) {
        long n_paired = count_paired(model);
        if (!n_paired) continue;
        int n_configurations = 0, n_ranked = 0;
        Mat rates = mat_new(MAX_COP, 2 * N_GRIDS + 1);
        for (int cop = 0; cop < MAX_COP; cop++) {
            n_ranked += noise_rank[model][cop] > 0;
            double own[N_GRIDS], off[N_GRIDS];
            if (isnan(noise[model][cop]) || !configuration_rates(model, cop, n_paired, own, off)) continue;
            AT(rates, n_configurations, 0) = noise[model][cop];
            for (int grid = 0; grid < N_GRIDS; grid++) {
                AT(rates, n_configurations, 1 + grid) = own[grid];
                AT(rates, n_configurations, 1 + N_GRIDS + grid) = off[grid];
            }
            n_configurations++;
        }
        Mat noise_column = mat_slice(rates, 0, n_configurations, 0, 1);
        for (int grid = 0; grid < N_GRIDS; grid++) {
            Mat own = mat_slice(rates, 0, n_configurations, 1 + grid, 2 + grid);
            Mat off = mat_slice(rates, 0, n_configurations, 1 + N_GRIDS + grid, 2 + N_GRIDS + grid);
            fprintf(out, "%s, %s, %d, %.3f, %.3f, %.1f, %.1f, %.1f\n", model_label[model], grid_label[grid],
                    n_configurations, (double)stats_spearman(noise_column, own),
                    (double)stats_spearman(noise_column, off), lowest_noise_rank_sum[grid][model] / n_paired,
                    member_noise_rank_sum[grid][model] / member_count[grid][model], (n_ranked + 1) / 2.0);
        }
        mat_free(rates);
    }
}

int main(void) {
    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) {
            cells[grid][model] = malloc((size_t)MAX_COP * MAX_RUN * sizeof(Cell));
            assert(cells[grid][model]);
            for (size_t at = 0; at < (size_t)MAX_COP * MAX_RUN; at++) cells[grid][model][at].in_set = -1;
        }
    read_noise();

    struct stat st;
    const char *uncorrected_path = stat(UNCORRECTED_RESULT_PATH, &st) == 0 ? UNCORRECTED_RESULT_PATH
                                                                           : UNCORRECTED_PROGRESS_PATH;
    read_grid(GRID_CORRECTED, CORRECTED_PATH, store_row);
    read_grid(GRID_UNCORRECTED, uncorrected_path, store_row);
    read_grid(GRID_CORRECTED, CORRECTED_PATH, count_members);
    read_grid(GRID_UNCORRECTED, uncorrected_path, count_members);

    mkdir("out", 0755);
    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "seri_corrected_against_uncorrected_grid: cannot open the report");
    fprintf(out, "Corrected: %s. Uncorrected: %s.\n", CORRECTED_PATH, uncorrected_path);
    for (int grid = 0; grid < N_GRIDS; grid++)
        fprintf(out, "%s: %ld rows that did not parse, %ld benchmarks found twice.\n", grid_label[grid],
                n_unparsed[grid], n_duplicated[grid]);
    fprintf(out, "\n");
    write_per_grid(out);
    write_paired(out);
    write_per_configuration(out);
    write_noise(out);
    fclose(out);

    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) free(cells[grid][model]);
    return 0;
}
