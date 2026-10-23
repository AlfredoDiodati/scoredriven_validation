/*
The per-configuration tables of the discarded naive-loss grid, computed for
the grid of montecarlo/seri_true_configuration_recovery.c and, from the same
code, for the naive-loss grid of old/naive_mae_loss/montecarlo/sweep_grid.c,
so the two can be read side by side. They are the tables the Python scripts
old/naive_mae_loss/montecarlo/sweep_grid_cop_0191.py,
sweep_grid_presence.py and sweep_grid_identifiability_plots.py wrote for the
naive grid alone.

- Seri: squared error minus each run's own deviation from its configuration's
  mean response, scaled by n/(n-1), then et_al's seri_mcs at alpha 0.05.
- Naive: mean absolute error between response vectors, no correction, then
  et_al's mcs, MCS_TR, alpha 0.05, block length 1, bootstrap variance, 2000
  resamples, seed 123 stream 0.

Writes, per grid:

1. out/montecarlo/seri_against_naive_by_configuration_cop_0191.txt and .tex: per
   model, the share of benchmarks whose set holds the configuration, holds it
   alone, and gives it the smallest mean loss, in three tables: cop_0191 as the
   truth (its own 1000 benchmarks), cop_0191 as a wrong answer (the other
   999,000), and every configuration as the truth (all 1,000,000, the
   configuration being the benchmark's own).
2. out/montecarlo/seri_against_naive_by_configuration_presence.txt: per model
   and configuration, the own rate (share of its own benchmarks whose set holds
   it) and the off-target rate (share of the other configurations' benchmarks
   whose set holds it), summarised as in the naive grid's
   sweep_grid_presence.txt.
3. out/montecarlo/seri_against_naive_by_configuration_identifiability.csv: per
   grid, model and configuration, the percentage of its own sets that hold it,
   and the other configuration found in most of its sets with the percentage of
   its sets that hold that one; a tie goes to the lower configuration number;
   and the percentage of the other configurations' sets that hold it.
   montecarlo/seri_against_naive_by_configuration_plots.py draws it.

Every row of both files is read; a row that does not parse, a benchmark found
twice, and a row whose member list disagrees with its set size are counted
and reported at the top of the text files.
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define SERI_PATH "out/montecarlo/seri_true_configuration_recovery.csv.gz"
#define NAIVE_PATH "out/old/naive_mae_loss/montecarlo/sweep_grid.csv.gz"
#define COP_0191_TEXT_PATH "out/montecarlo/seri_against_naive_by_configuration_cop_0191.txt"
#define COP_0191_LATEX_PATH "out/montecarlo/seri_against_naive_by_configuration_cop_0191.tex"
#define PRESENCE_PATH "out/montecarlo/seri_against_naive_by_configuration_presence.txt"
#define IDENTIFIABILITY_PATH "out/montecarlo/seri_against_naive_by_configuration_identifiability.csv"
#define MAX_COP 1001
#define MAX_RUN 1000
#define FOCUS_COP 191
#define TOP 10
#define N_GRIDS 2

enum { GRID_SERI, GRID_NAIVE };
static const char *const grid_label[N_GRIDS] = { "Seri", "naive" };
static const char *const model_name[N_MODELS] = { "qvarma", "linear", "state 1", "state 2", "both states" };

/* The three outcomes of the cop_0191 tables, for one configuration in one set. */
enum { OUTCOME_IN_SET, OUTCOME_ALONE, OUTCOME_SMALLEST_LOSS, N_OUTCOMES };
static const char *const outcome_name[N_OUTCOMES] = { "in the set", "alone in the set", "smallest mean loss" };

/* The three cop_0191 tables. */
enum { TABLE_TRUTH, TABLE_WRONG_ANSWER, TABLE_EVERY_TRUTH, N_TABLES };

typedef struct {
    long benchmarks, outcome[N_OUTCOMES];
} Tally;

typedef struct {
    Tally table[N_TABLES];
    long benchmarks[MAX_COP], own_hits[MAX_COP], off_hits[MAX_COP];
    /* appearances[benchmark_cop * MAX_COP + member], member != benchmark_cop */
    int *appearances;
} ModelCounts;

static ModelCounts counts[N_GRIDS][N_MODELS];
/* set_size_of[grid][model][cop * MAX_RUN + run], 0 for a benchmark not in the file. */
static short *set_size_of[N_GRIDS][N_MODELS];
static long n_unparsed[N_GRIDS], n_duplicated[N_GRIDS], n_size_mismatch[N_GRIDS];

static void tally(Tally *t, const int *members, int n_members, int lowest_loss_cop, int cop) {
    int in_set = 0;
    for (int i = 0; i < n_members; i++) in_set |= members[i] == cop;
    t->benchmarks++;
    t->outcome[OUTCOME_IN_SET] += in_set;
    t->outcome[OUTCOME_ALONE] += n_members == 1 && members[0] == cop;
    t->outcome[OUTCOME_SMALLEST_LOSS] += lowest_loss_cop == cop;
}

static void add_row(int grid, const char *line) {
    char label[32];
    int cop, run, in_set, rank, set_size, decided, dropped, lowest, offset;
    double pvalue, final_pvalue;
    if (sscanf(line, "%31[^,],%d,%d,%d,%d,%lf,%d,%d,%lf,%d,%d,%n", label, &cop, &run, &in_set, &rank, &pvalue,
               &set_size, &decided, &final_pvalue, &dropped, &lowest, &offset) != 11
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

    int members[MAX_COP], n_members = 0;
    const char *at = line + offset;
    char *end;
    for (long member = strtol(at, &end, 10); end != at; member = strtol(at, &end, 10)) {
        if (member < 0 || member >= MAX_COP || n_members == MAX_COP) {
            n_unparsed[grid]++;
            return;
        }
        members[n_members++] = (int)member;
        at = end;
    }
    if (n_members != set_size) n_size_mismatch[grid]++;

    short *stored_size = &set_size_of[grid][model][(size_t)cop * MAX_RUN + run];
    if (*stored_size) n_duplicated[grid]++;
    *stored_size = (short)n_members;

    ModelCounts *c = &counts[grid][model];
    tally(&c->table[TABLE_EVERY_TRUTH], members, n_members, lowest, cop);
    tally(&c->table[cop == FOCUS_COP ? TABLE_TRUTH : TABLE_WRONG_ANSWER], members, n_members, lowest, FOCUS_COP);
    c->benchmarks[cop]++;
    c->own_hits[cop] += in_set;
    for (int i = 0; i < n_members; i++)
        if (members[i] != cop) {
            c->off_hits[members[i]]++;
            c->appearances[(size_t)cop * MAX_COP + members[i]]++;
        }
}

static void read_grid(int grid, const char *path) {
    FILE *f = fopen(path, "rb");
    assert(f && "seri_against_naive_by_configuration: cannot open a result file");
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

static void write_row_checks(FILE *out) {
    fprintf(out, "%s against %s.\n", SERI_PATH, NAIVE_PATH);
    for (int grid = 0; grid < N_GRIDS; grid++)
        fprintf(out, "%s: %ld rows that did not parse, %ld benchmarks found twice, %ld rows whose member list "
                     "disagrees with their set size.\n", grid_label[grid], n_unparsed[grid], n_duplicated[grid],
                n_size_mismatch[grid]);
    fprintf(out, "\n");
}

static double percent(const Tally *t, int outcome) {
    return t->benchmarks ? 100.0 * t->outcome[outcome] / t->benchmarks : 0.0;
}

/* The median number of configurations in the set over the benchmarks of one
   cop_0191 table: cop_0191's own, the others, or all. */
static double median_set_size(int grid, int model, int table) {
    const short *sizes = set_size_of[grid][model];
    long n = 0;
    for (int cop = 0; cop < MAX_COP; cop++)
        if (table == TABLE_EVERY_TRUTH || (table == TABLE_TRUTH) == (cop == FOCUS_COP))
            for (int run = 0; run < MAX_RUN; run++) n += sizes[(size_t)cop * MAX_RUN + run] > 0;
    if (!n) return NAN;
    Mat selected = mat_new((int)n, 1);
    long at = 0;
    for (int cop = 0; cop < MAX_COP; cop++)
        if (table == TABLE_EVERY_TRUTH || (table == TABLE_TRUTH) == (cop == FOCUS_COP))
            for (int run = 0; run < MAX_RUN; run++) {
                short size = sizes[(size_t)cop * MAX_RUN + run];
                if (size > 0) AT(selected, at++, 0) = size;
            }
    double median = stats_median(selected);
    mat_free(selected);
    return median;
}

/* The title of a cop_0191 table, with the benchmark count of the t-QVARMA as
   the old script had it; every model holds the same benchmarks. */
static void table_title(char *title, size_t size, int table) {
    long n = counts[GRID_SERI][MODEL_QVARMA].table[table].benchmarks;
    char count[32];
    if (n >= 1000000) snprintf(count, sizeof count, "%ld,%03ld,%03ld", n / 1000000, n / 1000 % 1000, n % 1000);
    else if (n >= 1000) snprintf(count, sizeof count, "%ld,%03ld", n / 1000, n % 1000);
    else snprintf(count, sizeof count, "%ld", n);
    if (table == TABLE_TRUTH)
        snprintf(title, size, "cop_%04d is the true configuration: share of its own %s benchmarks", FOCUS_COP, count);
    else if (table == TABLE_WRONG_ANSWER)
        snprintf(title, size, "cop_%04d is not the true configuration: share of the other %s benchmarks where "
                              "cop_%04d is", FOCUS_COP, count, FOCUS_COP);
    else
        snprintf(title, size, "every true configuration: share of all %s benchmarks where the true configuration is",
                 count);
}

static void write_cop_0191_text(void) {
    FILE *out = fopen(COP_0191_TEXT_PATH, "w");
    assert(out && "seri_against_naive_by_configuration: cannot open the cop_0191 text file");
    write_row_checks(out);
    for (int table = 0; table < N_TABLES; table++) {
        char title[256];
        table_title(title, sizeof title, table);
        fprintf(out, "%s\n", title);
        fprintf(out, "order: %s, %s, %s, median configurations in the set; Seri grid, then naive grid\n",
                outcome_name[0], outcome_name[1], outcome_name[2]);
        for (int model = 0; model < N_MODELS; model++) {
            fprintf(out, "%s:", model_name[model]);
            for (int grid = 0; grid < N_GRIDS; grid++) {
                const Tally *t = &counts[grid][model].table[table];
                fprintf(out, "%s %s %.2f%%, %.2f%%, %.2f%%, %g", grid ? ";" : "", grid_label[grid],
                        percent(t, OUTCOME_IN_SET), percent(t, OUTCOME_ALONE), percent(t, OUTCOME_SMALLEST_LOSS),
                        median_set_size(grid, model, table));
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }
    fclose(out);
}

static void write_cop_0191_latex(void) {
    static const char *const label[N_TABLES] = { "truth", "wrong_answer", "every_truth" };
    FILE *out = fopen(COP_0191_LATEX_PATH, "w");
    assert(out && "seri_against_naive_by_configuration: cannot open the cop_0191 LaTeX file");
    for (int table = 0; table < N_TABLES; table++) {
        char title[256], caption[512];
        table_title(title, sizeof title, table);
        size_t length = 0;
        for (const char *c = title; *c && length + 2 < sizeof caption; c++) {
            if (*c == '_' || *c == '%') caption[length++] = '\\';
            caption[length++] = *c;
        }
        caption[length] = '\0';
        fprintf(out, "\\begin{table}[htbp]\n\\centering\n");
        fprintf(out, "\\caption{%s}\n\\label{tab:seri_against_naive_cop_%04d_%s}\n", caption, FOCUS_COP,
                label[table]);
        fprintf(out, "\\begin{tabular}{lrrrrrrrr}\n\\hline\n");
        fprintf(out, " & \\multicolumn{4}{c}{Seri} & \\multicolumn{4}{c}{naive} \\\\\n");
        fprintf(out, "model & %s & %s & %s & median set size & %s & %s & %s & median set size \\\\\n\\hline\n",
                outcome_name[0], outcome_name[1], outcome_name[2], outcome_name[0], outcome_name[1], outcome_name[2]);
        for (int model = 0; model < N_MODELS; model++) {
            fprintf(out, "%s", model_name[model]);
            for (int grid = 0; grid < N_GRIDS; grid++) {
                for (int outcome = 0; outcome < N_OUTCOMES; outcome++)
                    fprintf(out, " & %.2f\\%%", percent(&counts[grid][model].table[table], outcome));
                fprintf(out, " & %g", median_set_size(grid, model, table));
            }
            fprintf(out, " \\\\\n");
        }
        fprintf(out, "\\hline\n\\end{tabular}\n\\end{table}\n\n");
    }
    fclose(out);
}

static const double *off_rate_for_sort;

/* Highest off-target rate first, a tie to the lower configuration number. */
static int by_off_rate(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    if (off_rate_for_sort[x] != off_rate_for_sort[y]) return off_rate_for_sort[x] < off_rate_for_sort[y] ? 1 : -1;
    return x - y;
}

static void write_presence_block(FILE *out, int grid, int model) {
    const ModelCounts *c = &counts[grid][model];
    int cops[MAX_COP], n_cops = 0;
    long total = 0, off_slots = 0;
    for (int cop = 0; cop < MAX_COP; cop++) {
        if (c->benchmarks[cop]) cops[n_cops++] = cop;
        total += c->benchmarks[cop];
        off_slots += c->off_hits[cop];
    }
    if (!n_cops) return;
    double own[MAX_COP], off[MAX_COP], mean_own = 0, mean_off = 0;
    for (int i = 0; i < n_cops; i++) {
        int cop = cops[i];
        own[cop] = (double)c->own_hits[cop] / c->benchmarks[cop];
        off[cop] = (double)c->off_hits[cop] / (total - c->benchmarks[cop]);
        mean_own += own[cop];
        mean_off += off[cop];
    }
    mean_own /= n_cops;
    mean_off /= n_cops;
    int off_above_own = 0, off_above_mean_own = 0;
    for (int i = 0; i < n_cops; i++) {
        off_above_own += off[cops[i]] > own[cops[i]];
        off_above_mean_own += off[cops[i]] > mean_own;
    }
    off_rate_for_sort = off;
    qsort(cops, (size_t)n_cops, sizeof cops[0], by_off_rate);
    int n_top = n_cops < TOP ? n_cops : TOP;
    long top_slots = 0;
    for (int i = 0; i < n_top; i++) top_slots += c->off_hits[cops[i]];

    fprintf(out, "%s, %s grid: mean own-CoP rate %.4f; mean off-target presence per CoP %.4f; median configurations "
                 "in the set %g\n", model_label[model], grid_label[grid], mean_own, mean_off,
            median_set_size(grid, model, TABLE_EVERY_TRUTH));
    fprintf(out, "  CoPs whose off-target presence rate exceeds their own recovery rate: %d of %d\n", off_above_own,
            n_cops);
    fprintf(out, "  CoPs whose off-target presence exceeds the mean own rate %.4f: %d\n", mean_own,
            off_above_mean_own);
    fprintf(out, "  top %d by off-target presence (cop, off rate, own rate):\n", TOP);
    for (int i = 0; i < n_top; i++) fprintf(out, "    cop_%04d %.4f %.4f\n", cops[i], off[cops[i]], own[cops[i]]);
    fprintf(out, "  share of all off-target set slots taken by the top %d: %.3f\n", TOP,
            off_slots ? (double)top_slots / off_slots : 0.0);
}

static void write_presence(void) {
    FILE *out = fopen(PRESENCE_PATH, "w");
    assert(out && "seri_against_naive_by_configuration: cannot open the presence file");
    write_row_checks(out);
    fprintf(out, "Own rate: the share of a configuration's own benchmarks whose set holds it. Off-target rate: the "
                 "share of the other configurations' benchmarks whose set holds it.\n\n");
    for (int model = 0; model < N_MODELS; model++)
        for (int grid = 0; grid < N_GRIDS; grid++) {
            write_presence_block(out, grid, model);
            fprintf(out, "\n");
        }
    fclose(out);
}

static void write_identifiability(void) {
    FILE *out = fopen(IDENTIFIABILITY_PATH, "w");
    assert(out && "seri_against_naive_by_configuration: cannot open the identifiability file");
    fprintf(out, "grid,model,benchmark_cop,runs,true_cop_percent,top_other_cop,top_other_percent,off_target_percent\n");
    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) {
            const ModelCounts *c = &counts[grid][model];
            long total = 0;
            for (int cop = 0; cop < MAX_COP; cop++) total += c->benchmarks[cop];
            for (int cop = 0; cop < MAX_COP; cop++) {
                long runs = c->benchmarks[cop];
                if (!runs) continue;
                const int *row = &c->appearances[(size_t)cop * MAX_COP];
                int top_other = -1;
                for (int member = 0; member < MAX_COP; member++)
                    if (row[member] && (top_other < 0 || row[member] > row[top_other])) top_other = member;
                fprintf(out, "%s,%s,%d,%ld,%.3f,", grid_label[grid], model_label[model], cop, runs,
                        100.0 * c->own_hits[cop] / runs);
                if (top_other >= 0) fprintf(out, "%d,%.3f,", top_other, 100.0 * row[top_other] / runs);
                else fprintf(out, ",0.000,");
                fprintf(out, "%.4f\n", 100.0 * c->off_hits[cop] / (total - runs));
            }
        }
    fclose(out);
}

int main(void) {
    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) {
            counts[grid][model].appearances = calloc((size_t)MAX_COP * MAX_COP, sizeof(int));
            set_size_of[grid][model] = calloc((size_t)MAX_COP * MAX_RUN, sizeof(short));
            assert(counts[grid][model].appearances && set_size_of[grid][model]);
        }
    read_grid(GRID_SERI, SERI_PATH);
    read_grid(GRID_NAIVE, NAIVE_PATH);

    mkdir("out", 0755);
    mkdir("out/montecarlo", 0755);
    write_cop_0191_text();
    write_cop_0191_latex();
    write_presence();
    write_identifiability();

    for (int grid = 0; grid < N_GRIDS; grid++)
        for (int model = 0; model < N_MODELS; model++) {
            free(counts[grid][model].appearances);
            free(set_size_of[grid][model]);
        }
    return 0;
}
