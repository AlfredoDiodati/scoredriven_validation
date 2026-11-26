/*
One LaTeX table for the true configuration recovery experiment, the four grids
side by side: the confidence set of Hansen et al. (MCS_TR, 2000 bootstrap
resamples) or of Seri et al. (seri_mcs), each with the naive loss or the
adjusted one, both at alpha 0.05.

- Naive: squared error between response vectors averaged over the entries.
- Adjusted: the same squared error minus the variance correction of the grid,
  montecarlo/recovery_loss.h for Hansen et al. and
  montecarlo/own_deviation_loss.h for Seri et al.

Per model and grid, over every benchmark (every run of every configuration):
the share whose set holds the benchmark's configuration, the share where that
configuration has the smallest mean loss, and the median size of the final
set. Rows are both states in one vector under the three state weights of
applications/lp_system.h other than the causal HP cycle at the estimated
lambda: the rows montecarlo/identification_table.c leaves out.
docs/OLD_EXPERIMENTS.md says why they are left out. A row label carries C for
the causal HP cycle (none: two-sided) and E for the lambda estimated on each
run's switching series (none: lambda 1600). A model's cell is filled only
once the grid holds all 1,000,000 of its benchmarks, and is "-" otherwise,
so the table can be written while a grid is still running. A grid is read
from its result, or from its progress file when the result does not exist yet.

Writes, in out/old/montecarlo_other_state_weights/identification_table/,
identification_table.tex, identification_table_counts.txt (the file each grid
was read from and the benchmarks it holds per model) and
identification_table.log (when each grid was read and how long it took).
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <et_al./stats.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define OUTPUT_DIR "out/old/montecarlo_other_state_weights/identification_table"
#define TABLE_PATH OUTPUT_DIR "/identification_table.tex"
#define COUNTS_PATH OUTPUT_DIR "/identification_table_counts.txt"
#define LOG_PATH OUTPUT_DIR "/identification_table.log"
#define MAX_COP 1024
#define MAX_RUN 1000
#define FULL_GRID 1000000L

enum { HANSEN_NAIVE, HANSEN_ADJUSTED, SERI_NAIVE, SERI_ADJUSTED, N_GRIDS };

static const char *const grid_label[N_GRIDS] = { "Hansen et al., naive", "Hansen et al., adjusted",
                                                 "Seri et al., naive", "Seri et al., adjusted" };
static const char *const result_path[N_GRIDS] = {
    "out/old/montecarlo_other_state_weights/uncorrected_true_configuration_recovery.csv.gz",
    "out/old/montecarlo_other_state_weights/true_configuration_recovery.csv.gz",
    "out/old/montecarlo_other_state_weights/seri_uncorrected_true_configuration_recovery.csv.gz",
    "out/old/montecarlo_other_state_weights/seri_true_configuration_recovery.csv.gz",
};
static const char *const progress_path[N_GRIDS] = {
    "out/old/montecarlo_other_state_weights/uncorrected_true_configuration_recovery_progress.csv",
    "out/old/montecarlo_other_state_weights/true_configuration_recovery_progress.csv",
    "out/old/montecarlo_other_state_weights/seri_uncorrected_true_configuration_recovery_progress.csv",
    "out/old/montecarlo_other_state_weights/seri_true_configuration_recovery_progress.csv",
};

/* The table's rows, in the order they are printed. */
#define N_ROWS 3
static const int row_model[N_ROWS] = { MODEL_LP_NL, MODEL_LP_NL_CAUSAL_HP, MODEL_LP_NL_ESTIMATED_LAMBDA };
static const char *const row_label[N_ROWS] = { "NLP", "NLP C", "NLP E" };

typedef struct {
    long n, in_set, lowest, unparsed, duplicated;
    unsigned char *seen;
    Mat set_sizes;
} Tally;

typedef struct {
    const char *source;
    Tally tally[N_MODELS];
} Grid;

static Grid grids[N_GRIDS];

static int model_of(const char *label) {
    for (int model = 0; model < N_MODELS; model++)
        if (strcmp(label, model_label[model]) == 0) return model;
    return -1;
}

static void add_row(Grid *grid, const char *line) {
    char label[MODEL_LABEL_CAPACITY];
    int cop, run, in_set, rank, set_size, model = -1;
    double pvalue;
    if (sscanf(line, "%63[^,],%d,%d,%d,%d,%lf,%d,", label, &cop, &run, &in_set, &rank, &pvalue, &set_size) == 7)
        model = model_of(label);
    /* A progress file being appended can end in half a row, which fails here
       unless it reaches the eleventh comma. */
    int commas = 0;
    for (const char *c = line; *c; c++) commas += *c == ',';
    if (model < 0 || commas < 11 || cop < 0 || cop >= MAX_COP || run < 0 || run >= MAX_RUN) {
        grid->tally[0].unparsed++;
        return;
    }
    Tally *t = &grid->tally[model];
    size_t at = (size_t)cop * MAX_RUN + run;
    if (t->seen[at]) {
        t->duplicated++;
        return;
    }
    t->seen[at] = 1;
    AT(t->set_sizes, t->n, 0) = set_size;
    t->n++;
    t->in_set += in_set == 1;
    t->lowest += rank == 1;
}

static void read_lines(FILE *rows, Grid *grid) {
    char *line = NULL;
    size_t capacity = 0;
    int header = 1;
    while (getline(&line, &capacity, rows) > 0) {
        if (header) { header = 0; continue; }
        add_row(grid, line);
    }
    free(line);
}

static int exists(const char *path) {
    struct stat info;
    return stat(path, &info) == 0;
}

static void read_grid(int index) {
    Grid *grid = &grids[index];
    for (int model = 0; model < N_MODELS; model++) {
        grid->tally[model].seen = calloc((size_t)MAX_COP * MAX_RUN, 1);
        assert(grid->tally[model].seen);
        grid->tally[model].set_sizes = mat_new(MAX_COP * MAX_RUN, 1);
    }
    if (exists(result_path[index])) {
        grid->source = result_path[index];
        FILE *f = fopen(grid->source, "rb");
        assert(f && "other_state_weights_identification_table: cannot open a result file");
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
        assert(text && "other_state_weights_identification_table: a result did not inflate, out of memory or not gzip");
        free(compressed);
        FILE *rows = fmemopen(text, text_length, "r");
        assert(rows);
        read_lines(rows, grid);
        fclose(rows);
        free(text);
    } else if (exists(progress_path[index])) {
        grid->source = progress_path[index];
        FILE *rows = fopen(grid->source, "r");
        assert(rows && "other_state_weights_identification_table: cannot open a progress file");
        read_lines(rows, grid);
        fclose(rows);
    }
    for (int model = 0; model < N_MODELS; model++) {
        free(grid->tally[model].seen);
        grid->tally[model].seen = NULL;
    }
}

static void write_cells(FILE *out, int index, int model) {
    const Tally *t = &grids[index].tally[model];
    if (t->n != FULL_GRID) {
        fprintf(out, " & - & - & -");
        return;
    }
    double median_size = stats_median(mat_slice(t->set_sizes, 0, (int)t->n, 0, 1));
    fprintf(out, " & %.2f\\%% & %.2f\\%% & %g", 100.0 * t->in_set / t->n, 100.0 * t->lowest / t->n, median_size);
}

static void write_table(void) {
    FILE *out = fopen(TABLE_PATH, "w");
    assert(out && "other_state_weights_identification_table: cannot open the table");
    fprintf(out, "\\begin{table}[t]\n\\centering\n");
    fprintf(out, "\\caption{Recovery Study for Simulated Validation Data}\n");
    fprintf(out, "\\label{tab:other_state_weights_identification_table}\n\n");
    fprintf(out, "\\begin{tabular}{l rrr rrr}\n\\toprule\n");
    fprintf(out, " & \\multicolumn{3}{c}{Naive} & \\multicolumn{3}{c}{Adjusted} \\\\\n");
    fprintf(out, "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n");
    /* Block m is the procedure whose set is M^(m) in the paper: 1 for Hansen et
       al., 2 for Seri et al. */
    const int block_naive[2] = { HANSEN_NAIVE, SERI_NAIVE };
    const int block_adjusted[2] = { HANSEN_ADJUSTED, SERI_ADJUSTED };
    for (int block = 0; block < 2; block++) {
        int m = block + 1;
        if (block > 0) fprintf(out, "&&&&&&\\\\\n");
        fprintf(out, "\\rule{0pt}{3ex}");
        for (int column_group = 0; column_group < 2; column_group++)
            fprintf(out, " & $\\params^*_%d \\in \\mathcal{M}\\explab{%d}$ & $\\params^*_%d=\\params\\udagg_%d$ & "
                         "$\\text{Med}|\\mathcal{M}\\explab{%d}|$", m, m, m, m, m);
        fprintf(out, " \\\\\n");
        for (int column = 2; column <= 7; column++) fprintf(out, "\\cmidrule(lr){%d-%d}", column, column);
        fprintf(out, "\n");
        for (int row = 0; row < N_ROWS; row++) {
            fprintf(out, "%s", row_label[row]);
            write_cells(out, block_naive[block], row_model[row]);
            write_cells(out, block_adjusted[block], row_model[row]);
            fprintf(out, " \\\\\n");
        }
    }
    fprintf(out, "\\bottomrule\n\\end{tabular}\n");
    fprintf(out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    fprintf(out, "\\begin{tablenotes}\n");
    fprintf(out, "\\item \\footnotesize \\textit{Note:} share of the 1,000,000 Monte Carlo runs whose \\mcs holds "
                 "the true configuration $\\params^*_m \\in \\mathcal{M}\\explab{m}$; share where the true "
                 "configuration has the smallest mean loss $\\params^*_m=\\params\\udagg_m$; and median size of "
                 "the final MCS $\\text{Med}|\\mathcal{M}\\explab{m}|$; at a 5\\%% confidence level for both "
                 "procedures of \\cite{Hansen-Lunde-Nason-Eco-2011} and \\cite{Seri-et-al-ES-2021} and all "
                 "auxiliary models. Naive column uses squared error as loss; adjusted column uses squared error "
                 "minus the variance corrections of Equations \\eqref{eq:hansen_mcs_loss} and "
                 "\\eqref{eq:seri_mcs_loss}. NLP compares the responses of both states in one vector; C "
                 "is the causal HP cycle of the switching series, at $t$ from periods $1$ to $t$ only, otherwise "
                 "two-sided; E is the $\\lambda$ found by maximum likelihood on each run's own switching series, "
                 "otherwise $\\lambda = 1600$.\n");
    fprintf(out, "\\end{tablenotes}\n\\end{table}\n");
    fclose(out);
}

static void write_counts(void) {
    FILE *out = fopen(COUNTS_PATH, "w");
    assert(out && "other_state_weights_identification_table: cannot open the counts");
    fprintf(out, "Benchmarks per model and grid; a cell of %s is filled only at %ld.\n\n", TABLE_PATH, FULL_GRID);
    for (int index = 0; index < N_GRIDS; index++) {
        const Grid *grid = &grids[index];
        fprintf(out, "%s: %s\n", grid_label[index], grid->source ? grid->source : "no result and no progress file");
        if (!grid->source) continue;
        fprintf(out, "  rows that did not parse %ld\n", grid->tally[0].unparsed);
        for (int model = 0; model < N_MODELS; model++)
            fprintf(out, "  %s: %ld benchmarks, %ld found twice\n", model_label[model], grid->tally[model].n,
                    grid->tally[model].duplicated);
    }
    fclose(out);
}

static void log_line(FILE *log, const char *what) {
    char stamp[32];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", localtime(&now));
    fprintf(log, "%s %s\n", stamp, what);
    fflush(log);
}

int main(void) {
    mkdir("out", 0755);
    mkdir("out/old", 0755);
    mkdir("out/old/montecarlo_other_state_weights", 0755);
    mkdir(OUTPUT_DIR, 0755);
    FILE *log = fopen(LOG_PATH, "a");
    assert(log && "other_state_weights_identification_table: cannot open the log");
    log_line(log, "started");
    for (int index = 0; index < N_GRIDS; index++) {
        time_t start = time(NULL);
        read_grid(index);
        char what[400];
        snprintf(what, sizeof what, "%s read from %s in %.0f s", grid_label[index],
                 grids[index].source ? grids[index].source : "nothing", difftime(time(NULL), start));
        log_line(log, what);
    }
    write_table();
    write_counts();
    log_line(log, "wrote " TABLE_PATH " and " COUNTS_PATH);
    fclose(log);
    for (int index = 0; index < N_GRIDS; index++)
        for (int model = 0; model < N_MODELS; model++) mat_free(grids[index].tally[model].set_sizes);
    return 0;
}
