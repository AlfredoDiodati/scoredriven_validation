/*
Summaries of the fit metric (docs/MONTECARLO_FIT_METRIC.md) from the rows
bin/fit_metric_oracle and bin/fit_metric_bootstrap write.

Tables, one row per model:
  oracle            every benchmark of the grid, 1,000,000 per model: medians of
                    d and sigma2, quantiles of s and v, and the mean of v
  oracle by winner  the median v where the configuration with the smallest mean
                    loss is the benchmark's own, and where it is another
  bootstrap         run 0 of every configuration, 1000 per model: the oracle and
                    the bootstrap sigma2, quantiles of their ratio, their
                    Spearman rank correlation across configurations, and v
                    against v_star
  resamples         how many resamples gave a usable response, with how the
                    t-QVARMA refits went underneath
  block length      the median ratio sigma2_star / sigma2 under each block
                    length run on the first configurations, one row per file
                    fit_metric_bootstrap_run0_B<B>_block<L>_iterations<I>_first<n>.csv
The bootstrap tables are left out while montecarlo/out/fit_metric_bootstrap.csv
does not exist.

Quantiles, medians and rank correlations are et_al's stats_quantile,
stats_median and stats_spearman. Reads montecarlo/out/fit_metric_oracle.csv.gz
and the bootstrap files beside it. Writes montecarlo/out/fit_metric_report.txt
and the same tables as LaTeX in montecarlo/out/fit_metric_report.tex. Nothing
printed.
*/

#include "montecarlo/fit_metric.h"
#include <et_al./stats.h>
#include <et_al./frame/gzip.h>
#include <dirent.h>

#define OUT_DIR "montecarlo/out"
#define ORACLE_PATH OUT_DIR "/fit_metric_oracle.csv.gz"
#define BOOTSTRAP_PATH OUT_DIR "/fit_metric_bootstrap.csv"
#define REPORT_PATH OUT_DIR "/fit_metric_report.txt"
#define LATEX_PATH OUT_DIR "/fit_metric_report.tex"
#define MAX_COLUMNS 12
#define MAX_ROWS 16
#define CELL 48

static const char *const model_name[N_MODELS] = { "qvarma", "linear", "state 1", "state 2", "both states" };

typedef struct {
    const char *label, *title;
    int n_columns, n_rows;
    char header[MAX_COLUMNS][CELL];
    char cell[MAX_ROWS][MAX_COLUMNS][CELL];
    char note[256];
} Table;

static void set_header(Table *t, int n_columns, const char *const *names) {
    t->n_columns = n_columns;
    for (int i = 0; i < n_columns; i++) snprintf(t->header[i], CELL, "%s", names[i]);
}

/* A growable column of doubles. */
typedef struct { double *value; int n, cap; } Column;

static void push(Column *c, double value) {
    if (c->n == c->cap) {
        c->cap = c->cap ? 2 * c->cap : 1024;
        c->value = realloc(c->value, (size_t)c->cap * sizeof(double));
        assert(c->value && "fit_metric_report: out of memory");
    }
    c->value[c->n++] = value;
}

static Mat as_mat(const Column *c) { return (Mat){ c->n, 1, 1, c->value }; }
static double quantile(const Column *c, double p) { return (double)stats_quantile(as_mat(c), (mreal)p); }
static double median(const Column *c) { return (double)stats_median(as_mat(c)); }
static double spearman(const Column *x, const Column *y) { return (double)stats_spearman(as_mat(x), as_mat(y)); }

static int model_of(const char *label) {
    for (int model = 0; model < N_MODELS; model++)
        if (strcmp(label, model_label[model]) == 0) return model;
    assert(0 && "fit_metric_report: unknown model label");
    return -1;
}

/* The oracle file's columns by model: d, sigma2, s, v, and v split by whether
   the configuration with the smallest mean loss is the benchmark's own. */
typedef struct { Column d, sigma2, s, v, v_own, v_other; } OracleColumns;

static void read_oracle(OracleColumns *oracle) {
    FILE *f = fopen(ORACLE_PATH, "rb");
    assert(f && "fit_metric_report: montecarlo/out/fit_metric_oracle.csv.gz is missing");
    fseek(f, 0, SEEK_END);
    size_t length = (size_t)ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc(length);
    size_t n_read = fread(compressed, 1, length, f);
    assert(n_read == length && "fit_metric_report: short read of the oracle file");
    (void)n_read;
    fclose(f);
    size_t text_length;
    char *text = (char *)gzip_inflate(compressed, length, &text_length);
    free(compressed);

    char *line = strchr(text, '\n') + 1;
    while (line < text + text_length && *line) {
        char *end = strchr(line, '\n');
        if (end) *end = 0;
        char label[32];
        int cop, run, lowest;
        double d, sigma2, s, v;
        int matched = sscanf(line, "%31[^,],%d,%d,%d,%lf,%lf,%lf,%lf", label, &cop, &run, &lowest, &d, &sigma2, &s, &v);
        assert(matched == 8 && "fit_metric_report: a malformed row in the oracle file");
        (void)matched;
        OracleColumns *o = &oracle[model_of(label)];
        push(&o->d, d);
        push(&o->sigma2, sigma2);
        push(&o->s, s);
        push(&o->v, v);
        push(lowest == cop ? &o->v_own : &o->v_other, v);
        if (!end) break;
        line = end + 1;
    }
    free(text);
}

/* One bootstrap file's columns by model, rows whose sigma2_star is NaN left
   out, with the resample counts. */
typedef struct {
    Column sigma2, sigma2_star, ratio, v, v_star;
    int n_rows, n_bootstrap, fewest_usable, all_usable;
    long qvarma_converged;
    double qvarma_iterations;
} BootstrapColumns;

static int read_bootstrap(const char *path, BootstrapColumns *boot) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    for (int model = 0; model < N_MODELS; model++) {
        memset(&boot[model], 0, sizeof boot[model]);
        boot[model].fewest_usable = 1 << 30;
    }
    char line[1024];
    int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        char label[32];
        int cop, run, n_bootstrap, block, usable, converged, lowest;
        double iterations, d, sigma2, sigma2_star, s, v, s_star, v_star;
        int matched = sscanf(line, "%31[^,],%d,%d,%d,%d,%d,%d,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf", label, &cop, &run,
                             &n_bootstrap, &block, &usable, &converged, &iterations, &lowest, &d, &sigma2,
                             &sigma2_star, &s, &v, &s_star, &v_star);
        assert(matched == 16 && "fit_metric_report: a malformed row in a bootstrap file");
        (void)matched;
        BootstrapColumns *b = &boot[model_of(label)];
        b->n_rows++;
        b->n_bootstrap = n_bootstrap;
        if (usable < b->fewest_usable) b->fewest_usable = usable;
        b->all_usable += usable == n_bootstrap;
        b->qvarma_converged += converged;
        b->qvarma_iterations += iterations;
        if (MISNAN(sigma2_star)) continue;
        push(&b->sigma2, sigma2);
        push(&b->sigma2_star, sigma2_star);
        push(&b->ratio, sigma2_star / sigma2);
        push(&b->v, v);
        push(&b->v_star, v_star);
    }
    fclose(f);
    return 1;
}

static void free_bootstrap(BootstrapColumns *boot) {
    for (int model = 0; model < N_MODELS; model++) {
        free(boot[model].sigma2.value);
        free(boot[model].sigma2_star.value);
        free(boot[model].ratio.value);
        free(boot[model].v.value);
        free(boot[model].v_star.value);
    }
}

static void oracle_table(Table *t, const OracleColumns *oracle) {
    static const char *const names[] = { "model", "benchmarks", "median d", "median sigma2", "s 10%", "s 50%",
                                         "s 90%", "v 10%", "v 50%", "v 90%", "mean v" };
    set_header(t, 11, names);
    t->n_rows = N_MODELS;
    for (int m = 0; m < N_MODELS; m++) {
        const OracleColumns *o = &oracle[m];
        char (*row)[CELL] = t->cell[m];
        snprintf(row[0], CELL, "%s", model_name[m]);
        snprintf(row[1], CELL, "%d", o->v.n);
        snprintf(row[2], CELL, "%.4g", median(&o->d));
        snprintf(row[3], CELL, "%.4g", median(&o->sigma2));
        snprintf(row[4], CELL, "%.3f", quantile(&o->s, 0.1));
        snprintf(row[5], CELL, "%.3f", quantile(&o->s, 0.5));
        snprintf(row[6], CELL, "%.3f", quantile(&o->s, 0.9));
        snprintf(row[7], CELL, "%.3f", quantile(&o->v, 0.1));
        snprintf(row[8], CELL, "%.3f", quantile(&o->v, 0.5));
        snprintf(row[9], CELL, "%.3f", quantile(&o->v, 0.9));
        snprintf(row[10], CELL, "%.3f", (double)stats_mean(as_mat(&o->v)));
    }
}

static void winner_table(Table *t, const OracleColumns *oracle) {
    static const char *const names[] = { "model", "winner is the true configuration", "median v there",
                                         "winner is another", "median v there" };
    set_header(t, 5, names);
    t->n_rows = N_MODELS;
    for (int m = 0; m < N_MODELS; m++) {
        const OracleColumns *o = &oracle[m];
        char (*row)[CELL] = t->cell[m];
        snprintf(row[0], CELL, "%s", model_name[m]);
        snprintf(row[1], CELL, "%d", o->v_own.n);
        snprintf(row[2], CELL, "%.3f", o->v_own.n ? median(&o->v_own) : NAN);
        snprintf(row[3], CELL, "%d", o->v_other.n);
        snprintf(row[4], CELL, "%.3f", o->v_other.n ? median(&o->v_other) : NAN);
    }
}

static void bootstrap_table(Table *t, const BootstrapColumns *boot) {
    static const char *const names[] = { "model", "benchmarks", "median sigma2", "median sigma2_star", "ratio 10%",
                                         "ratio 50%", "ratio 90%", "rank correlation sigma2", "median v",
                                         "median v_star", "rank correlation v" };
    set_header(t, 11, names);
    t->n_rows = N_MODELS;
    for (int m = 0; m < N_MODELS; m++) {
        const BootstrapColumns *b = &boot[m];
        char (*row)[CELL] = t->cell[m];
        snprintf(row[0], CELL, "%s", model_name[m]);
        snprintf(row[1], CELL, "%d", b->ratio.n);
        snprintf(row[2], CELL, "%.4g", median(&b->sigma2));
        snprintf(row[3], CELL, "%.4g", median(&b->sigma2_star));
        snprintf(row[4], CELL, "%.4f", quantile(&b->ratio, 0.1));
        snprintf(row[5], CELL, "%.4f", quantile(&b->ratio, 0.5));
        snprintf(row[6], CELL, "%.4f", quantile(&b->ratio, 0.9));
        snprintf(row[7], CELL, "%.3f", spearman(&b->sigma2, &b->sigma2_star));
        snprintf(row[8], CELL, "%.3f", median(&b->v));
        snprintf(row[9], CELL, "%.3f", median(&b->v_star));
        snprintf(row[10], CELL, "%.3f", spearman(&b->v, &b->v_star));
    }
}

static void resample_table(Table *t, const BootstrapColumns *boot) {
    static const char *const names[] = { "model", "benchmarks", "resamples each", "fewest usable",
                                         "benchmarks with every resample usable" };
    set_header(t, 5, names);
    t->n_rows = N_MODELS;
    for (int m = 0; m < N_MODELS; m++) {
        const BootstrapColumns *b = &boot[m];
        char (*row)[CELL] = t->cell[m];
        snprintf(row[0], CELL, "%s", model_name[m]);
        snprintf(row[1], CELL, "%d", b->n_rows);
        snprintf(row[2], CELL, "%d", b->n_bootstrap);
        snprintf(row[3], CELL, "%d", b->fewest_usable);
        snprintf(row[4], CELL, "%d", b->all_usable);
    }
    const BootstrapColumns *q = &boot[MODEL_QVARMA];
    long fits = (long)q->n_rows * q->n_bootstrap;
    snprintf(t->note, sizeof t->note, "t-QVARMA refits: %ld of %ld converged (%.1f%%), %.0f iterations each on average",
             q->qvarma_converged, fits, 100.0 * q->qvarma_converged / fits, q->qvarma_iterations / q->n_rows);
}

/* One row per block-length file, in order of block length. */
static int block_length_table(Table *t) {
    static const char *const names[] = { "block length", "resamples", "configurations", "qvarma ratio 50%",
                                         "linear ratio 50%", "state 1 ratio 50%", "state 2 ratio 50%",
                                         "both states ratio 50%" };
    set_header(t, 8, names);
    t->n_rows = 0;
    int block[MAX_ROWS];
    DIR *dir = opendir(OUT_DIR);
    assert(dir && "fit_metric_report: cannot open montecarlo/out");
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL && t->n_rows < MAX_ROWS) {
        int run, n_bootstrap, length, iterations, first, consumed = 0;
        if (sscanf(entry->d_name, "fit_metric_bootstrap_run%d_B%d_block%d_iterations%d_first%d.csv%n", &run,
                   &n_bootstrap, &length, &iterations, &first, &consumed) != 5
            || entry->d_name[consumed] != 0 || run != 0) continue;
        char path[600];
        snprintf(path, sizeof path, "%s/%s", OUT_DIR, entry->d_name);
        BootstrapColumns boot[N_MODELS];
        read_bootstrap(path, boot);
        int at = t->n_rows++;
        block[at] = length;
        snprintf(t->cell[at][0], CELL, "%d", length);
        snprintf(t->cell[at][1], CELL, "%d", n_bootstrap);
        snprintf(t->cell[at][2], CELL, "%d", first);
        for (int m = 0; m < N_MODELS; m++) snprintf(t->cell[at][3 + m], CELL, "%.4f", median(&boot[m].ratio));
        free_bootstrap(boot);
    }
    closedir(dir);
    for (int i = 1; i < t->n_rows; i++)
        for (int j = i; j > 0 && block[j] < block[j - 1]; j--) {
            int swap = block[j]; block[j] = block[j - 1]; block[j - 1] = swap;
            char row[MAX_COLUMNS][CELL];
            memcpy(row, t->cell[j], sizeof row);
            memcpy(t->cell[j], t->cell[j - 1], sizeof row);
            memcpy(t->cell[j - 1], row, sizeof row);
        }
    return t->n_rows;
}

static void write_text(FILE *f, const Table *t) {
    int width[MAX_COLUMNS];
    for (int c = 0; c < t->n_columns; c++) {
        width[c] = (int)strlen(t->header[c]);
        for (int r = 0; r < t->n_rows; r++)
            if ((int)strlen(t->cell[r][c]) > width[c]) width[c] = (int)strlen(t->cell[r][c]);
    }
    fprintf(f, "%s\n", t->title);
    for (int c = 0; c < t->n_columns; c++) fprintf(f, c ? "  %*s" : "%*s", width[c], t->header[c]);
    fprintf(f, "\n");
    for (int r = 0; r < t->n_rows; r++) {
        for (int c = 0; c < t->n_columns; c++) fprintf(f, c ? "  %*s" : "%*s", width[c], t->cell[r][c]);
        fprintf(f, "\n");
    }
    if (t->note[0]) fprintf(f, "%s\n", t->note);
    fprintf(f, "\n");
}

/* text with _ and % escaped for LaTeX */
static void latex_text(FILE *f, const char *text) {
    for (; *text; text++) {
        if (*text == '_' || *text == '%') fputc('\\', f);
        fputc(*text, f);
    }
}

static void write_latex(FILE *f, const Table *t) {
    fprintf(f, "\\begin{table}[htbp]\n\\centering\n\\caption{");
    latex_text(f, t->title);
    fprintf(f, "}\n\\label{tab:fit_metric_%s}\n\\begin{tabular}{l", t->label);
    for (int c = 1; c < t->n_columns; c++) fputc('r', f);
    fprintf(f, "}\n\\hline\n");
    for (int c = 0; c < t->n_columns; c++) {
        if (c) fprintf(f, " & ");
        latex_text(f, t->header[c]);
    }
    fprintf(f, " \\\\\n\\hline\n");
    for (int r = 0; r < t->n_rows; r++) {
        for (int c = 0; c < t->n_columns; c++) {
            if (c) fprintf(f, " & ");
            latex_text(f, t->cell[r][c]);
        }
        fprintf(f, " \\\\\n");
    }
    fprintf(f, "\\hline\n\\end{tabular}\n");
    if (t->note[0]) {
        fprintf(f, "\\par\\smallskip ");
        latex_text(f, t->note);
        fprintf(f, "\n");
    }
    fprintf(f, "\\end{table}\n\n");
}

int main(void) {
    static OracleColumns oracle[N_MODELS];
    read_oracle(oracle);

    static Table table[5];
    int n_tables = 0;
    table[n_tables] = (Table){ .label = "oracle", .title = "Oracle fit metric over every benchmark of the grid" };
    oracle_table(&table[n_tables++], oracle);
    table[n_tables] = (Table){ .label = "winner",
                               .title = "Oracle v by whether the configuration with the smallest mean loss is the true one" };
    winner_table(&table[n_tables++], oracle);

    BootstrapColumns boot[N_MODELS];
    if (read_bootstrap(BOOTSTRAP_PATH, boot)) {
        table[n_tables] = (Table){ .label = "bootstrap",
                                   .title = "Bootstrap against oracle sigma2 on run 0 of every configuration" };
        bootstrap_table(&table[n_tables++], boot);
        table[n_tables] = (Table){ .label = "resamples", .title = "Bootstrap resamples with a usable response" };
        resample_table(&table[n_tables++], boot);
        free_bootstrap(boot);
    }
    table[n_tables] = (Table){ .label = "block_length",
                               .title = "Median ratio of bootstrap to oracle sigma2 by block length" };
    if (block_length_table(&table[n_tables])) n_tables++;

    FILE *report = fopen(REPORT_PATH, "w"), *latex = fopen(LATEX_PATH, "w");
    assert(report && latex && "fit_metric_report: cannot open the outputs");
    for (int i = 0; i < n_tables; i++) {
        write_text(report, &table[i]);
        write_latex(latex, &table[i]);
    }
    fclose(report);
    fclose(latex);

    for (int m = 0; m < N_MODELS; m++) {
        free(oracle[m].d.value);
        free(oracle[m].sigma2.value);
        free(oracle[m].s.value);
        free(oracle[m].v.value);
        free(oracle[m].v_own.value);
        free(oracle[m].v_other.value);
    }
    return 0;
}
