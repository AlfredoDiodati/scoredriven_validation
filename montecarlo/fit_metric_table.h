#ifndef MONTECARLO_FIT_METRIC_TABLE_H
#define MONTECARLO_FIT_METRIC_TABLE_H

/*
The tables bin/fit_metric_pairwise writes: cells
held as text, written once as a plain-text table and once as a LaTeX table
labelled tab:fit_metric_<label>.
*/

#include "montecarlo/response_cache.h"

#define MAX_COLUMNS 12
#define MAX_ROWS 16
#define CELL 48

static const char *const model_name[N_MODELS] = { "qvarma", "linear", "state 1", "state 2", "both states" };

typedef struct {
    const char *label, *title;
    int n_columns, n_rows;
    char header[MAX_COLUMNS][CELL];
    char cell[MAX_ROWS][MAX_COLUMNS][CELL];
    char note[512];
} Table;

static void set_header(Table *t, int n_columns, const char *const *names) {
    t->n_columns = n_columns;
    for (int i = 0; i < n_columns; i++) snprintf(t->header[i], CELL, "%s", names[i]);
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

#endif /* MONTECARLO_FIT_METRIC_TABLE_H */
