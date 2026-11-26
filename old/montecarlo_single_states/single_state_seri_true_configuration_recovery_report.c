/*
What the rows of
old/montecarlo_single_states/single_state_seri_true_configuration_recovery.c,
state 1 and state 2 of the state-dependent local projection on their own,
under the Seri et al. confidence set and the corrected loss, say so far. Reads
its progress file, or its compressed result when the progress file is gone,
and writes
out/old/montecarlo_single_states/seri_true_configuration_recovery_report.txt:

1. Per model: benchmarks done, how many run indices that covers, how often the
   benchmark's configuration is in the set, alone in it, and first by mean
   loss, the mean and largest set size, and how often the set was decided by
   an accepted test.
2. Per model, per configuration: the share of its own benchmarks whose set
   holds it, summarised as how many configurations reach at least half, at
   least 90%, and none.

montecarlo/seri_true_configuration_recovery_report.c reports the t-QVARMA,
the linear projection and both states under the causal HP cycle at the
estimated lambda. Rows are read as written; a benchmark appearing twice is
counted twice.
*/

/* getline and fmemopen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#define PROGRESS_PATH "out/old/montecarlo_single_states/seri_true_configuration_recovery_progress.csv"
#define RESULT_PATH "out/old/montecarlo_single_states/seri_true_configuration_recovery.csv.gz"
#define REPORT_PATH "out/old/montecarlo_single_states/seri_true_configuration_recovery_report.txt"
#define MAX_COP 4096
#define MAX_RUN 1000

typedef struct {
    long rows, in_set, alone, first, decided, set_size_sum, largest_set;
    int highest_run;
    long own_benchmarks[MAX_COP], own_in_set[MAX_COP];
} ModelSummary;

static ModelSummary summary[N_MODELS];

/* The models reported: state 1 and state 2 on their own. */
static int reported(int model) {
    return model_is_single_state(model);
}

/* A row of a model not reported is skipped like a row that does not parse. */
static int model_of(const char *label) {
    for (int model = 0; model < N_MODELS; model++)
        if (strcmp(label, model_label[model]) == 0) return reported(model) ? model : -1;
    return -1;
}

/* One data row: the fields before the member list. Returns 0 on a row that
   does not parse. */
static int parse_row(const char *line, int *model, int *cop, int *run, int *in_set, int *rank, int *set_size,
                     int *decided) {
    char label[MODEL_LABEL_CAPACITY];
    double pvalue;
    if (sscanf(line, "%63[^,],%d,%d,%d,%d,%lf,%d,%d,", label, cop, run, in_set, rank, &pvalue, set_size, decided)
        != 8)
        return 0;
    *model = model_of(label);
    return *model >= 0 && *cop >= 0 && *cop < MAX_COP && *run >= 0 && *run < MAX_RUN;
}

static void add_seri_row(const char *line) {
    int model, cop, run, in_set, rank, set_size, decided;
    if (!parse_row(line, &model, &cop, &run, &in_set, &rank, &set_size, &decided)) return;
    ModelSummary *s = &summary[model];
    s->rows++;
    s->in_set += in_set;
    s->alone += in_set && set_size == 1;
    s->first += rank == 1;
    s->decided += decided;
    s->set_size_sum += set_size;
    if (set_size > s->largest_set) s->largest_set = set_size;
    if (run > s->highest_run) s->highest_run = run;
    s->own_benchmarks[cop]++;
    s->own_in_set[cop] += in_set;
}

static void read_lines(FILE *f, void (*add)(const char *)) {
    char *line = NULL;
    size_t capacity = 0;
    int first = 1;
    while (getline(&line, &capacity, f) > 0) {
        if (first) { first = 0; continue; }
        add(line);
    }
    free(line);
}

/* The seri rows from the progress file, or from the result when the progress
   file is gone. */
static const char *read_seri_rows(void) {
    FILE *f = fopen(PROGRESS_PATH, "r");
    if (f) {
        read_lines(f, add_seri_row);
        fclose(f);
        return PROGRESS_PATH;
    }
    f = fopen(RESULT_PATH, "rb");
    assert(f && "single_state_seri_true_configuration_recovery_report: neither the progress file nor the result exists");
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
    FILE *memory = fmemopen(text, text_length, "r");
    assert(memory);
    read_lines(memory, add_seri_row);
    fclose(memory);
    free(text);
    return RESULT_PATH;
}

static void write_per_model(FILE *out) {
    fprintf(out, "1. Per model, over the benchmarks done so far. Percentages are of the benchmarks done.\n\n");
    fprintf(out, "model, benchmarks, run indices reached, own configuration in the set, alone in the set, "
                 "first by mean loss, mean set size, largest set, decided by an accepted test\n");
    for (int model = 0; model < N_MODELS; model++) {
        const ModelSummary *s = &summary[model];
        if (!reported(model)) continue;
        if (!s->rows) {
            fprintf(out, "%s, 0 benchmarks so far\n", model_label[model]);
            continue;
        }
        double n = (double)s->rows;
        fprintf(out, "%s, %ld, 0 to %d, %ld (%.1f%%), %ld (%.1f%%), %ld (%.1f%%), %.2f, %ld, %ld (%.1f%%)\n",
                model_label[model], s->rows, s->highest_run, s->in_set, 100 * s->in_set / n, s->alone,
                100 * s->alone / n, s->first, 100 * s->first / n, s->set_size_sum / n, s->largest_set, s->decided,
                100 * s->decided / n);
    }
    fprintf(out, "\n");
}

static void write_per_configuration(FILE *out) {
    fprintf(out, "2. Per configuration: the share of its own benchmarks whose set holds it. Over the "
                 "configurations with at least one benchmark done.\n\n");
    fprintf(out, "model, configurations, at least half, at least 90%%, none, largest share\n");
    for (int model = 0; model < N_MODELS; model++) {
        const ModelSummary *s = &summary[model];
        if (!reported(model)) continue;
        int n_configurations = 0, half = 0, ninety = 0, none = 0;
        double largest = 0;
        for (int cop = 0; cop < MAX_COP; cop++) {
            if (!s->own_benchmarks[cop]) continue;
            n_configurations++;
            double share = (double)s->own_in_set[cop] / s->own_benchmarks[cop];
            half += share >= 0.5;
            ninety += share >= 0.9;
            none += s->own_in_set[cop] == 0;
            if (share > largest) largest = share;
        }
        if (!n_configurations) continue;
        fprintf(out, "%s, %d, %d, %d, %d, %.3f\n", model_label[model], n_configurations, half, ninety, none,
                largest);
    }
    fprintf(out, "\n");
}

int main(void) {
    for (int model = 0; model < N_MODELS; model++) summary[model].highest_run = -1;
    const char *source = read_seri_rows();

    mkdir("out", 0755);
    mkdir("out/old", 0755);
    mkdir("out/old/montecarlo_single_states", 0755);
    FILE *out = fopen(REPORT_PATH, "w");
    assert(out && "single_state_seri_true_configuration_recovery_report: cannot open the report");
    fprintf(out, "Rows of old/montecarlo_single_states/single_state_seri_true_configuration_recovery.c read from "
                 "%s.\n\n", source);
    write_per_model(out);
    write_per_configuration(out);
    fclose(out);

    return 0;
}
