/*
Whether the tails studies/abm_system_tail_origin.c measures are the model's or
this project's.

That study reads every replicate under dataset/abm_system and reports that the
simulated series are thinner tailed than the US ones: over 400 periods the
median configuration has excess kurtosis near 0.61 in GDP growth against 2.87
in the US series, and 0.17 in the first difference of the interest rate against
16.09. The archives it reads were written by this project's build of the DSK
simulator, which is upstream's code with the rewrites docs/DSK_MODEL_CHANGES.md
records. If any of those rewrites thinned a tail, the study would be reporting
on the rewrites rather than on the model.

tests/dsk_dataset_reproduction.c already puts both builds on the same
configuration and seed and requires the results file, the error log and the
stored series to agree exactly. It samples configurations evenly across the
design and takes replications 0 to 4 of each. Tails are not spread evenly over
replicates: they are carried by the few runs that went furthest from the
ordinary path, and an even sample of low replication numbers is not where those
are. The rewrites that could act on those runs alone rather than on all of them
are the long-horizon ones in docs/DSK_LONG_HORIZON.md, which change nothing
until a threshold is crossed, so the extreme replicates are both where a
difference would first appear and where the even sample does not look.

So the sample here is chosen by how heavy a replicate's tails are, on two
criteria, because a replicate can carry a large fourth moment with no single
large observation and the other way round:

    heaviest kurtosis   the largest excess kurtosis m4 / m2^2 - 3 over the six
                        series, with m_k the k-th central sample moment divided
                        by n, the definition the study uses
    largest deviation   the largest |x_t - mean| / sd over the six series and
                        all stored periods

The six series are the five stored ones and the first difference of the
interest rate, again the study's own set: the interest rate is a persistent
level, so its own moments describe where the level wandered rather than how
large its shocks were.

Both criteria are computed for every replicate of every configuration, a
read-only sweep of all the archives, and the top N by each are taken, N = 12
unless given as the first argument. The two sets overlap, so the number of runs
compared is between N and 2N.

Each selected replicate is then run in both builds at the configuration it
belongs to, 600 periods, seed = replicate + 1 as the experiment was run, and
four things are required of it:

    the 83-column results file equals upstream's byte for byte
    the error log equals upstream's, so a run that ends early ends early in
        both for the same stated reason
    the five series rebuilt from the run equal the stored archive exactly, not
        to a tolerance: the stored values came through the same arithmetic from
        the same printed file, so anything but equality means the model now
        computes something else
    no long-horizon log was written, since none of those mechanisms may fire
        within 600 periods

A pass says the tail-carrying replicates are upstream's own runs, so the
thinness the study reports is the model's. It does not say the study's reading
of those runs is right, which is a separate question, and it covers only the
replicates selected here.

The upstream build is unoptimised and takes about 35 seconds a run, so the
comparisons are spread over the cores with OpenMP, as is the sweep.

    ./bin/dsk_tail_replicate_reproduction [SELECTED]

Needs `make model`, `make model-upstream` and dataset/abm_system. Nothing is
written anywhere but a scratch directory of this process's own and
out/dsk_tail_replicate_reproduction.txt; the archives are only read.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>
#include <omp.h>

#include "applications/abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./json.h>
#include "tests/dsk_scratch.h"

#define MODEL "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define DESIGN "dataset/abm_system_design.csv"
#define DATASET "dataset/abm_system"
#define REPORT "out/dsk_tail_replicate_reproduction.txt"
#define RUN_NAME "r"
#define N_STEPS 600
#define N_MODEL_COLUMNS 83
#define DEFAULT_SELECTED 12

/* The five stored series and the first difference of the interest rate. */
#define N_SERIES (ABM_SYSTEM_K + 1)

/* The model's own columns the five series are built from, one per row and in
   the same order, as applications/abm_system_simulate.c reads them. */
static const int model_column[ABM_SYSTEM_K] = {
    [LEVEL_GDP] = 2,
    [LEVEL_ENERGY] = 8,
    [LEVEL_EMPLOYMENT] = 5,
    [LEVEL_PRICE] = 34,
    [LEVEL_INTEREST] = 50
};

/* What the sweep keeps about one replicate. The ranks are over every replicate
   of every configuration, 1 being the most extreme on that criterion. */
typedef struct {
    int configuration;
    int replicate;
    double heaviest_kurtosis;
    double largest_deviation;
    int kurtosis_rank;
    int deviation_rank;
} Tail;

typedef struct {
    Tail tail;
    int ran;                 /* both builds produced a results file */
    int results_differ;      /* against the upstream build, byte for byte */
    int errors_differ;
    int long_horizon_logs;   /* a mechanism that must not fire in 600 periods */
    Mat series;              /* the five series this run produces */
} Comparison;

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0)
        assert(errno == EEXIST && "dsk_tail_replicate_reproduction: mkdir failed");
}

static int file_exists(const char *path) { return access(path, F_OK) == 0; }

static void link_build(const char *dir, const char *binary) {
    make_directory(dir);
    char link[512];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(binary, link) == 0 &&
           "dsk_tail_replicate_reproduction: cannot link a build into place");
}

/* Excess kurtosis m4 / m2^2 - 3 and the largest absolute deviation from the
   mean in standard deviations, of one series of n values. Both are NAN when
   the series has fewer than two values or no spread, where neither ratio is
   defined. */
static void series_tails(const double *x, int n, double *excess_kurtosis,
                         double *largest_deviation) {
    *excess_kurtosis = NAN;
    *largest_deviation = NAN;
    if (n < 2) return;

    double mean = 0;
    for (int t = 0; t < n; t++) mean += x[t];
    mean /= (double)n;

    double m2 = 0, m4 = 0, widest = 0;
    for (int t = 0; t < n; t++) {
        double d = x[t] - mean, d2 = d * d;
        m2 += d2;
        m4 += d2 * d2;
        if (fabs(d) > widest) widest = fabs(d);
    }
    m2 /= (double)n;
    m4 /= (double)n;
    if (!(m2 > 0)) return;

    *excess_kurtosis = m4 / (m2 * m2) - 3.0;
    *largest_deviation = widest / sqrt(m2);
}

/* The two criteria for one ABM_SYSTEM_K x periods block, each the largest of
   the six series' own value. scratch holds at least periods doubles. A series
   with no spread contributes nothing, and a block where none of the six has
   any spread comes back at minus infinity, which ranks it last. */
static void block_tails(Mat y, double *scratch, double *heaviest_kurtosis,
                        double *largest_deviation) {
    int periods = y.c;
    *heaviest_kurtosis = -INFINITY;
    *largest_deviation = -INFINITY;

    for (int k = 0; k < N_SERIES; k++) {
        int n;
        if (k < ABM_SYSTEM_K) {
            n = periods;
            for (int t = 0; t < periods; t++) scratch[t] = (double)AT(y, k, t);
        } else {
            n = periods - 1;
            for (int t = 1; t < periods; t++)
                scratch[t - 1] = (double)AT(y, ROW_INTEREST_RATE, t) -
                                 (double)AT(y, ROW_INTEREST_RATE, t - 1);
        }

        double kurtosis, deviation;
        series_tails(scratch, n, &kurtosis, &deviation);
        if (!isnan(kurtosis) && kurtosis > *heaviest_kurtosis) *heaviest_kurtosis = kurtosis;
        if (!isnan(deviation) && deviation > *largest_deviation) *largest_deviation = deviation;
    }
}

/* The configurations stored under DATASET, as the numbers in their names,
   ascending. Caller must free. */
static int *list_configurations(int *count) {
    DIR *handle = opendir(DATASET);
    assert(handle && "dsk_tail_replicate_reproduction: cannot open dataset/abm_system");

    int *cop = NULL, n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int here;
        if (sscanf(entry->d_name, "cop_%d", &here) != 1) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 1024;
            cop = (int*)realloc(cop, (size_t)cap * sizeof(int));
            assert(cop && "dsk_tail_replicate_reproduction: out of memory");
        }
        cop[n++] = here;
    }
    closedir(handle);

    for (int i = 1; i < n; i++) {
        int current = cop[i], j = i;
        while (j > 0 && cop[j - 1] > current) { cop[j] = cop[j - 1]; j--; }
        cop[j] = current;
    }
    *count = n;
    return cop;
}

/* The largest batch index an archive exists for, or -1. Batches are walked by
   index rather than in the order readdir returns them, so a batch that failed
   entirely is simply absent rather than a hole to stop at. */
static int highest_batch(const char *dir) {
    DIR *handle = opendir(dir);
    assert(handle && "dsk_tail_replicate_reproduction: cannot open a configuration directory");

    int highest = -1;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int batch;
        if (sscanf(entry->d_name, "batch_%d.npz", &batch) == 1 && batch > highest) highest = batch;
    }
    closedir(handle);
    return highest;
}

static int by_kurtosis(const void *a, const void *b) {
    double left = ((const Tail*)a)->heaviest_kurtosis, right = ((const Tail*)b)->heaviest_kurtosis;
    if (left > right) return -1;
    if (left < right) return 1;
    return 0;
}

static int by_deviation(const void *a, const void *b) {
    double left = ((const Tail*)a)->largest_deviation, right = ((const Tail*)b)->largest_deviation;
    if (left > right) return -1;
    if (left < right) return 1;
    return 0;
}

static void write_inputs(const char *base_json, const char *path, const Mat *parameter,
                         int design_row) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 &&
           "dsk_tail_replicate_reproduction: no params block");

    JsonValue *block = json_array_get(params, 0);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        json_object_set(block, abm_system_parameter_names()[p],
                        json_number((double)AT(parameter[p], design_row, 0)));
    json_object_set(block, "T", json_number(N_STEPS));

    json_write_file(inputs, path);
    json_free(inputs);
}

static int run_build(const char *dir, const char *json_path, int seed) {
    char command[8192];
    snprintf(command, sizeof command,
             "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, json_path, RUN_NAME, seed);
    return system(command);
}

/* Whole file, or NULL. */
static char *read_file(const char *path, long *size) {
    FILE *f = fopen(path, "rb");
    if (!f) { *size = -1; return NULL; }
    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *text = malloc((size_t)*size + 1);
    assert(text && "dsk_tail_replicate_reproduction: out of memory");
    if (*size && fread(text, 1, (size_t)*size, f) != (size_t)*size) {
        free(text); fclose(f); *size = -1; return NULL;
    }
    text[*size] = '\0';
    fclose(f);
    return text;
}

/* The five series, rebuilt from a results file the way the driver rebuilds
   them. Returns an empty matrix if the file is not what it should be. */
static Mat series_from_results(const char *path) {
    Mat empty = {0};
    FILE *f = fopen(path, "r");
    if (!f) return empty;

    Mat levels = mat_new(ABM_SYSTEM_K, N_STEPS);
    for (int t = 0; t < N_STEPS; t++) {
        for (int c = 1; c <= N_MODEL_COLUMNS; c++) {
            double value;
            if (fscanf(f, "%lf", &value) != 1) { fclose(f); mat_free(levels); return empty; }
            for (int series = 0; series < ABM_SYSTEM_K; series++)
                if (model_column[series] == c) AT(levels, series, t) = (mreal)value;
        }
    }
    fclose(f);

    Mat y = abm_system_transform(levels, ABM_SYSTEM_BURN_IN);
    mat_free(levels);
    return y;
}

static int count_long_horizon_logs(const char *dir, int seed) {
    static const char *const suffix[] = {"_redenominations", "_machine_lots", "_good_units"};
    int found = 0;
    for (int k = 0; k < 3; k++) {
        char path[1024];
        snprintf(path, sizeof path, "%s/output/results_%s_%d%s.txt", dir, RUN_NAME, seed, suffix[k]);
        found += file_exists(path);
    }
    return found;
}

int main(int argc, char **argv) {
    int selected_per_criterion = argc > 1 ? atoi(argv[1]) : DEFAULT_SELECTED;
    assert(selected_per_criterion >= 1 &&
           "dsk_tail_replicate_reproduction: nothing to select");

    char model[PATH_MAX], upstream[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "dsk_tail_replicate_reproduction: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(UPSTREAM, upstream) &&
           "dsk_tail_replicate_reproduction: bin/dsk_SFC_upstream is not built - run make model-upstream");
    assert(realpath(INPUTS, base_json) &&
           "dsk_tail_replicate_reproduction: the parameter file is missing");

    DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
    Mat parameter[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        parameter[p] = df_col_numeric(&design, abm_system_parameter_names()[p]);

    int n_configurations;
    int *cop = list_configurations(&n_configurations);
    assert(n_configurations > 0 && "dsk_tail_replicate_reproduction: no configurations stored");

    /* The sweep. One configuration per thread, each archive read once, and the
       two criteria of every replicate it holds appended to one list. */
    Tail *tail = NULL;
    long n_tails = 0, capacity = 0;
    const double sweep_started = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_configurations; i++) {
        char dir[512];
        snprintf(dir, sizeof dir, "%s/cop_%04d", DATASET, cop[i]);

        Tail *local = NULL;
        int n_local = 0, cap_local = 0;
        double *scratch = NULL;
        int scratch_periods = 0;

        int highest = highest_batch(dir);
        for (int batch = 0; batch <= highest; batch++) {
            Mat block[ABM_SYSTEM_BATCH];
            int replicate[ABM_SYSTEM_BATCH];
            int count = abm_system_read_batch(dir, batch, block, replicate);

            for (int b = 0; b < count; b++) {
                if (block[b].c > scratch_periods) {
                    scratch_periods = block[b].c;
                    scratch = (double*)realloc(scratch, (size_t)scratch_periods * sizeof(double));
                    assert(scratch && "dsk_tail_replicate_reproduction: out of memory");
                }
                if (n_local == cap_local) {
                    cap_local = cap_local ? 2 * cap_local : 128;
                    local = (Tail*)realloc(local, (size_t)cap_local * sizeof(Tail));
                    assert(local && "dsk_tail_replicate_reproduction: out of memory");
                }

                Tail here = {0};
                here.configuration = cop[i];
                here.replicate = replicate[b];
                block_tails(block[b], scratch, &here.heaviest_kurtosis, &here.largest_deviation);
                local[n_local++] = here;
                mat_free(block[b]);
            }
        }
        free(scratch);

        #pragma omp critical(tail_append)
        if (n_local) {
            if (n_tails + n_local > capacity) {
                while (n_tails + n_local > capacity) capacity = capacity ? 2 * capacity : 4096;
                tail = (Tail*)realloc(tail, (size_t)capacity * sizeof(Tail));
                assert(tail && "dsk_tail_replicate_reproduction: out of memory");
            }
            memcpy(tail + n_tails, local, (size_t)n_local * sizeof(Tail));
            n_tails += n_local;
        }
        free(local);
    }

    const double sweep_elapsed = omp_get_wtime() - sweep_started;
    assert(n_tails > 0 && "dsk_tail_replicate_reproduction: the archives hold no replicates");

    /* Ranked once on each criterion, so both ranks travel with the record and
       the selection is one pass over the list in its final order. */
    qsort(tail, (size_t)n_tails, sizeof(Tail), by_kurtosis);
    for (long k = 0; k < n_tails; k++) tail[k].kurtosis_rank = (int)k + 1;
    qsort(tail, (size_t)n_tails, sizeof(Tail), by_deviation);
    for (long k = 0; k < n_tails; k++) tail[k].deviation_rank = (int)k + 1;

    int n_chosen = 0;
    for (long k = 0; k < n_tails; k++)
        if (tail[k].kurtosis_rank <= selected_per_criterion ||
            tail[k].deviation_rank <= selected_per_criterion) n_chosen++;

    Comparison *chosen = (Comparison*)calloc((size_t)n_chosen, sizeof(Comparison));
    assert(chosen && "dsk_tail_replicate_reproduction: out of memory");
    int filled = 0;
    for (long k = 0; k < n_tails; k++)
        if (tail[k].kurtosis_rank <= selected_per_criterion ||
            tail[k].deviation_rank <= selected_per_criterion)
            chosen[filled++].tail = tail[k];

    /* Where the selection sits in the two distributions, for the report. The
       list is ordered by deviation now and by kurtosis before that, so the
       quantiles are read off the sorted order rather than recomputed. */
    double deviation_quantile[3], kurtosis_quantile[3];
    static const double probability[3] = {0.5, 0.99, 1.0};
    for (int q = 0; q < 3; q++) {
        long index = (long)((1.0 - probability[q]) * (double)(n_tails - 1));
        deviation_quantile[q] = tail[index].largest_deviation;
    }
    qsort(tail, (size_t)n_tails, sizeof(Tail), by_kurtosis);
    for (int q = 0; q < 3; q++) {
        long index = (long)((1.0 - probability[q]) * (double)(n_tails - 1));
        kurtosis_quantile[q] = tail[index].heaviest_kurtosis;
    }

    const char *scratch_root = getenv("TMPDIR");
    if (!scratch_root) scratch_root = "/tmp";

    /* Short directory names: upstream builds every output filename in a
       64-byte buffer, which a long path overruns. */
    char root[128];
    snprintf(root, sizeof root, "%s/t%d", scratch_root, (int)getpid());
    make_directory(root);

    const int threads = omp_get_max_threads();
    for (int k = 0; k < threads; k++) {
        char dir[192];
        snprintf(dir, sizeof dir, "%s/m%d", root, k);
        link_build(dir, model);
        snprintf(dir, sizeof dir, "%s/u%d", root, k);
        link_build(dir, upstream);
    }

    /* One parameter file per configuration the selection reaches, written
       once even when several of its replicates were chosen. */
    for (int k = 0; k < n_chosen; k++) {
        int configuration = chosen[k].tail.configuration;
        assert(configuration >= 1 && configuration <= design.r &&
               "dsk_tail_replicate_reproduction: a stored configuration is not in the design");
        char json_path[256];
        snprintf(json_path, sizeof json_path, "%s/c%d.json", root, configuration);
        if (!file_exists(json_path))
            write_inputs(base_json, json_path, parameter, configuration - 1);
    }

    const double run_started = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < n_chosen; k++) {
        const int thread = omp_get_thread_num();
        const int seed = chosen[k].tail.replicate + 1;

        char json_path[256], mine[192], theirs[192];
        snprintf(json_path, sizeof json_path, "%s/c%d.json", root, chosen[k].tail.configuration);
        snprintf(mine, sizeof mine, "%s/m%d", root, thread);
        snprintf(theirs, sizeof theirs, "%s/u%d", root, thread);

        int my_status = run_build(mine, json_path, seed);
        int their_status = run_build(theirs, json_path, seed);

        char my_results[320], their_results[320], my_errors[320], their_errors[320];
        snprintf(my_results, sizeof my_results, "%s/output/results_%s_%d.txt", mine, RUN_NAME, seed);
        snprintf(their_results, sizeof their_results, "%s/output/results_%s_%d.txt", theirs, RUN_NAME, seed);
        snprintf(my_errors, sizeof my_errors, "%s/output/errors/Errors_%s_%d.txt", mine, RUN_NAME, seed);
        snprintf(their_errors, sizeof their_errors, "%s/output/errors/Errors_%s_%d.txt", theirs, RUN_NAME, seed);

        long my_size, their_size;
        char *my_text = read_file(my_results, &my_size);
        char *their_text = read_file(their_results, &their_size);

        chosen[k].ran = my_status == 0 && their_status == 0 && my_text && their_text;
        chosen[k].results_differ = !chosen[k].ran || my_size != their_size ||
                                   memcmp(my_text, their_text, (size_t)my_size) != 0;

        long my_error_size, their_error_size;
        char *my_error_text = read_file(my_errors, &my_error_size);
        char *their_error_text = read_file(their_errors, &their_error_size);
        chosen[k].errors_differ = my_error_size != their_error_size ||
                                  (my_error_size > 0 && memcmp(my_error_text, their_error_text,
                                                               (size_t)my_error_size) != 0);
        free(my_error_text);
        free(their_error_text);

        if (chosen[k].ran) chosen[k].series = series_from_results(my_results);
        chosen[k].long_horizon_logs = count_long_horizon_logs(mine, seed);

        free(my_text);
        free(their_text);
        unlink(my_results);
        unlink(their_results);
        unlink(my_errors);
        unlink(their_errors);
    }

    const double run_elapsed = omp_get_wtime() - run_started;

    int compared = 0, missing = 0, differing_runs = 0, failed_runs = 0, fired = 0;
    int results_differing = 0, errors_differing = 0;
    long differing_values = 0;
    double worst_gap = 0;
    int worst_configuration = 0, worst_replicate = 0;

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_tail_replicate_reproduction: cannot open the report path");
    fprintf(report, "The replicates that carry the tails, against the authors' build and against\n"
                    "the archives studies/abm_system_tail_origin.c reads.\n\n");
    fprintf(report, "%d configurations, %ld replicates swept in %.1f minutes.\n",
            n_configurations, n_tails, sweep_elapsed / 60.0);
    fprintf(report, "Top %d by each criterion, %d distinct runs, %d threads, %.1f minutes of\n"
                    "wall clock on the runs.\n",
            selected_per_criterion, n_chosen, threads, run_elapsed / 60.0);
    fprintf(report, "Seeds are the replicate number plus one, as the experiment was run.\n\n");

    fprintf(report, "the two criteria over everything swept, context rather than a verdict\n");
    fprintf(report, "  criterion           median       p99       max\n");
    fprintf(report, "  heaviest kurtosis %9.4g %9.4g %9.4g\n",
            kurtosis_quantile[0], kurtosis_quantile[1], kurtosis_quantile[2]);
    fprintf(report, "  largest deviation %9.4g %9.4g %9.4g\n\n",
            deviation_quantile[0], deviation_quantile[1], deviation_quantile[2]);

    fprintf(report, "the runs compared, with their rank on each criterion\n");
    fprintf(report, "  configuration replicate seed  kurtosis  rank  deviation  rank  verdict\n");

    for (int k = 0; k < n_chosen; k++) {
        const Tail *here = &chosen[k].tail;
        char verdict[256];
        verdict[0] = '\0';

        if (!chosen[k].ran) {
            failed_runs++;
            snprintf(verdict, sizeof verdict, "a build did not produce a results file  FAILED");
        } else {
            if (chosen[k].results_differ) {
                results_differing++;
                strcat(verdict, "results file differs  FAILED  ");
            }
            if (chosen[k].errors_differ) {
                errors_differing++;
                strcat(verdict, "error log differs  FAILED  ");
            }
            fired += chosen[k].long_horizon_logs;
            if (chosen[k].long_horizon_logs)
                strcat(verdict, "a long-horizon mechanism fired  FAILED  ");

            char cop_dir[512];
            snprintf(cop_dir, sizeof cop_dir, "%s/cop_%04d", DATASET, here->configuration);
            char batch_path[640];
            abm_system_batch_path(batch_path, sizeof batch_path, cop_dir, here->replicate);

            if (!file_exists(batch_path)) {
                missing++;
                strcat(verdict, "the archive holding it is gone");
            } else {
                Mat stored = abm_system_read_replicate(cop_dir, here->replicate);
                Mat produced = chosen[k].series;
                if (produced.r != stored.r || produced.c != stored.c) {
                    differing_runs++;
                    char shape[128];
                    snprintf(shape, sizeof shape, "%d x %d produced against %d x %d stored  FAILED",
                             produced.r, produced.c, stored.r, stored.c);
                    strcat(verdict, shape);
                } else {
                    long differing_here = 0;
                    for (int series = 0; series < stored.r; series++) {
                        for (int t = 0; t < stored.c; t++) {
                            double a = (double)AT(stored, series, t);
                            double b = (double)AT(produced, series, t);
                            if (a != b) {
                                differing_here++;
                                double gap = a != 0 ? fabs(b - a) / fabs(a) : fabs(b);
                                if (gap > worst_gap) {
                                    worst_gap = gap;
                                    worst_configuration = here->configuration;
                                    worst_replicate = here->replicate;
                                }
                            }
                        }
                    }
                    if (differing_here) {
                        differing_runs++;
                        differing_values += differing_here;
                        char counted[128];
                        snprintf(counted, sizeof counted, "%ld of %d stored values differ  FAILED",
                                 differing_here, stored.r * stored.c);
                        strcat(verdict, counted);
                    } else {
                        strcat(verdict, "reproduced");
                    }
                    compared++;
                }
                mat_free(stored);
            }
            mat_free(chosen[k].series);
        }

        fprintf(report, "  cop_%04d      %9d %4d %9.4g %5d %10.4g %5d  %s\n",
                here->configuration, here->replicate, here->replicate + 1,
                here->heaviest_kurtosis, here->kurtosis_rank,
                here->largest_deviation, here->deviation_rank, verdict);
    }

    fprintf(report, "\nruns compared                   %d\n", n_chosen - failed_runs);
    fprintf(report, "results files against upstream  %d compared, %d differing\n",
            n_chosen - failed_runs, results_differing);
    fprintf(report, "error logs against upstream     %d compared, %d differing\n",
            n_chosen - failed_runs, errors_differing);
    fprintf(report, "runs against the archives       %d compared, %d differing\n",
            compared, differing_runs);
    if (differing_values)
        fprintf(report, "values differing                %ld, worst relative gap %.3e at cop_%04d replicate %d\n",
                differing_values, worst_gap, worst_configuration, worst_replicate);
    fprintf(report, "archives missing                %d\n", missing);
    fprintf(report, "long-horizon logs written       %d (must be 0)\n", fired);
    fprintf(report, "runs that did not complete      %d\n", failed_runs);

    int failed = failed_runs || results_differing || errors_differing || differing_runs ||
                 fired || !compared;
    fprintf(report, "\n%s\n", failed ? "FAILED" : "PASSED");
    fclose(report);

    printf("tail replicate reproduction: %ld replicates swept, %d runs compared, "
           "%d differing from upstream, %d differing from the archives, %.1f minutes\n",
           n_tails, n_chosen, results_differing, differing_runs,
           (sweep_elapsed + run_elapsed) / 60.0);
    printf("%s\n", failed ? "FAILED" : "PASSED, 0 failures");

    free(chosen);
    free(tail);
    free(cop);
    df_free(&design);
    dsk_scratch_finish(root, !failed);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
