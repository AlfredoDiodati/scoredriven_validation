/*
Whether the model as it stands still computes what the authors' model computes,
and still produces the archives the 1000 x 1000 experiment was run on.

Everything done to model/dsk_sfc is meant to leave a 600-period run exactly
where it was and to change what happens only past thresholds no 600-period run
reaches; docs/DSK_LONG_HORIZON.md lists each one and what it takes to fire. The
other equivalence tests check that against upstream at the model's own parameter
settings and at four points of the design. This one checks it where the
experiment actually ran, and against two references at once:

  - the authors' code, built unchanged by `make model-upstream`, run on the same
    configuration and seed. The 83-column results file and the error log must
    match byte for byte.
  - the archives in dataset/abm_system, which applications/abm_system_simulate.c
    wrote from the simulator as it was then. The five series are rebuilt from
    this run exactly as the driver rebuilds them, and every one of the 5 x 400
    numbers must equal the stored one exactly. Not to a tolerance: the stored
    values came through the same arithmetic from the same printed file, so
    anything but equality means the model now computes something else.

It also requires that no long-horizon mechanism fired: a run that redenominated,
changed the machine lot or changed the good's unit writes a log beside its
results file, and a 600-period run must write none.

The sample is configurations spread evenly across the design, and the first
replications of each. Replications carry the same seeds in every configuration -
seed = replication + 1 - which is how the experiment was run, so a replication
number means the same draw everywhere. Both builds are run on every sampled
pair, in parallel: the upstream build is unoptimised and takes about 35 seconds
a run, so the work is divided over the threads OpenMP gives it.

    ./bin/dsk_dataset_reproduction [CONFIGURATIONS [REPLICATIONS]]

Default 20 and 5, which is 100 pairs. Nothing is written anywhere but a scratch
directory of this process's own and out/dsk_dataset_reproduction.txt; the
archives are only read.

Needs `make model`, `make model-upstream`, and dataset/abm_system.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>
#include <omp.h>

#include "applications/abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./json.h>

#define MODEL "model/dsk_sfc/dsk_SFC"
#define UPSTREAM "bin/dsk_SFC_upstream"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define DESIGN "dataset/abm_system_design.csv"
#define DATASET "dataset/abm_system"
#define REPORT "out/dsk_dataset_reproduction.txt"
#define RUN_NAME "r"
#define N_STEPS 600
#define N_MODEL_COLUMNS 83
#define DEFAULT_CONFIGURATIONS 20
#define DEFAULT_REPLICATIONS 5

/* The model's own columns the five series are built from, as
   applications/abm_system_simulate.c reads them. */
static const int model_column[ABM_SYSTEM_K] = {
    [LEVEL_GDP] = 2,
    [LEVEL_ENERGY] = 8,
    [LEVEL_EMPLOYMENT] = 5,
    [LEVEL_PRICE] = 34,
    [LEVEL_INTEREST] = 50
};

typedef struct {
    int cop;
    int replicate;
    int ran;                 /* both builds produced a results file */
    int results_differ;      /* against the upstream build, byte for byte */
    int errors_differ;
    int long_horizon_logs;   /* a mechanism that must not fire in 600 periods */
    Mat series;              /* the five series this run produces */
} Pair;

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_dataset_reproduction: mkdir failed");
}

static int file_exists(const char *path) { return access(path, F_OK) == 0; }

static void link_build(const char *dir, const char *binary) {
    make_directory(dir);
    char link[512];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(binary, link) == 0 && "dsk_dataset_reproduction: cannot link a build into place");
}

static void write_inputs(const char *base_json, const char *path, const Mat *parameter, int design_row) {
    JsonValue *inputs = json_parse_file(base_json);
    JsonValue *params = json_object_get(inputs, "params");
    assert(params && json_array_len(params) >= 1 && "dsk_dataset_reproduction: no params block");

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
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
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
    assert(text && "dsk_dataset_reproduction: out of memory");
    if (*size && fread(text, 1, (size_t)*size, f) != (size_t)*size) { free(text); fclose(f); *size = -1; return NULL; }
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
    int configurations = argc > 1 ? atoi(argv[1]) : DEFAULT_CONFIGURATIONS;
    int replications = argc > 2 ? atoi(argv[2]) : DEFAULT_REPLICATIONS;
    assert(configurations >= 1 && replications >= 1 && "dsk_dataset_reproduction: nothing to compare");

    char model[PATH_MAX], upstream[PATH_MAX], base_json[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "dsk_dataset_reproduction: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(UPSTREAM, upstream) &&
           "dsk_dataset_reproduction: bin/dsk_SFC_upstream is not built - run make model-upstream");
    assert(realpath(INPUTS, base_json) && "dsk_dataset_reproduction: the parameter file is missing");

    DataFrame design = df_read_csv(DESIGN, csv_read_options_default());
    Mat parameter[ABM_SYSTEM_N_PARAMETERS];
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++)
        parameter[p] = df_col_numeric(&design, abm_system_parameter_names()[p]);

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    /* Short directory names: upstream builds every output filename in a
       64-byte buffer, which a long path overruns. */
    char root[128];
    snprintf(root, sizeof root, "%s/d%d", scratch, (int)getpid());
    make_directory(root);

    const int threads = omp_get_max_threads();
    for (int k = 0; k < threads; k++) {
        char dir[192];
        snprintf(dir, sizeof dir, "%s/m%d", root, k);
        link_build(dir, model);
        snprintf(dir, sizeof dir, "%s/u%d", root, k);
        link_build(dir, upstream);
    }

    const int stride = design.r / configurations > 0 ? design.r / configurations : 1;
    const int n_pairs = configurations * replications;
    Pair *pair = calloc((size_t)n_pairs, sizeof(Pair));
    assert(pair && "dsk_dataset_reproduction: out of memory");
    for (int k = 0; k < configurations; k++) {
        int cop = 1 + k * stride;
        if (cop > design.r) cop = design.r;
        for (int r = 0; r < replications; r++) {
            pair[k * replications + r].cop = cop;
            pair[k * replications + r].replicate = r;
        }
    }

    /* One parameter file per configuration, written once. */
    for (int k = 0; k < configurations; k++) {
        char json_path[256];
        snprintf(json_path, sizeof json_path, "%s/c%d.json", root, pair[k * replications].cop);
        write_inputs(base_json, json_path, parameter, pair[k * replications].cop - 1);
    }

    const double started = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < n_pairs; k++) {
        const int thread = omp_get_thread_num();
        const int seed = pair[k].replicate + 1;

        char json_path[256], mine[192], theirs[192];
        snprintf(json_path, sizeof json_path, "%s/c%d.json", root, pair[k].cop);
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

        pair[k].ran = my_status == 0 && their_status == 0 && my_text && their_text;
        pair[k].results_differ = !pair[k].ran || my_size != their_size ||
                                 memcmp(my_text, their_text, (size_t)my_size) != 0;

        long my_error_size, their_error_size;
        char *my_error_text = read_file(my_errors, &my_error_size);
        char *their_error_text = read_file(their_errors, &their_error_size);
        pair[k].errors_differ = my_error_size != their_error_size ||
                                (my_error_size > 0 && memcmp(my_error_text, their_error_text,
                                                             (size_t)my_error_size) != 0);
        free(my_error_text);
        free(their_error_text);

        if (pair[k].ran) pair[k].series = series_from_results(my_results);
        pair[k].long_horizon_logs = count_long_horizon_logs(mine, seed);

        free(my_text);
        free(their_text);
        unlink(my_results);
        unlink(their_results);
        unlink(my_errors);
        unlink(their_errors);
    }

    const double elapsed = omp_get_wtime() - started;

    /* The archives, read one at a time: the comparison is arithmetic, the
       reading is not what takes the time. */
    int compared = 0, missing = 0, differing_runs = 0, failed_runs = 0, fired = 0;
    int results_differing = 0, errors_differing = 0;
    long differing_values = 0;
    double worst_gap = 0;
    int worst_cop = 0, worst_replicate = 0;

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_dataset_reproduction: cannot open the report path");
    fprintf(report, "The model as it stands against the authors' build and against the archives\n"
                    "of the 1000 x 1000 experiment.\n\n");
    fprintf(report, "%d configurations spread across the design, replications 0 to %d of each,\n"
                    "%d pairs, 600 periods, %d threads, %.1f minutes of wall clock.\n",
            configurations, replications - 1, n_pairs, threads, elapsed / 60.0);
    fprintf(report, "Seeds are the replication number plus one, the same in every configuration,\n"
                    "as the experiment was run.\n\n");

    for (int k = 0; k < n_pairs; k++) {
        if (!pair[k].ran) {
            failed_runs++;
            fprintf(report, "cop_%04d replicate %d: a build did not produce a results file  FAILED\n",
                    pair[k].cop, pair[k].replicate);
            continue;
        }
        if (pair[k].results_differ) {
            results_differing++;
            fprintf(report, "cop_%04d replicate %d: the results file differs from the authors' build  FAILED\n",
                    pair[k].cop, pair[k].replicate);
        }
        if (pair[k].errors_differ) {
            errors_differing++;
            fprintf(report, "cop_%04d replicate %d: the error log differs from the authors' build  FAILED\n",
                    pair[k].cop, pair[k].replicate);
        }
        fired += pair[k].long_horizon_logs;

        char cop_dir[512];
        snprintf(cop_dir, sizeof cop_dir, "%s/cop_%04d", DATASET, pair[k].cop);
        char batch_path[640];
        abm_system_batch_path(batch_path, sizeof batch_path, cop_dir, pair[k].replicate);
        if (!file_exists(batch_path)) { missing++; mat_free(pair[k].series); continue; }

        Mat stored = abm_system_read_replicate(cop_dir, pair[k].replicate);
        Mat produced = pair[k].series;
        if (produced.r != stored.r || produced.c != stored.c) {
            differing_runs++;
            fprintf(report, "cop_%04d replicate %d: %d x %d produced against %d x %d stored  FAILED\n",
                    pair[k].cop, pair[k].replicate, produced.r, produced.c, stored.r, stored.c);
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
                            worst_gap = gap; worst_cop = pair[k].cop; worst_replicate = pair[k].replicate;
                        }
                    }
                }
            }
            if (differing_here) {
                differing_runs++;
                differing_values += differing_here;
                fprintf(report, "cop_%04d replicate %d: %ld of %d stored values differ  FAILED\n",
                        pair[k].cop, pair[k].replicate, differing_here, stored.r * stored.c);
            }
            compared++;
        }
        mat_free(stored);
        mat_free(produced);
    }

    fprintf(report, "pairs run                       %d\n", n_pairs - failed_runs);
    fprintf(report, "results files against upstream  %d compared, %d differing\n",
            n_pairs - failed_runs, results_differing);
    fprintf(report, "error logs against upstream     %d compared, %d differing\n",
            n_pairs - failed_runs, errors_differing);
    fprintf(report, "runs against the archives       %d compared, %d differing\n", compared, differing_runs);
    if (differing_values)
        fprintf(report, "values differing                %ld, worst relative gap %.3e at cop_%04d replicate %d\n",
                differing_values, worst_gap, worst_cop, worst_replicate);
    fprintf(report, "archives missing                %d\n", missing);
    fprintf(report, "long-horizon logs written       %d (must be 0)\n", fired);
    fprintf(report, "runs that did not complete      %d\n", failed_runs);

    int failed = failed_runs || results_differing || errors_differing || differing_runs || fired || !compared;
    fprintf(report, "\n%s\n", failed ? "FAILED" : "PASSED");
    fclose(report);

    printf("dataset reproduction: %d pairs, %d differing from upstream, %d differing from the archives, "
           "%d long-horizon logs, %.1f minutes\n",
           n_pairs, results_differing, differing_runs, fired, elapsed / 60.0);
    printf("%s\n", failed ? "FAILED" : "PASSED, 0 failures");

    free(pair);
    df_free(&design);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
