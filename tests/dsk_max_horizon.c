/*
How long a single run of the simulator can be before it breaks.

docs/DSK_LONG_HORIZON.md lists what stops a long run: two upstream bugs near
periods 1,640 and 2,788, the cost of rationing near 8,889, machine counts
passing 2^53 near 12,908, and the levels of the economy growing towards the
largest double. model/dsk_sfc fixes or rescales all of them except world
emissions, which overflow near period 127,672. The longest run any other test
makes is dsk_machine_lot_rebase's 14,000 periods. This one runs further and
records where, if anywhere, the run stops being usable.

One run, seed 1, the model's own parameter file with only the horizon
overwritten, -f 0 -c 0 -v 0. The horizon is HORIZON unless the first argument
names another. What limits it on a given machine is memory: four arrays are
sized by the full horizon, g, gtemp and g_price one double per period,
capital-good firm and consumption-good firm, and acquired one int, so with
N1 = 20 and N2 = 200 a run holds about 112 KB per period of horizon. The test
refuses to start a run that would not fit in the memory the kernel reports as
available, because the out-of-memory killer would end it with nothing learned.

While the run goes, the results file is watched once a second and its rows
counted. A run whose results file gains no row for STALL_SECONDS is killed and
counted as stalled at the period after its last row: every hang the docs record
is a loop that never exits, and the model's own alarm would only end it after
twice the horizon in seconds.

Afterwards every value of the results file is read. What is required:

  - the run exits with status 0 and writes all HORIZON rows;
  - its error log is empty, which also covers the model's own accounting check
    and its carbon-cycle convergence message;
  - every row splits into 83 separate numbers. SAVE() in dsk_sfc_main.cpp
    writes each value right-aligned in FIELD_WIDTH characters, fixed point,
    ten decimals, so a value of 1e49 or more fills its field and is written
    against the previous value with no space between them. Read back, the two
    become one malformed number, or two wrong ones;
  - every value in every one of the 83 columns is finite.

What is reported besides: the first period at which any column, and any of the
five columns the pipeline reads, stops being finite; the first line of the
error log; how many times money was redenominated, machines rebased and the
unit of the good changed, from the logs the model writes for each; the elapsed
time every PROGRESS_STEP periods; peak resident memory.

Writes out/tests/dsk_max_horizon.txt.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <errno.h>
#include <assert.h>
#include "tests/dsk_scratch.h"

#define MODEL "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/tests/dsk_max_horizon.txt"
#define RUN_NAME "horizon"
#define MODEL_COLUMNS 83
#define SEED 1
#define HORIZON 40000
#define BYTES_PER_PERIOD 112000.0
#define STALL_SECONDS 600
#define FIELD_WIDTH 60
#define PROGRESS_STEP 2000
#define MAX_PROGRESS_MARKS 1024

/* The five columns applications/abm_system_simulate.c keeps: real GDP, energy
   demand, employment, the price level and the policy rate. */
static const int pipeline_column[] = {2, 8, 5, 34, 50};
#define PIPELINE_COLUMNS (int)(sizeof pipeline_column / sizeof pipeline_column[0])

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_max_horizon: mkdir failed");
}

static void write_inputs(const char *source, const char *target, int periods) {
    FILE *in = fopen(source, "rb");
    assert(in && "dsk_max_horizon: cannot read the model's parameter file");
    FILE *out = fopen(target, "wb");
    assert(out && "dsk_max_horizon: cannot write a parameter file");

    char line[4096];
    while (fgets(line, sizeof line, in)) {
        if (strstr(line, "\"T\"")) fprintf(out, "      \"T\": %d,\n", periods);
        else fputs(line, out);
    }
    fclose(out);
    fclose(in);
}

/* MemAvailable from /proc/meminfo, in bytes, or -1 when it cannot be read. */
static double available_memory(void) {
    FILE *f = fopen("/proc/meminfo", "r");
    if (!f) return -1;
    char line[256];
    double kilobytes = -1;
    while (fgets(line, sizeof line, f))
        if (sscanf(line, "MemAvailable: %lf kB", &kilobytes) == 1) break;
    fclose(f);
    return kilobytes < 0 ? -1 : kilobytes * 1024.0;
}

static double seconds_since(const struct timespec *start) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double)(now.tv_sec - start->tv_sec) + 1e-9 * (double)(now.tv_nsec - start->tv_nsec);
}

static pid_t start_model(const char *executable, const char *inputs) {
    char seed[32];
    snprintf(seed, sizeof seed, "%d", SEED);
    pid_t pid = fork();
    assert(pid >= 0 && "dsk_max_horizon: fork failed");
    if (pid == 0) {
        int null = open("/dev/null", O_WRONLY);
        if (null >= 0) {
            dup2(null, STDOUT_FILENO);
            dup2(null, STDERR_FILENO);
            close(null);
        }
        char *const argv[] = {(char *)executable, (char *)inputs, "-r", RUN_NAME, "-s", seed,
                              "-f", "0", "-c", "0", "-v", "0", NULL};
        execv(executable, argv);
        _exit(127);
    }
    return pid;
}

/* Counts the newlines added to the results file since the last call. The file
   does not exist until the model has set up, so a missing file counts nothing. */
static long count_new_rows(const char *path, FILE **results) {
    if (!*results) *results = fopen(path, "r");
    if (!*results) return 0;
    char chunk[1 << 16];
    long rows = 0;
    size_t got;
    while ((got = fread(chunk, 1, sizeof chunk, *results)) > 0)
        for (size_t i = 0; i < got; i++) rows += chunk[i] == '\n';
    clearerr(*results);
    return rows;
}

static long file_size(const char *path) {
    struct stat info;
    if (stat(path, &info) != 0) return -1;
    return (long)info.st_size;
}

static long count_lines(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    long lines = 0;
    int c;
    while ((c = fgetc(f)) != EOF) lines += c == '\n';
    fclose(f);
    return lines;
}

static void first_nonblank_line(const char *path, char *line, size_t room) {
    line[0] = '\0';
    FILE *f = fopen(path, "r");
    if (!f) return;
    while (fgets(line, (int)room, f)) {
        line[strcspn(line, "\n")] = '\0';
        if (line[strspn(line, " \t")] != '\0') break;
        line[0] = '\0';
    }
    fclose(f);
}

int main(int argc, char **argv) {
    int horizon = argc > 1 ? atoi(argv[1]) : HORIZON;
    assert(horizon > 0 && "dsk_max_horizon: the horizon must be a positive number of periods");

    char model[PATH_MAX], inputs[PATH_MAX];
    assert(realpath(MODEL, model) && "dsk_max_horizon: model/dsk_sfc/dsk_SFC is not built - run make model");
    assert(realpath(INPUTS, inputs) && "dsk_max_horizon: the parameter file is missing");

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_max_horizon: cannot open the report path");
    fprintf(report, "One run of %d periods, seed %d, the model's own parameter file.\n\n", horizon, SEED);

    double needed = BYTES_PER_PERIOD * horizon;
    double available = available_memory();
    fprintf(report, "memory the horizon-sized arrays need: %.2f GB\n", needed / 1e9);
    fprintf(report, "memory available at the start: %.2f GB\n\n", available / 1e9);
    if (available > 0 && needed > available) {
        fprintf(report, "the run would not fit in available memory and was not started\n\nFAILED\n");
        fclose(report);
        printf("max horizon: %d periods need %.2f GB, %.2f GB available\nFAILED\n",
               horizon, needed / 1e9, available / 1e9);
        return EXIT_FAILURE;
    }

    const char *scratch = dsk_scratch_root();
    char root[512], dir[640], link[1024], run_inputs[1024];
    snprintf(root, sizeof root, "%s/dsk_max_horizon_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(dir, sizeof dir, "%s/run", root);
    make_directory(dir);
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    assert(symlink(model, link) == 0 && "dsk_max_horizon: cannot link the model");
    snprintf(run_inputs, sizeof run_inputs, "%s/horizon.json", root);
    write_inputs(inputs, run_inputs, horizon);

    char results_path[1024], errors_path[1024];
    snprintf(results_path, sizeof results_path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, SEED);
    snprintf(errors_path, sizeof errors_path, "%s/output/errors/Errors_%s_%d.txt", dir, RUN_NAME, SEED);

    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
    pid_t pid = start_model(link, run_inputs);

    FILE *watched = NULL;
    long rows = 0;
    double last_growth = 0;
    int marks = 0;
    double mark_seconds[MAX_PROGRESS_MARKS];
    int status = 0, stalled = 0;
    const struct timespec poll_interval = {1, 0};

    for (;;) {
        pid_t done = waitpid(pid, &status, WNOHANG);
        long added = count_new_rows(results_path, &watched);
        double elapsed = seconds_since(&start);
        if (added > 0) last_growth = elapsed;
        rows += added;
        while (marks < MAX_PROGRESS_MARKS && rows >= (long)(marks + 1) * PROGRESS_STEP)
            mark_seconds[marks++] = elapsed;
        if (done == pid) break;
        if (elapsed - last_growth > STALL_SECONDS) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            rows += count_new_rows(results_path, &watched);
            stalled = 1;
            break;
        }
        nanosleep(&poll_interval, NULL);
    }
    double wall = seconds_since(&start);
    if (watched) fclose(watched);

    struct rusage usage;
    getrusage(RUSAGE_CHILDREN, &usage);

    /* Every row, for the first period whose fields run together and the first
       period anything stops being a number. A row that does not split into 83
       clean fields cannot be read value by value, so only clean rows are
       checked for finiteness. */
    int written_rows = 0, clean_rows = 0, first_nonfinite = 0, first_nonfinite_column = 0, first_pipeline_nonfinite = 0;
    int first_pipeline_column = 0, first_merged = 0, first_merged_column = 0, first_merged_fields = 0;
    FILE *results = fopen(results_path, "r");
    if (results) {
        double value[MODEL_COLUMNS];
        char *line = NULL;
        size_t room = 0;
        while (getline(&line, &room, results) > 0) {
            written_rows++;
            int fields = 0, merged_column = 0;
            for (char *field = strtok(line, " \t\n"); field; field = strtok(NULL, " \t\n")) {
                char *end;
                double parsed = strtod(field, &end);
                if ((*end != '\0' || strlen(field) >= FIELD_WIDTH) && !merged_column) merged_column = fields + 1;
                if (fields < MODEL_COLUMNS) value[fields] = parsed;
                fields++;
            }
            if (fields != MODEL_COLUMNS || merged_column) {
                if (!first_merged) {
                    first_merged = written_rows;
                    first_merged_column = merged_column;
                    first_merged_fields = fields;
                }
                continue;
            }
            clean_rows++;
            for (int c = 0; c < MODEL_COLUMNS && !first_nonfinite; c++)
                if (!isfinite(value[c])) {
                    first_nonfinite = written_rows;
                    first_nonfinite_column = c + 1;
                }
            for (int k = 0; k < PIPELINE_COLUMNS && !first_pipeline_nonfinite; k++)
                if (!isfinite(value[pipeline_column[k] - 1])) {
                    first_pipeline_nonfinite = written_rows;
                    first_pipeline_column = pipeline_column[k];
                }
        }
        free(line);
        fclose(results);
    }

    char log_path[1200];
    snprintf(log_path, sizeof log_path, "%s/output/results_%s_%d_redenominations.txt", dir, RUN_NAME, SEED);
    long redenominations = count_lines(log_path);
    snprintf(log_path, sizeof log_path, "%s/output/results_%s_%d_machine_lots.txt", dir, RUN_NAME, SEED);
    long machine_rebases = count_lines(log_path);
    snprintf(log_path, sizeof log_path, "%s/output/results_%s_%d_good_units.txt", dir, RUN_NAME, SEED);
    long good_unit_changes = count_lines(log_path);

    long error_bytes = file_size(errors_path);
    char first_error[512];
    first_nonblank_line(errors_path, first_error, sizeof first_error);

    int exited = WIFEXITED(status);
    int exit_code = exited ? WEXITSTATUS(status) : -1;
    int signal_number = WIFSIGNALED(status) ? WTERMSIG(status) : 0;

    fprintf(report, "wall time: %.1f s\n", wall);
    fprintf(report, "peak resident memory: %.1f MB\n", usage.ru_maxrss / 1024.0);
    if (stalled)
        fprintf(report, "how it ended: killed after %d s without a new row\n", STALL_SECONDS);
    else if (exited)
        fprintf(report, "how it ended: exit status %d\n", exit_code);
    else
        fprintf(report, "how it ended: killed by signal %d\n", signal_number);
    fprintf(report, "rows written: %d of %d\n", written_rows, horizon);
    fprintf(report, "rows that split into %d clean numbers: %d\n", MODEL_COLUMNS, clean_rows);
    if (first_merged)
        fprintf(report, "first row whose fields run together: period %d, at field %d, %d fields on the row\n",
                first_merged, first_merged_column, first_merged_fields);
    else
        fprintf(report, "first row whose fields run together: none\n");
    fprintf(report, "error log: %ld bytes\n", error_bytes < 0 ? 0 : error_bytes);
    if (first_error[0]) fprintf(report, "first error logged: %s\n", first_error);
    if (first_nonfinite)
        fprintf(report, "first non-finite value: period %d, column %d\n", first_nonfinite, first_nonfinite_column);
    else
        fprintf(report, "first non-finite value: none\n");
    if (first_pipeline_nonfinite)
        fprintf(report, "first non-finite value in a pipeline column: period %d, column %d\n",
                first_pipeline_nonfinite, first_pipeline_column);
    else
        fprintf(report, "first non-finite value in a pipeline column: none\n");
    fprintf(report, "money redenominated: %ld times\n", redenominations);
    fprintf(report, "machines rebased: %ld times\n", machine_rebases);
    fprintf(report, "unit of the good: changed %ld times\n\n", good_unit_changes);

    fprintf(report, "elapsed time by period reached, read from the results file once a second\n");
    for (int m = 0; m < marks; m++) {
        double previous = m > 0 ? mark_seconds[m - 1] : 0;
        fprintf(report, "  period %d: %.1f s, %.2f ms a period over the last %d\n", (m + 1) * PROGRESS_STEP,
                mark_seconds[m], 1000.0 * (mark_seconds[m] - previous) / PROGRESS_STEP, PROGRESS_STEP);
    }
    fprintf(report, "\n");

    int failures = 0;
    if (stalled) {
        fprintf(report, "the run stalled at period %d  FAILED\n", written_rows + 1);
        failures++;
    }
    if (!stalled && (!exited || exit_code != 0)) {
        fprintf(report, "the run did not exit cleanly, last row written was period %d  FAILED\n", written_rows);
        failures++;
    }
    if (written_rows != horizon) {
        fprintf(report, "the run wrote %d of %d rows  FAILED\n", written_rows, horizon);
        failures++;
    }
    if (first_merged) {
        fprintf(report, "from period %d the results file cannot be read back  FAILED\n", first_merged);
        failures++;
    }
    if (error_bytes > 0) {
        fprintf(report, "the error log is not empty  FAILED\n");
        failures++;
    }
    if (first_nonfinite) {
        fprintf(report, "a value stops being finite at period %d  FAILED\n", first_nonfinite);
        failures++;
    }
    if (failures == 0) fprintf(report, "all %d periods written, finite and without errors\n", horizon);

    fprintf(report, "\n%s\n", failures == 0 ? "PASSED" : "FAILED");
    fclose(report);

    printf("max horizon: %d of %d periods written, %d failures\n", written_rows, horizon, failures);
    printf("%s\n", failures == 0 ? "PASSED, 0 failures" : "FAILED");

    dsk_scratch_finish(root, failures == 0);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
