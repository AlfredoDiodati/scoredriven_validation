/*
Moves the rows of the local projections fitted with inflation and the interest
rate in first differences out of the results of the two Seri grids of
montecarlo/ (adjusted, with the mean from the last 100 runs, and unadjusted)
and into the results of their copies in this folder.

For each grid, the rows of out/montecarlo/<grid>.csv.gz are parted by model:
the rows of the differenced models (model_is_differenced of
montecarlo/response_cache.h) go to
out/old/montecarlo_differenced_lp/<grid>.csv.gz, every other row back to
out/montecarlo/<grid>.csv.gz, both under the same header and in the order they
had. The new main result is written beside the old one and renamed over it only
after the check below.

Check, before anything is renamed: both new files are read back, and walking
the original rows in order, each one equals byte for byte the next row of the
file its model went to, with no row left over in either file. The two files are
then exactly a partition of the original rows.

A result is up to a few GB of text and et_al's gzip works on whole buffers, so
the rows are streamed through the system's zcat and gzip -9 instead.

A grid whose result holds no rows to move is left alone, a grid whose old
result already exists is refused, so rows cannot be moved twice, and a grid
with a progress file is refused, since it is still running. The log goes to
out/old/montecarlo_differenced_lp/differenced_lp_rows_split.log.
*/

/* getline and popen are POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define OUTPUT_DIR "out/old/montecarlo_differenced_lp"
#define LOG_PATH OUTPUT_DIR "/differenced_lp_rows_split.log"
#define N_GRIDS 2

static const char *const grid_name[N_GRIDS] = {
    "seri_held_out_mean_true_configuration_recovery_100_mean_runs", "seri_uncorrected_true_configuration_recovery",
};

/* The rows that move: those of the models this folder runs. */
static int moves(int model) {
    return model_is_differenced(model);
}

static FILE *log_file;

static void log_line(const char *what) {
    char stamp[32];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", localtime(&now));
    fprintf(log_file, "%s %s\n", stamp, what);
    fflush(log_file);
}

static int exists(const char *path) {
    struct stat info;
    return stat(path, &info) == 0;
}

static FILE *open_pipe(const char *format, const char *path, const char *mode) {
    char command[700];
    snprintf(command, sizeof command, format, path);
    FILE *pipe = popen(command, mode);
    assert(pipe && "rows_split: cannot start zcat or gzip");
    return pipe;
}

/* Whether a row's model moves. A row whose label is no model's stays. */
static int row_moves(const char *line) {
    const char *comma = strchr(line, ',');
    if (!comma) return 0;
    size_t length = (size_t)(comma - line);
    for (int model = 0; model < N_MODELS; model++)
        if (strlen(model_label[model]) == length && strncmp(line, model_label[model], length) == 0)
            return moves(model);
    return 0;
}

/* Whether the next line of stream equals line, length bytes. */
static int next_line_equals(FILE *stream, char **buffer, size_t *capacity, const char *line, ssize_t length) {
    ssize_t got = getline(buffer, capacity, stream);
    return got == length && memcmp(*buffer, line, (size_t)length) == 0;
}

static int split_grid(int index) {
    char main_path[256], old_path[256], progress_path[256], staged_path[300], what[1200];
    snprintf(main_path, sizeof main_path, "out/montecarlo/%s.csv.gz", grid_name[index]);
    snprintf(progress_path, sizeof progress_path, "out/montecarlo/%s_progress.csv", grid_name[index]);
    snprintf(old_path, sizeof old_path, "%s/%s.csv.gz", OUTPUT_DIR, grid_name[index]);
    snprintf(staged_path, sizeof staged_path, "%s.staged", main_path);
    if (exists(progress_path)) {
        snprintf(what, sizeof what, "%s: refused, %s exists, the grid is still running", grid_name[index],
                 progress_path);
        log_line(what);
        return 1;
    }
    if (exists(old_path)) {
        snprintf(what, sizeof what, "%s: refused, %s exists already", grid_name[index], old_path);
        log_line(what);
        return 1;
    }

    FILE *in = open_pipe("zcat %s", main_path, "r");
    FILE *main_out = open_pipe("gzip -9 > %s", staged_path, "w");
    FILE *old_out = open_pipe("gzip -9 > %s", old_path, "w");
    char *line = NULL;
    size_t capacity = 0;
    long kept = 0, moved = 0;
    ssize_t length = getline(&line, &capacity, in);
    assert(length > 0 && "rows_split: a result has no header");
    fwrite(line, 1, (size_t)length, main_out);
    fwrite(line, 1, (size_t)length, old_out);
    while ((length = getline(&line, &capacity, in)) > 0) {
        if (row_moves(line)) {
            fwrite(line, 1, (size_t)length, old_out);
            moved++;
        } else {
            fwrite(line, 1, (size_t)length, main_out);
            kept++;
        }
    }
    int read_ok = pclose(in) == 0;
    int written_ok = pclose(main_out) == 0;
    written_ok &= pclose(old_out) == 0;
    if (!read_ok || !written_ok) {
        snprintf(what, sizeof what, "%s: FAILED, zcat or gzip exited with an error; %s unchanged, %s and %s left "
                                    "for inspection", grid_name[index], main_path, staged_path, old_path);
        log_line(what);
        free(line);
        return 1;
    }
    if (!moved) {
        remove(staged_path);
        remove(old_path);
        snprintf(what, sizeof what, "%s: no rows to move, left alone (%ld rows)", grid_name[index], kept);
        log_line(what);
        free(line);
        return 0;
    }

    in = open_pipe("zcat %s", main_path, "r");
    FILE *main_back = open_pipe("zcat %s", staged_path, "r");
    FILE *old_back = open_pipe("zcat %s", old_path, "r");
    char *back = NULL;
    size_t back_capacity = 0;
    length = getline(&line, &capacity, in);
    int ok = next_line_equals(main_back, &back, &back_capacity, line, length)
          && next_line_equals(old_back, &back, &back_capacity, line, length);
    while (ok && (length = getline(&line, &capacity, in)) > 0)
        ok = next_line_equals(row_moves(line) ? old_back : main_back, &back, &back_capacity, line, length);
    ok = ok && getline(&back, &back_capacity, main_back) < 0 && getline(&back, &back_capacity, old_back) < 0;
    pclose(in);
    pclose(main_back);
    pclose(old_back);
    free(back);
    free(line);

    if (!ok) {
        snprintf(what, sizeof what, "%s: FAILED, the parted files are not a partition of the result; %s "
                                    "unchanged, %s and %s left for inspection", grid_name[index], main_path,
                 staged_path, old_path);
        log_line(what);
        return 1;
    }
    int renamed = rename(staged_path, main_path) == 0;
    assert(renamed && "rows_split: cannot rename the new main result");
    (void)renamed;
    snprintf(what, sizeof what, "%s: %ld rows kept in %s, %ld rows moved to %s; every original row found in "
                                "order in one of the two", grid_name[index], kept, main_path, moved, old_path);
    log_line(what);
    return 0;
}

int main(void) {
    mkdir("out", 0755);
    mkdir("out/old", 0755);
    mkdir(OUTPUT_DIR, 0755);
    log_file = fopen(LOG_PATH, "a");
    assert(log_file && "rows_split: cannot open the log");
    log_line("started");
    int failures = 0;
    for (int index = 0; index < N_GRIDS; index++) failures += split_grid(index);
    log_line(failures ? "ended with failures" : "ended");
    fclose(log_file);
    return failures ? 1 : 0;
}
