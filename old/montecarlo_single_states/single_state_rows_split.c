/*
Moves the rows of state 1 and state 2 on their own out of the results of the
four grids of montecarlo/ and into the results of their copies in this folder,
so each grid's result holds the models it runs and nothing else.

For each grid, out/montecarlo/<grid>.csv.gz is read whole and its rows are
parted by model: the rows of the eight single-state models
(model_is_single_state of montecarlo/response_cache.h) go to
out/old/montecarlo_single_states/<grid>.csv.gz, every other row back to
out/montecarlo/<grid>.csv.gz, both under the same header and in the order they
had. The new main result is written beside the old one and renamed over it only
after the check below.

Check, before anything is renamed: both new files are read back, and walking
the original rows in order, each one equals byte for byte the next row of the
file its model went to, with no row left over in either file. The two files are
then exactly a partition of the original rows.

A grid whose result holds no single-state rows is left alone, and a grid whose
old result already exists is refused, so the program cannot move rows twice.
A grid with a progress file is refused too: it is still running. The log goes
to out/old/montecarlo_single_states/single_state_rows_split.log.
*/

/* fmemopen is POSIX, hidden by -std=c11 otherwise. */
#define _POSIX_C_SOURCE 200809L

#include "montecarlo/response_cache.h"
#include <et_al./frame/gzip.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define OUTPUT_DIR "out/old/montecarlo_single_states"
#define LOG_PATH OUTPUT_DIR "/single_state_rows_split.log"
#define N_GRIDS 4

static const char *const grid_name[N_GRIDS] = {
    "true_configuration_recovery", "uncorrected_true_configuration_recovery",
    "seri_true_configuration_recovery", "seri_uncorrected_true_configuration_recovery",
};

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

/* The file inflated whole, NULL-terminated; its length in *length. */
static char *read_inflated(const char *path, size_t *length) {
    FILE *f = fopen(path, "rb");
    assert(f && "single_state_rows_split: cannot open a result");
    fseek(f, 0, SEEK_END);
    long compressed_length = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *compressed = malloc((size_t)compressed_length);
    assert(compressed && "single_state_rows_split: out of memory reading a result");
    size_t n_read = fread(compressed, 1, (size_t)compressed_length, f);
    assert(n_read == (size_t)compressed_length && "single_state_rows_split: short read of a result");
    (void)n_read;
    fclose(f);
    size_t text_length;
    unsigned char *text = gzip_inflate(compressed, (size_t)compressed_length, &text_length);
    assert(text && "single_state_rows_split: a result did not inflate");
    free(compressed);
    char *terminated = realloc(text, text_length + 1);
    assert(terminated);
    terminated[text_length] = '\0';
    *length = text_length;
    return terminated;
}

static void write_deflated(const char *path, const char *text, size_t length) {
    size_t compressed_length;
    unsigned char *compressed = gzip_deflate_level((const unsigned char *)text, length, 9, &compressed_length);
    assert(compressed && "single_state_rows_split: deflate failed");
    FILE *out = fopen(path, "wb");
    assert(out && "single_state_rows_split: cannot open an output for writing");
    size_t n_written = fwrite(compressed, 1, compressed_length, out);
    assert(n_written == compressed_length && "single_state_rows_split: short write");
    (void)n_written;
    fclose(out);
    free(compressed);
}

/* Whether the row starting at line belongs to a single-state model. A row whose
   label is no model's stays with the main result. */
static int is_single_state_row(const char *line) {
    const char *comma = strchr(line, ',');
    if (!comma) return 0;
    size_t length = (size_t)(comma - line);
    for (int model = 0; model < N_MODELS; model++)
        if (strlen(model_label[model]) == length && strncmp(line, model_label[model], length) == 0)
            return model_is_single_state(model);
    return 0;
}

/* The line starting at text, newline included; its length. */
static size_t line_length(const char *text, const char *end) {
    const char *newline = memchr(text, '\n', (size_t)(end - text));
    return newline ? (size_t)(newline - text) + 1 : (size_t)(end - text);
}

static int split_grid(int index) {
    char main_path[256], old_path[256], progress_path[256], staged_path[300], what[1200];
    snprintf(main_path, sizeof main_path, "out/montecarlo/%s.csv.gz", grid_name[index]);
    snprintf(progress_path, sizeof progress_path, "out/montecarlo/%s_progress.csv", grid_name[index]);
    snprintf(old_path, sizeof old_path, "%s/%s.csv.gz", OUTPUT_DIR, grid_name[index]);
    snprintf(staged_path, sizeof staged_path, "%s.without_single_states", main_path);
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

    size_t length;
    char *text = read_inflated(main_path, &length);
    const char *end = text + length;
    size_t header_length = line_length(text, end);
    char *main_text = malloc(length + 1), *old_text = malloc(length + 1);
    assert(main_text && old_text && "single_state_rows_split: out of memory for the parted rows");
    memcpy(main_text, text, header_length);
    memcpy(old_text, text, header_length);
    size_t main_length = header_length, old_length = header_length;
    long main_rows = 0, old_rows = 0;
    for (const char *line = text + header_length; line < end;) {
        size_t n = line_length(line, end);
        if (is_single_state_row(line)) {
            memcpy(old_text + old_length, line, n);
            old_length += n;
            old_rows++;
        } else {
            memcpy(main_text + main_length, line, n);
            main_length += n;
            main_rows++;
        }
        line += n;
    }
    if (!old_rows) {
        snprintf(what, sizeof what, "%s: no single-state rows, left alone (%ld rows)", grid_name[index], main_rows);
        log_line(what);
        free(text);
        free(main_text);
        free(old_text);
        return 0;
    }

    write_deflated(staged_path, main_text, main_length);
    write_deflated(old_path, old_text, old_length);
    free(main_text);
    free(old_text);

    size_t main_back_length, old_back_length;
    char *main_back = read_inflated(staged_path, &main_back_length);
    char *old_back = read_inflated(old_path, &old_back_length);
    const char *main_end = main_back + main_back_length, *old_end = old_back + old_back_length;
    int ok = main_back_length >= header_length && old_back_length >= header_length
          && memcmp(main_back, text, header_length) == 0 && memcmp(old_back, text, header_length) == 0;
    const char *main_at = main_back + header_length, *old_at = old_back + header_length;
    for (const char *line = text + header_length; ok && line < end;) {
        size_t n = line_length(line, end);
        const char **at = is_single_state_row(line) ? &old_at : &main_at;
        const char *at_end = at == &old_at ? old_end : main_end;
        ok = (size_t)(at_end - *at) >= n && memcmp(*at, line, n) == 0;
        *at += n;
        line += n;
    }
    ok = ok && main_at == main_end && old_at == old_end;
    free(main_back);
    free(old_back);
    free(text);

    if (!ok) {
        snprintf(what, sizeof what, "%s: FAILED, the parted files are not a partition of the result; %s left "
                                    "unchanged, %s and %s left for inspection", grid_name[index], main_path,
                 staged_path, old_path);
        log_line(what);
        return 1;
    }
    int renamed = rename(staged_path, main_path) == 0;
    assert(renamed && "single_state_rows_split: cannot rename the new main result");
    (void)renamed;
    snprintf(what, sizeof what, "%s: %ld rows kept in %s, %ld single-state rows moved to %s; every original row "
                                "found in order in one of the two", grid_name[index], main_rows, main_path,
             old_rows, old_path);
    log_line(what);
    return 0;
}

int main(void) {
    mkdir("out", 0755);
    mkdir("out/old", 0755);
    mkdir(OUTPUT_DIR, 0755);
    log_file = fopen(LOG_PATH, "a");
    assert(log_file && "single_state_rows_split: cannot open the log");
    log_line("started");
    int failures = 0;
    for (int index = 0; index < N_GRIDS; index++) failures += split_grid(index);
    log_line(failures ? "ended with failures" : "ended");
    fclose(log_file);
    return failures ? 1 : 0;
}
