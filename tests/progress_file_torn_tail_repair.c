/*
Whether progress_file_drop_torn_tail (montecarlo/progress_file.h) leaves a
grid's progress file holding only whole rows, and leaves a sound one alone.

1. A file of whole rows is unchanged byte for byte.
2. A file whose last row lacks its newline loses that row only.
3. A file ending in a cut row followed by NUL bytes loses the cut row and the
   NUL bytes only.
4. A file of NUL bytes only is removed.
5. A progress file cut by an unclean shutdown,
   out/montecarlo/state_weight_trial/seri_uncorrected.csv, repaired in a copy:
   it keeps its bytes up to the line holding its first NUL byte, ends in a
   newline and holds no NUL byte.

Files are written to
out/tests/progress_file_torn_tail_repair/<YYYYmmdd_HHMMSS>/ with the report.
Exits 1 on any failure.
*/

#include "montecarlo/progress_file.h"
#include <stdlib.h>
#include <time.h>

#define OUTPUT_ROOT "out/tests/progress_file_torn_tail_repair"
#define TRIAL_PROGRESS "out/montecarlo/state_weight_trial/seri_uncorrected.csv"

static int failures = 0;
static FILE *report;

static void check(int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    fflush(report);
    if (!ok) failures++;
}

static void write_bytes(const char *path, const char *bytes, size_t length) {
    FILE *f = fopen(path, "wb");
    assert(f);
    fwrite(bytes, 1, length, f);
    fclose(f);
}

/* The whole file in a fresh buffer, NULL when it does not exist. */
static char *read_bytes(const char *path, size_t *length) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    *length = (size_t)ftell(f);
    fseek(f, 0, SEEK_SET);
    char *bytes = malloc(*length + 1);
    assert(bytes);
    size_t got = fread(bytes, 1, *length, f);
    assert(got == *length);
    fclose(f);
    return bytes;
}

static int file_equals(const char *path, const char *expected, size_t expected_length) {
    size_t length = 0;
    char *bytes = read_bytes(path, &length);
    int equal = bytes && length == expected_length && memcmp(bytes, expected, length) == 0;
    free(bytes);
    return equal;
}

int main(void) {
    char stamp[32], dir[256], path[400];
    time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", localtime(&now));
    mkdir("out", 0755);
    mkdir("out/tests", 0755);
    mkdir(OUTPUT_ROOT, 0755);
    snprintf(dir, sizeof dir, "%s/%s", OUTPUT_ROOT, stamp);
    int created = mkdir(dir, 0755) == 0;
    assert(created && "progress_file_torn_tail_repair: the output directory exists already");
    (void)created;
    snprintf(path, sizeof path, "%s/progress_file_torn_tail_repair_report.txt", dir);
    report = fopen(path, "w");
    assert(report);

    const char whole[] = "model,cop,replicate\nlp_s1,1,0,0.5\nlp_s1,2,0,0.25\n";
    size_t whole_length = sizeof whole - 1;

    snprintf(path, sizeof path, "%s/whole_rows.csv", dir);
    write_bytes(path, whole, whole_length);
    long dropped = progress_file_drop_torn_tail(path);
    check(dropped == 0 && file_equals(path, whole, whole_length), "whole rows: unchanged byte for byte");

    char cut[256];
    size_t cut_length = (size_t)snprintf(cut, sizeof cut, "%slp_s1,3,0,0.1", whole);
    snprintf(path, sizeof path, "%s/cut_row.csv", dir);
    write_bytes(path, cut, cut_length);
    dropped = progress_file_drop_torn_tail(path);
    check(dropped == (long)(cut_length - whole_length) && file_equals(path, whole, whole_length),
          "a last row without its newline: that row dropped, the rest unchanged");

    char nul_tail[512];
    memcpy(nul_tail, cut, cut_length);
    memset(nul_tail + cut_length, 0, 200);
    snprintf(path, sizeof path, "%s/cut_row_then_nul.csv", dir);
    write_bytes(path, nul_tail, cut_length + 200);
    dropped = progress_file_drop_torn_tail(path);
    check(dropped == (long)(cut_length + 200 - whole_length) && file_equals(path, whole, whole_length),
          "a cut row then NUL bytes: both dropped, the rest unchanged");

    char nul_only[300];
    memset(nul_only, 0, sizeof nul_only);
    snprintf(path, sizeof path, "%s/nul_only.csv", dir);
    write_bytes(path, nul_only, sizeof nul_only);
    dropped = progress_file_drop_torn_tail(path);
    struct stat st;
    check(dropped == (long)sizeof nul_only && stat(path, &st) != 0, "NUL bytes only: the file removed");

    size_t trial_length = 0;
    char *trial = read_bytes(TRIAL_PROGRESS, &trial_length);
    check(trial != NULL, TRIAL_PROGRESS " readable");
    if (trial) {
        size_t first_nul = trial_length;
        for (size_t i = 0; i < trial_length; i++)
            if (trial[i] == '\0') { first_nul = i; break; }
        size_t keep = 0;
        for (size_t i = 0; i < first_nul; i++)
            if (trial[i] == '\n') keep = i + 1;
        snprintf(path, sizeof path, "%s/seri_uncorrected_trial_copy.csv", dir);
        write_bytes(path, trial, trial_length);
        dropped = progress_file_drop_torn_tail(path);
        size_t length = 0;
        char *repaired = read_bytes(path, &length);
        int no_nul = repaired && memchr(repaired, '\0', length) == NULL;
        char what[600];
        snprintf(what, sizeof what, "trial file of %zu bytes, first NUL at byte %zu: %zu bytes kept, %ld dropped, "
                                    "ends in a newline, no NUL left", trial_length, first_nul, length, dropped);
        check(repaired && length == keep && memcmp(repaired, trial, keep) == 0 && length > 0
                  && repaired[length - 1] == '\n' && no_nul,
              what);
        free(repaired);
        free(trial);
    }

    fprintf(report, "\n%s, %d failures\n", failures ? "FAILED" : "passed", failures);
    fclose(report);
    printf("%s, %d failures, report in %s/progress_file_torn_tail_repair_report.txt\n", failures ? "FAILED" : "passed",
           failures, dir);
    return failures ? 1 : 0;
}
