/*
Whether the simulator still works when it is invoked through a long path.

It did not. The model builds every output file name by strcpy and strcat into
fixed 64-byte buffers, starting from the directory of the path it was invoked
as, so the length of that path decides whether the name fits. The error file's
name is the longest one: the executable's own directory, plus "output", plus
"/errors/Errors", plus the run name and seed, plus ".txt". Past about 26
characters of directory the writes ran off the end of those buffers and into
the globals that follow them. Unoptimised the damage was silent; optimised it
was a segmentation fault, and both builds crashed once the path was long
enough - measured at 101 bytes of error-file path for the unoptimised build
and 66 for the optimised one.

This is not a hypothetical. applications/abm_system_simulate.c reaches the
model through a symlink in a scratch directory, and on a cluster that
directory is $SLURM_TMPDIR/dsk_<pid>, which is exactly the length that
triggers it. A million runs would have failed as nonzero exit statuses.

model/dsk_sfc widens those buffers to PATH_MAX. This checks the fix from both
directions: the simulator must complete from a path far past the old limit,
and it must produce the same output it produces from a short one, since a name
that silently overflowed into a neighbouring global could just as easily have
changed a number as crashed.

Needs model/dsk_sfc/dsk_SFC built - `make model` - and writes
out/dsk_long_path.txt.
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#define MODEL "model/dsk_sfc/dsk_SFC"
#define INPUTS "model/dsk_sfc/dsk_sfc_inputs.json"
#define REPORT "out/dsk_long_path.txt"
#define RUN_NAME "longpath"
#define SEED 1

/* Long enough that the error file's name passes the old 64-byte buffers by a
   wide margin, short enough to stay well inside PATH_MAX. */
#define LONG_DIRECTORY_NAME \
    "a_directory_deep_enough_that_the_old_sixtyfour_byte_name_buffers_overflowed"

static void make_directory(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk_long_path: mkdir failed");
}

/* A directory holding a symlink to the model, so the model's output lands in
   it: the model writes beside the path it was invoked as, which is the whole
   reason the invocation path's length matters here. */
static void prepare(const char *dir, const char *model) {
    make_directory(dir);

    char link[1024];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(model, link) == 0 && "dsk_long_path: cannot link the model into place");
}

static int run(const char *dir, const char *inputs) {
    char command[8192];
    snprintf(command, sizeof command, "\"%s/dsk_SFC\" \"%s\" -r %s -s %d -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, RUN_NAME, SEED);
    return system(command);
}

static char *read_results(const char *dir, long *size) {
    char path[1024];
    snprintf(path, sizeof path, "%s/output/results_%s_%d.txt", dir, RUN_NAME, SEED);

    FILE *f = fopen(path, "rb");
    if (!f) { *size = -1; return NULL; }

    fseek(f, 0, SEEK_END);
    *size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *text = malloc((size_t)*size);
    assert(text && fread(text, 1, (size_t)*size, f) == (size_t)*size &&
           "dsk_long_path: cannot read a results file");
    fclose(f);
    return text;
}

int main(void) {
    char model[PATH_MAX];
    assert(realpath(MODEL, model) &&
           "dsk_long_path: model/dsk_sfc/dsk_SFC is not built - run make model");
    char inputs[PATH_MAX];
    assert(realpath(INPUTS, inputs) && "dsk_long_path: the model's own inputs JSON is missing");

    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    char root[512], short_dir[640], long_dir[1024];
    snprintf(root, sizeof root, "%s/dsk_long_path_%d", scratch, (int)getpid());
    make_directory(root);
    snprintf(short_dir, sizeof short_dir, "%s/s", root);
    snprintf(long_dir, sizeof long_dir, "%s/%s", root, LONG_DIRECTORY_NAME);

    prepare(short_dir, model);
    prepare(long_dir, model);

    int short_status = run(short_dir, inputs);
    int long_status = run(long_dir, inputs);

    long short_size, long_size;
    char *short_text = read_results(short_dir, &short_size);
    char *long_text = read_results(long_dir, &long_size);

    int completed = (short_status == 0 && long_status == 0 && short_text && long_text);
    int identical = (completed && short_size == long_size &&
                     memcmp(short_text, long_text, (size_t)short_size) == 0);

    /* The name whose length the old buffers could not hold. */
    size_t error_name_length = strlen(long_dir) + strlen("/output/errors/Errors_") +
                               strlen(RUN_NAME) + strlen("_1.txt");

    FILE *report = fopen(REPORT, "w");
    assert(report && "dsk_long_path: cannot open the report path");
    fprintf(report, "Invoking model/dsk_sfc/dsk_SFC through a short and a long path.\n\n");
    fprintf(report, "short path  %s\n", short_dir);
    fprintf(report, "long path   %s\n", long_dir);
    fprintf(report, "the long run's error file name is %zu bytes, against the 64 the\n"
                    "buffers held before the fix\n\n", error_name_length);
    fprintf(report, "short path: exit %d, %ld bytes of results\n", short_status, short_size);
    fprintf(report, "long path:  exit %d, %ld bytes of results\n", long_status, long_size);
    fprintf(report, "outputs identical: %s\n\n", identical ? "yes" : "NO");
    fprintf(report, "%s\n", identical ? "PASSED" : "FAILED");
    fclose(report);

    printf("long-path invocation: error file name %zu bytes\n", error_name_length);
    printf("  short path exit %d, long path exit %d\n", short_status, long_status);
    printf("  outputs identical: %s\n", identical ? "yes" : "NO");
    printf("%s\n", identical ? "PASSED, 0 failures" : "FAILED");

    free(short_text);
    free(long_text);
    return identical ? EXIT_SUCCESS : EXIT_FAILURE;
}
