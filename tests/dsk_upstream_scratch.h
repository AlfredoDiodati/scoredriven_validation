/*
A scratch directory the upstream DSK build can actually be run in.

Upstream copies the output path into a buffer it sizes at exactly the path's
length, so every run writes the terminating byte one past the end of it
(dsk_sfc_main.cpp, the pathname array beside outstr). Whether that byte matters
depends on what the compiler placed next to the array, which depends on the
length of the path: at some lengths the run writes its files, and at others it
writes nothing at all, reports nothing and still exits 0. Upstream also declares
every output filename as char[64] and builds the names with strcat, so a
directory deep enough overruns those as well.

This project's build fixed both (docs/DSK_MODEL_CHANGES.md), but the reference
build the equivalence tests compare against is upstream's code and has to be run
somewhere it works. A test that did not check would compare its own output
against a file that was never written, or against one left over from an earlier
run.

So the directory name is lengthened one character at a time until upstream
writes a results file there, and the length that worked is reported. The
comparison then runs both builds in sibling directories of that same length.
*/

#ifndef DSK_UPSTREAM_SCRATCH_H
#define DSK_UPSTREAM_SCRATCH_H

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

/* The size of every filename buffer in model/dsk_sfc/upstream. */
#define DSK_UPSTREAM_FILENAME_BYTES 64

/* How many characters of padding to try before giving up. Two lengths in a row
   have never both failed here, but the search is cheap and a failure to find
   one is a clear message rather than a silent wrong answer. */
#define DSK_SCRATCH_ATTEMPTS 8

typedef struct {
    char root[256];
    char upstream[288];
    char modified[288];
    int padding;
} DskScratch;

static void dsk_scratch_mkdir(const char *path) {
    if (mkdir(path, 0755) != 0) assert(errno == EEXIST && "dsk scratch: mkdir failed");
}

static void dsk_scratch_link(const char *dir, const char *executable) {
    dsk_scratch_mkdir(dir);

    char link[320];
    snprintf(link, sizeof link, "%s/dsk_SFC", dir);
    unlink(link);
    assert(symlink(executable, link) == 0 && "dsk scratch: cannot link a build into place");
}

/* One run of the build in dir, with whatever the caller's tests use, and
   whether it left a results file behind. The file is removed either way. */
static int dsk_scratch_writes(const char *dir, const char *inputs, const char *run_name) {
    char command[8192];
    snprintf(command, sizeof command,
             "\"%s/dsk_SFC\" \"%s\" -r %s -s 1 -f 0 -c 0 -v 0 >/dev/null 2>&1",
             dir, inputs, run_name);
    if (system(command) != 0) return 0;

    char path[512];
    snprintf(path, sizeof path, "%s/output/results_%s_1.txt", dir, run_name);
    if (access(path, R_OK) != 0) return 0;

    unlink(path);
    snprintf(path, sizeof path, "%s/output/errors/Errors_%s_1.txt", dir, run_name);
    unlink(path);
    return 1;
}

/* Two sibling directories with the two builds linked into them, at a name
   length upstream tolerates. */
static DskScratch dsk_scratch_open(const char *tag, const char *run_name, const char *inputs,
                                   const char *upstream_binary, const char *modified_binary,
                                   int n_seeds) {
    const char *scratch = getenv("TMPDIR");
    if (!scratch) scratch = "/tmp";

    int digits = 1;
    for (int seed = n_seeds; seed >= 10; seed /= 10) digits++;

    DskScratch found;
    found.padding = -1;

    for (int pad = 0; pad < DSK_SCRATCH_ATTEMPTS; pad++) {
        char suffix[DSK_SCRATCH_ATTEMPTS + 1];
        for (int k = 0; k < pad; k++) suffix[k] = 'x';
        suffix[pad] = '\0';

        DskScratch tried;
        tried.padding = pad;
        snprintf(tried.root, sizeof tried.root, "%s/dsk_%s%s_%d", scratch, tag, suffix, (int)getpid());
        snprintf(tried.upstream, sizeof tried.upstream, "%s/u", tried.root);
        snprintf(tried.modified, sizeof tried.modified, "%s/m", tried.root);

        /* The longest name the model builds is the error log. */
        size_t longest = strlen(tried.upstream) + strlen("/output/errors/Errors_") +
                         strlen(run_name) + 1 + (size_t)digits + strlen(".txt") + 1;
        if (longest > DSK_UPSTREAM_FILENAME_BYTES) continue;

        dsk_scratch_mkdir(tried.root);
        dsk_scratch_link(tried.upstream, upstream_binary);
        dsk_scratch_link(tried.modified, modified_binary);

        if (dsk_scratch_writes(tried.upstream, inputs, run_name)) {
            found = tried;
            break;
        }
    }

    assert(found.padding >= 0 &&
           "dsk scratch: no scratch path of any length let the upstream build write its output");
    return found;
}

#endif
