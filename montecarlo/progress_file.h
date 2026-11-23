#ifndef MONTECARLO_PROGRESS_FILE_H
#define MONTECARLO_PROGRESS_FILE_H

/*
The progress file of a grid after the machine went down without a clean
shutdown. Rows still in the page cache are lost, and the file system can keep
the file's new length while the blocks behind it read as zero bytes, so the
file ends in a run of NUL bytes, possibly with a row cut in half before it.
A cut row read as done would keep its benchmark from being redone and leave a
short row in the result, and rows appended after the NUL bytes would carry them
into the result.

progress_file_drop_torn_tail cuts the file at the start of the line holding its
first NUL byte, or after its last newline when it holds none, so every row left
is whole and the benchmarks of the dropped rows are run again. A file left
empty is removed, so the grid treats it as absent. Returns the bytes dropped.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static inline long progress_file_drop_torn_tail(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) return 0;
    FILE *f = fopen(path, "rb");
    assert(f && "progress_file_drop_torn_tail: the progress file exists and cannot be read");
    static char chunk[1 << 16];
    long offset = 0, after_last_newline = 0;
    size_t got;
    int nul_found = 0;
    while (!nul_found && (got = fread(chunk, 1, sizeof chunk, f)) > 0) {
        for (size_t i = 0; i < got; i++) {
            if (chunk[i] == '\0') { nul_found = 1; break; }
            if (chunk[i] == '\n') after_last_newline = offset + (long)i + 1;
        }
        offset += (long)got;
    }
    fclose(f);
    long dropped = (long)st.st_size - after_last_newline;
    if (dropped == 0) return 0;
    if (after_last_newline == 0) {
        int removed = remove(path) == 0;
        assert(removed && "progress_file_drop_torn_tail: cannot remove the empty progress file");
        (void)removed;
    } else {
        FILE *g = fopen(path, "r+b");
        int cut = g && ftruncate(fileno(g), (off_t)after_last_newline) == 0;
        if (g) fclose(g);
        assert(cut && "progress_file_drop_torn_tail: cannot truncate the progress file");
        (void)cut;
    }
    fprintf(stderr, "%s: dropped its last %ld bytes, a row cut short or NUL bytes after an unclean shutdown\n", path,
            dropped);
    return dropped;
}

#endif
