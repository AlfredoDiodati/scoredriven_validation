/*
Where a DSK test runs the simulator, and removing that directory only when the
test passed.

Every test that runs the simulator gives it a directory of its own under
DSK_SCRATCH_ROOT, inside the repository, because the model writes its output
beside the executable it was invoked as and
two tests sharing a directory would read each other's files. A 600-period run
writes about 3 MB, a -f 1 run twelve files and 28 MB, and several of these tests
run the model tens of times, so the directories are not small and one per test
per run accumulates quickly.

They are removed on success and kept on failure. What a failed comparison leaves
behind is the two output files that differ, and deleting those would mean
rerunning the test by hand to see anything. So a test that fails says where its
directory is and stops touching it; the next passing run of that test does not
clean it up either, since the path carries the process id and belongs to the run
that made it.

The root is a path relative to the repository, which every test is run from, and
not TMPDIR or /tmp: nothing a test writes leaves the repository. Relative also
keeps the names short. The upstream build stores every output filename in a
64-byte buffer, and the repository's absolute path is 56 characters on its own,
so an absolute root would leave upstream no room to write anything.

Nothing here runs a shell, so no path is ever interpreted as a command.
*/

#ifndef DSK_SCRATCH_H
#define DSK_SCRATCH_H

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <ftw.h>
#include <sys/stat.h>

#define DSK_SCRATCH_ROOT "out/tests"

/* The directory every test's scratch directory goes in, created if missing. */
static const char *dsk_scratch_root(void) {
    if (mkdir("out", 0755) != 0) assert(errno == EEXIST && "dsk scratch: cannot create out/");
    if (mkdir(DSK_SCRATCH_ROOT, 0755) != 0) assert(errno == EEXIST && "dsk scratch: cannot create the root");
    return DSK_SCRATCH_ROOT;
}

/* Enough open descriptors for nftw to walk without reopening, few enough that
   a test running several of these at once cannot run out. */
#define DSK_SCRATCH_WALK_DEPTH 16

static int _dsk_scratch_remove_entry(const char *path, const struct stat *info,
                                     int type, struct FTW *ftw) {
    (void)info; (void)type; (void)ftw;
    return remove(path);
}

/* Depth first, so a directory is removed after what is inside it, and physical,
   so a symlink to a build is unlinked rather than followed into the repository.
   A path that is already gone is not an error. */
static void dsk_scratch_remove(const char *root) {
    if (!root || !*root) return;
    if (access(root, F_OK) != 0) return;
    if (nftw(root, _dsk_scratch_remove_entry, DSK_SCRATCH_WALK_DEPTH,
             FTW_DEPTH | FTW_PHYS) != 0)
        fprintf(stderr, "dsk scratch: could not remove %s\n", root);
}

/* Call once, on the way out, with whatever the test is about to return as its
   verdict. Keeps the directory and names it when the test failed, so the files
   that disagree are still there to look at. */
static void dsk_scratch_finish(const char *root, int passed) {
    if (passed) dsk_scratch_remove(root);
    else printf("  scratch kept for inspection: %s\n", root);
}

#endif
