/*
Which shared resource the taped evaluation contends on, tested away from et_al
so the answer does not depend on anything the tape does: a tight loop of small
cblas_dgemm calls, a tight loop of malloc and free at the tape's block size,
and a loop of pure arithmetic, each at one and four threads. Not part of the
pipeline; run explicitly. Writes out/small_call_scaling.txt.

The shapes are the ones the qvarma filter issues: 5x5 times 5x1, and 5x5 times
5x5, several thousand per evaluation.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <omp.h>
#include <malloc.h>

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* The tape frees its blocks on every objective evaluation, and handing those
   back to the kernel costs more than the arithmetic does: 22 million minor
   page faults over eight fits, against 2,331 once the arena is kept resident.
   Set here rather than left to MALLOC_TRIM_THRESHOLD_ in the environment, so
   the numbers this writes do not depend on how it was launched. The pipeline
   itself does not set these yet. */
static void keep_the_arena_resident(void) {
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

int main(void) {
    keep_the_arena_resident();
    openblas_set_num_threads(1);
    FILE *report = fopen("out/small_call_scaling.txt", "w");
    /* Per job rather than shared, so each loop runs long enough to time: a
       dgemm call is three hundred times a malloc and five hundred times the
       arithmetic, and the same count for all three leaves two of them at a
       millisecond. */
    const long calls[3] = { 2000000, 20000000, 400000000 };
    fprintf(report, "%-26s %8s %10s %12s %10s %14s\n", "operation", "threads", "wall_s",
            "ns_per_op", "speedup", "ops_per_worker");

    for (int job = 0; job < 3; job++) {
        const char *name = job == 0 ? "cblas_dgemm 5x5 by 5x1"
                         : job == 1 ? "malloc+free 64 KiB" : "pure arithmetic";
        double baseline = 0;
        for (int c = 0; c < 2; c++) {
            int threads = c == 0 ? 1 : 4;
            double t0 = now_seconds();
            #pragma omp parallel num_threads(threads)
            {
                double a[25], b[5], out[5];
                for (int i = 0; i < 25; i++) a[i] = 0.01 * i;
                for (int i = 0; i < 5; i++) b[i] = 0.1 * i;
                volatile double sink = 0;
                long n = calls[job];
                if (job == 0) {
                    for (long i = 0; i < n; i++) {
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 5, 1, 5,
                                    1.0, a, 5, b, 1, 0.0, out, 1);
                        sink += out[0];
                    }
                } else if (job == 1) {
                    for (long i = 0; i < n; i++) {
                        void *p = malloc(1 << 16);
                        ((char*)p)[0] = (char)i;
                        sink += ((char*)p)[0];
                        free(p);
                    }
                } else {
                    double x = 1.000001;
                    for (long i = 0; i < n; i++) { x = x * 1.0000001 + 1e-9; sink += x; }
                }
                (void)sink;
            }
            double wall = now_seconds() - t0;
            if (c == 0) baseline = wall;
            fprintf(report, "%-26s %8d %10.3f %12.2f %10.2f %14ld\n", name, threads, wall,
                    1e9 * wall / (double)calls[job], threads * baseline / wall, calls[job]);
            fflush(report);
        }
    }
    fclose(report);
    return 0;
}
