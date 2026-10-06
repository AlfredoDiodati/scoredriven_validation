/*
How often the confidence set returns the benchmark's own configuration, for
any per-run output vector stored in the layout of montecarlo/response_cache.h:
four ints (model, configurations, runs, entries per run), then float32
[configuration][run][entry]. The raw response caches and the compressions
montecarlo/compressed_response_learn.py writes are both in that layout.

Each run in [first_benchmark, last_benchmark] of each configuration stands in
for the data in turn. The loss of a configuration on a row is the mean
absolute difference between that row's vector and the benchmark's, rows are
the runs [row_from, row_to) minus the benchmark's run, and the confidence set
is MCS_TR at alpha 0.05 with 2000 resamples, block length 1, bootstrap
variance, seed 123 stream 0: the protocol of montecarlo/sweep_grid.c, run on
a different vector.

Benchmarks run side by side, one per thread, each with its own loss matrix.
The file is mapped rather than read, so pass a file that fits in memory: a
mapped file larger than the free memory is reread from disk for every
benchmark.

usage: compressed_response_recovery vectors.f32 first_benchmark last_benchmark row_from row_to out.csv [rank_only]

With rank_only the confidence set is skipped and only the rank of the
benchmark's configuration by mean loss is recorded; in_set and set_size are
then -1. That takes seconds where the confidence set takes minutes.

Columns: configuration, benchmark run, whether the configuration is in the
set, its rank by mean loss (1 is the smallest), the set size, the
configuration with the smallest mean loss, the final p-value.
*/

#include <et_al./inference/mcs.h>
#include <et_al./frame/frame.h>
#include <cblas.h>
#include <fcntl.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char **argv) {
    assert((argc == 7 || argc == 8) && "usage: compressed_response_recovery vectors.f32 first_benchmark "
                                       "last_benchmark row_from row_to out.csv [rank_only]");
    int first_benchmark = atoi(argv[2]), last_benchmark = atoi(argv[3]);
    int row_from = atoi(argv[4]), row_to = atoi(argv[5]);
    int rank_only = argc == 8;
    openblas_set_num_threads(1);
    /* A confidence set inside a worker runs on that worker's thread alone. */
    omp_set_max_active_levels(1);

    int fd = open(argv[1], O_RDONLY);
    assert(fd >= 0 && "compressed_response_recovery: cannot open the vectors");
    struct stat st;
    fstat(fd, &st);
    const char *map = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    assert(map != MAP_FAILED);
    const int *header = (const int *)map;
    int n_cops = header[1], n_runs = header[2], dim = header[3];
    assert(st.st_size == 16 + (off_t)n_cops * n_runs * dim * (off_t)sizeof(float) &&
           "compressed_response_recovery: the file is shorter or longer than its header says");
    assert(row_to <= n_runs && last_benchmark < n_runs);
    const float *data = (const float *)(map + 16);

    const char **names = malloc((size_t)n_cops * sizeof(char *));
    char (*name_buffer)[16] = malloc((size_t)n_cops * sizeof *name_buffer);
    for (int c = 0; c < n_cops; c++) {
        snprintf(name_buffer[c], sizeof name_buffer[c], "cop_%04d", c + 1);
        names[c] = name_buffer[c];
    }

    int n_benchmarks = n_cops * (last_benchmark - first_benchmark + 1);
    int *in_set = malloc((size_t)n_benchmarks * sizeof(int));
    int *rank = malloc((size_t)n_benchmarks * sizeof(int));
    int *set_size = malloc((size_t)n_benchmarks * sizeof(int));
    int *lowest = malloc((size_t)n_benchmarks * sizeof(int));
    double *final_pvalue = malloc((size_t)n_benchmarks * sizeof(double));

    #pragma omp parallel
    {
        double *reference = malloc((size_t)dim * sizeof(double));
        double *column_mean = malloc((size_t)n_cops * sizeof(double));
        int *kept = malloc((size_t)(row_to - row_from) * sizeof(int));

        #pragma omp for schedule(dynamic, 1)
        for (int k = 0; k < n_benchmarks; k++) {
            int cop = k % n_cops, run = first_benchmark + k / n_cops;
            int n_kept = 0;
            for (int r = row_from; r < row_to; r++) if (r != run) kept[n_kept++] = r;
            const float *own = data + ((size_t)cop * n_runs + run) * dim;
            for (int i = 0; i < dim; i++) reference[i] = own[i];

            Mat losses = mat_new(n_kept, n_cops);
            for (int c = 0; c < n_cops; c++) {
                double total = 0;
                for (int row = 0; row < n_kept; row++) {
                    const float *cell = data + ((size_t)c * n_runs + kept[row]) * dim;
                    double sum = 0;
                    for (int i = 0; i < dim; i++) sum += fabs((double)cell[i] - reference[i]);
                    AT(losses, row, c) = (mreal)(sum / dim);
                    total += sum / dim;
                }
                column_mean[c] = total / n_kept;
            }

            int r = 1, low = 0;
            for (int c = 0; c < n_cops; c++) {
                if (column_mean[c] < column_mean[low]) low = c;
                if (c != cop && column_mean[c] < column_mean[cop]) r++;
            }
            rank[k] = r;
            lowest[k] = low;
            if (rank_only) {
                mat_free(losses);
                in_set[k] = -1;
                set_size[k] = -1;
                final_pvalue[k] = NAN;
                continue;
            }

            DataFrame frame = df_from_matrix(losses, names);
            mat_free(losses);
            MCSOptions opt = mcs_options_default();
            opt.bootstrap = 2000;
            opt.block_length = 1;
            opt.variance = MCS_VARIANCE_BOOTSTRAP;
            opt.stat = MCS_TR;
            MCSResult res = mcs(&frame, opt);
            in_set[k] = mcs_in_set(&res, cop);
            set_size[k] = res.n_surviving;
            final_pvalue[k] = res.final_pvalue;
            mcs_free(&res);
            df_free(&frame);
        }
        free(reference);
        free(column_mean);
        free(kept);
    }

    FILE *out = fopen(argv[6], "w");
    assert(out && "compressed_response_recovery: cannot open the output");
    fprintf(out, "cop,run,in_set,rank,set_size,lowest_loss_cop,final_pvalue\n");
    for (int k = 0; k < n_benchmarks; k++)
        fprintf(out, "%d,%d,%d,%d,%d,%d,%.6f\n", k % n_cops + 1, first_benchmark + k / n_cops, in_set[k], rank[k],
                set_size[k], lowest[k] + 1, final_pvalue[k]);
    fclose(out);

    free(in_set);
    free(rank);
    free(set_size);
    free(lowest);
    free(final_pvalue);
    free(name_buffer);
    free(names);
    munmap((void *)map, (size_t)st.st_size);
    close(fd);
    return 0;
}
