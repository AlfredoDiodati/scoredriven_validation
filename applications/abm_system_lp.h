#ifndef ABM_SYSTEM_LP_H
#define ABM_SYSTEM_LP_H

#include "applications/abm_system.h"
#include "applications/lp_system.h"
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/*
One configuration's impulse responses under the three local projections, read
back from the archive applications/abm_system_fit_lp.c writes, for the loss
tables of applications/abm_system_lp_irf_loss.c against the US data and of
montecarlo/lp_irf_loss.c and montecarlo/lp_sweep.c against a simulated run.

The three responses of a replicate are irf_lin_mean, irf_s1_mean and
irf_s2_mean, each K x (hor + 1) x K, taken whole, horizon 0 included, as the R
pipeline takes them with as.vector. The order of the entries does not matter to
the loss, which is a mean over entry-by-entry differences between two vectors
laid out the same way.

Every stored row is checked against the replicate's own series before it is
used: the dataset is read and transformed again here and the row's fingerprints
compared, so a row fitted on other data counts as missing rather than being
read as if it belonged to this replicate.
*/

enum { LP_LOSS_LIN, LP_LOSS_S1, LP_LOSS_S2, LP_N_LOSSES };

static inline const char *lp_loss_name(int loss) {
    static const char *const name[LP_N_LOSSES] = { "lin", "s1", "s2" };
    return name[loss];
}

#define LP_RESPONSE_DIM LP_SYSTEM_RESPONSE_DIM

#define LP_FIT_DIR_DEFAULT "out/abm_system_fit_lp"
#define LP_INPUT_DIR_DEFAULT "dataset/abm_system"

/* One configuration, indexed by replicate. row[replicate] is the stored row, or
   NULL when the archive has none for it or the one it has was fitted on other
   data. ok[replicate * LP_N_LOSSES + loss] says whether that model's response
   exists and is finite; et_al's mcs refuses a loss matrix with a hole. */
typedef struct {
    int n_replicates;
    LpSystemRecord *records;
    const LpSystemRecord **row;
    unsigned char *ok;
} LpConfiguration;

static inline const mreal *lp_response(const LpSystemRecord *r, int loss) {
    return loss == LP_LOSS_LIN ? r->lin : loss == LP_LOSS_S1 ? r->s1 : r->s2;
}

static inline int _lp_all_finite(const mreal *values) {
    for (int i = 0; i < LP_RESPONSE_DIM; i++)
        if (MISNAN(values[i]) || MISINF(values[i])) return 0;
    return 1;
}

/* Replicates 0 to n_replicates - 1 of sample, their rows checked against the
   series of the given layout. A configuration with no usable archive comes
   back with every row NULL. */
static inline LpConfiguration lp_configuration_load(const char *fit_dir, const char *input_dir,
                                                    const char *sample, int n_replicates, LpLayout layout) {
    LpConfiguration c;
    c.n_replicates = n_replicates;
    c.row = calloc((size_t)n_replicates, sizeof(LpSystemRecord *));
    c.ok = calloc((size_t)n_replicates * LP_N_LOSSES, 1);
    assert(c.row && c.ok && "abm_system_lp: out of memory for a configuration");

    char path[640];
    snprintf(path, sizeof path, "%s/%s.npz", fit_dir, sample);
    int n_records = 0;
    c.records = lp_system_read_archive(path, lp_system_spec(), &n_records);

    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", input_dir, sample);
    int n_batches = (n_replicates + ABM_SYSTEM_BATCH - 1) / ABM_SYSTEM_BATCH;
    for (int batch = 0; batch < n_batches && c.records; batch++) {
        Mat block[ABM_SYSTEM_BATCH];
        int replicate[ABM_SYSTEM_BATCH];
        int count = abm_system_read_batch(dir, batch, block, replicate);
        for (int b = 0; b < count; b++) {
            int index = replicate[b];
            const LpSystemRecord *found = NULL;
            /* Rows are written in the dataset's order, so the row at the
               replicate's own index is almost always the one; the scan is the
               fallback for a configuration with a replicate missing. */
            if (index < n_records && c.records[index].replicate == index) found = &c.records[index];
            for (int i = 0; !found && i < n_records; i++)
                if (c.records[i].replicate == index) found = &c.records[i];
            if (found && index < n_replicates && found->has_responses) {
                Mat y, switching;
                lp_system_series(block[b], layout, &y, &switching);
                if (lp_system_record_matches(found, y, switching)) {
                    c.row[index] = found;
                    for (int loss = 0; loss < LP_N_LOSSES; loss++)
                        c.ok[(size_t)index * LP_N_LOSSES + loss] = (unsigned char)_lp_all_finite(lp_response(found, loss));
                }
                mat_free(y);
                mat_free(switching);
            }
            mat_free(block[b]);
        }
    }
    return c;
}

static inline int lp_configuration_ok(const LpConfiguration *c, int replicate, int loss) {
    return c->ok[(size_t)replicate * LP_N_LOSSES + loss];
}

static inline void lp_configuration_free(LpConfiguration *c) {
    free(c->records);
    free(c->row);
    free(c->ok);
}

typedef struct { char *name; int index; } LpSample;

static inline int _lp_compare_samples(const void *a, const void *b) {
    return ((const LpSample *)a)->index - ((const LpSample *)b)->index;
}

/* The configurations stored under fit_dir, one archive each, ordered by their
   trailing number, so column j of every table is the same configuration. */
static inline LpSample *lp_list_samples(const char *fit_dir, int *count) {
    DIR *handle = opendir(fit_dir);
    assert(handle && "abm_system_lp: cannot open the LP fits - run abm_system_fit_lp first");

    LpSample *samples = NULL;
    int n = 0, cap = 0;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        size_t length = strlen(entry->d_name);
        if (length < 5 || strcmp(entry->d_name + length - 4, ".npz") != 0) continue;
        const char *underscore = strrchr(entry->d_name, '_');
        assert(underscore && "abm_system_lp: an archive name has no trailing _<N>");
        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            LpSample *grown = realloc(samples, (size_t)cap * sizeof(LpSample));
            assert(grown && "abm_system_lp: out of memory listing configurations");
            samples = grown;
        }
        samples[n].name = malloc(length - 3);
        assert(samples[n].name && "abm_system_lp: out of memory copying a name");
        memcpy(samples[n].name, entry->d_name, length - 4);
        samples[n].name[length - 4] = 0;
        samples[n].index = atoi(underscore + 1);
        n++;
    }
    closedir(handle);
    qsort(samples, (size_t)n, sizeof(LpSample), _lp_compare_samples);
    *count = n;
    return samples;
}

static inline void lp_free_samples(LpSample *samples, int count) {
    for (int i = 0; i < count; i++) free(samples[i].name);
    free(samples);
}

/* One past the highest replicate index the dataset holds for sample, read off
   its batch file names and the last archive, so the tables' row count comes
   from the data rather than from a constant. */
static inline int lp_count_replicates(const char *input_dir, const char *sample) {
    char dir[560];
    snprintf(dir, sizeof dir, "%s/%s", input_dir, sample);
    DIR *handle = opendir(dir);
    assert(handle && "abm_system_lp: cannot open a configuration's dataset directory");
    int highest = -1;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int batch;
        if (sscanf(entry->d_name, "batch_%d.npz", &batch) == 1 && batch > highest) highest = batch;
    }
    closedir(handle);
    assert(highest >= 0 && "abm_system_lp: a configuration holds no archive");

    Mat block[ABM_SYSTEM_BATCH];
    int replicate[ABM_SYSTEM_BATCH];
    int count = abm_system_read_batch(dir, highest, block, replicate);
    int last = -1;
    for (int b = 0; b < count; b++) {
        if (replicate[b] > last) last = replicate[b];
        mat_free(block[b]);
    }
    return last + 1;
}

#endif /* ABM_SYSTEM_LP_H */
