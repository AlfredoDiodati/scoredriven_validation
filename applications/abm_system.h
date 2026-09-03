#ifndef ABM_SYSTEM_H
#define ABM_SYSTEM_H

#include <et_al./linalg/mat.h>
#include <et_al./frame/npz.h>
#include <et_al./frame/rdata.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

/*
Converts one replicate of one ABM-simulated dataset (an .Rdata file under
dataset/simulated/) into the same K x T layout
applications/us_qvarmad_employment_change.c's own build_block() produces
from the real data - GDP_growth, EN_growth, Employment_change, Inflation,
InterestRate - so the t-QVARMAd(1,1,2)/(1,1,4) auxiliary models this
project settled on can be fit to a simulated series with no change to
qvarma_d.h or to how a fit is called. Named for what the data is (the ABM's
own simulated economic system, converted into this project's own "system"
layout - us_data.h's own load_us_system is the real-data counterpart), not
for the model that eventually consumes it - this is a data-conversion
primitive, not a QVARMA of any kind, and does not fit anything or loop over
replicates or files itself. It exists to be called from inside whichever
script needs one (file, replicate) pair at a time -
applications/abm_system_extract.c converts every one of them up front.

Each .Rdata file holds one saved object ("estimation"), an 8-variable x
n_periods x n_replicates double array (dataset/simulated's own
EstimationSeriesSample1_<i>.Rdata, n_replicates = 108 as of this writing,
n_periods = 600). Named dim[0] entries: GDP, Consumption, Investment,
"Employment rate", "Gross inflation", "Interest rate", Emissions, Energy.
Consumption, Investment and Emissions are not used - the auxiliary models
this converts for do not have a slot for them. Variables are looked up by
name, not position, so this does not break if a future save orders them
differently.

Every transformation below was checked against the real data's own
equivalent series (out/us_system.csv, the same five, computed the same way
applications/us_qvarmad_employment_change.c's own build_block() does) before
being adopted - see the session that built this file for the actual
comparison (mean and sd of each of the five, real against simulated,
replicate 0 of EstimationSeriesSample1_1.Rdata): the same order of magnitude
throughout, which a wrong scale (log10 instead of ln, or a missed
proportion-to-percentage-point rescale) would not have produced. None of
this is guaranteed by any documentation of the ABM's own output - it is
inferred from the numbers, and stated here as exactly that, an inference:

  - GDP, Energy: raw levels (GDP grows monotonically from about 2.2e5 to
    1.4e6 over 600 periods in the one file checked; Energy does not trend
    the same way but is on a comparable raw scale, about 1.8e5 to 2.5e5) -
    100 ln(level), first-differenced, the identical transform
    applications/us_prepare_data.c's own header comment documents for the
    real data's LogGDP and LogEnergyDemand (100 ln(x_t), Fabiano's own
    convention, not this project's invention).
  - "Employment rate": a 0-1 proportion (checked range 0.86-1.00), not the
    real data's own 0-100 scale (Employment = 100 - Unemployment,
    us_prepare_data.c's own construction) - rescaled by 100 before
    differencing, so a one-point change in the simulated series means the
    same thing a one-point change in the real Employment series does.
  - "Gross inflation": NOT read as a per-period rate despite the name -
    checked range 1.26 to 24.6, rising essentially monotonically over 600
    periods in the one file checked, which a rate hovering near 1 could
    not do. Read instead as a price *level* (like the real data's own
    LogCPI), which a sustained compounding inflation process would produce
    exactly this shape from: 100 ln(level), first-differenced, the same
    transform as LogCPI.
  - "Interest rate": a decimal proportion (checked range 0-0.0221, i.e. 0
    to 2.21 per cent), not the real data's own percentage-point units
    (values like 8.49) - rescaled by 100, kept as a level, not
    differenced, matching InterestRate's own treatment in the real data.

Nothing here handles missing data, non-finite values, or a file whose
"estimation" object has a different shape than described above - asserts
fire instead. Extending this to tolerate those is a decision for whoever
writes the actual loop, once it is known whether any of the 100 files or
108 replicates per file actually need it.
*/

enum { ROW_GDP_GROWTH, ROW_EN_GROWTH, ROW_EMPLOYMENT_CHANGE, ROW_INFLATION, ROW_INTEREST_RATE,
      ABM_SYSTEM_K };

/* The model's own untransformed series each output row is built from, one
   per row and in the same order, so a caller fills this matrix with the
   indices it reads the result with. */
enum { LEVEL_GDP, LEVEL_ENERGY, LEVEL_EMPLOYMENT, LEVEL_PRICE, LEVEL_INTEREST };

/* levels is ABM_SYSTEM_K x n_obs of the model's own series in the units it
   writes them; start is the index of the first period reported, and period
   start - 1 anchors its first difference. Caller must mat_free the result.
   The .Rdata reader below and applications/abm_system_simulate.c both come
   through here. */
static inline Mat abm_system_transform(Mat levels, int start) {
    int n_obs = levels.c;
    assert(levels.r == ABM_SYSTEM_K && "abm_system: levels has the wrong number of rows");
    assert(start >= 1 && "abm_system: no period before the first one to anchor a difference on");
    assert(n_obs > start && "abm_system: fewer periods than the burn-in alone needs");

    Mat y = mat_new(ABM_SYSTEM_K, n_obs - start);
    for (int t = start; t < n_obs; t++) {
        int c = t - start;
        AT(y, ROW_GDP_GROWTH, c) = (mreal)(100.0 * (log((double)AT(levels, LEVEL_GDP, t))
                                                   - log((double)AT(levels, LEVEL_GDP, t - 1))));
        AT(y, ROW_EN_GROWTH, c) = (mreal)(100.0 * (log((double)AT(levels, LEVEL_ENERGY, t))
                                                  - log((double)AT(levels, LEVEL_ENERGY, t - 1))));
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = (mreal)(100.0 * ((double)AT(levels, LEVEL_EMPLOYMENT, t)
                                                          - (double)AT(levels, LEVEL_EMPLOYMENT, t - 1)));
        AT(y, ROW_INFLATION, c) = (mreal)(100.0 * (log((double)AT(levels, LEVEL_PRICE, t))
                                                  - log((double)AT(levels, LEVEL_PRICE, t - 1))));
        AT(y, ROW_INTEREST_RATE, c) = (mreal)(100.0 * (double)AT(levels, LEVEL_INTEREST, t));
    }
    return y;
}

/* Reads path once. Caller must abm_system_free when done slicing replicates
   out of it - the primitive half of the "read once, slice many" split
   rdata.h's own df_read_rdata/rdata_read pair already establishes, needed
   here because a caller looping over a file's own 108 replicates should not
   decompress and reparse the same file 108 times. */
static inline RData abm_system_read(const char *path) {
    RData d = rdata_read(path);
    assert(d.n == 1 && "abm_system: expected exactly one saved object per file");
    return d;
}

static inline void abm_system_free(RData *d) { rdata_free(d); }

static inline RValue *_abm_system_array(const RData *d) {
    RValue *v = rdata_get(d, "estimation");
    assert(v->type == RVALUE_REAL && v->ndim == 3 &&
           "abm_system: expected a 3D double array named 'estimation'");
    return v;
}

static inline int abm_system_n_replicates(const RData *d) { return _abm_system_array(d)->dim[2]; }

static inline int _abm_system_variable_index(const RValue *v, const char *name) {
    for (int i = 0; i < v->dim[0]; i++)
        if (strcmp(v->dimnames[0][i], name) == 0) return i;
    assert(0 && "abm_system: variable name not found in this file's own dimnames");
    return -1;
}

/* Periods 1..ABM_SYSTEM_BURN_IN (1-indexed, matching R's own convention)
   are discarded as transient before anything else - not this project's
   own choice, matched from a collaborator's own R pipeline for the same
   simulated data (`tau <- 600`, `vt <- 201:tau`, `iT <- length(vt)`,
   i.e. periods 201:600 kept, 400 of them). Period 200 (1-indexed), the
   one immediately before the kept window, is still used as the anchor
   for period 201's own first difference, the same way this project's own
   real-data series anchor their own first differenced period on the
   observation immediately before the reported sample rather than losing
   an extra one - so the output below has exactly 400 periods, not 399,
   matching the collaborator's own iT exactly.

   A file does not have to carry the transient for this to work. The files
   dataset/simulated/ was populated with hold every period from the first
   and carry no period labels; a file written with the transient already
   dropped labels its second dimension with the model periods it does hold,
   and both give the same 400 periods here. */
#define ABM_SYSTEM_BURN_IN 200

/* The model period this file's first stored period is, 1-indexed. Absent
   labels mean the file starts at period 1. */
static inline int _abm_system_first_period(const RValue *v) {
    if (v->dimnames && v->dimnames[1]) return atoi(v->dimnames[1][0]);
    return 1;
}

/* The index of the first period this header reports, which is also the
   first one past the anchor. */
static inline int _abm_system_slice_start(const RValue *v) {
    int start = ABM_SYSTEM_BURN_IN - (_abm_system_first_period(v) - 1);
    assert(start >= 1 && "abm_system: this file starts too late to anchor its own first difference");
    return start;
}

/* Periods one slice of this file returns, after the burn-in and the
   differencing. */
static inline int abm_system_n_periods(const RData *d) {
    const RValue *v = _abm_system_array(d);
    return v->dim[1] - _abm_system_slice_start(v);
}

/* One replicate of an already-read file, converted to ABM_SYSTEM_K x
   abm_system_n_periods(d) in the layout this
   header's own comment above describes. Caller must mat_free. */
static inline Mat abm_system_slice(const RData *d, int replicate) {
    RValue *v = _abm_system_array(d);
    int n_var = v->dim[0], n_obs = v->dim[1], n_rep = v->dim[2];
    int start = _abm_system_slice_start(v);
    assert(replicate >= 0 && replicate < n_rep && "abm_system: replicate out of range");
    assert(n_obs > start &&
           "abm_system: fewer periods in this file than the burn-in alone needs");

    int index[ABM_SYSTEM_K];
    index[LEVEL_GDP] = _abm_system_variable_index(v, "GDP");
    index[LEVEL_ENERGY] = _abm_system_variable_index(v, "Energy");
    index[LEVEL_EMPLOYMENT] = _abm_system_variable_index(v, "Employment rate");
    index[LEVEL_PRICE] = _abm_system_variable_index(v, "Gross inflation");
    index[LEVEL_INTEREST] = _abm_system_variable_index(v, "Interest rate");

    Mat levels = mat_new(ABM_SYSTEM_K, n_obs);
    for (int series = 0; series < ABM_SYSTEM_K; series++)
        for (int t = 0; t < n_obs; t++)
            AT(levels, series, t) =
                v->v.real[index[series] + t * n_var + replicate * n_var * n_obs];

    Mat y = abm_system_transform(levels, start);
    mat_free(levels);
    return y;
}

/* One-shot convenience for a caller that wants exactly one replicate from
   exactly one file: read, slice, free, done - the same shape
   rdata.h's own df_read_rdata gives a caller that does not want to manage
   RData's own lifetime. A caller processing many replicates from the same
   file should use abm_system_read/abm_system_slice/abm_system_free
   directly instead, to pay the decompress-and-parse cost once per file
   rather than once per replicate. */
static inline Mat abm_system_block(const char *path, int replicate) {
    RData d = abm_system_read(path);
    Mat y = abm_system_slice(&d, replicate);
    abm_system_free(&d);
    return y;
}


/* How the transformed replicates are stored, which both writers -
   applications/abm_system_extract.c from the .Rdata files and
   applications/abm_system_simulate.c from the model itself - and every reader
   go through, so the layout is written down once.

   Replicates are grouped ABM_SYSTEM_BATCH to an archive. Grouping is purely
   about size: a deflate stream that sees ten replicates compresses better
   than ten streams that each see one, and the gain is nearly all of what
   compressing a whole configuration together would give. The cost is that a
   crash loses the batch being filled rather than the single run in flight.

   Every archive carries a "replicate" column beside the five series, holding
   the index of the replicate each row belongs to, so an archive says what is
   in it rather than leaving that to be inferred from its name. A batch whose
   runs did not all succeed is short, and reading it still resolves each
   replicate correctly. */
#define ABM_SYSTEM_BATCH 10

/* The nine structural parameters of the ABM this project varies, under the
   names the model's own JSON uses for them. applications/abm_system_design.c
   samples them and applications/abm_system_simulate.c writes them back into
   that JSON, and both reach the list here so a design column and the field it
   ends up in cannot drift apart. The bounds are not here: only the design
   needs them. */
#define ABM_SYSTEM_N_PARAMETERS 9

static inline const char *const *abm_system_parameter_names(void) {
    static const char *const name[ABM_SYSTEM_N_PARAMETERS] = {
        "Gamma", "chi", "psi1", "psi3", "alfa", "taylor1", "taylor2", "taylor", "kappa"
    };
    return name;
}

static inline const char *const *abm_system_column_names(void) {
    static const char *const name[ABM_SYSTEM_K + 1] = {
        "GDP_growth", "EN_growth", "Employment_change", "Inflation", "InterestRate",
        "replicate"
    };
    return name;
}

static inline void abm_system_batch_path(char *out, size_t n, const char *dir, int replicate) {
    int written = snprintf(out, n, "%s/batch_%03d.npz", dir, replicate / ABM_SYSTEM_BATCH);
    assert(written > 0 && (size_t)written < n && "abm_system: batch path does not fit");
}

/* One archive from count blocks, each ABM_SYSTEM_K x periods in the layout
   abm_system_transform returns, stacked in the order given. Blocks are laid
   out one after another down the rows, which is what puts a whole column of
   one series into one deflate stream. */
static inline void abm_system_write_batch(const char *dir, const Mat *block,
                                          const int *replicate, int count) {
    assert(count >= 1 && "abm_system: an archive with no blocks records nothing");
    int periods = block[0].c;

    Mat table = mat_new(periods * count, ABM_SYSTEM_K + 1);
    for (int b = 0; b < count; b++) {
        assert(block[b].r == ABM_SYSTEM_K && block[b].c == periods &&
               "abm_system: blocks in one archive must have the same shape");
        for (int t = 0; t < periods; t++) {
            for (int k = 0; k < ABM_SYSTEM_K; k++)
                AT(table, b * periods + t, k) = AT(block[b], k, t);
            AT(table, b * periods + t, ABM_SYSTEM_K) = (mreal)replicate[b];
        }
    }

    char path[512];
    abm_system_batch_path(path, sizeof path, dir, replicate[0]);

    DataFrame df = df_from_matrix(table, abm_system_column_names());
    df_write_npz_compressed(&df, path);

    df_free(&df);
    mat_free(table);
}

/* One replicate, ABM_SYSTEM_K x periods, out of the archive holding it.
   Caller must mat_free. */
static inline Mat abm_system_read_replicate(const char *dir, int replicate) {
    char path[512];
    abm_system_batch_path(path, sizeof path, dir, replicate);

    DataFrame df = df_read_npz(path);
    Mat index = df_col_numeric(&df, abm_system_column_names()[ABM_SYSTEM_K]);

    int start = -1, periods = 0;
    for (int t = 0; t < df.r; t++)
        if ((int)AT(index, t, 0) == replicate) {
            if (start < 0) start = t;
            periods++;
        }
    assert(start >= 0 && "abm_system: this archive does not hold that replicate");

    Mat y = mat_new(ABM_SYSTEM_K, periods);
    for (int k = 0; k < ABM_SYSTEM_K; k++) {
        Mat column = df_col_numeric(&df, abm_system_column_names()[k]);
        for (int t = 0; t < periods; t++) AT(y, k, t) = AT(column, start + t, 0);
    }

    df_free(&df);
    return y;
}

/* Every replicate index stored under dir, ascending. Archives are visited in
   batch order rather than in the order readdir returns them, so a caller can
   rely on the order without sorting, and a batch that failed entirely is
   simply absent rather than a hole to step over. Caller must free. */
static inline int *abm_system_list_replicates(const char *dir, int *count) {
    DIR *handle = opendir(dir);
    assert(handle && "abm_system: cannot open a sample directory");

    int highest = -1;
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        int batch;
        if (sscanf(entry->d_name, "batch_%d.npz", &batch) == 1 && batch > highest) highest = batch;
    }
    closedir(handle);

    int *replicate = NULL, n = 0, cap = 0;
    for (int batch = 0; batch <= highest; batch++) {
        char path[512];
        int written = snprintf(path, sizeof path, "%s/batch_%03d.npz", dir, batch);
        assert(written > 0 && (size_t)written < sizeof path && "abm_system: batch path does not fit");

        struct stat info;
        if (stat(path, &info) != 0) continue;

        DataFrame df = df_read_npz(path);
        Mat index = df_col_numeric(&df, abm_system_column_names()[ABM_SYSTEM_K]);
        int previous = -1;
        for (int t = 0; t < df.r; t++) {
            int here = (int)AT(index, t, 0);
            if (here == previous) continue;
            if (n == cap) { cap = cap ? cap * 2 : 32; replicate = (int*)realloc(replicate, (size_t)cap * sizeof(int)); }
            replicate[n++] = here;
            previous = here;
        }
        df_free(&df);
    }

    *count = n;
    return replicate;
}

#endif /* ABM_SYSTEM_H */
