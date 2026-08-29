#ifndef ABM_SYSTEM_H
#define ABM_SYSTEM_H

#include <et_al./linalg/mat.h>
#include <et_al./frame/rdata.h>
#include <math.h>

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
/* Periods in the raw file, before this header's own burn-in and differencing. */
static inline int abm_system_n_periods(const RData *d) { return _abm_system_array(d)->dim[1]; }

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
   an extra one - so the output below has exactly
   abm_system_n_periods(d) - ABM_SYSTEM_BURN_IN periods (400 of 600), not
   399, matching the collaborator's own iT exactly. */
#define ABM_SYSTEM_BURN_IN 200

/* One replicate of an already-read file, converted to ABM_SYSTEM_K x
   (abm_system_n_periods(d) - ABM_SYSTEM_BURN_IN) in the layout this
   header's own comment above describes. Caller must mat_free. */
static inline Mat abm_system_slice(const RData *d, int replicate) {
    RValue *v = _abm_system_array(d);
    int n_var = v->dim[0], n_obs = v->dim[1], n_rep = v->dim[2];
    assert(replicate >= 0 && replicate < n_rep && "abm_system: replicate out of range");
    assert(n_obs > ABM_SYSTEM_BURN_IN &&
           "abm_system: fewer periods in this file than the burn-in alone needs");

    int idx_gdp = _abm_system_variable_index(v, "GDP");
    int idx_emp = _abm_system_variable_index(v, "Employment rate");
    int idx_infl = _abm_system_variable_index(v, "Gross inflation");
    int idx_ir = _abm_system_variable_index(v, "Interest rate");
    int idx_energy = _abm_system_variable_index(v, "Energy");

    #define ABM_VAL(var, t) \
        ((double)v->v.real[(var) + (t) * n_var + replicate * n_var * n_obs])

    Mat y = mat_new(ABM_SYSTEM_K, n_obs - ABM_SYSTEM_BURN_IN);
    for (int t = ABM_SYSTEM_BURN_IN; t < n_obs; t++) {
        int c = t - ABM_SYSTEM_BURN_IN;
        AT(y, ROW_GDP_GROWTH, c) = (mreal)(100.0 * (log(ABM_VAL(idx_gdp, t))
                                                   - log(ABM_VAL(idx_gdp, t - 1))));
        AT(y, ROW_EN_GROWTH, c) = (mreal)(100.0 * (log(ABM_VAL(idx_energy, t))
                                                  - log(ABM_VAL(idx_energy, t - 1))));
        AT(y, ROW_EMPLOYMENT_CHANGE, c) = (mreal)(100.0 * (ABM_VAL(idx_emp, t)
                                                          - ABM_VAL(idx_emp, t - 1)));
        AT(y, ROW_INFLATION, c) = (mreal)(100.0 * (log(ABM_VAL(idx_infl, t))
                                                  - log(ABM_VAL(idx_infl, t - 1))));
        AT(y, ROW_INTEREST_RATE, c) = (mreal)(100.0 * ABM_VAL(idx_ir, t));
    }
    #undef ABM_VAL
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

#endif /* ABM_SYSTEM_H */
