#ifndef LP_SYSTEM_H
#define LP_SYSTEM_H

#include <et_al./linalg/mat.h>
#include <et_al./lp/lp.h>
#include <et_al./frame/npz.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
The series the local projections are fitted on, built from the five stored
series every dataset of this project holds (applications/abm_system.h's rows:
GDP_growth, EN_growth, Employment_change, Inflation, InterestRate).

lp_system_from_stored gives the fitted layout: the five stored series
themselves, the ones the t-QVARMA is fitted on, in the same order, and as the
series that decides the state the 4-period moving average of log GDP that the
collaborator's R pipeline builds. The first three periods are dropped, where
that moving average does not exist.

lp_system_r_transform is transform_data() of that R pipeline
(_temp/functions.R), the same arithmetic in the same order:

    GDP_growth        = log(init_GDP + cumsum(GDP_growth)) * 100
    Employment_change = init_Employment + cumsum(Employment_change)
    Inflation         = Inflation
    InterestRate      = InterestRate
    EN_growth         = log(init_Energy + cumsum(EN_growth)) * 100
    GDP_MA4           = rollmean(GDP_growth, k = 4, align = "right")

with the three periods where GDP_MA4 is not defined dropped from every series,
and the variables in R's column order: GDP, employment, inflation, interest
rate, energy. R keeps the input column names on the transformed series, so its
"GDP_growth" is a log level by then. Its input is this project's stored series
in their own units: the collaborator's data are this project's dataset.

The fits do not use R's levels. Fitted on them, the local projections
recovered the configuration a simulated run came from in 0 of 100 runs of the
Monte Carlo experiment, and fitted on the stored series in 90 of 100
(docs/MONTECARLO_LP_VALIDATION.md). The R layout is kept because the state
series comes from it and because it is what the collaborator's pipeline fits.
*/

/* The stored rows this header reads, abm_system.h's order, spelled out here so
   the header does not need abm_system.h and can be reached from the real-data
   side as well. They are also the rows of the fitted layout. */
enum { STORED_GDP_GROWTH, STORED_EN_GROWTH, STORED_EMPLOYMENT_CHANGE, STORED_INFLATION,
       STORED_INTEREST_RATE, STORED_K };
#define LP_SYSTEM_K STORED_K

/* The rows of R's layout. */
enum { LP_R_ROW_GDP, LP_R_ROW_EMPLOYMENT, LP_R_ROW_INFLATION, LP_R_ROW_INTEREST_RATE, LP_R_ROW_ENERGY };

#define LP_SYSTEM_INIT_GDP 100.0
#define LP_SYSTEM_INIT_ENERGY 100.0
#define LP_SYSTEM_INIT_EMPLOYMENT 0.90
/* zoo::rollmean's k in transform_data. */
#define LP_SYSTEM_SWITCHING_WINDOW 4

/* Periods the local projections see out of a stored block of T periods. */
static inline int lp_system_n_periods(int T) { return T - (LP_SYSTEM_SWITCHING_WINDOW - 1); }

/*
R's transform_data. stored is STORED_K x T. Fills y with the 5 x
lp_system_n_periods(T) block in R's layout and switching with the matching
lp_system_n_periods(T) x 1 GDP_MA4. Caller must mat_free both.

A cumulated growth that falls to -init or below has no logarithm and gives a
NaN, as it does in R. It is left in place rather than asserted on: lp.h reports
a NaN in y or in the switching series as var_ols_status -1 and computes
nothing, which a batch over simulated data can count.
*/
static inline void lp_system_r_transform(Mat stored, Mat *y, Mat *switching) {
    assert(stored.r == STORED_K && "lp_system: stored block has the wrong number of rows");
    int T = stored.c;
    int first = LP_SYSTEM_SWITCHING_WINDOW - 1;
    assert(T > first && "lp_system: fewer periods than the moving average needs");

    Mat cumulated = mat_cumsum(stored, 1);

    Mat indexed = mat_new(2, T);
    for (int t = 0; t < T; t++) {
        AT(indexed, 0, t) = (mreal)LP_SYSTEM_INIT_GDP + AT(cumulated, STORED_GDP_GROWTH, t);
        AT(indexed, 1, t) = (mreal)LP_SYSTEM_INIT_ENERGY + AT(cumulated, STORED_EN_GROWTH, t);
    }
    Mat logs = mat_log(indexed);
    Mat log_levels = mat_scale(logs, (mreal)100);

    Mat gdp = mat_slice(log_levels, 0, 1, 0, T);
    Mat gdp_ma4 = mat_rolling_mean(gdp, LP_SYSTEM_SWITCHING_WINDOW, 1);

    int n = lp_system_n_periods(T);
    *y = mat_new(5, n);
    *switching = mat_new(n, 1);
    for (int c = 0; c < n; c++) {
        int t = c + first;
        AT(*y, LP_R_ROW_GDP, c) = AT(log_levels, 0, t);
        AT(*y, LP_R_ROW_EMPLOYMENT, c) = (mreal)LP_SYSTEM_INIT_EMPLOYMENT
                                        + AT(cumulated, STORED_EMPLOYMENT_CHANGE, t);
        AT(*y, LP_R_ROW_INFLATION, c) = AT(stored, STORED_INFLATION, t);
        AT(*y, LP_R_ROW_INTEREST_RATE, c) = AT(stored, STORED_INTEREST_RATE, t);
        AT(*y, LP_R_ROW_ENERGY, c) = AT(log_levels, 1, t);
        AT(*switching, c, 0) = AT(gdp_ma4, 0, t);
    }

    mat_free(gdp_ma4);
    mat_free(log_levels);
    mat_free(logs);
    mat_free(indexed);
    mat_free(cumulated);
}

/*
The fitted layout. stored is STORED_K x T. Fills y with the stored series from
period 3 on, LP_SYSTEM_K x lp_system_n_periods(T) in the stored order, and
switching with R's GDP_MA4 over the same periods. Caller must mat_free both.
*/
static inline void lp_system_from_stored(Mat stored, Mat *y, Mat *switching) {
    Mat r_layout;
    lp_system_r_transform(stored, &r_layout, switching);
    mat_free(r_layout);
    int first = LP_SYSTEM_SWITCHING_WINDOW - 1;
    *y = mat_new(LP_SYSTEM_K, lp_system_n_periods(stored.c));
    for (int k = 0; k < LP_SYSTEM_K; k++)
        for (int t = first; t < stored.c; t++) AT(*y, k, t - first) = AT(stored, k, t);
}

/* Which series a fit uses: the stored ones (the default, and what the
   Monte Carlo experiment validated), or R's levels, the collaborator's own
   transformation. Everything a layout writes carries its suffix, so the two
   never share a file: out/abm_system_fit_lp and out/abm_system_fit_lp_r_levels,
   out/us_lp_fit.npz and out/us_lp_fit_r_levels.npz, and so on. */
typedef enum { LP_LAYOUT_GROWTH, LP_LAYOUT_R_LEVELS } LpLayout;

static inline const char *lp_layout_suffix(LpLayout layout) {
    return layout == LP_LAYOUT_R_LEVELS ? "_r_levels" : "";
}

static inline const char *lp_layout_description(LpLayout layout) {
    return layout == LP_LAYOUT_R_LEVELS
        ? "R's transform_data levels: GDP, Employment, Inflation, InterestRate, Energy"
        : "the stored series: GDP_growth, EN_growth, Employment_change, Inflation, InterestRate";
}

static inline void lp_system_series(Mat stored, LpLayout layout, Mat *y, Mat *switching) {
    if (layout == LP_LAYOUT_R_LEVELS) lp_system_r_transform(stored, y, switching);
    else lp_system_from_stored(stored, y, switching);
}

/* The command-line flags every local-projection script takes: --refit where
   it fits, --r-levels for R's layout. Anything else is a caller's mistake. */
static inline LpLayout lp_parse_flags(int argc, char **argv, int *force_refit) {
    LpLayout layout = LP_LAYOUT_GROWTH;
    if (force_refit) *force_refit = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--r-levels") == 0) layout = LP_LAYOUT_R_LEVELS;
        else if (force_refit && strcmp(argv[i], "--refit") == 0) *force_refit = 1;
        else assert(0 && "lp_system: unknown flag; the flags are --r-levels and, where fitting, --refit");
    }
    return layout;
}

/* The specification the R pipeline's estimate() calls lpirfs with: 4 lags, 15
   horizons, unit shocks; lp_nl with the switching series' HP cycle at lambda
   1600 through a logistic of slope 2, lagged one period (lpirfs' default
   lag_switching = TRUE). The estimator of Sigma_u does not change a unit
   shock. */
#define LP_SYSTEM_LAGS 4
#define LP_SYSTEM_HOR 15

static inline LpNlSpec lp_system_spec(void) {
    LpNlSpec spec;
    spec.lin = (LpSpec){ LP_SYSTEM_K, LP_SYSTEM_LAGS, LP_SYSTEM_HOR, LP_SHOCK_UNIT, VAR_SIGMA_ML };
    spec.lags_endog_nl = LP_SYSTEM_LAGS;
    spec.use_logistic = 1;
    spec.use_hp = 1;
    spec.lambda = 1600;
    spec.gamma = 2;
    spec.lag_switching = 1;
    return spec;
}

/*
How the fits are stored: one compressed .npz archive per configuration, one row
per replicate, written and read through et_al's df_write_npz_compressed and
df_read_npz. et_al's own LP cache (lp_lin_save_fit, lp_nl_save_fit) writes one
plain JSON file per fit, which for a million replicates is two million files
and about 33 GB of text; this layout holds the same numbers in binary.

A row holds what a later comparison reads and what is needed to trust it:

  replicate                   the replicate index
  has_responses               0 when the fit has no d, and then no responses
  var_ols_status, var_rank, var_chol_status
                              the VAR behind d, as LpNotes holds them
  lin_ols_status_<h>, nl_ols_status_<h>
                              ols's status at horizon h = 1..hor of each model,
                              above 0 where the rank-deficient fallback was used
  switching_is_constant       lp_nl's flag
  data_fingerprint, switching_fingerprint
                              lp_data_fingerprint of the fitted y and switching
                              series, so a reader can refuse a row computed
                              from other data
  spec_*                      the LpNlSpec the row was fitted under
  d_<i>                       d, K x K, row-major
  lin_<i>, s1_<i>, s2_<i>     irf_lin_mean, irf_s1_mean, irf_s2_mean, each
                              K x (hor + 1) x K flattened in the tensor's own
                              order: i = (k (hor + 1) + h) K + j for the
                              response of variable k at horizon h to shock j

fz, the state weight at every period, is not stored: it is a function of the
data and of lambda, gamma and lag_switching alone, and no comparison reads it.

The fingerprints are doubles, so the archive needs the float64 build, which is
the only one applications/ and montecarlo/ are compiled in.
*/

#define LP_SYSTEM_RESPONSE_DIM (LP_SYSTEM_K * (LP_SYSTEM_HOR + 1) * LP_SYSTEM_K)
#define LP_SYSTEM_N_SPEC 10

typedef struct {
    int replicate;
    int has_responses;
    int var_ols_status, var_rank, var_chol_status;
    int lin_ols_status[LP_SYSTEM_HOR], nl_ols_status[LP_SYSTEM_HOR];
    int switching_is_constant;
    double data_fingerprint, switching_fingerprint;
    mreal d[LP_SYSTEM_K * LP_SYSTEM_K];
    mreal lin[LP_SYSTEM_RESPONSE_DIM], s1[LP_SYSTEM_RESPONSE_DIM], s2[LP_SYSTEM_RESPONSE_DIM];
} LpSystemRecord;

static inline void _lp_system_spec_values(LpNlSpec spec, double *value) {
    const double values[LP_SYSTEM_N_SPEC] = {
        spec.lin.K, spec.lin.lags_endog_lin, spec.lags_endog_nl, spec.lin.hor, spec.lin.shock_type,
        spec.use_logistic, spec.use_hp, spec.lambda, spec.gamma, spec.lag_switching
    };
    memcpy(value, values, sizeof values);
}

/* The archive's column names, in the order a row is laid out. Caller frees
   each name and the array. */
static inline char **_lp_system_column_names(int *count) {
    static const char *const fixed[] = {
        "replicate", "has_responses", "var_ols_status", "var_rank", "var_chol_status",
        "switching_is_constant", "data_fingerprint", "switching_fingerprint"
    };
    static const char *const spec_name[LP_SYSTEM_N_SPEC] = {
        "spec_K", "spec_lags_endog_lin", "spec_lags_endog_nl", "spec_hor", "spec_shock_type",
        "spec_use_logistic", "spec_use_hp", "spec_lambda", "spec_gamma", "spec_lag_switching"
    };
    int n_fixed = (int)(sizeof fixed / sizeof fixed[0]);
    int n = n_fixed + 2 * LP_SYSTEM_HOR + LP_SYSTEM_N_SPEC + LP_SYSTEM_K * LP_SYSTEM_K
          + 3 * LP_SYSTEM_RESPONSE_DIM;
    char **names = (char **)malloc((size_t)n * sizeof(char *));
    assert(names && "lp_system: out of memory naming the archive's columns");
    char buffer[64];
    int at = 0;
    for (int i = 0; i < n_fixed; i++) names[at++] = frame_strdup(fixed[i]);
    for (int h = 1; h <= LP_SYSTEM_HOR; h++) {
        snprintf(buffer, sizeof buffer, "lin_ols_status_%d", h);
        names[at++] = frame_strdup(buffer);
    }
    for (int h = 1; h <= LP_SYSTEM_HOR; h++) {
        snprintf(buffer, sizeof buffer, "nl_ols_status_%d", h);
        names[at++] = frame_strdup(buffer);
    }
    for (int i = 0; i < LP_SYSTEM_N_SPEC; i++) names[at++] = frame_strdup(spec_name[i]);
    for (int i = 0; i < LP_SYSTEM_K * LP_SYSTEM_K; i++) {
        snprintf(buffer, sizeof buffer, "d_%d", i);
        names[at++] = frame_strdup(buffer);
    }
    static const char *const response_prefix[3] = { "lin", "s1", "s2" };
    for (int part = 0; part < 3; part++)
        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) {
            snprintf(buffer, sizeof buffer, "%s_%d", response_prefix[part], i);
            names[at++] = frame_strdup(buffer);
        }
    assert(at == n && "lp_system: the archive's column count is off");
    *count = n;
    return names;
}

static inline void _lp_system_free_names(char **names, int count) {
    for (int i = 0; i < count; i++) free(names[i]);
    free(names);
}

/* One replicate's row, from the fit lp_lin_and_nl returned on y and switching. */
static inline LpSystemRecord lp_system_record(const LpLinNlFit *fit, Mat y, Mat switching, int replicate) {
    LpSystemRecord r;
    memset(&r, 0, sizeof r);
    r.replicate = replicate;
    r.has_responses = fit->lin.d.d != NULL && fit->nl.d.d != NULL;
    r.var_ols_status = fit->lin.notes.var_ols_status;
    r.var_rank = fit->lin.notes.var_rank;
    r.var_chol_status = fit->lin.notes.var_chol_status;
    for (int h = 0; h < LP_SYSTEM_HOR; h++) {
        r.lin_ols_status[h] = fit->lin.notes.ols_status ? fit->lin.notes.ols_status[h] : 0;
        r.nl_ols_status[h] = fit->nl.notes.ols_status ? fit->nl.notes.ols_status[h] : 0;
    }
    r.switching_is_constant = fit->nl.switching_is_constant;
    r.data_fingerprint = lp_data_fingerprint(y);
    r.switching_fingerprint = lp_data_fingerprint(switching);
    if (r.has_responses) {
        memcpy(r.d, fit->lin.d.d, sizeof r.d);
        memcpy(r.lin, fit->lin.irf_lin_mean.d, sizeof r.lin);
        memcpy(r.s1, fit->nl.irf_s1_mean.d, sizeof r.s1);
        memcpy(r.s2, fit->nl.irf_s2_mean.d, sizeof r.s2);
    } else {
        for (int i = 0; i < LP_SYSTEM_K * LP_SYSTEM_K; i++) r.d[i] = (mreal)NAN;
        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) r.lin[i] = r.s1[i] = r.s2[i] = (mreal)NAN;
    }
    return r;
}

/* Whether a stored row was computed from exactly these series. */
static inline int lp_system_record_matches(const LpSystemRecord *r, Mat y, Mat switching) {
    return r->data_fingerprint == lp_data_fingerprint(y)
        && r->switching_fingerprint == lp_data_fingerprint(switching);
}

/*
Writes count rows fitted under spec to path, through a temporary file renamed
into place once complete, so an interrupted run never leaves a truncated
archive for the next one to read.
*/
static inline void lp_system_write_archive(const char *path, const LpSystemRecord *records, int count,
                                           LpNlSpec spec) {
    _Static_assert(sizeof(mreal) == sizeof(double), "lp_system: the archive needs the float64 build");
    assert(count >= 1 && "lp_system: an archive with no rows records nothing");
    int n_columns;
    char **names = _lp_system_column_names(&n_columns);
    double spec_value[LP_SYSTEM_N_SPEC];
    _lp_system_spec_values(spec, spec_value);

    Mat table = mat_new(count, n_columns);
    for (int row = 0; row < count; row++) {
        const LpSystemRecord *r = &records[row];
        mreal *out = &AT(table, row, 0);
        int at = 0;
        out[at++] = r->replicate;
        out[at++] = r->has_responses;
        out[at++] = r->var_ols_status;
        out[at++] = r->var_rank;
        out[at++] = r->var_chol_status;
        out[at++] = r->switching_is_constant;
        out[at++] = r->data_fingerprint;
        out[at++] = r->switching_fingerprint;
        for (int h = 0; h < LP_SYSTEM_HOR; h++) out[at++] = r->lin_ols_status[h];
        for (int h = 0; h < LP_SYSTEM_HOR; h++) out[at++] = r->nl_ols_status[h];
        for (int i = 0; i < LP_SYSTEM_N_SPEC; i++) out[at++] = spec_value[i];
        for (int i = 0; i < LP_SYSTEM_K * LP_SYSTEM_K; i++) out[at++] = r->d[i];
        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) out[at++] = r->lin[i];
        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) out[at++] = r->s1[i];
        for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) out[at++] = r->s2[i];
        assert(at == n_columns && "lp_system: a row does not fill the archive's columns");
    }

    DataFrame frame = df_from_matrix(table, (const char *const *)names);
    char partial[1024];
    int written = snprintf(partial, sizeof partial, "%s.partial", path);
    assert(written > 0 && (size_t)written < sizeof partial && "lp_system: archive path does not fit");
    df_write_npz_compressed(&frame, partial);
    int renamed = rename(partial, path);
    assert(renamed == 0 && "lp_system: could not move a finished archive into place");

    df_free(&frame);
    mat_free(table);
    _lp_system_free_names(names, n_columns);
}

/* The value of row in column *at, advancing *at to the next column. */
static inline mreal _lp_system_next(const Mat *column, int *at, int row) {
    mreal value = AT(column[*at], row, 0);
    (*at)++;
    return value;
}

/*
The rows of the archive at path, in the order stored, or NULL when there is no
archive there or it was written for a specification other than spec or in
another layout. A NULL is a cache miss, never a fault. Caller frees the result.
*/
static inline LpSystemRecord *lp_system_read_archive(const char *path, LpNlSpec spec, int *count) {
    *count = 0;
    FILE *probe = fopen(path, "rb");
    if (!probe) return NULL;
    fclose(probe);

    DataFrame frame = df_read_npz(path);
    int n_columns;
    char **names = _lp_system_column_names(&n_columns);
    Mat *column = (Mat *)malloc((size_t)n_columns * sizeof(Mat));
    assert(column && "lp_system: out of memory reading an archive");
    int layout_matches = frame.n_cols == n_columns && frame.r >= 1;
    for (int j = 0; j < n_columns && layout_matches; j++) {
        if (frame.columns[j].type != COL_NUMERIC || strcmp(frame.columns[j].name, names[j]) != 0)
            layout_matches = 0;
        else
            column[j] = df_col_numeric(&frame, names[j]);
    }

    double spec_value[LP_SYSTEM_N_SPEC];
    _lp_system_spec_values(spec, spec_value);
    int spec_start = 8 + 2 * LP_SYSTEM_HOR;
    for (int row = 0; row < frame.r && layout_matches; row++)
        for (int i = 0; i < LP_SYSTEM_N_SPEC; i++)
            if (AT(column[spec_start + i], row, 0) != spec_value[i]) layout_matches = 0;

    LpSystemRecord *records = NULL;
    if (layout_matches) {
        records = (LpSystemRecord *)malloc((size_t)frame.r * sizeof(LpSystemRecord));
        assert(records && "lp_system: out of memory reading an archive");
        for (int row = 0; row < frame.r; row++) {
            LpSystemRecord *r = &records[row];
            int at = 0;
            r->replicate = (int)_lp_system_next(column, &at, row);
            r->has_responses = (int)_lp_system_next(column, &at, row);
            r->var_ols_status = (int)_lp_system_next(column, &at, row);
            r->var_rank = (int)_lp_system_next(column, &at, row);
            r->var_chol_status = (int)_lp_system_next(column, &at, row);
            r->switching_is_constant = (int)_lp_system_next(column, &at, row);
            r->data_fingerprint = (double)_lp_system_next(column, &at, row);
            r->switching_fingerprint = (double)_lp_system_next(column, &at, row);
            for (int h = 0; h < LP_SYSTEM_HOR; h++) r->lin_ols_status[h] = (int)_lp_system_next(column, &at, row);
            for (int h = 0; h < LP_SYSTEM_HOR; h++) r->nl_ols_status[h] = (int)_lp_system_next(column, &at, row);
            at += LP_SYSTEM_N_SPEC;
            for (int i = 0; i < LP_SYSTEM_K * LP_SYSTEM_K; i++) r->d[i] = _lp_system_next(column, &at, row);
            for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) r->lin[i] = _lp_system_next(column, &at, row);
            for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) r->s1[i] = _lp_system_next(column, &at, row);
            for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) r->s2[i] = _lp_system_next(column, &at, row);
        }
        *count = frame.r;
    }

    free(column);
    _lp_system_free_names(names, n_columns);
    df_free(&frame);
    return records;
}

#endif /* LP_SYSTEM_H */
