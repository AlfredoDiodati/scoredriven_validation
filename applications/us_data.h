#ifndef US_DATA_H
#define US_DATA_H

#include <et_al./linalg/mat.h>
#include <et_al./frame/csv.h>
#include <stdio.h>

/*
The five variables and the sample, shared by every script under applications/
that needs them: log GDP, employment, log CPI, the interest rate and log
energy demand, over 1973Q1 to 2019Q4. GDP, CPI and energy demand are Fabiano's
non-stationary, co-integrated three and are 100 ln(level); employment and the
interest rate are his stationary two and are their own untransformed units.

This header does not read us_real.csv itself. applications/us_prepare_data.c
does that once, following papers/fabiano.txt's transformation, and writes
out/us_system.csv; this header only loads that CSV back into the K x T matrix
layout the rest of the project uses. Run `make app-us_prepare_data` before
anything that includes this header, which `make applications` and the
per-script `app-<name>` targets already do for you - see the Makefile's
ordering on APPLICATION_STEMS.

qvarma_data.txt (Blazsek, Escribano and Licht's own series) plays no part
here or anywhere downstream of it: it is a diagnostic file, read only by
applications/us_transformation_search.c, which checks a candidate
transformation of us_real.csv against a published one and nothing more.
*/

#define ESTIMATION_PERIODS 188

enum { LOG_GDP, EMPLOYMENT, LOG_CPI, INTEREST_RATE, LOG_ENERGY_DEMAND, N_VARIABLES };

static const char *variable_name[N_VARIABLES] = {
    "LogGDP", "Employment", "LogCPI", "InterestRate", "LogEnergyDemand"
};

/*
The five variables as an N_VARIABLES x ESTIMATION_PERIODS matrix, one column per
quarter, the layout the rest of the project uses, read from out/us_system.csv.
Caller must mat_free.
*/
static inline Mat load_us_system(void) {
    FILE *probe = fopen("out/us_system.csv", "r");
    assert(probe && "out/us_system.csv is missing; run make app-us_prepare_data first");
    fclose(probe);
    DataFrame system = df_read_csv("out/us_system.csv", csv_read_options_default());
    assert(system.r == ESTIMATION_PERIODS
           && "out/us_system.csv does not have ESTIMATION_PERIODS rows; rerun "
              "make app-us_prepare_data");

    Mat y = mat_new(N_VARIABLES, ESTIMATION_PERIODS);
    for (int k = 0; k < N_VARIABLES; k++) {
        Mat column = df_col_numeric(&system, variable_name[k]);
        for (int t = 0; t < ESTIMATION_PERIODS; t++) AT(y, k, t) = AT(column, t, 0);
    }
    df_free(&system);
    return y;
}

#endif /* US_DATA_H */
