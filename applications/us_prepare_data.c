/*
The one place dataset/us_real.csv is read and turned into the five-variable
system the rest of applications/ uses: GDP, employment, CPI, the interest rate
and energy demand.

GDP, CPI and energy demand trend upwards and are put in log-levels times 100,

    x_var_t = 100 ln(x_t)

to stabilise their variance. Employment and the interest rate are already
rates and enter in their own units. The models are fitted on first
differences of these levels, built later (build_block in
applications/us_qvarma_spec_choice.c and us_system_growth_block in
applications/us_data.h); this script writes the levels only.

Employment is not a column us_real.csv has; it is 100 minus Unemployment.
Energy demand uses the seasonally adjusted column, Des_Energy_demand: the raw
column has a pronounced quarterly pattern that a level transformation does
nothing to remove.

The sample is us_real.csv's rows 0 to 187, 1973Q1 to 2019Q4, 188 quarters.
Post-2019 is excluded to avoid the COVID quarters. Differencing then leaves 187
quarters, 1973Q2 to 2019Q4, the sample the t-QVARMA is estimated on.

Writes out/applications/us_system.csv, the five variables over that sample, with
a Quarter column so a reader does not have to count rows against this comment.
Nothing printed. docs/DATA_DOCUMENTATION.md describes both steps.
*/

#include <et_al./linalg/mat.h>
#include <et_al./frame/csv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define RAW_PERIODS 193

/* 1973Q1 to 2019Q4, us_real.csv rows 0 to 187. */
#define ESTIMATION_PERIODS 188

enum { LOG_GDP, EMPLOYMENT, LOG_CPI, INTEREST_RATE, LOG_ENERGY_DEMAND, N_VARIABLES };

/* No spaces: these are column names in out/applications/us_system.csv, and a name with a
   space in it is one more thing a reader has to quote correctly. */
static const char *variable_name[N_VARIABLES] = {
    "LogGDP", "Employment", "LogCPI", "InterestRate", "LogEnergyDemand"
};

/* Which of the five are logged, in the same order as the enum: GDP, CPI and
   energy demand trend upwards; employment and the interest rate are rates and
   stay in their own units. */
static const int is_logged[N_VARIABLES] = { 1, 0, 1, 0, 1 };

#define QUARTER_LABEL_SIZE 8

static void quarter_label(int year0, int q0, int index, char *out) {
    int total = (q0 - 1) + index;
    snprintf(out, QUARTER_LABEL_SIZE, "%dQ%d", year0 + total / 4, 1 + total % 4);
}

static DataFrame new_dated_dataframe(int rows, int year0, int q0) {
    DataFrame df = df_new(rows);
    char **labels = (char**)malloc((size_t)rows * sizeof(char*));
    for (int t = 0; t < rows; t++) {
        labels[t] = (char*)malloc(QUARTER_LABEL_SIZE);
        quarter_label(year0, q0, t, labels[t]);
    }
    df_add_string_col(&df, "Quarter", (const char *const *)labels);
    for (int t = 0; t < rows; t++) free(labels[t]);
    free(labels);
    return df;
}

int main(void) {
    DataFrame real = df_read_csv("dataset/us_real.csv", csv_read_options_default());
    assert(real.r == RAW_PERIODS && "us_real.csv is 193 rows");
    assert(ESTIMATION_PERIODS <= RAW_PERIODS && "the estimation sample runs past the file");

    Mat gdp = df_col_numeric(&real, "GDP");
    Mat cpi = df_col_numeric(&real, "Cpi");
    Mat unemployment = df_col_numeric(&real, "Unemployment");
    Mat fed_rate = df_col_numeric(&real, "Fed_rate");
    Mat energy = df_col_numeric(&real, "Des_Energy_demand");

    DataFrame system = new_dated_dataframe(ESTIMATION_PERIODS, 1973, 1);
    for (int k = 0; k < N_VARIABLES; k++) {
        Vec column = mat_new(ESTIMATION_PERIODS, 1);
        for (int t = 0; t < ESTIMATION_PERIODS; t++) {
            double level;
            switch (k) {
            case LOG_GDP: level = (double)AT(gdp, t, 0); break;
            case EMPLOYMENT: level = 100.0 - (double)AT(unemployment, t, 0); break;
            case LOG_CPI: level = (double)AT(cpi, t, 0); break;
            case INTEREST_RATE: level = (double)AT(fed_rate, t, 0); break;
            case LOG_ENERGY_DEMAND: level = (double)AT(energy, t, 0); break;
            default: level = 0; break;
            }
            column.d[t] = (mreal)(is_logged[k] ? 100.0 * log(level) : level);
        }
        df_add_numeric_col(&system, variable_name[k], column);
        mat_free(column);
    }
    df_write_csv(&system, "out/applications/us_system.csv", csv_write_options_default());

    df_free(&system);
    df_free(&real);
    return 0;
}
