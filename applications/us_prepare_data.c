/*
The one place us_real.csv is read and turned into the five-variable system the
rest of applications/ uses. qvarma_data.txt (Blazsek, Escribano and Licht's own
three-series file) is not read here or by anything downstream of this script: it
is a diagnostic file, used only by applications/us_transformation_search.c to
check a candidate transformation against a published one, and is never a source
for the variables the model is actually estimated on.

The five variables follow papers/fabiano.txt (Fabiano, "Evaluating Nonlinear
Simulation Models with Model Confidence Sets"), section 4.2, which uses this
same FRED-QD-sourced US data: GDP, employment, CPI, the interest rate and
energy demand. Only the three Fabiano's own ADF and Johansen results find
trend non-stationary and co-integrated - GDP, CPI and energy demand - are put
in log-levels times 100,

    x_var_t = 100 ln(x_t)

which is what keeps a co-integrating relation from being differenced away
(Sims et al. 1990; Cochrane 1997). Employment and the interest rate are
Fabiano's own two stationary variables and enter in their natural units,
untransformed: logging a variable that is already stationary is not needed to
preserve anything, and static_model.md's own state vector, which this system
descends from, keeps Emp_t and IR_t unlogged for exactly that reason.

Employment is not a column us_real.csv has; it is 100 minus Unemployment, the
same identity docs/DATA_DOCUMENTATION.md records for `_other/modified.csv`'s
Employement column. Energy demand uses the deseasonalised column,
Des_Energy_demand: the raw column has a pronounced quarterly pattern that a
level transformation does nothing to remove, and FRED-QD's own series are
seasonally adjusted to begin with.

Fabiano's own sample is 1973:Q2 to 2019:Q4, T = 188, post-2019 excluded to
avoid the COVID quarters. us_real.csv's row 0 is 1973Q1; 188 quarters from
there lands on 2019Q4 exactly, so this takes rows 0 to 187 rather than trying
to reproduce a Q2 start that does not fit that count - either the thesis
summary's quarter is off by one or its own source starts a quarter later than
this file, and the row count is the fact that is checkable here.

Writes out/us_system.csv, the five variables over that sample, with a Quarter
column so a reader does not have to count rows against this comment. Nothing
printed.
*/

#include <et_al./linalg/mat.h>
#include <et_al./frame/csv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RAW_PERIODS 193

/* 1973Q1 to 2019Q4, us_real.csv rows 0 to 187. */
#define ESTIMATION_PERIODS 188

enum { LOG_GDP, EMPLOYMENT, LOG_CPI, INTEREST_RATE, LOG_ENERGY_DEMAND, N_VARIABLES };

/* No spaces: these are column names in out/us_system.csv, and a name with a
   space in it is one more thing a reader has to quote correctly. */
static const char *variable_name[N_VARIABLES] = {
    "LogGDP", "Employment", "LogCPI", "InterestRate", "LogEnergyDemand"
};

/* Which of the five are logged, in the same order as the enum: GDP, CPI and
   energy demand are Fabiano's non-stationary, co-integrated three; employment
   and the interest rate are his stationary two and stay in their own units. */
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
    df_write_csv(&system, "out/us_system.csv", csv_write_options_default());

    df_free(&system);
    df_free(&real);
    return 0;
}
