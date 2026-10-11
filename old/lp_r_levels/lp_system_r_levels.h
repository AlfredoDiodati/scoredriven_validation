#ifndef LP_SYSTEM_R_LEVELS_H
#define LP_SYSTEM_R_LEVELS_H

/*
Discarded from both pipelines: the local projections fitted on R's levels
instead of on the five series the t-QVARMA is fitted on. The main
transformation is kept because it loses less information; on the Monte Carlo
experiment the levels recovered the configuration a simulated run came from in
0 of 100 runs of cop_0191 under the linear model, against 90 of 100 on the
main transformation. docs/OLD_EXPERIMENTS.md.

transform_data() of the collaborator's R pipeline (_temp/functions.R), the same
arithmetic in the same order:

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

Built on applications/lp_system.h, which keeps the state series GDP_MA4 for the
fitted series; this header adds the levels.
*/

#include "applications/lp_system.h"

/* The rows of R's layout. */
enum { LP_R_ROW_GDP, LP_R_ROW_EMPLOYMENT, LP_R_ROW_INFLATION, LP_R_ROW_INTEREST_RATE, LP_R_ROW_ENERGY };

#define LP_SYSTEM_INIT_ENERGY 100.0
#define LP_SYSTEM_INIT_EMPLOYMENT 0.90

#define LP_R_LEVELS_SERIES_DESCRIPTION "R's transform_data levels: GDP, Employment, Inflation, InterestRate, Energy"

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

#endif /* LP_SYSTEM_R_LEVELS_H */
