/*
The parameter design every simulation run is indexed by: 1000 configurations
of the nine DSK parameters this project varies, drawn as a Latin hypercube
over the box below and written to dataset/abm_system_design.csv, one row per
configuration and one column per parameter.

Latin hypercube rather than independent uniform sampling because every
one-dimensional projection is then perfectly stratified: with 1000 points each
parameter takes 1000 distinct values, one in each 1/1000 interval of its
range, where an independent sample would leave some intervals empty. Rather
than a grid because three levels in nine dimensions is 3^9 = 19683 points,
twenty times this budget, and still only three distinct values per parameter.

The design is a deterministic function of the seed, the bounds and et_al's
own sampler, so it is reproducible rather than something to be carried
between machines. It is not reproducible against R's lhs package: both draw
from the same distribution, but the generators differ, so the same seed does
not give the same points.

An existing design is kept rather than overwritten. The row index is the
identity every stored replication, every fit and every confidence set entry
is named by, so regenerating after any simulation has run would renumber all
of them against results already on disk. Pass --force to overwrite anyway.
*/

#include "abm_system.h"

#include <et_al./frame/csv.h>
#include <et_al./random/lhs.h>
#include <et_al./random/random.h>

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define DESIGN_PATH "dataset/abm_system_design.csv"
#define N_COP 1000
#define DESIGN_SEED 1

/* Same order as abm_system_parameter_names(). taylor1 starting at 1 is where
   the Taylor principle binds and taylor inside (0,1) is an interest-rate
   smoothing coefficient; the rest are the model's own and the bounds are all
   this file asserts about them. */
static const double lower_bound[ABM_SYSTEM_N_PARAMETERS] = {
    0.05, -1.50, 0.10, 0.10, 0.10, 1.00, 0.00, 0.50, 0.50
};
static const double upper_bound[ABM_SYSTEM_N_PARAMETERS] = {
    0.25, -1.25, 0.50, 0.50, 0.50, 1.50, 0.50, 0.95, 0.95
};

int main(int argc, char **argv) {
    int force = argc > 1 && strcmp(argv[1], "--force") == 0;

    if (!force && access(DESIGN_PATH, F_OK) == 0) {
        FILE *note = fopen("out/abm_system_design.txt", "w");
        assert(note && "abm_system_design: cannot open the report path");
        fprintf(note, "%s already exists and was kept.\n", DESIGN_PATH);
        fprintf(note, "Run with --force to draw a new design, which renumbers every "
                      "configuration against results already on disk.\n");
        fclose(note);
        return 0;
    }

    Rng rng = rng_new(DESIGN_SEED, 0);
    Mat unit_design = lhs_random(&rng, N_COP, ABM_SYSTEM_N_PARAMETERS);

    Mat lower = mat_new(1, ABM_SYSTEM_N_PARAMETERS), upper = mat_new(1, ABM_SYSTEM_N_PARAMETERS);
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
        AT(lower, 0, p) = (mreal)lower_bound[p];
        AT(upper, 0, p) = (mreal)upper_bound[p];
    }

    Mat design = lhs_scale(unit_design, lower, upper);
    DataFrame df = df_from_matrix(design, abm_system_parameter_names());
    df_write_csv(&df, DESIGN_PATH, csv_write_options_default());

    /* What the design is, recorded next to it rather than left to be
       recovered from the numbers. */
    FILE *report = fopen("out/abm_system_design.txt", "w");
    assert(report && "abm_system_design: cannot open the report path");
    fprintf(report, "%s: %d configurations, Latin hypercube, seed %d\n\n",
            DESIGN_PATH, N_COP, DESIGN_SEED);
    fprintf(report, "%-10s %10s %10s %14s %14s\n", "parameter", "lower", "upper", "min drawn", "max drawn");
    for (int p = 0; p < ABM_SYSTEM_N_PARAMETERS; p++) {
        mreal smallest = AT(design, 0, p), largest = AT(design, 0, p);
        for (int i = 1; i < N_COP; i++) {
            if (AT(design, i, p) < smallest) smallest = AT(design, i, p);
            if (AT(design, i, p) > largest) largest = AT(design, i, p);
        }
        fprintf(report, "%-10s %10.4f %10.4f %14.6f %14.6f\n",
                abm_system_parameter_names()[p], lower_bound[p], upper_bound[p],
                (double)smallest, (double)largest);
    }
    fclose(report);

    df_free(&df);
    mat_free(design);
    mat_free(unit_design);
    mat_free(lower);
    mat_free(upper);
    return 0;
}
