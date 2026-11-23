/*
Whether applications/lp_system.h fits and stores the state-dependent local
projections under each of the four state weights (LpSystemStateWeight): the
two-sided HP cycle at lambda 1600 of lpirfs, the causal cycle at 1600, and
each of the two with lambda estimated on the replicate's switching series.

A weight that silently reached the wrong et_al function, or a lambda that was
estimated and then not used, would change every state-dependent response of
a whole tree of fits while every downstream program kept reading it as the
weight its directory names. So would an archive of one weight read back as
another.

Data: four replicates of 188 stored periods (the US sample length), each built
by stored_block below from rng_new(2027, i), the fourth with a NaN in GDP
growth at period 50, so that its switching series is NaN from there on.

1. Each weight's rows against et_al called directly: lambda is 1600, or
   exactly hp_estimate's lambda on the same switching series under the default
   search; the response tensors and d are bit for bit those of lp_lin_and_nl
   or lp_lin_and_nl_causal at that lambda; causal_hp and hp_estimate's two
   flags are recorded as they came back.
2. The NaN replicate: no responses, and under an estimated lambda the row
   records lambda NaN and both flags 0 rather than reaching hp_estimate,
   which asserts on a NaN.
3. The weights differ where they should: the causal and the two-sided state 1
   responses are not equal, nor those at 1600 and at the estimate.
4. The archive: every field of every row read back bit for bit under its own
   weight, and a miss under each of the three others. An archive of lpirfs'
   weight has exactly the columns the archives stored before the other weights
   existed. An archive of any other weight stores no linear projection and no
   d: its rows read back with has_linear 0, lin_ols_status 0 and d and lin NaN,
   every other field bit for bit, and it has exactly lpirfs' columns less
   lin_ols_status_<h>, d_<i> and lin_<i>, plus the weight's own. An
   estimated-lambda archive reads whatever spec.lambda is and is a miss under
   another search tolerance; a given-lambda archive is a miss under another
   lambda.

Writes out/tests/lp_system_state_weight_correctness_report.txt, and a scratch
archive under out/tests/ that it removes on the way out. Exits 1 on any
failure.
*/

#include "applications/lp_system.h"
#include <et_al./random/random.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define REPORT "out/tests/lp_system_state_weight_correctness_report.txt"
#define SCRATCH_ARCHIVE "out/tests/lp_system_state_weight_correctness_scratch.npz"
#define T_STORED 188
#define N_REPLICATES 4
#define NAN_REPLICATE 3
#define N_WEIGHTS 4

static int failures = 0;

static void check(FILE *report, int ok, const char *what) {
    fprintf(report, "%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) failures++;
}

/* The stored scales of tests/lp_system_transform.c: GDP and energy growth near
   0.5 and 0.1 with unit noise, employment changes near 0, inflation near 0.5,
   the rate near 1. */
static Mat stored_block(Rng *rng, int T) {
    Mat stored = mat_new(STORED_K, T);
    for (int t = 0; t < T; t++) {
        AT(stored, STORED_GDP_GROWTH, t) = (mreal)(0.5 + rng_normal(rng));
        AT(stored, STORED_EN_GROWTH, t) = (mreal)(0.1 + rng_normal(rng));
        AT(stored, STORED_EMPLOYMENT_CHANGE, t) = (mreal)(0.3 * rng_normal(rng));
        AT(stored, STORED_INFLATION, t) = (mreal)(0.5 + 0.2 * rng_normal(rng));
        AT(stored, STORED_INTEREST_RATE, t) = (mreal)(1.0 + 0.3 * rng_normal(rng));
    }
    return stored;
}

static int same_bits(double a, double b) { return memcmp(&a, &b, sizeof a) == 0; }

/* Every field but the linear projection's, the doubles and the response
   tensors by bit pattern, so a NaN compares equal to itself and finite values
   to exact equality. Padding is not compared, since a row read back is not
   zeroed first. */
static int state_fields_identical(const LpSystemRecord *a, const LpSystemRecord *b) {
    return a->replicate == b->replicate && a->has_responses == b->has_responses
        && a->var_ols_status == b->var_ols_status && a->var_rank == b->var_rank
        && a->var_chol_status == b->var_chol_status && a->switching_is_constant == b->switching_is_constant
        && memcmp(a->nl_ols_status, b->nl_ols_status, sizeof a->nl_ols_status) == 0
        && same_bits(a->data_fingerprint, b->data_fingerprint)
        && same_bits(a->switching_fingerprint, b->switching_fingerprint)
        && same_bits(a->lambda, b->lambda) && a->causal_hp == b->causal_hp
        && a->lambda_is_converged == b->lambda_is_converged && a->lambda_at_bound == b->lambda_at_bound
        && memcmp(a->s1, b->s1, sizeof a->s1) == 0 && memcmp(a->s2, b->s2, sizeof a->s2) == 0;
}

static int linear_fields_identical(const LpSystemRecord *a, const LpSystemRecord *b) {
    return a->has_linear == b->has_linear
        && memcmp(a->lin_ols_status, b->lin_ols_status, sizeof a->lin_ols_status) == 0
        && memcmp(a->d, b->d, sizeof a->d) == 0 && memcmp(a->lin, b->lin, sizeof a->lin) == 0;
}

/* What a row read from an archive without the linear projection holds there. */
static int linear_fields_absent(const LpSystemRecord *r) {
    int absent = r->has_linear == 0;
    for (int h = 0; h < LP_SYSTEM_HOR; h++) absent = absent && r->lin_ols_status[h] == 0;
    for (int i = 0; i < LP_SYSTEM_K * LP_SYSTEM_K; i++) absent = absent && MISNAN(r->d[i]);
    for (int i = 0; i < LP_SYSTEM_RESPONSE_DIM; i++) absent = absent && MISNAN(r->lin[i]);
    return absent;
}

/* Whether the archive has a column whose name starts with prefix. */
static int has_column_prefix(const DataFrame *frame, const char *prefix) {
    for (int j = 0; j < frame->n_cols; j++)
        if (strncmp(frame->columns[j].name, prefix, strlen(prefix)) == 0) return 1;
    return 0;
}

int main(void) {
    FILE *report = fopen(REPORT, "w");
    if (!report) { fprintf(stderr, "cannot open " REPORT "\n"); return 1; }

    Mat y[N_REPLICATES], switching[N_REPLICATES];
    for (int i = 0; i < N_REPLICATES; i++) {
        Rng rng = rng_new(2027, (uint64_t)i);
        Mat block = stored_block(&rng, T_STORED);
        if (i == NAN_REPLICATE) AT(block, STORED_GDP_GROWTH, 50) = (mreal)NAN;
        lp_system_from_stored(block, &y[i], &switching[i]);
        mat_free(block);
    }

    LpNlSpec spec = lp_system_spec();
    LpSystemStateWeight weights[N_WEIGHTS] = {
        lp_system_state_weight(0, 0), lp_system_state_weight(1, 0),
        lp_system_state_weight(0, 1), lp_system_state_weight(1, 1)
    };
    LpSystemRecord rows[N_WEIGHTS][N_REPLICATES];
    char what[256];

    for (int w = 0; w < N_WEIGHTS; w++) {
        LpSystemStateWeight weight = weights[w];
        const char *name = lp_system_state_weight_description(weight);
        int lambda_right = 1, flags_right = 1, responses_right = 1;
        for (int i = 0; i < N_REPLICATES; i++) {
            rows[w][i] = lp_system_fit_record(y[i], switching[i], spec, weight, 100 + i);
            if (i == NAN_REPLICATE) continue;
            LpNlSpec direct = spec;
            HpEstimate estimate;
            memset(&estimate, 0, sizeof estimate);
            if (weight.estimate_lambda) {
                estimate = hp_estimate(switching[i], hp_estimate_default_options());
                direct.lambda = estimate.lambda;
            }
            LpLinNlFit fit = weight.causal_hp ? lp_lin_and_nl_causal(y[i], switching[i], direct, (LpBands){0})
                                              : lp_lin_and_nl(y[i], switching[i], direct, (LpBands){0});
            const LpSystemRecord *r = &rows[w][i];
            lambda_right = lambda_right && same_bits(r->lambda, direct.lambda) && r->causal_hp == weight.causal_hp;
            flags_right = flags_right && r->lambda_is_converged == (weight.estimate_lambda && estimate.is_converged)
                       && r->lambda_at_bound == (weight.estimate_lambda && estimate.at_bound);
            responses_right = responses_right && r->has_responses && fit.nl.causal_hp == weight.causal_hp
                           && memcmp(r->d, fit.lin.d.d, sizeof r->d) == 0
                           && memcmp(r->lin, fit.lin.irf_lin_mean.d, sizeof r->lin) == 0
                           && memcmp(r->s1, fit.nl.irf_s1_mean.d, sizeof r->s1) == 0
                           && memcmp(r->s2, fit.nl.irf_s2_mean.d, sizeof r->s2) == 0;
            lp_lin_nl_fit_free(&fit);
        }
        snprintf(what, sizeof what, "1. %s: lambda %s and the cycle recorded", name,
                 weight.estimate_lambda ? "hp_estimate's, bit for bit," : "1600");
        check(report, lambda_right, what);
        snprintf(what, sizeof what, "1. %s: hp_estimate's flags recorded as it returned them, 0 under 1600", name);
        check(report, flags_right, what);
        snprintf(what, sizeof what, "1. %s: d and the three response tensors bit for bit those of et_al", name);
        check(report, responses_right, what);

        const LpSystemRecord *hole = &rows[w][NAN_REPLICATE];
        int hole_right = !hole->has_responses && hole->causal_hp == weight.causal_hp
                      && (weight.estimate_lambda ? MISNAN(hole->lambda) && !hole->lambda_is_converged
                                                       && !hole->lambda_at_bound
                                                 : hole->lambda == spec.lambda);
        snprintf(what, sizeof what, "2. %s: a NaN switching series gives no responses and %s", name,
                 weight.estimate_lambda ? "lambda NaN, flags 0" : "lambda 1600");
        check(report, hole_right, what);
    }

    /* The estimate is a number away from 1600 on these series, so any of the
       three non-lpirfs weights that fell back on lpirfs' would show here. */
    int causal_differs = 1, estimate_differs = 1;
    for (int i = 0; i < NAN_REPLICATE; i++) {
        causal_differs = causal_differs && memcmp(rows[0][i].s1, rows[1][i].s1, sizeof rows[0][i].s1) != 0
                      && memcmp(rows[2][i].s1, rows[3][i].s1, sizeof rows[2][i].s1) != 0;
        estimate_differs = estimate_differs && rows[2][i].lambda != spec.lambda
                        && memcmp(rows[0][i].s1, rows[2][i].s1, sizeof rows[0][i].s1) != 0
                        && memcmp(rows[1][i].s1, rows[3][i].s1, sizeof rows[1][i].s1) != 0;
    }
    check(report, causal_differs, "3. the causal and the two-sided cycle give different state 1 responses");
    check(report, estimate_differs, "3. the estimated lambda differs from 1600 and gives different state 1 responses");
    fprintf(report, "     estimated lambdas, two-sided cycle:");
    for (int i = 0; i < NAN_REPLICATE; i++) fprintf(report, " %.6g", rows[2][i].lambda);
    fprintf(report, "\n");

    for (int w = 0; w < N_WEIGHTS; w++) {
        const char *name = lp_system_state_weight_description(weights[w]);
        lp_system_write_weighted_archive(SCRATCH_ARCHIVE, rows[w], N_REPLICATES, spec, weights[w]);
        int n_read = 0;
        LpSystemRecord *read_back = lp_system_read_weighted_archive(SCRATCH_ARCHIVE, spec, weights[w], &n_read);
        int lpirfs_weight = lp_system_state_weight_is_lpirfs(weights[w]);
        int identical = read_back && n_read == N_REPLICATES;
        for (int i = 0; identical && i < N_REPLICATES; i++)
            identical = state_fields_identical(&rows[w][i], &read_back[i])
                     && (lpirfs_weight ? linear_fields_identical(&rows[w][i], &read_back[i])
                                       : linear_fields_absent(&read_back[i]));
        free(read_back);
        snprintf(what, sizeof what, "4. %s: %s", name,
                 lpirfs_weight ? "every field of every row read back bit for bit"
                               : "every field of every row read back bit for bit, the linear projection and d "
                                 "absent");
        check(report, identical, what);

        if (!lpirfs_weight) {
            DataFrame frame = df_read_npz(SCRATCH_ARCHIVE);
            int n_weight_columns = 2 + (weights[w].estimate_lambda ? LP_SYSTEM_N_LAMBDA_SEARCH + 2 : 0);
            int expected = 8 + LP_SYSTEM_HOR + LP_SYSTEM_N_SPEC + 2 * LP_SYSTEM_RESPONSE_DIM + n_weight_columns;
            int no_linear = !has_column_prefix(&frame, "lin_") && !has_column_prefix(&frame, "d_");
            snprintf(what, sizeof what, "4. %s: the archive has no lin_ or d_ column and %d columns in all",
                     name, expected);
            check(report, no_linear && frame.n_cols == expected, what);
            df_free(&frame);
        }

        int others_miss = 1;
        for (int v = 0; v < N_WEIGHTS; v++) {
            if (v == w) continue;
            int n_other = -1;
            LpSystemRecord *other = lp_system_read_weighted_archive(SCRATCH_ARCHIVE, spec, weights[v], &n_other);
            others_miss = others_miss && other == NULL && n_other == 0;
            free(other);
        }
        snprintf(what, sizeof what, "4. %s: a miss under each of the three other weights", name);
        check(report, others_miss, what);

        if (w == 0) {
            int n_plain = -1;
            LpSystemRecord *plain = lp_system_read_archive(SCRATCH_ARCHIVE, spec, &n_plain);
            check(report, plain != NULL && n_plain == N_REPLICATES,
                  "4. lpirfs' weight: lp_system_read_archive reads it as before");
            free(plain);
            DataFrame frame = df_read_npz(SCRATCH_ARCHIVE);
            int expected = 8 + 2 * LP_SYSTEM_HOR + LP_SYSTEM_N_SPEC + LP_SYSTEM_K * LP_SYSTEM_K
                         + 3 * LP_SYSTEM_RESPONSE_DIM;
            check(report, frame.n_cols == expected,
                  "4. lpirfs' weight: the archive has the columns of the archives stored before the other weights");
            df_free(&frame);

            LpNlSpec other_lambda = spec;
            other_lambda.lambda = 100;
            int n_other = -1;
            LpSystemRecord *other = lp_system_read_archive(SCRATCH_ARCHIVE, other_lambda, &n_other);
            check(report, other == NULL && n_other == 0, "4. lpirfs' weight: a miss under lambda 100");
            free(other);
        }
        if (weights[w].estimate_lambda) {
            LpNlSpec other_lambda = spec;
            other_lambda.lambda = 100;
            int n_any = -1;
            LpSystemRecord *any = lp_system_read_weighted_archive(SCRATCH_ARCHIVE, other_lambda, weights[w], &n_any);
            snprintf(what, sizeof what, "4. %s: read whatever spec.lambda is", name);
            check(report, any != NULL && n_any == N_REPLICATES, what);
            free(any);

            LpSystemStateWeight other_search = weights[w];
            other_search.lambda_options.xatol = 1e-5;
            int n_other = -1;
            LpSystemRecord *other = lp_system_read_weighted_archive(SCRATCH_ARCHIVE, spec, other_search, &n_other);
            snprintf(what, sizeof what, "4. %s: a miss under another search tolerance", name);
            check(report, other == NULL && n_other == 0, what);
            free(other);
        }
        remove(SCRATCH_ARCHIVE);
    }

    fprintf(report, "\n%s\n", failures ? "FAILED" : "all passed");
    fclose(report);
    for (int i = 0; i < N_REPLICATES; i++) {
        mat_free(y[i]);
        mat_free(switching[i]);
    }
    if (failures) fprintf(stderr, "lp_system_state_weight_correctness: %d check(s) failed, see " REPORT "\n", failures);
    return failures ? 1 : 0;
}
