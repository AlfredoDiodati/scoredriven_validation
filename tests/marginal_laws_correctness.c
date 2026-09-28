/*
Whether studies/marginal_laws.h computes what it claims to compute.

Five things are checked, in an order that puts the pieces the rest depend on
first. Nothing here reads the experiment; every number is either a published
constant, a closed form, or data this test simulated itself.

  the Bessel function   log K_nu(x) against the two orders that have closed
                        forms, K_{1/2}(x) = sqrt(pi/2x) exp(-x) and
                        K_{3/2}(x) = sqrt(pi/2x) exp(-x) (1 + 1/x); against
                        published values at orders that have none; and against
                        the recurrence K_{nu+1} = K_{nu-1} + (2 nu / x) K_nu,
                        which ties orders together that the closed forms do not
                        reach.
  the densities         each one integrates to 1 over a range wide enough that
                        what is outside it is below the tolerance. A density
                        that is off by a constant factor fits as well as one
                        that is not, and the likelihood comparison the study
                        makes would be meaningless.
  one law inside        the generalised hyperbolic at lambda = -1/2 is the
  another               normal inverse Gaussian. The two are written
                        separately, so this is what says the five-parameter
                        form is the family it claims to contain.
  against et_al         the scalar Gaussian and Student t against
                        dist/gauss.h and dist/student.h, which are the same
                        laws in their Mat form.
  the fit               maximum likelihood recovers the parameters it was
                        given, on samples drawn from each law at a size where
                        the estimate should be close. This is what stands in
                        for the analytic gradient the fit does not have.

Writes out/marginal_laws_correctness.txt.
*/

#include "studies/marginal_laws.h"

#include <et_al./dist/gauss.h>
#include <et_al./dist/student.h>
#include <et_al./random/random.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define REPORT "out/marginal_laws_correctness.txt"
#define SEED 20260923
#define SAMPLE 20000

static int failures = 0;
static FILE *report = NULL;

static void check(const char *what, double got, double want, double tolerance) {
    double gap = fabs(got - want);
    double relative = fabs(want) > 1e-12 ? gap / fabs(want) : gap;
    int ok = relative <= tolerance;
    if (!ok) failures++;
    fprintf(report, "  %-52s %16.10g %16.10g %10.2e  %s\n",
            what, got, want, relative, ok ? "ok" : "FAILED");
}

/* An expected infinity cannot go through either tolerance above: the
   difference of two infinities is not a number and every comparison against it
   is false, so an exactly right answer reads as a failure. */
static void check_is_minus_infinity(const char *what, double got) {
    int ok = got == -INFINITY;
    if (!ok) failures++;
    fprintf(report, "  %-52s %16.10g %16s %10s  %s\n",
            what, got, "-inf", "", ok ? "ok" : "FAILED");
}

/* A location is the one parameter a relative tolerance is wrong for. Its true
   value can sit anywhere including at zero, and what an estimate of it should
   be held to is the spread of the data it came from, not its own size. */
static void check_absolute(const char *what, double got, double want, double tolerance) {
    double gap = fabs(got - want);
    int ok = gap <= tolerance;
    if (!ok) failures++;
    fprintf(report, "  %-52s %16.10g %16.10g %10.2e  %s\n",
            what, got, want, gap, ok ? "ok" : "FAILED");
}

/* Closed forms of the two half-integer orders. */
static double bessel_k_half(double x) { return sqrt(MARGINAL_LAWS_PI / (2.0 * x)) * exp(-x); }
static double bessel_k_three_half(double x) { return bessel_k_half(x) * (1.0 + 1.0 / x); }

/* Trapezoid over a fixed grid, used only to ask whether a density integrates
   to one. The grid is fine enough and wide enough that the answer is decided
   by the density rather than by the rule. */
static double integrates_to(const MarginalLaw *law, const double *theta,
                            double from, double to, int steps) {
    double h = (to - from) / (double)steps;
    double total = 0;
    for (int i = 0; i <= steps; i++) {
        double x = from + h * (double)i;
        double density = exp(law->log_density(x, theta));
        if (!isfinite(density)) continue;
        total += (i == 0 || i == steps ? 0.5 : 1.0) * density;
    }
    return total * h;
}

static double laplace_draw(Rng *rng, double location, double scale) {
    double u = rng_uniform(rng) - 0.5;
    return location - scale * (u < 0 ? -1.0 : 1.0) * log(1.0 - 2.0 * fabs(u));
}

/* |X|^shape is Gamma(1/shape, 1), which is where a Subbotin draw comes from. */
static double subbotin_draw(Rng *rng, double location, double scale, double shape) {
    double g = rng_gamma(rng, 1.0 / shape);
    double magnitude = scale * pow(g, 1.0 / shape);
    return location + (rng_uniform(rng) < 0.5 ? -magnitude : magnitude);
}

static double student_t_draw(Rng *rng, double location, double scale, double nu) {
    double chi_squared = 2.0 * rng_gamma(rng, 0.5 * nu);
    return location + scale * rng_normal(rng) / sqrt(chi_squared / nu);
}

/* A normal whose variance is a gamma draw with mean 1 and variance kappa. */
static double variance_gamma_draw(Rng *rng, double mu, double sigma, double kappa, double skew) {
    double g = kappa * rng_gamma(rng, 1.0 / kappa);
    return mu + skew * g + sigma * sqrt(g) * rng_normal(rng);
}

/* Michael, Schucany and Haas for the inverse Gaussian, which is the mixing law
   the normal inverse Gaussian is a normal over. */
static double inverse_gaussian_draw(Rng *rng, double mean, double shape) {
    double z = rng_normal(rng);
    double v = z * z;
    double w = mean + mean * mean * v / (2.0 * shape)
               - mean / (2.0 * shape) * sqrt(4.0 * mean * shape * v + mean * mean * v * v);
    return rng_uniform(rng) <= mean / (mean + w) ? w : mean * mean / w;
}

static double nig_draw(Rng *rng, double mu, double delta, double alpha, double beta) {
    double gamma = sqrt(alpha * alpha - beta * beta);
    double mixing = inverse_gaussian_draw(rng, delta / gamma, delta * delta);
    return mu + beta * mixing + sqrt(mixing) * rng_normal(rng);
}

static double lognormal_mixture_draw(Rng *rng, double location, double median, double spread) {
    double sd = median * exp(spread * rng_normal(rng));
    return location + sd * rng_normal(rng);
}

/* Fits one law to a sample drawn from it and checks every parameter came
   back. */
static void recovers(const char *name, int index, const double *x, int n,
                     const double *truth, const double *tolerance) {
    const MarginalLaw *law = &marginal_laws_all()[index];
    MarginalFit fit = marginal_laws_fit(law, x, n);
    if (!fit.is_converged) {
        failures++;
        fprintf(report, "  %-52s %16s %16s %10s  FAILED\n",
                name, "did not", "converge", "");
        return;
    }
    for (int p = 0; p < law->parameters; p++) {
        char what[128];
        snprintf(what, sizeof what, "%s, parameter %d", name, p + 1);
        if (p == 0) check_absolute(what, fit.theta[p], truth[p], tolerance[p]);
        else check(what, fit.theta[p], truth[p], tolerance[p]);
    }
}

int main(void) {
    report = fopen(REPORT, "w");
    if (!report) { fprintf(stderr, "cannot open %s\n", REPORT); return EXIT_FAILURE; }

    fprintf(report, "Whether studies/marginal_laws.h computes what it claims to compute.\n\n");
    fprintf(report, "  %-52s %16s %16s %10s\n", "what", "got", "want", "relative");

    fprintf(report, "\nThe Bessel function, against the closed forms at orders 1/2 and 3/2.\n");
    static const double argument[] = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 500.0};
    for (size_t i = 0; i < sizeof argument / sizeof argument[0]; i++) {
        char what[128];
        snprintf(what, sizeof what, "log K_{1/2}(%g)", argument[i]);
        check(what, marginal_laws_log_bessel_k(0.5, argument[i]),
              log(bessel_k_half(argument[i])), 1e-10);
        snprintf(what, sizeof what, "log K_{3/2}(%g)", argument[i]);
        check(what, marginal_laws_log_bessel_k(1.5, argument[i]),
              log(bessel_k_three_half(argument[i])), 1e-10);
    }

    fprintf(report, "\nAgainst published values at orders with no closed form.\n");
    check("K_0(1)", exp(marginal_laws_log_bessel_k(0.0, 1.0)), 0.4210244382407083, 1e-10);
    check("K_1(1)", exp(marginal_laws_log_bessel_k(1.0, 1.0)), 0.6019072301972346, 1e-10);
    check("K_0(2)", exp(marginal_laws_log_bessel_k(0.0, 2.0)), 0.1138938727495334, 1e-10);
    check("K_1(2)", exp(marginal_laws_log_bessel_k(1.0, 2.0)), 0.1398658818165224, 1e-10);
    check("K_2(1)", exp(marginal_laws_log_bessel_k(2.0, 1.0)), 1.6248388986351774, 1e-10);
    check("K_0(0.1)", exp(marginal_laws_log_bessel_k(0.0, 0.1)), 2.4270690247020164, 1e-10);
    check("K_1(10)", exp(marginal_laws_log_bessel_k(1.0, 10.0)), 1.8648773453825585e-05, 1e-10);

    fprintf(report, "\nAgainst the recurrence K_{nu+1} = K_{nu-1} + (2 nu / x) K_nu.\n");
    static const double order[] = {0.3, 1.0, 1.7, 3.2};
    for (size_t i = 0; i < sizeof order / sizeof order[0]; i++) {
        for (size_t j = 0; j < sizeof argument / sizeof argument[0]; j++) {
            double nu = order[i], x = argument[j];
            if (x < 0.05) continue;
            double below = exp(marginal_laws_log_bessel_k(nu - 1.0 < 0 ? 1.0 - nu : nu - 1.0, x));
            double middle = exp(marginal_laws_log_bessel_k(nu, x));
            double above = exp(marginal_laws_log_bessel_k(nu + 1.0, x));
            if (!isfinite(below) || below == 0 || !isfinite(above) || above == 0) continue;
            char what[128];
            snprintf(what, sizeof what, "recurrence at nu = %g, x = %g", nu, x);
            check(what, above, below + 2.0 * nu / x * middle, 1e-8);
        }
    }

    fprintf(report, "\nK is even in its order, which the variance gamma reaches whenever its\n"
                    "shape goes above 2.\n");
    for (size_t i = 0; i < sizeof order / sizeof order[0]; i++) {
        for (size_t j = 0; j < sizeof argument / sizeof argument[0]; j++) {
            char what[128];
            snprintf(what, sizeof what, "K_{-%g}(%g) against K_{%g}(%g)",
                     order[i], argument[j], order[i], argument[j]);
            check(what, marginal_laws_log_bessel_k(-order[i], argument[j]),
                  marginal_laws_log_bessel_k(order[i], argument[j]), 1e-15);
        }
    }

    fprintf(report, "\nThe variance gamma at a shape above 2, where its order turns negative and\n"
                    "its density diverges at the location. The fit does not go there, because\n"
                    "the likelihood is unbounded once it can, so this is the density on its own.\n");
    {
        double steep[4] = {0.0, 1.0, 3.0, 0.2};
        check_is_minus_infinity("minus infinity at the location",
                                marginal_laws_variance_gamma(0.0, steep));
        check("finite just beside the location",
              isfinite(marginal_laws_variance_gamma(1e-6, steep)) ? 1.0 : 0.0, 1.0, 0.0);
        check("larger closer in, since the divergence is real",
              marginal_laws_variance_gamma(1e-8, steep) >
              marginal_laws_variance_gamma(1e-4, steep) ? 1.0 : 0.0, 1.0, 0.0);
        double bounded[4];
        static const double reach[] = {40.0, 1e3, 1e30, -1e30};
        for (size_t r = 0; r < sizeof reach / sizeof reach[0]; r++) {
            double unconstrained[4] = {0.0, 0.0, reach[r], 0.0};
            marginal_laws_unpack_variance_gamma(unconstrained, bounded);
            char what[128];
            snprintf(what, sizeof what, "shape below 2 at an unconstrained %g", reach[r]);
            check(what, bounded[2] < 2.0 && bounded[2] > 0.0 ? 1.0 : 0.0, 1.0, 0.0);
            snprintf(what, sizeof what, "order above 0 at an unconstrained %g", reach[r]);
            check(what, 1.0 / bounded[2] - 0.5 > 0.0 ? 1.0 : 0.0, 1.0, 0.0);
        }
    }

    fprintf(report, "\nParameter values an optimiser reaches but no density admits come back as\n"
                    "minus infinity rather than stopping the program.\n");
    {
        /* An asymmetry so close to the tail parameter that alpha^2 - beta^2
           underflows, and a scale small enough to underflow on its own. */
        double edge_nig[4] = {0.0, 1e-320, 1.0, 0.0};
        check_is_minus_infinity("normal inverse Gaussian, scale underflowed",
                                marginal_laws_nig(0.0, edge_nig));
        double edge_gh[5] = {0.0, 1.0, 1e-160, 1e-160 * (1.0 - 1e-17), -0.5};
        check_is_minus_infinity("generalised hyperbolic, alpha and beta together",
                                marginal_laws_hyperbolic(0.5, edge_gh));
        double edge_vg[4] = {0.0, 1e150, 0.5, 0.0};
        check_is_minus_infinity("variance gamma, argument underflowed",
                                marginal_laws_variance_gamma(1e-300, edge_vg));
    }

    fprintf(report, "\nEvery density integrates to 1.\n");
    {
        const MarginalLaw *law = marginal_laws_all();
        double gaussian[2] = {0.4, 1.3};
        check("Gaussian", integrates_to(&law[0], gaussian, -40, 40, 400000), 1.0, 1e-9);
        double laplace[2] = {-0.2, 0.9};
        check("Laplace", integrates_to(&law[1], laplace, -60, 60, 400000), 1.0, 1e-7);
        double subbotin[3] = {0.1, 1.1, 1.4};
        check("Subbotin", integrates_to(&law[2], subbotin, -60, 60, 400000), 1.0, 1e-9);
        double student[3] = {0.3, 1.2, 6.0};
        check("Student t", integrates_to(&law[3], student, -4000, 4000, 4000000), 1.0, 1e-5);
        double mixture[3] = {0.2, 1.0, 0.5};
        check("lognormal mixture", integrates_to(&law[4], mixture, -80, 80, 400000), 1.0, 1e-8);
        double nig[4] = {0.1, 1.0, 1.2, 0.3};
        check("normal inverse Gaussian", integrates_to(&law[5], nig, -80, 80, 400000), 1.0, 1e-8);
        double vg[4] = {0.0, 1.0, 0.5, 0.2};
        check("variance gamma", integrates_to(&law[6], vg, -80, 80, 800000), 1.0, 1e-6);
        double gh[5] = {0.1, 1.0, 1.2, 0.3, -0.5};
        check("generalised hyperbolic", integrates_to(&law[7], gh, -80, 80, 400000), 1.0, 1e-8);
    }

    fprintf(report, "\nThe generalised hyperbolic at lambda = -1/2 is the normal inverse Gaussian.\n");
    for (int i = -30; i <= 30; i += 10) {
        double x = 0.3 * i;
        double nig[4] = {0.1, 1.0, 1.2, 0.3};
        double gh[5] = {0.1, 1.0, 1.2, 0.3, -0.5};
        char what[128];
        snprintf(what, sizeof what, "log density at x = %g", x);
        check(what, marginal_laws_hyperbolic(x, gh), marginal_laws_nig(x, nig), 1e-9);
    }

    fprintf(report, "\nThe scalar Gaussian and Student t against et_al's Mat forms.\n");
    for (int i = -20; i <= 20; i += 5) {
        double x = 0.4 * i;
        Mat xm = mat_new(1, 1), loc = mat_new(1, 1), scale = mat_new(1, 1), nu = mat_new(1, 1);
        xm.d[0] = (mreal)x; loc.d[0] = (mreal)0.3; scale.d[0] = (mreal)1.2; nu.d[0] = (mreal)6.0;

        Mat gauss = gauss_logpdf(xm, loc, scale);
        double theta_gauss[2] = {0.3, 1.2};
        char what[128];
        snprintf(what, sizeof what, "Gaussian log density at x = %g", x);
        check(what, marginal_laws_gaussian(x, theta_gauss), (double)gauss.d[0], 1e-12);
        mat_free(gauss);

        Mat t = student_logpdf(xm, loc, scale, nu);
        double theta_student[3] = {0.3, 1.2, 6.0};
        snprintf(what, sizeof what, "Student t log density at x = %g", x);
        check(what, marginal_laws_student(x, theta_student), (double)t.d[0], 1e-12);
        mat_free(t);

        mat_free(xm); mat_free(loc); mat_free(scale); mat_free(nu);
    }

    fprintf(report, "\nMaximum likelihood recovers the parameters it was given, %d draws each,\n"
                    "seed %d.\n", SAMPLE, SEED);
    {
        Rng rng = rng_new(SEED, 1);
        double *x = (double*)malloc((size_t)SAMPLE * sizeof(double));
        if (!x) { fprintf(stderr, "out of memory\n"); return EXIT_FAILURE; }

        double gaussian[2] = {0.4, 1.3};
        for (int i = 0; i < SAMPLE; i++) x[i] = gaussian[0] + gaussian[1] * rng_normal(&rng);
        double gaussian_tolerance[2] = {0.05, 0.02};
        recovers("Gaussian", 0, x, SAMPLE, gaussian, gaussian_tolerance);

        double laplace[2] = {-0.2, 0.9};
        for (int i = 0; i < SAMPLE; i++) x[i] = laplace_draw(&rng, laplace[0], laplace[1]);
        double laplace_tolerance[2] = {0.05, 0.03};
        recovers("Laplace", 1, x, SAMPLE, laplace, laplace_tolerance);

        double subbotin[3] = {0.1, 1.1, 1.4};
        for (int i = 0; i < SAMPLE; i++)
            x[i] = subbotin_draw(&rng, subbotin[0], subbotin[1], subbotin[2]);
        double subbotin_tolerance[3] = {0.06, 0.04, 0.06};
        recovers("Subbotin", 2, x, SAMPLE, subbotin, subbotin_tolerance);

        double student[3] = {0.3, 1.2, 6.0};
        for (int i = 0; i < SAMPLE; i++)
            x[i] = student_t_draw(&rng, student[0], student[1], student[2]);
        double student_tolerance[3] = {0.06, 0.04, 0.15};
        recovers("Student t", 3, x, SAMPLE, student, student_tolerance);

        double mixture[3] = {0.2, 1.0, 0.5};
        for (int i = 0; i < SAMPLE; i++)
            x[i] = lognormal_mixture_draw(&rng, mixture[0], mixture[1], mixture[2]);
        double mixture_tolerance[3] = {0.08, 0.06, 0.10};
        recovers("lognormal mixture", 4, x, SAMPLE, mixture, mixture_tolerance);

        double nig[4] = {0.1, 1.0, 1.2, 0.3};
        for (int i = 0; i < SAMPLE; i++)
            x[i] = nig_draw(&rng, nig[0], nig[1], nig[2], nig[3]);
        double nig_tolerance[4] = {0.10, 0.12, 0.12, 0.25};
        recovers("normal inverse Gaussian", 5, x, SAMPLE, nig, nig_tolerance);

        double vg[4] = {0.0, 1.0, 0.5, 0.2};
        for (int i = 0; i < SAMPLE; i++)
            x[i] = variance_gamma_draw(&rng, vg[0], vg[1], vg[2], vg[3]);
        double vg_tolerance[4] = {0.10, 0.10, 0.15, 0.30};
        recovers("variance gamma", 6, x, SAMPLE, vg, vg_tolerance);

        free(x);
    }

    fprintf(report, "\n%s\n", failures ? "FAILED" : "PASSED");
    fclose(report);

    printf("marginal laws: %d checks failed\n", failures);
    printf("%s\n", failures ? "FAILED" : "PASSED, 0 failures");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
