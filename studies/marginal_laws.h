#ifndef MARGINAL_LAWS_H
#define MARGINAL_LAWS_H

/*
Scalar log-densities and maximum likelihood for the laws a marginal with finite
variance and a tail heavier than Gaussian can belong to, and the modified
Bessel function of the second kind three of them need.

None of this is specific to this project and none of it is in et_al yet.
et_al has dist/gauss.h and dist/student.h in their Mat forms and no Laplace,
no Subbotin, no normal inverse Gaussian, no variance gamma, no generalised
hyperbolic, and no Bessel function of any kind. All of that belongs there. It
is written here so the fit can happen now, in the shape et_al's dist/ headers
use, so that porting it is a move rather than a rewrite.

The Gaussian and the Student t are written here too rather than called from
et_al. Their Mat forms allocate a matrix per evaluation, which a likelihood
called tens of thousands of times inside a line search cannot afford, and a
comparison across laws is only worth reading if every law came through the same
arithmetic. The scalar forms here agree with et_al's to the last bit, and
tests/marginal_laws_correctness.c is what says so.

Every law is written in a location-scale parameterisation with the location
first and the scale second, so the same fitting routine drives all of them.

THE BESSEL FUNCTION

log K_nu(x) for real nu >= 0 and x > 0, from

    K_nu(x) = integral_0^inf exp(-x cosh t) cosh(nu t) dt.

The logarithm is what is computed, never K itself: at the arguments the normal
inverse Gaussian reaches, x is in the hundreds and K underflows a double while
its logarithm is an ordinary number.

The substitution t = u / sqrt(x) turns the exponent into -u^2/2 + O(u^4/x), so
the integrand is close to a Gaussian in u whatever x is, and a fixed
Gauss-Legendre rule on u in [0, U] is accurate across the whole range rather
than only near one scale. Small x is the exception: there the integrand decays
like exp(-x cosh t) rather than like a Gaussian and needs t out to roughly
log(2/x) plus the working precision, so the upper limit is set from x instead.

tests/marginal_laws_correctness.c checks it against the two half-integer orders
that have closed forms, K_{1/2}(x) = sqrt(pi/2x) exp(-x) and
K_{3/2}(x) = sqrt(pi/2x) exp(-x) (1 + 1/x), at orders where no closed form
exists against published values, and against the recurrence
K_{nu+1} = K_{nu-1} + (2 nu / x) K_nu.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <et_al./linalg/mat.h>
#include <et_al./solver/lbfgs.h>

#define MARGINAL_LAWS_PI 3.14159265358979323846

/* Nodes of a fixed Gauss-Legendre rule on [-1, 1], 64 points, enough that the
   integrand's remaining departure from a Gaussian is below the tolerance the
   correctness test holds it to. Only the non-negative half is stored, since the
   rule is symmetric. */
#define MARGINAL_LAWS_NODES 32

static const double marginal_laws_node[MARGINAL_LAWS_NODES] = {
    0.0243502926634244, 0.0729931217877990, 0.1214628192961206, 0.1696444204239928,
    0.2174236437400071, 0.2646871622087674, 0.3113228719902110, 0.3572201583376681,
    0.4022701579639916, 0.4463660172534641, 0.4894031457070530, 0.5312794640198946,
    0.5718956462026340, 0.6111553551723933, 0.6489654712546573, 0.6852363130542333,
    0.7198818501716109, 0.7528199072605319, 0.7839723589433414, 0.8132653151227975,
    0.8406292962525803, 0.8659993981540928, 0.8893154459951141, 0.9105221370785028,
    0.9295691721319396, 0.9464113748584028, 0.9610087996520538, 0.9733268277899110,
    0.9833362538846260, 0.9910133714767443, 0.9963401167719553, 0.9993050417357722
};

static const double marginal_laws_weight[MARGINAL_LAWS_NODES] = {
    0.0486909570091397, 0.0485754674415034, 0.0483447622348030, 0.0479993885964583,
    0.0475401657148303, 0.0469681828162100, 0.0462847965813144, 0.0454916279274181,
    0.0445905581637566, 0.0435837245293235, 0.0424735151236536, 0.0412625632426235,
    0.0399537411327203, 0.0385501531786156, 0.0370551285402400, 0.0354722132568824,
    0.0338051618371416, 0.0320579283548516, 0.0302346570724025, 0.0283396726142595,
    0.0263774697150547, 0.0243527025687109, 0.0222701738083833, 0.0201348231535302,
    0.0179517157756973, 0.0157260304760247, 0.0134630478967186, 0.0111681394601311,
    0.0088467598263639, 0.0065044579689784, 0.0041470332605625, 0.0017832807216964
};

/* log K_nu(x) for real nu and x > 0. K is even in its order, K_{-nu} = K_nu,
   which the integral above shows directly since cosh is even, so a negative
   order is turned round here rather than at every call site. The variance
   gamma reaches one whenever its shape goes above 2. */
static inline double marginal_laws_log_bessel_k(double nu, double x) {
    assert(x > 0 && "marginal_laws: log_bessel_k needs a positive argument");
    nu = fabs(nu);

    /* Where to stop, and in which variable. For x of order one and above the
       integrand is a Gaussian in u = t sqrt(x) and twelve standard deviations
       is past anything a double keeps; for smaller x the decay is in t itself.
       The order widens the integrand through cosh(nu t), so the limit carries
       nu as well. */
    double upper, scale;
    int panels;
    if (x >= 1.0) {
        /* In u the integrand is a Gaussian of unit width whatever x is, and
           one rule across twelve of those resolves it to the last few bits. */
        scale = sqrt(x);
        upper = 12.0 + 2.0 * nu / scale;
        panels = 1;
    } else {
        /* Below one the decay is in t rather than in u and the range is tens
           of units wide, which one rule does not resolve. */
        scale = 1.0;
        upper = log(2.0 / x) + 40.0 + nu;
        panels = (int)(upper / 1.5) + 1;
        if (panels > 64) panels = 64;
    }

    double width = upper / (double)panels;
    double half = 0.5 * width;

    /* The integrand is exp(-x (cosh t - 1) + log cosh(nu t)) and the -x is
       carried outside, so nothing underflows before the logarithm is taken. */
    double total = 0;
    for (int panel = 0; panel < panels; panel++) {
        double middle = ((double)panel + 0.5) * width;
        for (int i = 0; i < MARGINAL_LAWS_NODES; i++) {
            for (int side = 0; side < 2; side++) {
                double u = middle + (side ? -half : half) * marginal_laws_node[i];
                double t = u / scale;
                double exponent = -x * (cosh(t) - 1.0);
                double order = nu * t;
                /* log cosh, written so a large order does not overflow. */
                double log_cosh = order > 30.0 ? order - log(2.0) : log(cosh(order));
                double value = exponent + log_cosh;
                if (value > -745.0) total += marginal_laws_weight[i] * exp(value);
            }
        }
    }
    total *= half / scale;

    return total > 0 ? -x + log(total) : -INFINITY;
}

/* Whether a Bessel argument a density has built is one the function can take.
   The function asserts on a non-positive argument, which is right: for it,
   that is a caller's mistake. For a density it is not. An optimiser walks into
   parameter values that are merely infeasible, and the two that produce a
   non-positive argument here are both reachable by ordinary steps: a scale
   driven small enough to underflow, and an asymmetry driven close enough to
   the tail parameter that alpha^2 - beta^2 underflows to zero. Those return a
   sentinel and let the line search back out, per et_al's own rule that an
   assert is for programmer error and an infeasible parameter is not one. */
static inline int marginal_laws_bessel_argument_ok(double x) {
    return x > 0 && isfinite(x);
}

/*
THE LAWS

Each is a scalar log-density of one observation. theta is the law's own
parameter vector in the order its comment gives, already on the natural scale.
A parameter value the law does not admit returns -INFINITY rather than
asserting, because an optimiser will propose them.
*/

/* Gaussian: location, scale. */
static inline double marginal_laws_gaussian(double x, const double *theta) {
    double location = theta[0], scale = theta[1];
    if (!(scale > 0)) return -INFINITY;
    double z = (x - location) / scale;
    return -0.5 * z * z - log(scale) - 0.5 * log(2.0 * MARGINAL_LAWS_PI);
}

/* Laplace: location, scale. Excess kurtosis 3 whatever the parameters. */
static inline double marginal_laws_laplace(double x, const double *theta) {
    double location = theta[0], scale = theta[1];
    if (!(scale > 0)) return -INFINITY;
    return -fabs(x - location) / scale - log(2.0 * scale);
}

/* Subbotin, the generalised error: location, scale, shape b. Density
   proportional to exp(-|z|^b), so b = 2 is Gaussian and b = 1 Laplace. The
   scale here is the a of exp(-|(x-mu)/a|^b), not the standard deviation. */
static inline double marginal_laws_subbotin(double x, const double *theta) {
    double location = theta[0], scale = theta[1], shape = theta[2];
    if (!(scale > 0) || !(shape > 0)) return -INFINITY;
    double z = fabs(x - location) / scale;
    return -pow(z, shape) - log(2.0 * scale) + log(shape) - lgamma(1.0 / shape);
}

/* Student t: location, scale, degrees of freedom nu. Tail index nu, so the
   variance exists for nu > 2 and the fourth moment for nu > 4. */
static inline double marginal_laws_student(double x, const double *theta) {
    double location = theta[0], scale = theta[1], nu = theta[2];
    if (!(scale > 0) || !(nu > 0)) return -INFINITY;
    double z = (x - location) / scale;
    return lgamma(0.5 * (nu + 1.0)) - lgamma(0.5 * nu)
           - 0.5 * log(nu * MARGINAL_LAWS_PI) - log(scale)
           - 0.5 * (nu + 1.0) * log1p(z * z / nu);
}

/* Normal inverse Gaussian: location mu, scale delta, tail alpha, asymmetry
   beta, with |beta| < alpha. Semi-heavy tails, every moment finite. */
static inline double marginal_laws_nig(double x, const double *theta) {
    double mu = theta[0], delta = theta[1], alpha = theta[2], beta = theta[3];
    if (!(delta > 0) || !(alpha > 0) || !(fabs(beta) < alpha)) return -INFINITY;

    double gamma = sqrt(alpha * alpha - beta * beta);
    double centred = x - mu;
    double radius = sqrt(delta * delta + centred * centred);
    if (!marginal_laws_bessel_argument_ok(alpha * radius)) return -INFINITY;

    return log(alpha * delta / MARGINAL_LAWS_PI) - log(radius)
           + delta * gamma + beta * centred
           + marginal_laws_log_bessel_k(1.0, alpha * radius);
}

/* Variance gamma: location mu, scale sigma, shape kappa, asymmetry theta. A
   normal whose variance is a gamma draw, so every moment is finite; kappa is
   the variance of the mixing gamma and kappa -> 0 is the Gaussian. */
static inline double marginal_laws_variance_gamma(double x, const double *theta) {
    double mu = theta[0], sigma = theta[1], kappa = theta[2], skew = theta[3];
    if (!(sigma > 0) || !(kappa > 0)) return -INFINITY;

    double centred = x - mu;
    double variance = sigma * sigma;
    double root = sqrt(2.0 * variance / kappa + skew * skew);
    double order = 1.0 / kappa - 0.5;

    double common = log(2.0) + skew * centred / variance
                    - log(sigma) - 0.5 * log(2.0 * MARGINAL_LAWS_PI)
                    - log(kappa) / kappa - lgamma(1.0 / kappa);

    /* Exactly at the location the density is |x - mu|^order times
       K_order(c |x - mu|), which is zero times infinity rather than either.
       For a positive order the two cancel: K_order(z) -> Gamma(order) (2/z)^order
       / 2 as z falls to zero, the log |x - mu| goes with it, and what is left
       is the finite value below. For an order of zero or less the density does
       diverge there and there is nothing to return. An observation landing on
       the location is a measure-zero event in data, but an optimiser moving
       the location onto one is not. */
    double distance = fabs(centred);
    if (!(distance > 0)) {
        if (!(order > 0)) return -INFINITY;
        return common + log(0.5) + lgamma(order)
               + order * (log(2.0 * variance) - 2.0 * log(root));
    }

    if (!marginal_laws_bessel_argument_ok(distance * root / variance)) return -INFINITY;

    return common + order * (log(distance) - log(root))
           + marginal_laws_log_bessel_k(order, distance * root / variance);
}

/* Generalised hyperbolic: location mu, scale delta, tail alpha, asymmetry
   beta, class lambda. Contains the normal inverse Gaussian at lambda = -1/2
   and the hyperbolic at lambda = 1. Five parameters on a few hundred
   observations is not well identified, and the fit reports that rather than
   hiding it. */
static inline double marginal_laws_hyperbolic(double x, const double *theta) {
    double mu = theta[0], delta = theta[1], alpha = theta[2], beta = theta[3];
    double lambda = theta[4];
    if (!(delta > 0) || !(alpha > 0) || !(fabs(beta) < alpha)) return -INFINITY;

    double gamma = sqrt(alpha * alpha - beta * beta);
    double centred = x - mu;
    double radius = sqrt(delta * delta + centred * centred);
    if (!marginal_laws_bessel_argument_ok(delta * gamma)) return -INFINITY;
    if (!marginal_laws_bessel_argument_ok(alpha * radius)) return -INFINITY;

    double log_norm = lambda * (log(gamma) - log(delta))
                      - 0.5 * log(2.0 * MARGINAL_LAWS_PI)
                      - (lambda - 0.5) * log(alpha)
                      - marginal_laws_log_bessel_k(fabs(lambda), delta * gamma);
    return log_norm + (lambda - 0.5) * log(radius) + beta * centred
           + marginal_laws_log_bessel_k(fabs(lambda - 0.5), alpha * radius);
}

/* A Gaussian whose standard deviation is lognormal: location, the median of
   that standard deviation, and the standard deviation of its logarithm. The
   density has no closed form, so it is the mixture integral taken by
   Gauss-Hermite; s = 0 is the Gaussian. This is the member of the normal
   variance mixtures whose every moment is finite. */
#define MARGINAL_LAWS_MIXTURE_NODES 24

static inline double marginal_laws_lognormal_mixture(double x, const double *theta) {
    double location = theta[0], median = theta[1], spread = theta[2];
    if (!(median > 0) || !(spread >= 0)) return -INFINITY;
    if (spread == 0) {
        double flat[2] = {location, median};
        return marginal_laws_gaussian(x, flat);
    }

    /* Gauss-Legendre on the logarithm of the standard deviation, over the
       range the lognormal puts essentially all of its mass in. */
    double upper = 8.0 * spread;
    double half = upper;
    double total = 0;
    for (int i = 0; i < MARGINAL_LAWS_NODES; i++) {
        for (int side = 0; side < 2; side++) {
            double g = half * (side ? -marginal_laws_node[i] : marginal_laws_node[i]);
            double sd = median * exp(g);
            double z = (x - location) / sd;
            double log_weight = -0.5 * (g / spread) * (g / spread)
                                - log(spread) - 0.5 * log(2.0 * MARGINAL_LAWS_PI);
            double log_density = -0.5 * z * z - log(sd) - 0.5 * log(2.0 * MARGINAL_LAWS_PI);
            double value = log_weight + log_density;
            if (value > -745.0) total += marginal_laws_weight[i] * exp(value);
        }
    }
    total *= half;
    return total > 0 ? log(total) : -INFINITY;
}

/*
FITTING

One routine drives every law. A law is a log-density, a count of parameters, a
map from the unconstrained vector the optimiser steps to the natural one the
density reads, and a starting point built from the sample. Keeping the map
beside the density is what stops the two drifting: a law whose scale must be
positive says so once, here, rather than in every caller.

The optimiser is et_al's solver/lbfgs.h. Nothing is reimplemented. Its
objective wants a gradient, and none of these densities is differentiated in
closed form here, so the gradient is central differences on the unconstrained
vector. That is the one place this cuts a corner, and the reason is that the
normal inverse Gaussian, the variance gamma and the generalised hyperbolic all
differentiate through a Bessel function whose own derivative would need the
same quadrature again; the fit is checked by whether it recovers parameters it
was given, in tests/marginal_laws_correctness.c, rather than by the gradient
being exact.
*/

#define MARGINAL_LAWS_MAX_PARAMETERS 5

typedef double (*MarginalLogDensity)(double x, const double *theta);
typedef void (*MarginalUnpack)(const double *unconstrained, double *theta);
typedef void (*MarginalStart)(const double *x, int n, double *unconstrained);

typedef struct {
    const char *name;
    int parameters;
    MarginalLogDensity log_density;
    MarginalUnpack unpack;
    MarginalStart start;
} MarginalLaw;

/* Sample mean and standard deviation, which every starting point is built
   from. */
static inline void marginal_laws_moments(const double *x, int n, double *mean, double *sd) {
    double total = 0;
    for (int i = 0; i < n; i++) total += x[i];
    *mean = total / (double)n;
    double m2 = 0;
    for (int i = 0; i < n; i++) { double d = x[i] - *mean; m2 += d * d; }
    *sd = n > 1 ? sqrt(m2 / (double)(n - 1)) : 1.0;
    if (!(*sd > 0)) *sd = 1.0;
}

static inline void marginal_laws_unpack_gaussian(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
}
static inline void marginal_laws_start_gaussian(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd);
}

/* The Laplace scale that matches a given standard deviation is sd / sqrt(2). */
static inline void marginal_laws_start_laplace(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd / sqrt(2.0));
}

static inline void marginal_laws_unpack_subbotin(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = exp(u[2]);
}
/* Started at the Gaussian, shape 2, whose matching scale is sd sqrt(2). */
static inline void marginal_laws_start_subbotin(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd * sqrt(2.0));
    u[2] = log(2.0);
}

static inline void marginal_laws_unpack_student(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = exp(u[2]);
}
/* Started at nu = 8, where the variance and the fourth moment both exist, with
   the scale that matches the sample standard deviation there. */
static inline void marginal_laws_start_student(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd * sqrt((8.0 - 2.0) / 8.0));
    u[2] = log(8.0);
}

/* beta is held inside alpha by a tanh, which is what keeps |beta| < alpha
   without the optimiser meeting a wall. */
static inline void marginal_laws_unpack_nig(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = exp(u[2]);
    theta[3] = theta[2] * tanh(u[3]);
}
/* Symmetric, and at the alpha and delta whose variance is the sample's: for a
   symmetric normal inverse Gaussian the variance is delta / alpha, so setting
   both from the sample standard deviation puts the start in the right place
   without committing to a tail. */
static inline void marginal_laws_start_nig(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd);
    u[2] = log(1.0 / sd);
    u[3] = 0.0;
}

/* kappa is held below 2, which is not a convenience. The density carries
   |x - mu|^(1/kappa - 1/2), so above 2 that exponent turns negative and the
   density diverges at the location. The divergence is integrable and the law
   is perfectly good, but the likelihood is not: with mu free, putting it on an
   observation sends that observation's density and the whole likelihood to
   infinity, so there is no maximum to find. It is the pathology the
   three-parameter lognormal has. The search stays where a maximum exists and
   the density itself still handles the other side, which the negative-order
   checks in tests/marginal_laws_correctness.c exercise. */
static inline void marginal_laws_unpack_variance_gamma(const double *u, double *theta) {
    /* The logistic reaches exactly 2 in double precision once its argument is
       large enough, and exactly 2 is not safe either: the order is zero there
       and the density still diverges at the location, logarithmically rather
       than as a power, which is enough for the likelihood to be unbounded
       again. Clamping the argument leaves the result a representable distance
       below the bound, about 2e-13 at the value used here, which is closer to
       2 than any sample of a few hundred can tell apart. */
    const double reach = 30.0;
    double shape = u[2] > reach ? reach : (u[2] < -reach ? -reach : u[2]);

    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = 2.0 / (1.0 + exp(-shape));
    theta[3] = u[3];
}
/* Started symmetric at kappa = 0.5, away from the Gaussian limit kappa -> 0
   and from the bound at 2. */
static inline void marginal_laws_start_variance_gamma(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd);
    u[2] = log(0.5 / (2.0 - 0.5));
    u[3] = 0.0;
}

static inline void marginal_laws_unpack_hyperbolic(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = exp(u[2]);
    theta[3] = theta[2] * tanh(u[3]);
    theta[4] = u[4];
}
/* Started at the normal inverse Gaussian, lambda = -1/2, so the extra
   parameter begins where the four-parameter member sits. */
static inline void marginal_laws_start_hyperbolic(const double *x, int n, double *u) {
    marginal_laws_start_nig(x, n, u);
    u[4] = -0.5;
}

static inline void marginal_laws_unpack_lognormal_mixture(const double *u, double *theta) {
    theta[0] = u[0];
    theta[1] = exp(u[1]);
    theta[2] = exp(u[2]);
}
/* Started at a spread of 0.3, which is away from the Gaussian limit of 0 and
   gives an excess kurtosis of about 0.9. */
static inline void marginal_laws_start_lognormal_mixture(const double *x, int n, double *u) {
    double mean, sd;
    marginal_laws_moments(x, n, &mean, &sd);
    u[0] = mean;
    u[1] = log(sd);
    u[2] = log(0.3);
}

#define MARGINAL_LAWS_COUNT 8

static inline const MarginalLaw *marginal_laws_all(void) {
    static const MarginalLaw law[MARGINAL_LAWS_COUNT] = {
        {"Gaussian", 2, marginal_laws_gaussian, marginal_laws_unpack_gaussian,
         marginal_laws_start_gaussian},
        {"Laplace", 2, marginal_laws_laplace, marginal_laws_unpack_gaussian,
         marginal_laws_start_laplace},
        {"Subbotin", 3, marginal_laws_subbotin, marginal_laws_unpack_subbotin,
         marginal_laws_start_subbotin},
        {"Student t", 3, marginal_laws_student, marginal_laws_unpack_student,
         marginal_laws_start_student},
        {"lognormal mixture", 3, marginal_laws_lognormal_mixture,
         marginal_laws_unpack_lognormal_mixture, marginal_laws_start_lognormal_mixture},
        {"normal inverse Gaussian", 4, marginal_laws_nig, marginal_laws_unpack_nig,
         marginal_laws_start_nig},
        {"variance gamma", 4, marginal_laws_variance_gamma,
         marginal_laws_unpack_variance_gamma, marginal_laws_start_variance_gamma},
        {"generalised hyperbolic", 5, marginal_laws_hyperbolic,
         marginal_laws_unpack_hyperbolic, marginal_laws_start_hyperbolic}
    };
    return law;
}

typedef struct {
    const MarginalLaw *law;
    const double *x;
    int n;
} MarginalContext;

/* Minus the log-likelihood, which is what an optimiser descends. An
   unreachable parameter value returns a large finite number rather than an
   infinity, so the line search can back out of it instead of the search ending
   with LBFGS_NOT_FINITE. */
static inline mreal marginal_laws_objective_value(const MarginalContext *context,
                                                  const double *unconstrained) {
    double theta[MARGINAL_LAWS_MAX_PARAMETERS];
    context->law->unpack(unconstrained, theta);

    double total = 0;
    for (int i = 0; i < context->n; i++) {
        double here = context->law->log_density(context->x[i], theta);
        if (!isfinite(here)) return (mreal)1e100;
        total += here;
    }
    return (mreal)(-total);
}

static inline mreal marginal_laws_objective(Vec theta, Vec gradient, void *raw) {
    const MarginalContext *context = (const MarginalContext*)raw;
    int n = context->law->parameters;

    double point[MARGINAL_LAWS_MAX_PARAMETERS];
    for (int i = 0; i < n; i++) point[i] = (double)theta.d[i];
    mreal value = marginal_laws_objective_value(context, point);

    if (gradient.d) {
        for (int i = 0; i < n; i++) {
            double step = 1e-5 * (1.0 + fabs(point[i]));
            double keep = point[i];
            point[i] = keep + step;
            mreal up = marginal_laws_objective_value(context, point);
            point[i] = keep - step;
            mreal down = marginal_laws_objective_value(context, point);
            point[i] = keep;
            gradient.d[i] = (mreal)(((double)up - (double)down) / (2.0 * step));
        }
    }
    return value;
}

typedef struct {
    double log_likelihood;
    double theta[MARGINAL_LAWS_MAX_PARAMETERS];
    int parameters;
    int is_converged;
    double aic;
    double bic;
} MarginalFit;

/* Maximum likelihood for one law on one sample. */
static inline MarginalFit marginal_laws_fit(const MarginalLaw *law, const double *x, int n) {
    MarginalContext context = {law, x, n};

    double start[MARGINAL_LAWS_MAX_PARAMETERS];
    law->start(x, n, start);

    Vec begin = mat_new(law->parameters, 1);
    for (int i = 0; i < law->parameters; i++) begin.d[i] = (mreal)start[i];

    LbfgsOptions options = lbfgs_default_options();
    options.max_iterations = 500;
    options.gradient_tolerance = (mreal)1e-5;
    LbfgsResult result = lbfgs(marginal_laws_objective, &context, begin, options);

    MarginalFit fit;
    fit.parameters = law->parameters;
    fit.log_likelihood = -(double)result.value;
    fit.is_converged = result.is_converged;

    double unconstrained[MARGINAL_LAWS_MAX_PARAMETERS];
    for (int i = 0; i < law->parameters; i++) unconstrained[i] = (double)result.theta.d[i];
    law->unpack(unconstrained, fit.theta);
    for (int i = law->parameters; i < MARGINAL_LAWS_MAX_PARAMETERS; i++) fit.theta[i] = NAN;

    fit.aic = 2.0 * law->parameters - 2.0 * fit.log_likelihood;
    fit.bic = law->parameters * log((double)n) - 2.0 * fit.log_likelihood;

    mat_free(begin);
    mat_free(result.theta);
    return fit;
}

#endif /* MARGINAL_LAWS_H */
