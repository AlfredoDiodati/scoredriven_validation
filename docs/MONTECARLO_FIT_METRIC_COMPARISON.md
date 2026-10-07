# Comparing the fit metric across methods: test design

Status: design only. Nothing in this file is implemented yet; it is the
starting point for the code.

## The question

The oracle fit metric of `docs/MONTECARLO_FIT_METRIC.md` gives, for each of the
five methods, one value of $v$ per benchmark of the 1000 by 1000 Monte Carlo,
1,000,000 per method. Its medians are 0.774 for the t-QVARMA and 0.692, 0.706,
0.706 and 0.703 for the linear, state 1, state 2 and both-states local
projections. The question is whether the methods differ by more than chance,
answered by one test over all five methods at once rather than by pairwise
tests.

Two nulls are tested, separately from each other:

1. the five methods have the same mean $v$;
2. no method is systematically ahead: within a configuration, every ordering
   of the five methods' $v$ is equally likely.

What a difference means. $v$ is high when the best configuration in the final
confidence set is close to the benchmark compared with the true configuration's
sampling standard deviation. It is not a measure of finding the configuration
the data came from: within every method, $v$ is higher on the benchmarks where
another configuration wins than on those where the true one does
(`docs/MONTECARLO_FIT_METRIC.md`, "Results"). A test on $v$ therefore says which
method's best configuration fits more tightly. How often each method recovers
the true configuration is the "smallest mean loss" count of the grid, a
different question.

## The observations

Notation: configurations $j = 1, \dots, J$ with $J = 1000$; runs
$r = 1, \dots, S$ with $S = 1000$; methods $m = 1, \dots, M$ with $M = 5$.
$v_{j,r,m}$ is the oracle value of method $m$ when run $r$ of configuration $j$
is the benchmark, read from `montecarlo/out/fit_metric_oracle.csv.gz`.

The comparison is paired: every benchmark has a value under every method.

The runs of one configuration are not independent observations. The $d$ and
$\sigma^2$ of benchmark $(j, r)$ are both computed from the runs of
configuration $j$ other than $r$, so two benchmarks of the same configuration
share 998 of those runs. The unit of observation is therefore the
configuration, through the average over its runs,

$$
\bar v_{j,m} = \frac{1}{S} \sum_{r=1}^{S} v_{j,r,m},
$$

and every test below works on the $J \times M$ matrix of $\bar v_{j,m}$.

Two kinds of dependence between configurations remain.

- Local. Configurations with similar parameters have similar responses and,
  plausibly, similar $\bar v_{j,m}$. This dependence fades with the distance
  between the parameter vectors. The spatial versions of the tests below
  allow for it.
- The shared field. Every benchmark's best configuration is chosen among the
  same 1000 configurations, and its $d$ is a distance to the runs of one of
  them. This links every configuration to every other, near or far, and no
  variance estimator below removes it. Every result is conditional on this
  field of 1000 configurations.

## Null 1: equal means

$$
H_0: \; E\left[\bar v_{j,1}\right] = E\left[\bar v_{j,2}\right] = \dots = E\left[\bar v_{j,M}\right],
$$

the expectation taken over configurations drawn from the design box.

Take method 1 as reference and form the $M - 1 = 4$ contrasts

$$
\delta_j = \left( \bar v_{j,2} - \bar v_{j,1}, \; \dots, \; \bar v_{j,M} - \bar v_{j,1} \right)',
\qquad
\bar\delta = \frac{1}{J} \sum_{j=1}^{J} \delta_j .
$$

The Wald statistic is

$$
W_1 = J \, \bar\delta' \, \hat\Omega^{-1} \, \bar\delta \;\xrightarrow{d}\; \chi^2_{M-1}
\quad \text{under } H_0,
$$

computed by solving $\hat\Omega x = \bar\delta$ through a Cholesky factor, never
by inverting $\hat\Omega$. $W_1$ does not depend on which method is the
reference: another reference turns $\delta_j$ into $A \delta_j$ for an
invertible $A$, $\hat\Omega$ into $A \hat\Omega A'$, and the quadratic form is
unchanged.

### Independent configurations

$$
\hat\Omega = \frac{1}{J} \sum_{j=1}^{J} \left( \delta_j - \bar\delta \right) \left( \delta_j - \bar\delta \right)' .
$$

This is the multivariate analogue of the paired $t$ test. It assumes the
configurations are independent. They were drawn by a Latin hypercube
(`docs/ABM_SYSTEM_SIMULATION.md`), which is not independent sampling; for a
response additive in the parameters the Latin hypercube mean has smaller
variance than the independent one, which would make this test conservative,
but $\bar v_{j,m}$ is not known to be additive, so that is not guaranteed.

### Spatial HAC

Let $x_j \in [0, 1]^9$ be configuration $j$'s nine design parameters, each
rescaled by the lower and upper bounds of the design box, read from
`dataset/abm_system_design.csv`. With a kernel $K$ and a bandwidth $h$,

$$
\hat\Omega_h = \frac{1}{J} \sum_{j=1}^{J} \sum_{l=1}^{J} K\!\left( \frac{\lVert x_j - x_l \rVert}{h} \right) \left( \delta_j - \bar\delta \right) \left( \delta_l - \bar\delta \right)' ,
$$

the estimator of Conley (1999) with distance in parameter space in place of
geographic distance. As $h \to 0$ only the terms $j = l$ survive and
$\hat\Omega_h$ becomes the independent estimator above.

- Kernel. $\hat\Omega_h$ is positive semidefinite for every set of points only
  if $K(\lVert \cdot \rVert)$ is a positive definite function on
  $\mathbb{R}^9$. The Gaussian $K(u) = e^{-u^2/2}$ is. The triangular
  $K(u) = \max(0, 1 - u)$, the Newey-West choice on a time axis, is not
  guaranteed to be in nine dimensions. The plan is the Gaussian.
- Bandwidth. Not fixed yet. 1000 points in nine dimensions are sparse, so the
  distribution of nearest-neighbour distances between configurations is
  computed first, and $W_1$ is reported over a range of $h$ read against it,
  $h = 0$ included.
- Assumption. The correlation between $\delta_j$ and $\delta_l$ goes to zero as
  $\lVert x_j - x_l \rVert$ grows, and its range is small compared with the
  box, so that the box holds many groups of configurations that are nearly
  independent of each other. Conley's result is for a region that grows with
  the number of points; here the box is fixed and the points fill it, so the
  $\chi^2$ reference is an approximation that is good only under this
  assumption.

### Which methods

If $H_0$ is rejected, the Model Confidence Set over the five methods says which
ones are behind: loss $1 - \bar v_{j,m}$, the $J$ configurations as
observations, et_al's `MCS_TR` at $\alpha = 0.05$ with block length 1, seed 123.
It returns the methods whose mean $v$ cannot be told apart from the highest.
Its first equivalence test is a test of null 1. Block length 1 draws
configurations independently, so the set carries the independence assumption;
a version allowing for local dependence would need a bootstrap that resamples
neighbourhoods of the box, which is not planned.

## Null 2: no method systematically ahead

Let $R_{j,m}$ be the rank of $\bar v_{j,m}$ among the $M$ methods within
configuration $j$, 1 for the lowest $v$ and $M$ for the highest, ties given
their average rank, and $\bar R_m = \frac{1}{J} \sum_{j} R_{j,m}$.

$$
H_0: \; \text{within every configuration, all } M! \text{ orderings of } \bar v_{j,1}, \dots, \bar v_{j,M} \text{ are equally likely.}
$$

It implies $E[R_{j,m}] = \frac{M+1}{2} = 3$ for every method.

### Independent configurations: the Friedman test

$$
Q = \frac{12 J}{M (M + 1)} \sum_{m=1}^{M} \left( \bar R_m - \frac{M+1}{2} \right)^2 \;\xrightarrow{d}\; \chi^2_{M-1}
\quad \text{under } H_0 ,
$$

with the usual correction when there are ties. With $M = 2$ it reduces to the
sign test of the share of configurations in which one method has the higher
$\bar v$. It assumes the configurations are independent, as above. No normality
is needed.

### Spatial HAC

The same Wald construction as null 1, on ranks: contrasts
$\rho_j = (R_{j,2} - R_{j,1}, \dots, R_{j,M} - R_{j,1})'$, their mean
$\bar\rho$, $\hat\Omega_{\rho,h}$ from the spatial estimator with the same
kernel and bandwidths, and

$$
W_2 = J \, \bar\rho' \, \hat\Omega_{\rho,h}^{-1} \, \bar\rho \;\xrightarrow{d}\; \chi^2_{M-1}.
$$

The null this tests is equal mean ranks, which is weaker than every ordering
being equally likely: Friedman takes the variance of the ranks from that
stronger null, $W_2$ estimates it from the data.

### Which methods

The Model Confidence Set with loss $M - R_{j,m}$ (0 for the method with the
highest $v$ in configuration $j$) tests equal mean ranks at each step and
returns the methods whose mean rank cannot be told apart from the best. Same
settings and the same independence caveat as for null 1.

## What has to be written

In et_al, since these are general statistics and do not belong in this project:

- the Friedman test, with the tie correction; et_al has none;
- the spatial HAC covariance of an $n \times d$ sample given the points'
  coordinates, a kernel and a bandwidth; et_al has the Newey-West long-run
  covariance over a time index, not one over distances between points.

In this project:

- `montecarlo/fit_metric_comparison.c`: reads the oracle file and the design
  file, forms $\bar v_{j,m}$ and $R_{j,m}$, computes $W_1$ and $Q$ and $W_2$
  over the bandwidths, runs the two Model Confidence Sets, and writes
  `montecarlo/out/fit_metric_comparison.txt` and `.tex`;
- `tests/fit_metric_comparison_correctness.c`: $Q$ on a hand-computed case;
  $W_1$ the same under every reference method; $\hat\Omega_h$ equal to the
  independent estimator at $h = 0$; on simulated independent data with equal
  means and exchangeable orderings, rejection rates near $\alpha$ for $W_1$,
  $Q$ and $W_2$.

## Decisions left open

- The bandwidth range, once the nearest-neighbour distances are known.
- Whether to run the same tests on the bootstrap $v^*$ of
  `montecarlo/out/fit_metric_bootstrap.csv`, which already has one benchmark
  per configuration (run 0), so no averaging over runs is needed. Given how far
  the local projections' bootstrap variance is from the oracle one, its
  results would describe the bootstrap rather than the methods.
