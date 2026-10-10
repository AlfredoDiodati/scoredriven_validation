# Comparing the fit metric across auxiliary models: test design

Status: design. Nothing in this file is implemented yet.

## The question

The oracle fit metric of `docs/MONTECARLO_FIT_METRIC.md` gives, for each of the
five auxiliary models and each benchmark of the 1000 by 1000 Monte Carlo, one
value of $v$. Its medians over all benchmarks are 0.774 for the t-QVARMA and
0.692, 0.706, 0.706 and 0.703 for the linear, state 1, state 2 and both-states
local projections. The question is whether the auxiliary models differ by more
than chance, answered by one test over all five at once rather than by pairwise
tests.

Two nulls are tested, separately from each other:

1. the five auxiliary models have the same mean $v$;
2. the auxiliary model has no effect on $v$.

What a difference means. $v$ is high when the best CoP in the final confidence
set is close to the benchmark compared with the true CoP's sampling standard
deviation. It is not a measure of finding the CoP the data came from: within
every auxiliary model, $v$ is higher on the benchmarks where another CoP wins
than on those where the true one does (`docs/MONTECARLO_FIT_METRIC.md`,
"Results"). A test on $v$ therefore says which auxiliary model's best CoP fits
more tightly. How often each auxiliary model recovers the true CoP is the
"smallest mean loss" count of the grid, a different question.

## Notation and model

- $i = 1, \dots, I$, the auxiliary model, $I = 5$: t-QVARMA, linear, state 1,
  state 2, both states. State 1, state 2 and both states are outputs of one
  local-projection fit, so the treatment is which impulse responses are
  compared.
- $j = 1, \dots, J$, the CoP, $J = 1000$.
- $v_{ij}$, the $v$ statistic of CoP $j$ under auxiliary model $i$.

The Monte Carlo also has a third index, the run $k = 1, \dots, S$ with
$S = 1000$: $v_{ijk}$ is the oracle value when run $k$ of CoP $j$ is the
benchmark, read from `montecarlo/out/fit_metric_oracle.csv.gz`. It is left out
of the notation: the replicates $k$ are used to construct the p-values. The
$v_{ijk}$ of different runs of the same CoP are dependent, because the
denominator $\sigma$ of each is computed from the CoP's other runs.

The design is a randomized complete block design: the auxiliary models are the
treatments, the CoPs the blocks, one observation per cell. In general form,

$$
v_{ij} = f\left( \tau_i, \beta_j, \varepsilon_{ij} \right),
$$

- $\tau_i$, the effect on $v_{ij}$ of using auxiliary model $i$;
- $\beta_j$, the effect on $v_{ij}$ of being in CoP $j$;
- $\varepsilon_{ij}$, an arbitrary error term, potentially dependent across
  blocks. In this Monte Carlo, run $k$ of CoP $j$ is the same simulated data
  for every auxiliary model $i$, and run $k$ of every CoP is simulated from
  the same seed, $k$ (`applications/abm_system_simulate.c`). The shared seed
  induced no measurable correlation across CoPs: over all 499,500 CoP pairs the
  correlation of the per-run t-QVARMA losses against the US responses is
  0.0017 on average at the same seed against -0.0001 with the seeds shifted by
  one, because the model consumes a different number of random draws at
  different parameters and the streams desynchronise
  (`docs/ABM_SYSTEM_SIMULATION.md`, "Common random numbers"). That
  measurement is on those losses, not on $v$.

The comparison is paired: every CoP has a value under every auxiliary model.

## Null 1: equal means

$$
H_0: \; E\left[ v_{1j} \right] = E\left[ v_{2j} \right] = \dots = E\left[ v_{Ij} \right],
$$

the expectation taken over CoPs drawn from the design box.

Take auxiliary model 1 as reference and form the $I - 1 = 4$ contrasts

$$
\delta_j = \left( v_{2j} - v_{1j}, \; \dots, \; v_{Ij} - v_{1j} \right)',
\qquad
\bar\delta = \frac{1}{J} \sum_{j=1}^{J} \delta_j .
$$

The contrasts remove whatever is common to the auxiliary models within a CoP,
such as an additive $\beta_j$. The Wald statistic is

$$
W_1 = J \, \bar\delta' \, \hat\Omega^{-1} \, \bar\delta \;\xrightarrow{d}\; \chi^2_{I-1}
\quad \text{under } H_0,
$$

computed by solving $\hat\Omega x = \bar\delta$ through a Cholesky factor, never
by inverting $\hat\Omega$. $W_1$ does not depend on which auxiliary model is the
reference: another reference turns $\delta_j$ into $A \delta_j$ for an
invertible $A$, $\hat\Omega$ into $A \hat\Omega A'$, and the quadratic form is
unchanged.

### Independent contrasts

$$
\hat\Omega = \frac{1}{J} \sum_{j=1}^{J} \left( \delta_j - \bar\delta \right) \left( \delta_j - \bar\delta \right)' .
$$

This is Hotelling's $T^2$, the multivariate analogue of the paired $t$ test.
Its $\chi^2$ reference treats the contrasts $\delta_j$ as independent across
CoPs; the levels $v_{ij}$ may be dependent across CoPs.

The parametric test usually paired with the Friedman test, the F test of the
randomized complete block ANOVA, $MS_{\text{treatment}} / MS_{\text{error}}$ on
$I - 1$ and $(I - 1)(J - 1)$ degrees of freedom, tests the same null but also
needs sphericity: the difference between any two auxiliary models has the same
variance. That is unlikely here, since the t-QVARMA's $v$ is spread differently
from the local projections' and state 1, state 2 and both states come from one
fit. Hotelling's $T^2$ estimates the full covariance of the contrasts and does
not need it; the Greenhouse-Geisser correction of the F test's degrees of
freedom is the alternative.

### Spatial HAC

Let $x_j \in [0, 1]^9$ be CoP $j$'s nine design parameters, each rescaled by
the lower and upper bounds of the design box, read from
`dataset/abm_system_design.csv`. With a kernel $K$ and a bandwidth $h$,

$$
\hat\Omega_h = \frac{1}{J} \sum_{j=1}^{J} \sum_{j'=1}^{J} K\!\left( \frac{\lVert x_j - x_{j'} \rVert}{h} \right) \left( \delta_j - \bar\delta \right) \left( \delta_{j'} - \bar\delta \right)' ,
$$

the estimator of Conley (1999) with distance in parameter space in place of
geographic distance. As $h \to 0$ only the terms $j = j'$ survive and
$\hat\Omega_h$ becomes the independent estimator above. It allows the
contrasts of nearby CoPs to be dependent, for example a region of the
parameter space where the t-QVARMA leads and another where the local
projections do.

- Kernel. $\hat\Omega_h$ is positive semidefinite for every set of points only
  if $K(\lVert \cdot \rVert)$ is a positive definite function on
  $\mathbb{R}^9$. The Gaussian $K(u) = e^{-u^2/2}$ is. The triangular
  $K(u) = \max(0, 1 - u)$, the Newey-West choice on a time axis, is not
  guaranteed to be in nine dimensions. The plan is the Gaussian.
- Bandwidth. Not fixed yet. 1000 points in nine dimensions are sparse, so the
  distribution of nearest-neighbour distances between CoPs is computed first,
  and $W_1$ is reported over a range of $h$ read against it, $h = 0$ included.
- Assumption. The correlation between $\delta_j$ and $\delta_{j'}$ goes to zero
  as $\lVert x_j - x_{j'} \rVert$ grows, and its range is small compared with
  the box. Conley's result is for a region that grows with the number of
  points; here the box is fixed and the points fill it, so the $\chi^2$
  reference is an approximation that is good only under this assumption.

### Which auxiliary models

If $H_0$ is rejected, the Model Confidence Set over the five auxiliary models
says which ones are behind: loss $1 - v_{ij}$, the $J$ CoPs as observations,
et_al's `MCS_TR` at $\alpha = 0.05$ with block length 1, seed 123. It returns
the auxiliary models whose mean $v$ cannot be told apart from the highest. Its
first equivalence test is a test of null 1. Block length 1 draws CoPs
independently, so the set treats the contrasts as independent across CoPs; a
version allowing for local dependence would need a bootstrap that resamples
neighbourhoods of the box, which is not planned.

## Null 2: the auxiliary model has no effect

$$
H_0: \; \text{there exists a function } g_j \text{ such that } f\left( \tau_i, \beta_j, \varepsilon_{ij} \right) = g_j\left( \beta_j, \varepsilon_{ij} \right) \text{ for every admissible pair } \left( \beta_j, \varepsilon_{ij} \right),
$$

meaning the auxiliary model effect $\tau_i$ is not doing anything.

Let $R_{ij}$ be the rank of $v_{ij}$ among the $I$ auxiliary models within CoP
$j$, 1 for the lowest $v$ and $I$ for the highest, ties given their average
rank, and $\bar R_i = \frac{1}{J} \sum_{j} R_{ij}$.

Under $H_0$ the ordering within CoP $j$ depends on
$(\varepsilon_{1j}, \dots, \varepsilon_{Ij})$ only, through the same $g_j$ for
every auxiliary model. Every ordering is equally likely when that vector is
exchangeable across $i$; dependence across $i$ within a CoP is allowed. Then
$E[R_{ij}] = \frac{I+1}{2} = 3$ for every auxiliary model.

### The Friedman test

$$
Q = \frac{12 J}{I (I + 1)} \sum_{i=1}^{I} \left( \bar R_i - \frac{I+1}{2} \right)^2 ,
$$

with the usual correction when there are ties. $Q$ measures how far each
auxiliary model's mean rank is from the common mean rank $\frac{I+1}{2}$. With
$I = 2$ it reduces to the sign test. No normality and no additive form of $f$
are needed.

The reference distribution. The asymptotic $\chi^2_{I-1}$ treats the rankings
of different CoPs as independent: $\mathrm{Var}(\bar R_i)$ is then
$\frac{1}{J^2} \sum_j \mathrm{Var}(R_{ij})$, without the covariances between
CoPs. The levels $v_{ij}$, and the errors $\varepsilon_{ij}$, may be dependent
across CoPs. Here the p-value is to be constructed from the replicates $k$
instead (see "Decisions left open").

The Kruskal-Wallis test does not fit: it treats the auxiliary models as
independent samples, while here every auxiliary model is measured on the same
CoPs and runs, so it would ignore the pairing. Friedman is for matched blocks,
which is this design.

### Spatial HAC

The same Wald construction as null 1, on ranks: contrasts
$\rho_j = (R_{2j} - R_{1j}, \dots, R_{Ij} - R_{1j})'$, their mean $\bar\rho$,
$\hat\Omega_{\rho,h}$ from the spatial estimator with the same kernel and
bandwidths, and

$$
W_2 = J \, \bar\rho' \, \hat\Omega_{\rho,h}^{-1} \, \bar\rho \;\xrightarrow{d}\; \chi^2_{I-1}.
$$

The null this tests is equal mean ranks, which is weaker than $H_0$: Friedman
takes the variance of the ranks from $H_0$, $W_2$ estimates it from the data.

### Which auxiliary models

The Model Confidence Set with loss $I - R_{ij}$ (0 for the auxiliary model with
the highest $v$ in CoP $j$) tests equal mean ranks at each step and returns the
auxiliary models whose mean rank cannot be told apart from the best. Same
settings as for null 1.

## What has to be written

In et_al, since these are general statistics and do not belong in this project:

- the Friedman test, with the tie correction; et_al has none;
- the spatial HAC covariance of an $n \times d$ sample given the points'
  coordinates, a kernel and a bandwidth; et_al has the Newey-West long-run
  covariance over a time index, not one over distances between points.

In this project:

- `montecarlo/fit_metric_comparison.c`: reads the oracle file and the design
  file, forms $v_{ij}$ and $R_{ij}$, computes $W_1$, $Q$ and $W_2$ with their
  p-values, runs the two Model Confidence Sets, and writes
  `montecarlo/out/fit_metric_comparison.txt` and `.tex`;
- `tests/fit_metric_comparison_correctness.c`: $Q$ on a hand-computed case;
  $W_1$ the same under every reference auxiliary model; $\hat\Omega_h$ equal to
  the independent estimator at $h = 0$; on simulated data satisfying each
  null, rejection rates near $\alpha$.

## Decisions left open

- How the replicates $k$ construct the p-values, given that the $v_{ijk}$ of
  different runs of the same CoP are dependent through $\sigma$. One option is
  to average over runs, $v_{ij} = \frac{1}{S} \sum_k v_{ijk}$, and use the
  $\chi^2$ references above; the dependence across $k$ then stays inside each
  $v_{ij}$.
- The bandwidth range, once the nearest-neighbour distances are known.
- Whether to run the same tests on the bootstrap $v^*$ of
  `montecarlo/out/fit_metric_bootstrap.csv`, which has one benchmark per CoP
  (run 0). Given how far the local projections' bootstrap variance is from the
  oracle one, its results would describe the bootstrap rather than the
  auxiliary models.
