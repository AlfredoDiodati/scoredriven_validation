# Comparing the fit metric across auxiliary models: test design

## Status: discarded

This test design compares the fit metric of
`docs/OLD_MONTECARLO_FIT_METRIC.md`, which measures the winner's distance with
the naive loss, the mean absolute error between impulse-response vectors with
no variance correction, which biases the confidence set towards configurations
whose responses are less noisy (`docs/MCS_VALIDATION_PROCEDURE.md`, "Why the
variance correction"). Its code and outputs are under `old/naive_mae_loss/`,
and `docs/OLD_EXPERIMENTS.md` lists every discarded experiment. The text
below is the write-up as it stood when it was discarded.

Status: draft, in progress. The test for null 1 below is not designed yet, and
nothing in this file is implemented.
`old/naive_mae_loss/montecarlo/fit_metric_pairwise.c` (with its table writer
`old/naive_mae_loss/montecarlo/fit_metric_table.h`) is an earlier attempt that
tests the auxiliary models two at a time (et_al's `dm_test` on the differences
of $v$, across CoPs and per CoP with a Bonferroni correction); it builds with
`make old-fit_metric_pairwise` and has not produced results.

## The question

The oracle fit metric of `docs/OLD_MONTECARLO_FIT_METRIC.md` gives, for each of the
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
than on those where the true one does (`docs/OLD_MONTECARLO_FIT_METRIC.md`,
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
benchmark, read from `old/naive_mae_loss/montecarlo/out/fit_metric_oracle.csv.gz`. It is left out
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

## Null 1: the same mean $v$

Not designed yet.

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
across CoPs.

In the Monte Carlo a second p-value can be built without the $\chi^2_{I-1}$
approximation, by resampling the replicates $k$. Comparing the two p-values
says whether the $\chi^2_{I-1}$ p-value can be trusted in the empirical
setting, where this resampling is not possible and it is the only p-value
available.
