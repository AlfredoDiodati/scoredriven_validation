# Alternative calibration procedure

How a set of simulation-model configurations is calibrated against a
benchmark series by comparing impulse response functions. A configuration, a
"CoP" (Configuration of Parameters), is one point of the parameter design,
simulated for $n$ independent replicates.

## The steps

Notation: $d$ is the benchmark, $i$ a configuration, $s = 1, \dots, n$ the
replicates of a configuration, and $h = 1, \dots, K$ the entries of the
stacked response vector, every response at every horizon.
$\widehat{IR}_{d}$ is the benchmark's estimated response and
$\widehat{IR}_{is}$ the response estimated on replicate $s$ of
configuration $i$.

**Step 1.** For each configuration and entry, compute the mean response
across the configuration's $n$ replicates:

$$
\overline{IR}_{ih} = \frac{1}{n} \sum_{s=1}^{n} \widehat{IR}_{ish}.
$$

**Step 2: the corrected loss.** The loss of a replicate at entry $h$ is the
squared error minus the replicate's own squared deviation from its
configuration's mean, scaled by $n/(n-1)$:

$$
L^*_{ish} = \left( \widehat{IR}_{ish} - \widehat{IR}_{dh} \right)^2 - \frac{n}{n-1} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2
$$

The loss of the replicate is its average over the entries:

$$L^*_{is} = \frac{1}{K} \sum_{h=1}^{K} L^*_{ish}$$

**Step 3.** Run the Model Confidence Set over the configurations on the
per-replicate losses $L^*_{is}$, as described in the section on the Model
Confidence Set below. The losses are not averaged first: the confidence set
estimates the mean and variance of each configuration's loss from its
repeated observations, and an averaged number has none.

A loss $L^*_{is}$ can be negative. The confidence set reads only differences
between configurations' losses, so this changes nothing.

When the benchmark is itself a simulated replicate, it must not be among the
replicates used in steps 1 and 2.

## Why the correction

Let
$$
L_{ish} = \left( \widehat{IR}_{ish} - \widehat{IR}_{dh} \right)^2
$$
be the squared error between the benchmark's response and the response of
replicate $s$ of configuration $i$ at entry $h$. The simulated response is an
estimate: across replicates it is a random variable with mean $IR_{ih}$ and
variance $\sigma^2_{ih}$. With the benchmark fixed,
$$
E\left[ L_{ish} \mid \widehat{IR}_d \right] = \left( IR_{ih} - \widehat{IR}_{dh} \right)^2 + \sigma^2_{ih}.
$$
The Model Confidence Set compares two configurations $i$ and $j$ through the
difference of their losses, whose expectation is
$$
E\left[ L_{ish} - L_{jsh} \mid \widehat{IR}_d \right] = \left( IR_{ih} - \widehat{IR}_{dh} \right)^2 - \left( IR_{jh} - \widehat{IR}_{dh} \right)^2 + \sigma^2_{ih} - \sigma^2_{jh}.
$$
The last terms, $\sigma^2_{ih} - \sigma^2_{jh}$, relate only to the
estimation noise of the simulations, not to their distance from the
benchmark. The Model Confidence Set could therefore favour a configuration
whose mean response is far from the benchmark's but whose responses are less
noisy, giving a biased confidence set. Averaging the losses over the
replicates does not remove the term: the mean of $L_{ish}$ over $s$ still
contains the sample variance of the responses.

### The corrected loss removes the bias

Writing $\widehat{IR}_{ish} - \widehat{IR}_{dh}$ as
$(\widehat{IR}_{ish} - \overline{IR}_{ih}) + (\overline{IR}_{ih} - \widehat{IR}_{dh})$
and expanding the square,
$$
L^*_{ish} = \left( \overline{IR}_{ih} - \widehat{IR}_{dh} \right)^2
+ 2 \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right) \left( \overline{IR}_{ih} - \widehat{IR}_{dh} \right)
- \frac{1}{n-1} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2 .
$$
Averaged over the replicates, the middle term vanishes and
$$
\frac{1}{n} \sum_{s=1}^{n} L^*_{ish}
= \left( \overline{IR}_{ih} - \widehat{IR}_{dh} \right)^2 - \frac{\hat\sigma^2_{ih}}{n},
\qquad
\hat\sigma^2_{ih} = \frac{1}{n-1} \sum_{s=1}^{n} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2 .
$$
Since
$$
E\left[ \left( \overline{IR}_{ih} - \widehat{IR}_{dh} \right)^2 \mid \widehat{IR}_d \right]
= \left( IR_{ih} - \widehat{IR}_{dh} \right)^2 + \frac{\sigma^2_{ih}}{n}
$$
and $E[\hat\sigma^2_{ih}] = \sigma^2_{ih}$, the expected mean loss is
$$
E\left[ \frac{1}{n} \sum_{s=1}^{n} L^*_{ish} \mid \widehat{IR}_d \right]
= \left( IR_{ih} - \widehat{IR}_{dh} \right)^2,
$$
the distance term alone. The difference between two configurations'
corrected losses therefore has expectation
$$
\left( IR_{ih} - \widehat{IR}_{dh} \right)^2 - \left( IR_{jh} - \widehat{IR}_{dh} \right)^2 ,
$$
with no contribution from the estimation noise of the simulations.
Aggregating over the entries follows the same pattern.

The scaling $n/(n-1)$ is what makes the removal exact. Without it, that is
with $L_{ish} - (\widehat{IR}_{ish} - \overline{IR}_{ih})^2$, the mean loss
is $(\overline{IR}_{ih} - \widehat{IR}_{dh})^2$, whose expectation still
contains $\sigma^2_{ih}/n$.

### Why each replicate's own deviation

Subtracting the common estimate $\hat\sigma^2_{ih}$ from every replicate,
$L_{ish} - \hat\sigma^2_{ih}$, would remove the bias in the mean as well. It
would not give the right estimate of the mean's sampling variance. The
confidence set estimates that variance from the spread of the per-replicate
losses, and with a common term subtracted each replicate's loss still
contains its squared deviation $(\widehat{IR}_{ish} - \overline{IR}_{ih})^2$.
That term varies a great deal across replicates, but in the mean it cancels
against $\hat\sigma^2_{ih}$, which is computed from the same replicates. Its
spread would inflate the estimated variance of every configuration's mean
loss, and with it the confidence set.

With each replicate's own deviation subtracted, the part of $L^*_{ish}$ that
varies across replicates is the middle term of the expansion, which is linear
in the replicate's deviation, plus the last term, whose weight $1/(n-1)$ is
negligible. The spread of the losses then reflects only the uncertainty in
the configuration's mean response.

The replicates' losses are no longer exactly independent, since they share
$\overline{IR}_{ih}$. The dependence is of order $1/n$.

## The Model Confidence Set

The confidence set is built by the sequential testing and elimination of
Hansen, Lunde and Nason (2011), with the equivalence test and the
elimination rule of Seri et al. (2021, Section 4.2). Each configuration
contributes one model, and each replicate one observation.

### Mean and variance of each configuration's loss

For each configuration $i$, the mean and variance of its loss are estimated
from its $n$ replicates:

$$
\bar L_i = \frac{1}{n} \sum_{s=1}^{n} L^*_{is},
\qquad
\hat\sigma^2_i = \frac{1}{n} \sum_{s=1}^{n} \left( L^*_{is} \right)^2 - \bar L_i^2 .
$$

### The equivalence test

Let $M$ be a set of $m$ configurations, labelled $1, \dots, m$ in any order,
and $\bar{\mathbf L} = (\bar L_1, \dots, \bar L_m)'$. The null hypothesis is
that all configurations in $M$ have the same expected loss, and the
alternative is that at least two differ:

$$
H_{0,M}: E[L^*_{is}] = E[L^*_{js}] \;\; \forall i, j \in M,
\qquad
H_{A,M}: \exists\, i, j \in M : E[L^*_{is}] \neq E[L^*_{js}] .
$$

Write the null as $A \bar{\mathbf L} = \mathbf 0$ with the
$(m-1) \times m$ contrast matrix

$$
A = \begin{bmatrix} \mathbf u_{m-1} & -I_{m-1} \end{bmatrix},
$$

where $\mathbf u_{m-1}$ is a vector of $m-1$ ones and $I_{m-1}$ the identity,
so that the entries of $A \bar{\mathbf L}$ are $\bar L_1 - \bar L_k$ for
$k = 2, \dots, m$. The test statistic is

$$
W_M = n \left( A \bar{\mathbf L} \right)'
\left( \hat\sigma^2_1 U_{m-1} + \hat\Sigma_{-1} \right)^{-1}
\left( A \bar{\mathbf L} \right),
$$

where $U_{m-1}$ is the $(m-1) \times (m-1)$ matrix of ones and
$\hat\Sigma_{-1} = \mathrm{diag}(\hat\sigma^2_2, \dots, \hat\sigma^2_m)$.
Under the null, $W_M$ is asymptotically $\chi^2_{m-1}$, and the p-value of
the test is

$$
p_{H_0,M} = 1 - F_{\chi^2_{m-1}}\left( W_M \right).
$$

The same statistic can be computed without inverting a matrix. With the
precision-weighted mean loss

$$
\tilde L = \frac{\sum_{i \in M} \bar L_i / \hat\sigma^2_i}{\sum_{i \in M} 1 / \hat\sigma^2_i},
$$

it equals

$$
W_M = n \sum_{i \in M} \frac{\left( \bar L_i - \tilde L \right)^2}{\hat\sigma^2_i}.
$$

This form does not depend on which configuration is labelled first.

### Elimination and the confidence set

When the null is rejected, the configuration with the largest mean loss is
removed:

$$
e_M = \arg\max_{i \in M} \bar L_i .
$$

Since the elimination rule depends only on the mean losses, the order of
elimination is the order of the configurations by decreasing $\bar L_i$.
Sort the configurations in that order, $e_1, e_2, \dots, e_{m_0}$, where
$m_0$ is the total number of configurations, and let
$M_k = \{e_k, e_{k+1}, \dots, e_{m_0}\}$ be the set left before the $k$-th
elimination. For each $k$, compute the p-value $p_{H_0,M_k}$ of the test on
$M_k$, with the convention $p_{H_0,M_{m_0}} = 1$. The p-values of the test
are not necessarily decreasing in $k$, so the MCS p-value of configuration
$e_k$ is

$$
\hat p_{e_k} = \max_{j \le k} \, p_{H_0,M_j} .
$$

The Model Confidence Set at level $1-\alpha$ is the set of configurations
with $\hat p_{e_k} \ge \alpha$. It is the set obtained by starting from all
configurations, testing, removing $e_M$ while the null is rejected at level
$\alpha$, and stopping at the first set where it is not. The MCS p-values
give the confidence set for every $\alpha$ at once.

### Requirements

The test requires that every simulation run be independent of all the
others, across replicates and across configurations, that the runs of a
configuration be identically distributed, and that each configuration's
loss have finite mean and variance (Seri et al., 2021, assumptions A1 and
A5). The diagonal form of $\hat\Sigma_{-1}$ relies on the independence of
the runs of different configurations. The correction of step 2 adds a
dependence of order $1/n$ among the replicates of one configuration, through
$\overline{IR}_{ih}$.

Under these assumptions, as $n$ grows, the confidence set contains the
configurations with the smallest expected loss with probability at least
$1-\alpha$, and excludes every other configuration with probability
approaching one (Seri et al., 2021, Theorem 2).

## References

Hansen, P. R., Lunde, A., and Nason, J. M. (2011). The model confidence set.
*Econometrica*, 79(2), 453–497.

Seri, R., Martinoli, M., Secchi, D., and Centorrino, S. (2021). Model
calibration and validation via confidence sets. *Econometrics and
Statistics*, 20, 62–86.
