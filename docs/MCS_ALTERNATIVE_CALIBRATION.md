# Alternative calibration procedure

How a set of simulation-model configurations is calibrated against a
benchmark series by comparing impulse response functions. A configuration, a
"CoP" (Configuration of Parameters), is one point of the parameter design,
simulated for $n$ independent replicates.

In short: take the estimated impulse response function of the benchmark and
of every simulated replicate of every configuration, and measure each
replicate's distance from the benchmark's response by the squared error at
each entry of the response minus the replicate's own scaled squared deviation
from its configuration's mean response, averaged over the entries. Then run
the Model Confidence Set of Hansen, Lunde and Nason (2011) over the
configurations, one model per configuration and one observation per
replicate, using the equivalence test of Seri et al. (2021).

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
L^*_{ish} = \left( \widehat{IR}_{ish} - \widehat{IR}_{dh} \right)^2
- \frac{n}{n-1} \left( \widehat{IR}_{ish} - \overline{IR}_{ih} \right)^2 .
$$

The loss of the replicate is its average over the entries:

$$
L^*_{is} = \frac{1}{K} \sum_{h=1}^{K} L^*_{ish}.
$$

**Step 3.** Run the Model Confidence Set over the configurations on the
per-replicate losses $L^*_{is}$. The losses are not averaged first: the
confidence set estimates the mean and variance of each configuration's loss
from its repeated observations, and an averaged number has none.

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
E\left[ L_{ish} \mid \widehat{IR}_d \right]
= \left( IR_{ih} - \widehat{IR}_{dh} \right)^2 + \sigma^2_{ih}.
$$
The Model Confidence Set compares two configurations $i$ and $j$ through the
difference of their losses, whose expectation is
$$
E\left[ L_{ish} - L_{jsh} \mid \widehat{IR}_d \right]
= \left( IR_{ih} - \widehat{IR}_{dh} \right)^2
- \left( IR_{jh} - \widehat{IR}_{dh} \right)^2
+ \sigma^2_{ih} - \sigma^2_{jh}.
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

Each configuration contributes one model, and each replicate one
observation. The confidence set is built by the sequential procedure of
Hansen, Lunde and Nason (2011):

1. Start with the set of all configurations.
2. Test, at level $\alpha = 0.05$, the null hypothesis that every
   configuration in the set has the same expected loss, against the
   alternative that at least two differ.
3. If the null is rejected, remove the configuration with the largest mean
   loss and return to step 2.
4. When the null is no longer rejected, the configurations left form the
   confidence set.

The mean and variance of each configuration's loss, on which the test in
step 2 is based, are estimated by Gaussian quasi-likelihood, following Seri
et al. (2021). As set out by Martinoli et al. (2024), this estimation
requires that the vectors of losses across configurations be independent and
identically distributed over the replicates, and that each configuration's
loss have finite mean and variance. Independent replicates simulated from the
same configuration satisfy these, up to the $1/n$ dependence introduced by
the correction.

## References

Hansen, P. R., Lunde, A., and Nason, J. M. (2011). The model confidence set.
*Econometrica*, 79(2), 453–497.

Martinoli, M., Moneta, A., and Pallante, G. (2024). Calibration and
validation of macroeconomic simulation models by statistical causal search.
*Journal of Economic Behavior & Organization*, 228, 106786.

Seri, R., Martinoli, M., Secchi, D., and Centorrino, S. (2021). Model
calibration and validation via confidence sets. *Econometrics and
Statistics*, 20, 62–86.
