# A fit metric: the winner's distance against the sampling noise

## What it measures

The Model Confidence Set says which configurations cannot be told apart from
the best one. It does not say whether the best one is close to the data in any
absolute sense: a field of configurations that all miss the data by a wide
margin still has a winner, and still a confidence set around it.

The fit metric puts the winner's distance on a scale. The distance between a
configuration's impulse responses and the data's is compared with how much the
impulse responses themselves move from one sample to another. A winner whose
distance is small next to that sampling variation fits the data about as well
as any estimate could; a winner whose distance is large next to it does not.

It has two versions. The oracle version measures the sampling variation across
independent samples of the true configuration, which only a simulation has.
The feasible version measures it by resampling the one series the data
provides, which can be done on the US data as well. In the Monte Carlo
experiment both can be computed for the same benchmark, so comparing them
says whether the feasible version measures the uncertainty correctly.

## Definition

### Notation

- $j$ is the true configuration, the one the benchmark series was simulated
  from, and $m$ the method that turns a series into an impulse-response vector
  (the t-QVARMA or one of the local projections).
- $s = 1, \dots, S$ are the simulated runs of configuration $j$, and $s^*$ the
  run used as the validation data, the role the US data plays.
- $\widehat{ir}^{(s)}_{j,m} \in \mathbb{R}^K$ is the stacked impulse-response
  vector estimated with method $m$ on run $s$, with entries
  $\widehat{ir}^{(s)}_{j,m,k}$ for $k = 1, \dots, K$.

### Oracle version: the variance across the true configuration's runs

Use the runs of configuration $j$ other than the validation run. For each
entry $k$, the mean and the variance across those runs are

$$
\overline{ir}_{j,m,k} = \frac{1}{S-1} \sum_{s \neq s^*} \widehat{ir}^{(s)}_{j,m,k},
\qquad
\widehat{\mathrm{Var}}_{j,m,k} = \frac{1}{S-2} \sum_{s \neq s^*} \left( \widehat{ir}^{(s)}_{j,m,k} - \overline{ir}_{j,m,k} \right)^2 ,
$$

and the average over the $K$ entries is

$$
\sigma^2_{j,m} = \frac{1}{K} \sum_{k=1}^{K} \widehat{\mathrm{Var}}_{j,m,k}.
$$

This is the sampling variance of the impulse-response estimator under
configuration $j$. It is available only because the data-generating process is
known.

### Feasible version: the variance across bootstrap resamples of the validation run

Take only the validation series $y^{(s^*)}_j$ and draw $B$ bootstrap resamples
$y^{*(b)}_j$, $b = 1, \dots, B$, by a bootstrap that respects the time
dependence: a block bootstrap, or a residual bootstrap. Estimate the impulse
responses on each resample to get $\widehat{ir}^{*(b)}_{j,m}$ and apply the
same formulas over $b$:

$$
\overline{ir}^{*}_{j,m,k} = \frac{1}{B} \sum_{b=1}^{B} \widehat{ir}^{*(b)}_{j,m,k},
\qquad
\widehat{\mathrm{Var}}^{*}_{j,m,k} = \frac{1}{B-1} \sum_{b=1}^{B} \left( \widehat{ir}^{*(b)}_{j,m,k} - \overline{ir}^{*}_{j,m,k} \right)^2 ,
\qquad
\sigma^{2*}_{j,m} = \frac{1}{K} \sum_{k=1}^{K} \widehat{\mathrm{Var}}^{*}_{j,m,k}.
$$

This is what can also be computed on real data, where only one realisation
exists.

### The metric, for either version

$$
s_{j,m} = \frac{d_{j,m}}{\sigma_{j,m}}, \qquad \sigma_{j,m} = \sqrt{\sigma^2_{j,m}}, \qquad v_{j,m} = \frac{1}{1 + s_{j,m}},
$$

where $d_{j,m}$ is the average distance of the best configuration in the final
Model Confidence Set. Here the distance is the mean absolute difference
between two impulse-response vectors,

$$
\frac{1}{K} \sum_{k=1}^{K} \left| a_k - b_k \right|,
$$

which is in the units of the responses, so it is set against a standard
deviation, in the same units, rather than a variance. $s$ is then a pure
number: multiplying every response by a constant $c$ multiplies both $d_{j,m}$
and $\sigma_{j,m}$ by $c$ and leaves $s$ and $v$ as they were, so the values
of different methods can be compared. Had the distance been a mean squared
difference, the variance $\sigma^2_{j,m}$ itself would be its counterpart.

$v$ lies between 0 and 1: it approaches 1 when the winner's distance is small
next to the sampling variation and 0 when it is large. The feasible version
uses $\sigma^{*}_{j,m} = \sqrt{\sigma^{2*}_{j,m}}$ in place of
$\sigma_{j,m}$, giving $s^*_{j,m}$ and $v^*_{j,m}$, with the same $d_{j,m}$.

Comparing $\sigma^2_{j,m}$ with $\sigma^{2*}_{j,m}$ across $j$ shows whether
the bootstrap measures the uncertainty correctly.

## How it is computed here

### The ingredients from the Monte Carlo

Everything is computed on the 1000 by 1000 Monte Carlo of
`docs/MONTECARLO_VALIDATION.md`, "Every run of every configuration as the
benchmark": 1000 configurations, $S = 1000$ runs each, every run the benchmark
once, five methods. The methods and their vector lengths $K$ are those of
`montecarlo/response_cache.h`: the t-QVARMA (`qvarma`, total responses to
horizon 20, $K = 525$) and four local projections, linear, state 1, state 2
($K = 400$ each) and both states, state 1's and state 2's responses stacked
into one vector ($K = 800$).

**$d_{j,m}$** is the mean loss, over the benchmark's 999 rows, of the
configuration with the smallest mean loss. The loss is the one of the
experiment, the mean absolute difference between the two response vectors
over their $K$ entries, with the benchmark's run index held out of every
configuration. The configuration with the smallest mean loss is the one
`montecarlo/out/sweep_grid.csv.gz` records for each benchmark, and it is in the
final confidence set in all 5,000,000 rows, so it is the best configuration in
the set. $d$ is recomputed from the same response caches
(`out/sweep_grid_response_cache/`), over the same runs, with the same sums the
grid formed its loss matrices with.

**$\sigma^2_{j,m}$, the oracle**, is computed for every one of the 5,000,000
benchmarks. For each configuration the mean and the sum of squared deviations
of every entry over all 1000 runs are formed once; a benchmark's run $x$ is
then taken out of both, giving the mean over the other 999 runs,
$\overline{ir}' = (S\,\overline{ir} - x)/(S-1)$, and their sum of squares,
$Q - (x - \overline{ir})(x - \overline{ir}')$, which divided by $S - 2 = 998$
is $\widehat{\mathrm{Var}}_{j,m,k}$. That costs $K$ operations per benchmark
instead of $999K$. The responses are the float32 values of the caches, the
arithmetic is in double.

**$\sigma^{2*}_{j,m}$, the bootstrap**, is computed for run 0 of every
configuration, 1000 benchmarks per method: at $B = 200$ resamples each, all
1,000,000 benchmarks would take about 100 days of this machine, almost all of
it in t-QVARMA fits.

- The resample is a moving-block bootstrap of the benchmark's stored series,
  the five series the t-QVARMA is fitted on (GDP growth, energy growth,
  employment change, inflation, interest rate), 400 periods. Whole periods are
  drawn, so the five series keep their contemporaneous relation. Blocks are 8
  periods long, the cube root of 400 rounded up; block starts are uniform over
  the 393 possible ones, through et_al's `mcs_block_indices`.
- Every method is re-estimated on the same resamples. The local projections go
  through `applications/lp_system.h` and et_al's `lp_lin_and_nl` exactly as
  `applications/abm_system_fit_lp.c` fitted every simulated run, so the
  linear, state 1, state 2 and both-states vectors of one resample come from
  one fit. The t-QVARMA is fitted from the benchmark's own cached estimate
  (`out/abm_system_fit_qvarma/`) as its starting point, with et_al's default
  budget of 4000 iterations, and its responses are computed as the cache
  computes them.
- A resample whose response under a method is missing or not finite is left
  out of that method's $\sigma^{2*}$ and of nothing else.
- The resamples of a benchmark come from et_al's generator at seed 123 and a
  stream fixed by the configuration and the run, so they do not depend on the
  thread that draws them.

Only the block bootstrap is implemented. A residual bootstrap would fit a model
to the series (a VAR, say), resample its residuals and rebuild the series
through the model's recursion from the observed starting values. et_al's
`varima/var.h` can fit the VAR and simulate it with Gaussian shocks
(`var_simulate`), but has no recursion driven by given residuals from given
starting values, which is the general piece a residual bootstrap needs.

### Programs and outputs

| file | what it does |
| --- | --- |
| `montecarlo/fit_metric.h` | $\sigma^2$ over a set of vectors, the oracle shortcut, $d$, $s$ and $v$, the five response vectors of one series, and $\sigma^{2*}$ of one series |
| `montecarlo/fit_metric_oracle.c` | $d$, $\sigma^2$, $s$, $v$ for every benchmark of the grid |
| `montecarlo/fit_metric_bootstrap.c` | $\sigma^{2*}$, $s^*$, $v^*$ for one run of every configuration, beside the oracle values of the same benchmarks |
| `montecarlo/fit_metric_report.c` | the tables below |
| `montecarlo/fit_metric_run.sh` | runs the three in turn, the block-length comparison included, resuming where it stopped |
| `montecarlo/fit_metric_v_histograms.py` | histograms of the oracle $v$ per model, over every benchmark and over `cop_0191`'s runs, into `montecarlo/out/fit_metric_v_histograms/`; needs polars, plotly and kaleido, and is not run by the driver |
| `tests/fit_metric_correctness.c` | the checks under "Verification" |

    ./montecarlo/fit_metric_run.sh

| output | what it holds |
| --- | --- |
| `montecarlo/out/fit_metric_oracle.csv.gz` | one row per method and benchmark, 5,000,000 rows: method, benchmark configuration and run, the configuration with the smallest mean loss, $d$, $\sigma^2$, $s$, $v$; 130 MB, ignored by git |
| `montecarlo/out/fit_metric_bootstrap.csv` | one row per method and configuration, run 0: $B$, block length, usable resamples, converged t-QVARMA refits and their mean iterations, the configuration with the smallest mean loss, $d$, $\sigma^2$, $\sigma^{2*}$, $s$, $v$, $s^*$, $v^*$ |
| `montecarlo/out/fit_metric_bootstrap_run0_B50_block<L>_iterations4000_first100.csv` | the same for the first 100 configurations at 50 resamples and blocks of $L$ periods |
| `montecarlo/out/fit_metric_report.txt`, `.tex` | the tables, as text and as LaTeX tables |

`bin/fit_metric_bootstrap` takes its settings from the environment:
`FIT_METRIC_RUN` (which run of each configuration is the benchmark, default 0),
`FIT_METRIC_B` (200), `FIT_METRIC_BLOCK` (8), `FIT_METRIC_ITERATIONS` (4000) and
`FIT_METRIC_CONFIGURATIONS` (all). Under any setting other than the defaults its
files carry the settings in their names, so the default result is left alone.
It is resumable at benchmark granularity: each benchmark's five rows go to a
progress file when it finishes, and a rerun skips what that file holds.

Cost on this machine (AMD Ryzen 7 4800H, 16 threads): the oracle program took
4 minutes for all 5,000,000 benchmarks (28 to 61 s per method, the caches
already on disk); the bootstrap 9.1 s of wall-clock time per benchmark at
$B = 200$ with 16 benchmarks side by side, one per thread, almost all of it in
the t-QVARMA refits, so about 2.5 hours for 1000 benchmarks. Each block-length run (100 benchmarks at $B = 50$) took about 4 minutes.

### Verification

`tests/fit_metric_correctness.c` (`make test-fit_metric_correctness`) checks:

- $\sigma^2$ on a hand-computed case (vectors $(1,2)$, $(3,6)$, $(5,10)$ give
  variances 4 and 16, so $\sigma^2 = 10$), $s$ and $v$ at $d = 2$,
  $\sigma^2 = 4$ ($1$ and $1/2$), and the same $s$ at $d = 200$,
  $\sigma^2 = 40000$, as multiplying every response by 100 gives;
- the oracle shortcut against the variance computed directly over the other
  999 runs, on responses with an offset of 1000 and a spread of 0.01: largest
  relative gap $3.9 \times 10^{-15}$;
- that $d$ reads every run of the best configuration except the benchmark's;
- that the response vectors built for one real series, `cop_0191` run 706,
  reproduce the stored ones: the four local-projection vectors to $3.4 \times
  10^{-8}$ of the largest entry, the t-QVARMA vector from the cached
  parameters exactly;
- that blocks as long as the series give $\sigma^{2*} = 0$ under every method
  (every resample is the series itself; largest value $2.1 \times 10^{-32}$),
  and that the same seed gives the same $\sigma^{2*}$.

On the grid itself, $d$ for the single benchmark of
`docs/MONTECARLO_VALIDATION.md` (`cop_0191` run 706) reproduces the mean losses
published there: 0.01006 under the t-QVARMA, 0.689 linear, 1.336 state 1,
1.383 both states, and 1.408 under state 2, where the best configuration is
`cop_0437`. The benchmarks whose best configuration is the true one, counted
from `fit_metric_oracle.csv.gz`, are the "smallest mean loss" counts of the
grid's results table (53,926 for the t-QVARMA, 63,433 linear, 7,854, 7,305 and
6,746 for state 1, state 2 and both states).

## Results

### Oracle version, every benchmark of the grid

1,000,000 benchmarks per method. $d$ is in the units of the responses and
$\sigma^2$ in their square, which differ between methods: the t-QVARMA's
responses are about a hundredth of the local projections' in size. $s$ and $v$
do not depend on those units, and are given at their 10th, 50th and 90th
percentiles.

| method | median $d$ | median $\sigma^2$ | $s$ 10% | $s$ 50% | $s$ 90% | $v$ 10% | $v$ 50% | $v$ 90% | mean $v$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| qvarma | 0.01288 | 0.002071 | 0.150 | 0.292 | 0.532 | 0.653 | 0.774 | 0.870 | 0.765 |
| linear | 1.522 | 12.45 | 0.315 | 0.446 | 0.611 | 0.621 | 0.692 | 0.760 | 0.691 |
| state 1 | 2.771 | 48.9 | 0.288 | 0.416 | 0.577 | 0.634 | 0.706 | 0.776 | 0.706 |
| state 2 | 2.918 | 53.63 | 0.287 | 0.416 | 0.581 | 0.633 | 0.706 | 0.777 | 0.705 |
| both states | 2.905 | 51.2 | 0.303 | 0.422 | 0.575 | 0.635 | 0.703 | 0.767 | 0.702 |

Split by whether the best configuration is the true one:

| method | best is the true configuration | median $v$ there | best is another | median $v$ there |
| --- | --- | --- | --- | --- |
| qvarma | 53,926 | 0.720 | 946,074 | 0.777 |
| linear | 63,433 | 0.677 | 936,567 | 0.693 |
| state 1 | 7,854 | 0.635 | 992,146 | 0.707 |
| state 2 | 7,305 | 0.627 | 992,695 | 0.707 |
| both states | 6,746 | 0.631 | 993,254 | 0.704 |

Every benchmark here has its true configuration in the field, so these are the
values of $v$ when some configuration is exactly right. They sit around 0.7
under every method, not near 1: the benchmark is itself one draw, so even the
true configuration's runs differ from it by an amount of the order of the
sampling spread. A $v$ on the US data is read against these values rather
than against 1.

$v$ is lower where the best configuration is the true one than where another
configuration takes its place. Within a benchmark that is how the winner is
chosen: a configuration other than the true one is best only when its mean loss
is below the true configuration's, so its $d$ is smaller than the true
configuration's would be. A high $v$ therefore says that the best
configuration is close to the data, not that it is the configuration the data
came from.

### Feasible version against the oracle

Run 0 of every configuration is the benchmark, 1000 benchmarks per method.
$\sigma^{2*}$ comes from $B = 200$ moving-block resamples with blocks of 8
periods, seed 123; $\sigma^2$, $v$ and $d$ are the oracle values of the same
benchmarks. The ratio is $\sigma^{2*}_{j,m} / \sigma^2_{j,m}$, per benchmark,
at its 10th, 50th and 90th percentiles; a bootstrap that measured the sampling
variance correctly would put it near 1. The rank correlations are Spearman's,
across the 1000 configurations.

| method | median $\sigma^2$ | median $\sigma^{2*}$ | ratio 10% | ratio 50% | ratio 90% | rank correlation $\sigma^2$, $\sigma^{2*}$ | median $v$ | median $v^*$ | rank correlation $v$, $v^*$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| qvarma | 0.002072 | 0.0007799 | 0.116 | 0.411 | 1.113 | 0.519 | 0.771 | 0.671 | 0.487 |
| linear | 12.46 | 0.1887 | 0.0025 | 0.0152 | 0.0587 | 0.688 | 0.689 | 0.227 | -0.270 |
| state 1 | 48.96 | 0.7928 | 0.0032 | 0.0162 | 0.0567 | 0.702 | 0.704 | 0.243 | -0.151 |
| state 2 | 53.71 | 0.8368 | 0.0031 | 0.0159 | 0.0540 | 0.720 | 0.707 | 0.242 | -0.121 |
| both states | 51.2 | 0.8081 | 0.0032 | 0.0162 | 0.0554 | 0.714 | 0.704 | 0.240 | -0.232 |

- Under the four local projections the bootstrap variance is about 1.6% of the
  oracle one at the median, and below 6% in 90% of the benchmarks. $s^*$ is
  then about 8 times $s$ ($1/\sqrt{0.016} \approx 7.9$), and $v^*$ is about
  0.24 where $v$ is about 0.70.
- The bootstrap orders the configurations by their variance moderately well
  (rank correlation about 0.7), but it does not order them by $v$: the rank
  correlation between $v$ and $v^*$ is negative under every local projection.
- Under the t-QVARMA the bootstrap variance is 41% of the oracle one at the
  median, with a wide spread (12% to 111%), and $v^*$ (0.671) is closer to $v$
  (0.771).

Every resample gave a finite response under every method (200 usable of 200
in all 1000 benchmarks). The t-QVARMA refits mostly stop at the iteration cap:
26,822 of the 200,000 converged (13.4%), after 3714 iterations on average. The
t-QVARMA $\sigma^{2*}$ therefore mixes the variation across resamples with the
variation from stopping the optimizer early, in proportions this run does not
separate.

### Block length

The first 100 configurations, run 0, $B = 50$ resamples, seed 123, blocks of 8,
25, 50 and 100 periods out of 400; everything else as above. Median of
$\sigma^{2*}_{j,m} / \sigma^2_{j,m}$ over the 100 benchmarks:

| block length | qvarma | linear | state 1 | state 2 | both states |
| --- | --- | --- | --- | --- | --- |
| 8 | 0.390 | 0.0157 | 0.0163 | 0.0148 | 0.0155 |
| 25 | 0.547 | 0.0380 | 0.0537 | 0.0488 | 0.0527 |
| 50 | 0.529 | 0.0818 | 0.133 | 0.138 | 0.136 |
| 100 | 0.416 | 0.216 | 0.397 | 0.436 | 0.419 |

At blocks of 8 the 100-configuration, 50-resample run agrees with the full run
above (0.0157 against 0.0152 linear). Under the local projections the ratio
grows steadily with the block length and is still well below 1 at blocks of
100, a quarter of the series, where only 4 blocks make up a resample. Under the
t-QVARMA it stays between 0.39 and 0.55 with no trend. Why the moving-block
bootstrap of these series understates the local projections' variance by this
much is not established by these runs.

### Next

`docs/MONTECARLO_FIT_METRIC_COMPARISON.md` designs the tests of whether the
oracle $v$ differs across the five methods, one test over all five at once for
each of two nulls (equal means, and no method systematically ahead). Not
implemented yet.
