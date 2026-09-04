# Simulating the ABM over a Latin hypercube design

## What the experiment produces

A Monte Carlo experiment with four indices: variable, time, replication and
parameter configuration.

| index | size | what it is |
|---|---:|---|
| variable | 5 | the series the auxiliary model is fitted on |
| time | 400 | periods 201 to 600 of a 600-period simulation |
| replication | 1000 | independent seeds of the model's own randomness |
| configuration | 1000 | points of the structural parameter space |

One replication is one $400 \times 5$ table. One configuration - a CoP,
Configuration of Parameters - is 1000 of them, and the design is 1000
configurations, so the experiment is a conceptual
$5 \times 400 \times 1000 \times 1000$ object that exists nowhere as one
thing: it is a million small files, and the two outer indices are what the
work is split along.

## What is being sampled

The simulated economy is the DSK stock-flow-consistent agent-based model. It
has no likelihood that can be written down and no reduced form that can be
solved for; the only thing it offers is the ability to draw one path given a
parameter vector and a seed. Everything downstream in this project is built on
that restriction: an auxiliary model summarises simulated and real paths alike,
and configurations are compared in the auxiliary model's space rather than in
the structural model's own.

Nine structural parameters are varied. Every other element of the model's
calibration is held at its baseline value, so the experiment identifies nothing
about them; they are part of the maintained specification, not of the design.

| parameter | lower | upper |
|---|---:|---:|
| Gamma | 0.05 | 0.25 |
| chi | -1.50 | -1.25 |
| psi1 | 0.10 | 0.50 |
| psi3 | 0.10 | 0.50 |
| alfa | 0.10 | 0.50 |
| taylor1 | 1.00 | 1.50 |
| taylor2 | 0.00 | 0.50 |
| taylor | 0.50 | 0.95 |
| kappa | 0.50 | 0.95 |

The names are the model's own and the interpretation of each parameter belongs
to the model, not to this design; the supports in the table are the whole of
what the design asserts. Two readings are worth stating as inferences from the supports
alone, since they are what makes the bounds look deliberate rather than
arbitrary: `taylor1` starting exactly at 1 is where the Taylor principle binds,
and `taylor` confined to $(0,1)$ is the shape of an interest-rate smoothing
coefficient. Neither is verified against the model's source here.

## The design

The design is a file, `dataset/abm_system_design.csv`, drawn by
`applications/abm_system_design.c` and read by everything else. Splitting the
two is what lets a thousand simulation processes share one design without any
of them being able to redraw it. The sampler is et_al's `lhs_random`; the
design is a deterministic function of the seed, the bounds and that sampler,
so it is reproducible rather than something to carry between machines, and
regenerating it is refused unless forced, since the row index is the identity
every stored replication and every result downstream is named by.

The nine bounds define a box $\prod_{k=1}^{9}[a_k, b_k]$, and the design is a
sample of $n = 1000$ points from the uniform distribution on that box, drawn by
Latin hypercube sampling rather than independently.

Latin hypercube sampling partitions each margin into $n$ intervals of equal
probability, places exactly one point in each interval, and pairs the margins
by an independent random permutation. The construction is done on the unit
hypercube and mapped to the box by the componentwise inverse of the uniform
distribution function, $\theta_k = a_k + (b_k - a_k)u_k$, which is affine. The
map changes the scale of each margin and nothing else: no correlation is
imposed between parameters, and every point of the box has the same density as
every other.

Three properties are the reason for using it here.

- Every margin is exactly stratified. An independent uniform sample of size
  $n$ leaves some intervals of a margin empty and puts several points in
  others; the Latin hypercube leaves none empty. With 1000 points the design
  visits 1000 distinct values of each of the nine parameters.
- For a response that is additive across parameters, the sample mean over a
  Latin hypercube has strictly smaller variance than over an independent sample
  of the same size (McKay, Beckman and Conover 1979), and asymptotically the
  additive part of the response's variance decomposition is removed entirely
  (Stein 1987). The gain is in the additive component only: pure interaction is
  estimated no better than by independent sampling.
- A tensor grid is not an option at this dimension. Three levels per parameter
  is $3^9 = 19683$ points, six times the budget here, and it still visits only
  three distinct values of each parameter, so nothing about a response's shape
  along a margin can be resolved.

The sampler used is the plain randomised one: within-stratum positions are
uniform and the pairing across margins is a uniform permutation, with no
maximin-distance or orthogonality post-optimisation. It is not reproducible
against R's `lhs` package, whose API et_al's follows: both draw from the same
distribution, but the generators differ, so the same seed does not give the
same points. The design is therefore
space-filling in the marginal sense and only in that sense; a sample
correlation between two parameters is not zero by construction, it is merely
$O(n^{-1/2})$. Rows are checked for exact duplication, which a continuous
sampler should never produce and which would indicate the sampler was
misconfigured.

The design is drawn under a fixed seed, so the 1000 configurations are the same
object on every machine and the design table can be rebuilt rather than
transported. It is written once and then treated as read-only: regenerating it
after any simulation has run would silently renumber the configurations the
arrays refer to.

## The replication dimension

For a fixed configuration the model is a distribution over paths, not a path.
One seed draws one realisation of that distribution; the 1000 seeds within a
configuration are treated as an independent sample from it, which is a property
of the model's random number streams and is assumed here rather than tested.

The replications are not there to be averaged away at simulation time. The
validation step this dataset feeds computes a loss per replication and needs
the individual losses, because the Model Confidence Set estimates a sampling
variance from repeated observations of the same configuration and a
pre-averaged scalar supplies none. Any configuration-level quantity computed
from $R$ replications carries a sampling error of order $R^{-1/2}$, and $R$ is
a free parameter of the experiment in a way the model's own dynamics are not.

## Common random numbers

The same seeds, $1$ through $R$, are used for every configuration. This is the
common random numbers variance reduction: what the experiment is ultimately
about is differences across configurations, and for two random quantities

$$\operatorname{Var}(X - Y) = \operatorname{Var}(X) + \operatorname{Var}(Y) - 2\operatorname{Cov}(X, Y),$$

so inducing positive correlation between configurations shrinks the variance of
their comparison without touching either marginal distribution.

The size of that gain is not guaranteed. Common random numbers work through
synchronisation: the same seed must drive the same shocks in both runs. In an
agent-based model whose parameters change how many draws are consumed per
period, the two streams desynchronise after the first divergence, and the
induced correlation can be small. Nothing is biased when that happens - within
a configuration the replications remain a valid sample from that
configuration's path distribution - only the variance reduction is lost.

The seed a replication came from is its own index plus one, and the archive
records that index beside the series, so a stored replication can always be
traced back to the run that produced it.

The dependence across configurations that this creates is compatible with the
resampling scheme downstream: the Model Confidence Set bootstrap resamples the
replication index jointly for all configurations, so whatever cross-
configuration correlation the shared seeds induce is carried into the bootstrap
draws rather than destroyed by them.

## The observables

The model writes 83 whitespace-separated columns per period with no header.
Eight are kept:

| name stored | model expression | column |
|---|---|---:|
| GDP | `GDP_r(1)` | 2 |
| Employment rate | `1-U(1)` | 5 |
| Gross inflation | `cpi(1)` | 34 |
| Interest rate | `r` | 50 |
| Energy | `D_en_TOT` | 8 |

The names are the ones `applications/abm_system.h` resolves its variables by,
not descriptions chosen here. One of them is a misnomer inherited from the
existing files rather than a claim: the series stored as `Gross inflation` is
the price level $cpi_t$, and the header reads it as a level, differencing its
logarithm to get inflation.

The model also writes consumption, investment and emissions. They are not
stored. The auxiliary model has five equations and the real US data has no
counterpart series for any of the three, so nothing in this project can compare
them against anything; storing them costs 37.5% of the experiment's disk for a
comparison that would need real data this project does not have. This is
irreversible in the only sense that matters: bringing them back means running
the million simulations again.

Two of these positions are easy to get wrong and neither mistake announces
itself in the data.

- The price series is the level. Column 6 is not a level: it is
  $cpi_t / cpi_{t-4}$, a four-period gross inflation factor. Both are positive,
  smooth and plausible; a series that is a ratio where a level is expected
  survives every downstream transformation and produces wrong dynamics quietly.
- Column 50 is the central bank's policy rate. The average interest rate on
  commercial bank loans is column 82. They co-move and differ by a spread, so
  again nothing breaks, it only stops answering the question asked.

The five columns are not stored in the model's own units. Each run is
transformed as soon as it is read, by the same code the older dataset goes
through, into what the auxiliary model is actually fitted on: GDP and energy
as $100$ times the first difference of their logarithms, the price level the
same way, giving inflation, the employment rate as $100$ times its first
difference, and the interest rate as $100$ times its level, which turns the
model's decimal into the percentage points the real US data is measured in.

That transformation lives in one place, `applications/abm_system.h`, and both
this simulator and the reader for the older `.Rdata` files call it. Storing
the levels instead would keep the transformation revisable, at the price of
1.5 times the disk and a second pass over a million files; the transformation
has been fixed since the older dataset was built, and this makes the
simulator's output directly fittable with nothing in between.

The transient goes the same way. Each configuration is simulated for 600
periods and the first 200 are discarded as convergence to the model's own
attractor rather than the dynamics of interest; period 200 is the anchor for
period 201's first difference, so what is written is periods 201 to 600, 400
of them. The cut is 200 because that is where the auxiliary model's sample
starts.

The burn-in can still be lengthened after the fact by dropping rows from what
is stored, and never shortened. Asking whether 200 periods was enough needs a
re-simulation.

## Failed replications

A replication is accepted into the array only if all four of the following
hold: the executable exits with status zero, its error log is empty, the output
contains exactly $600 \times 83$ numbers, and every extracted value is finite.
A replication that fails any of them is stored as missing, and its seed and the
reason are written to a failure table beside the array.

This is a selection mechanism and not measurement noise. The probability that a
run fails depends on where in the parameter space it sits - explosive or
degenerate configurations are exactly the ones that produce non-finite output
or a non-empty error log - so the completed replications of a configuration are
a sample from its path distribution *conditional on completion*, and the number
of completions is itself information about the configuration. Two consequences
follow. Any statistic computed over completed runs only is conditional on
completion and should be reported with the completion count next to it. A
configuration that lost a substantial share of its replications is not a
configuration with a smaller sample; it is a different object, and the failure
table is what makes the difference visible instead of silently absorbed into an
average.

## What is written, and when

Simulation is the whole pipeline's only step that is not C, and it is not one
either: `applications/abm_system_simulate.c` reads the design, runs the model,
transforms each run and writes it. Nothing between the model executable and
the auxiliary model's fit is anything but C, the standard library and et_al.

Replications are written ten at a time into one `.npz` archive:

    dataset/abm_system/cop_0001/batch_000.npz

holding six members: the five series, `GDP_growth`, `EN_growth`,
`Employment_change`, `Inflation` and `InterestRate`, stacked one replication
after another down the rows, and `replicate`, the index each row belongs to.
That is exactly what `applications/abm_system_extract.c` produces from the
older `.Rdata` dataset and exactly what `applications/abm_system_fit_qvarma.c`
reads, so the fit does not know or care which of the two wrote a given file,
and there is no extraction pass after this one. Both writers and every reader
go through the same four functions in `applications/abm_system.h`, so the
layout is written down once.

Ten to an archive rather than one is what makes the compression work. A
deflate stream that sees ten replications of a series finds far more to reuse
than ten streams that each see one: 10,619 bytes per replication against
12,297, measured on the same real file. Ten is where that curve has
essentially flattened, so a larger batch would buy almost nothing more and
would put more work at risk.

The `replicate` member is what lets an archive be short. When runs inside a
batch fail, the archive holds only the ones that succeeded and still says
which they are, so nothing downstream has to infer identity from position or
from a file name.

The archive is a zip of NumPy arrays, so it carries the column names with it
and the fit still asks for its series by name rather than by position. It is
also `numpy.load`-able, which is what makes inspection outside this pipeline
possible at all.

When each archive is written is as much of the design as what is in it. An
archive is written the moment its tenth run finishes, and the model's own
83-column output is deleted as soon as it has been read, so at most ten
transformed series are ever held in memory and nothing is written at the end
of a configuration. An interruption costs the batch being filled, at most ten
runs, which is the price of the compression above. The failure log and the
per-configuration manifest line are appended and closed on every write, so
what is on disk is what has actually happened rather than what a buffer would
eventually hold.

Resuming is the same command again. A batch whose archive already exists is
skipped without running the model, so an interrupted run, a re-submitted job
array and a deliberate extension to more configurations are all the same
operation. A batch is skipped whole, its failures included: the model is
deterministic in its seed, so a run that failed once fails the same way again,
and the failure log is where that is recorded.

The model is `model/dsk_sfc`, vendored into this repository and built by
`make model`; `docs/DSK_MODEL_CHANGES.md` records what this project changed in
it, what those changes are tested to leave alone, and the bug they fixed. It
reads one JSON parameter file and one seed per execution and writes its output
into the directory holding the executable it was invoked as. Each
process therefore reaches the executable through a symlink inside its own
scratch directory: the baseline JSON is parsed, the nine design parameters and
the horizon are overwritten, and the result is written next to that symlink,
which is what lets many processes share one build without colliding. The
parameter file is rewritten in full rather than patched textually, and the
101 fields of the baseline that the design does not touch come through
unchanged.

## Cost, and what the cost estimate rests on

One replication is 10,619 bytes, the mean over the 108 replications of
`EstimationSeriesSample1_1` re-extracted through the current writer, so it is
measured on real dynamics rather than on synthetic test data. A configuration
is 10.6 MB and the design 10.6 GB. The same replications were 29,356 bytes
each as CSV, 29.4 GB, and 12,297 bytes each in one archive per replication,
12.3 GB.

What is left below this is not a format question. Storing the model's levels
rather than the transformed series would compress better, because levels are
smooth and growth rates are noise, but the fit would then have to transform
what it reads, which is not what the fit is for. The full experiment is
$1000 \times 1000 = 10^6$ model executions.

The only timing available is a single one: one 600-step run of the optimised
build took about 18 seconds on the machine where this design was prepared. It
is one run, not repeated, not on the hardware the experiment would use, and it
is quoted here because it is the entire basis of the arithmetic that follows:
$10^{6} \times 18\ \text{s} = 1.8 \times 10^{7}\ \text{s}$, about 5000
core-hours. That number should be treated as an order of magnitude and
re-measured on the target machine before anything is committed to a queue.

The experiment is embarrassingly parallel along the configuration index and
along nothing else that is worth exploiting: one task per configuration, its
replications run sequentially inside it. Which configurations a single
invocation simulates is an argument, so the whole design is one command on one
machine and one configuration per job on a cluster, with no separate driver
between them. An array already on disk is skipped rather than recomputed, which
is what makes an interrupted run resumable and a job array re-submittable. Each task copies the executable to
node-local scratch, because the model writes its output into the directory
holding the executable - concurrent tasks sharing one copy would overwrite each
other's results and put a million small writes on a shared filesystem.

Before the full design, run a handful of configurations at the corners of the
box with a few replications each. What that pilot checks is not correctness of
the pipeline but the failure rate: the completion count is the one quantity
whose behaviour across the parameter space cannot be predicted from the design
and determines whether the experiment as specified is worth its 5000 hours.
