# The simulator on long horizons

The experiment in `docs/ABM_SYSTEM_SIMULATION.md` runs the DSK model for 600
periods. A second pipeline wants one run of the winning configuration for
400,200 periods. The model cannot do that as it stands: every run longer than
about 1,640 periods is rejected, and every run longer than 2,787 periods never
finishes. This file records why, what each failure is, and which of them are
bugs rather than properties of the model.

`docs/DSK_MODEL_CHANGES.md` is the record of how `model/dsk_sfc` differs from
upstream and how that is tested. The three failures that are bugs or costs are
now changed in `model/dsk_sfc`; each section ends with what was changed and how
it was checked. The two limits that are properties of the model are not.

## How the failures were found

Every run below used the build `make model` produces, the parameter file
`model/dsk_sfc/dsk_sfc_inputs.json` with the nine design parameters of one row
of `dataset/abm_system_design.csv` written into it, the horizon `T` overwritten,
and `-f 0 -c 0 -v 0`. The row is `cop_0191`, the only configuration left in
`out/abm_system_mcs_joint.txt`, unless a table says otherwise. Runs were made
one or four at a time on the machine `docs/DSK_MODEL_CHANGES.md` describes, with
nothing else of note running.

A run counts as failed on the same rule `applications/abm_system_simulate.c`
applies: a non-empty error log, or a results file with fewer rows than `T`.

Horizon, `cop_0191`, seed 1:

| T | wall time | peak resident memory | rows written | first logged error |
|---:|---:|---:|---:|---|
| 600 | 0.51 s | 74 MB | 600 | none |
| 1,200 | 1.07 s | 139 MB | 1,200 | none |
| 2,400 | 2.35 s | 269 MB | 2,400 | period 1,640 |
| 4,800 | stopped by hand at 614 s | 529 MB | 2,787 | period 1,640 |

Seeds 1 to 16 at `T = 4,800`, each stopped after 120 s: all sixteen stop writing
between rows 2,767 and 2,787, and all sixteen log their first error between
periods 1,616 and 1,692.

Other configurations at `T = 4,800`, seed 1, stopped after 120 s: the base
parameter file, `cop_0001`, `cop_0500` and `cop_1000` stop writing between rows
2,785 and 2,816 and log their first error between periods 1,616 and 1,640. The
failures belong to the model, not to one configuration.

The unmodified upstream build, `bin/dsk_SFC_upstream`, run on `cop_0191`, seed
1, `T = 4,800`, stops at row 2,787 like this project's build, and its results
file and error log are byte-identical to this project's up to that point. Both
failures below that happen before period 2,788 are therefore upstream's.

Where a run stops was found with the model's own `-v 1` trace, which prints each
function as it returns, and then with counters and prints added to a scratch
copy of the source outside this repository. None of that instrumentation is in
`model/dsk_sfc`.

## The failures, in the order a long run meets them

| period | what is seen | kind |
|---:|---|---|
| 1,640 | "Carbon content loop did not converge", logged every simulated year from then on | bug, upstream's |
| 1,756 | every climate variable is NaN from then on | bug, upstream's, same cause |
| 2,788 | the run stops advancing and never exits | bug, upstream's |
| 8,889 | the run appears to stop; one loop has hundreds of millions of steps to do | cost that grows with the economy |
| 12,908 | the run stops advancing and never exits | machine counts past what a double counts exactly |
| about 127,700 | world emissions overflow to infinity | property of the model |
| about 142,800 | the price level overflows to infinity | property of the model |

The last two are computed, not observed; their sections say how.

### The carbon cycle's tolerance is below what a double can resolve

`CLIMATEBOX()` in `modules/module_climate_sfc.cpp` splits carbon between the
atmosphere and the upper ocean by iterating on the atmospheric stock until the
remaining imbalance `Cay` is below `1e-10` in absolute terms, for at most
`niterclim = 20` iterations. `dsk_sfc_main.cpp:302` logs the error when the
twentieth iteration is reached.

The stock it iterates on grows without limit (see "World emissions" below). At
period 1,640 the total carbon being split, `Ctot1`, is 1,206,784. The spacing
between adjacent doubles at that size is about `2.3e-10`, so an imbalance below
`1e-10` cannot be represented except as exactly zero, and the iteration settles
on an imbalance of `1.05e-10` and stays there for the remaining eighteen
iterations. The answer it has at that point is as accurate as the arithmetic
allows; the test simply asks for more.

The NaN at period 1,756 comes from the same place. Each iteration estimates the
slope of the imbalance from a second point `Caxx` placed a small step from the
current one, a step proportional to the current imbalance. Once the imbalance is
at the rounding floor that step is below the spacing of doubles too: at
`Ctot1 = 2,261,182`, iteration 3 puts `Caxx` on exactly the same double as
`Cax`, the slope is `0/0`, and every later value is NaN. `Cat(1)` takes the
NaN, the next year starts from it, and the climate block never recovers.

Measured on `cop_0191`, seed 1, by printing every iterate in periods 1,748 to
1,760.

What it affects: the error log, which on the rule above rejects the run, and
the climate columns of the results file (column 47, the temperature, is the
first to go NaN). It does not reach the economy in this experiment. The climate
block feeds the economy only through the shock functions, and all seven shock
flags are 0 in `dsk_sfc_inputs.json` and are never written by
`applications/abm_system_simulate.c`. The five columns the pipeline reads stay
finite in the same run to period 8,888.

**Changed.** The iteration now also stops when an iterate lands on the same
double as the one before it (`module_climate_sfc.cpp`, after the tolerance
test). From that point the upstream loop cannot move either, so where upstream
stays finite the answer is the one upstream reaches after twenty iterations, and
where upstream would divide 0 by 0 the answer is the last finite iterate. The
test cannot fire before the tolerance test does unless neighbouring doubles lie
further apart than the tolerance, which needs a carbon total of about a million,
far past anything a 600-period run reaches.

The same division by zero turns up a second way, much later. At period 6,288 of
`cop_0191`, seed 1, the carbon total is 9.8e16, where neighbouring doubles are 16
apart, and the gradient worked out before the first iteration already uses a
nearby point six units away, which rounds onto the guess itself. So the
iteration is also skipped when that first nearby point equals the guess, and the
guess is kept. Where the first gradient is not 0/0 nothing changes.

Checked on `cop_0191`, seed 1, `T = 2,400`, against the build without the change:
the results files are identical in every cell except column 47 from period
1,756 on, where the unchanged build holds NaN, and the error log is empty where
the unchanged build's has 9,180 bytes of the convergence message.

### Sold second-hand machines are marked with a price of 1,000,000

`ENTRYEXIT()` in `dsk_sfc_main.cpp` sells the machines of exiting
consumption-good firms on the second-hand market, cheapest vintage first. It
fills `C_secondhand` with the unit cost of each vintage on offer, and a vintage
qualifies when its unit cost is no larger than the smallest entry. After a
vintage is sold its entry is overwritten with `1000000` so that it drops out of
the minimum (`dsk_sfc_main.cpp:5286`, upstream's `dsk_sfc_main.cpp:4954`). The
loop repeats until every machine to be sold, `n_mach_exit2`, is gone.

`1000000` stands in for "larger than any unit cost". Unit costs are nominal and
the price level grows without limit, so that stops being true. At period 2,788
the cheapest machine held by an exiting firm costs 1,015,660. The first pass
sells it and writes `1000000` in its place, which is now the minimum; no other
machine costs 1,000,000 or less, none qualifies, and 114,755 machines remain to
be sold by a loop that can no longer sell any. The run does not exit, and
`alarm(T*2)` at the top of the program would end it only after `2T` seconds.

Measured on `cop_0191`, seed 1, by printing `n_mach_exit2` on every pass and
the smallest unit cost among the machines of exiting firms. Replacing `1000000`
with infinity in a scratch copy lets the same run continue past period 2,788 to
period 8,888, where the next section begins.

The rest of the matrix is already set to infinity every period
(`dsk_sfc_main.cpp:2109`), so infinity is also the value the code uses
elsewhere for the same purpose.

**Changed.** A sold vintage is now marked with infinity. While every unit cost
is below 1,000,000, which covers every 600-period run, the comparisons come out
the same as before. Checked on `cop_0191`, seed 1: `T = 3,200` finishes with all
3,200 rows and an empty error log.

### Rationed labour cancels machine orders one machine at a time

`LABOR()` in `modules/module_macro_sfc.cpp` scales production back when labour
demand exceeds supply. For a capital-good firm that means cancelling
`reduction` machines of its customers' orders, and `module_macro_sfc.cpp:46`
does it one machine per pass: draw one of the 200 consumption-good firms at
random, and cancel one machine if that firm is a customer with an order
outstanding. With about one customer in twenty that is roughly twenty random
draws per machine.

The number of machines scales with real output, which grows without limit. At
period 8,889 of `cop_0191`, seed 1, labour is rationed and capital-good firm 1
has to cancel 132,512,000 machines, which completes, and firm 3 873,124,000,
which did not complete in the roughly 300 seconds the run was left for. Nothing
is wrong with the result; the cost is the problem, and it grows with the
economy.

This file is identical upstream (`model/dsk_sfc/upstream/` holds no copy of it),
but upstream cannot reach this period, because the previous section stops it
first.

`ENTRYEXIT()` has three loops of the same shape, and none of them scales with
the economy. Two match a firm to `step` new customers, a handful. The third
hands entrants the second-hand machines left over after each has been given its
rounded-down share, which is fewer machines than there are entrants.

How the count grows, from a scratch build that printed every cancellation count
on `cop_0191`, seed 1: the largest in periods 1 to 1,000 was 94, in 2,001 to
3,000 1,093, in 6,001 to 7,000 77 million, in 8,001 to 9,000 576 million. At
600 periods, over 192 runs - the base parameter file, the lowest and the highest
value of every design parameter, `cop_0191` and 60 design rows drawn at random
with seed 7, each at seeds 1, 2 and 3 - the largest was 86.

**Changed.** Past 10,000 machines, `LABOR()` hands the cancellation to
`dsk_sfc_bulk_cancellation.h`, which draws how many machines each customer loses
as binomial counts from the model's own random stream, with the distribution the
loop has. At or below 10,000 the loop runs exactly as upstream wrote it, which is
more than a hundred times the largest count any 600-period run above produced.
Past the threshold a run is no longer the run the loop would have produced, only
one with the same distribution, because the random draws differ.

`tests/dsk_bulk_cancellation_distribution.cpp` checks that distribution against
exact probabilities worked out by enumerating the loop's own process, and runs
the loop itself and a deliberately wrong sampler through the same comparison.
The first version of the bulk draw failed it: `bnldev`, the model's binomial
generator, returns a Poisson draw rather than a binomial one when there are at
least 25 trials and the mean is below one, and at 30 trials and probability 1/40
that puts 12% too much weight on a count of four. The bulk draw inverts the
binomial distribution directly in that range and leaves `bnldev` as upstream
wrote it.

### Machine counts pass the largest whole number a double holds exactly

Machines are counted in doubles, and every whole number is exact in a double only
up to 2^53, about 9.007e15. The number of machines grows with real output. At
period 12,908 of `cop_0191`, seed 1, the exiting consumption-good firms hold
11,917,600,000,000,000 machines between them, above that limit, and the counts
added into the second-hand pool and the counts handed out of it no longer come to
the same total. Eleven firms enter; the last of them, firm 194, is owed one
machine when the pool is empty, and the loop in `ENTRYEXIT()` that hands out
second-hand machines draws a capital-good firm at random, looking for stock, for
ever.

Measured with a scratch build printing the machines to sell, the pool left and
the machines owed, on every pass of that loop whose number is a power of two.

This is not a bug in one loop. From here on every count in the model is rounded
to the nearest representable double, and any logic that needs two totals of
machines to agree exactly can fail the same way. It is the first of the limits
that come from the model's levels growing without bound, and it arrives much
earlier than the overflows below.

**Changed.** Nothing in the model depends on how big one machine is, so
`model/dsk_sfc/dsk_sfc_machine_lots.h` makes machines bigger: once a firm holds
more than 2^40 of them it doubles the output a machine makes and halves every
count, which buys back the headroom the counts had been growing into. It can be
done as often as needed, and on `cop_0191`, seed 1, it fires six times between
period 9,374 and period 11,000.

What moves with the lot: counts halve; the output one machine makes and the
payback threshold, which is in units of output, double; money for one machine -
the price of a machine, its production cost, the capital goods price index and
the price each machine was bought at - doubles; machines per worker and per unit
of energy, in use and in the innovation and imitation draws, halve. Capacity is
a count times the output a machine makes, so it and everything derived from it
are unchanged.

Two things make it work. It happens at one point in the period, straight after
`MACH()`, because that is where the machine state is self-contained: the working
copy and the holdings agree entry by entry, last period's investment has already
been added to capacity, and nothing is yet marked for scrapping. Elsewhere in
the period the capacity carried between periods and the holdings would be left
disagreeing by a machine, which shows up as a firm with negative investment in
the model's own error log. And the stocks the model keeps in step with the
holdings - a firm's machine count, its capacity, and the value of its capital at
the prices it paid - are re-derived from the rounded holdings rather than left
to drift.

Unlike redenominating money, this cannot be exact: halving an odd count leaves
half a machine. At the ceiling a holding is around 1e12 machines, so the rounding
moves it by about one part in 1e12. `tests/dsk_machine_lot_rebase.c` measures
what that does: in the period of the first rebase, real GDP, consumption,
investment, employment and the capital stock are within 6e-12 of a run that
never rebased, and every period before it is identical. After it the two runs
are different paths, which is what the model does with any perturbation - its
own `dsk_ulp_sensitivity` test shows a change of one part in 1e15 reaching the
output within a few periods - so what the test compares from there on is the
process: mean growth of real GDP and of the price level over the following 1,600
periods agree within half a standard error.

With this, a run reaches period 14,000, past the 12,908 where it used to stop.

### World emissions grow by a fixed factor every year, forever

`CLIMATEBOX()` multiplies the rest of the world's emissions, `Emiss_global`, by
`g_emiss_global = 1.02185` once per simulated year (every `freqclim = 4`
periods) and nothing else ever changes it (`module_climate_sfc.cpp:61`). It was
5.7963e12 at period 1,756. Multiplying by 1.02185 overflows the largest double,
about 1.8e308, after 31,479 more years, which is period 127,672.

This is not measured: it is the arithmetic of a constant factor, from the value
printed at period 1,756. What the climate block does after the overflow has not
been run.

The same growth is why the carbon stocks in the first section reach the
millions: atmospheric carbon at period 1,640 is 1,205,648 against a
pre-industrial 590 in the model's own units, and the temperature anomaly is
28 degrees.

### Redenominating the money side

Nothing in the model depends on the unit money is counted in, so the answer to
prices and money stocks that grow without limit is to change the unit:
`model/dsk_sfc/dsk_sfc_redenomination.h` divides every money quantity by a power
of two once the wage passes a ceiling, the way a currency reform drops zeros. A
power of two is exact in binary floating point and exactly reversible, so only
the exponents change. The ceiling is 2^512 and the step 2^256 unless the
parameter file names others, and a 600-period run never comes near either, so
the equivalence tests against upstream still compare identical files.

What counts as money is a list of 214 variables in that header, and it is the
whole difficulty. `tests/dsk_redenomination_invariance.c` is what makes the list
honest: it runs the model with the ceiling low enough to fire several times and
requires every value in the results file to come back either unchanged or scaled
by exactly the factor applied, which fails if a money variable is left out or a
real one is included. Getting it to pass took six rounds, and what it caught is
worth recording, because none of it was visible by reading the declarations:

- Four of my classifications were wrong. Consumption demand `D2`, its
  temporaries and `l2` are quantities, not money: `Qd = De + Ne - N` adds demand
  to inventories. `RDin` and `RDim` are R&D staff, not spending: they are
  `Ld1rd`, the money `RD` divided by the wage. `dw` and `dw2` are wage growth
  rates. Each wrong entry broke the model as thoroughly as a missing one.
- Three absolute amounts of money were written as bare numbers in the code, and
  each one silently changes what the model does once the unit changes. Sales
  carry `+ tolerance` to keep a ratio finite (`module_finance_sfc.cpp`), the run
  aborts if the price index falls below `0.01`, and consumption is shared out
  `while (Cres >= 1)`. They are now `sales_tolerance`, `cpi_floor` and
  `consumption_residual_floor`, set to the same values, and redenominated with
  everything else.
- One coefficient is per unit of money: the energy sector's chance of a
  successful innovation is `1-exp(-o1_en*RD_en)` with `RD_en` a share of
  revenue, so `o1_en` moves the other way. Its two manufacturing counterparts
  multiply R&D staff instead and are left alone.

This does not touch the machine counts, which pass what a double counts exactly
around period 11,500, nor the emissions.

### The price level grows without limit

Nominal quantities in the model are not rescaled. In the `cop_0191`, seed 1 run
the price level (results column 34) is 25.08 at period 600, 2.754e10 at period
4,800 and 1.816e19 at period 8,888; between the last two it grew at 0.497% per
period. Continued at that rate it overflows at period 142,850. GDP (column 2)
grew at 0.259% per period over the same stretch and would overflow at period
269,829.

These are extrapolations of a rate measured over 4,088 periods and should be
read as an order of magnitude. Long before either overflow, the gap between
neighbouring doubles at these magnitudes becomes larger than some of the flows
the model adds to them; where that first changes a result has not been
measured.

This limit and the machine counts above are not bugs. The model is written in
levels that grow, and no local repair removes them.

## What the long-horizon work costs

None of it is meant to make the model slower, and the experiment is billed in
runs a minute, so both were measured.

Setup: three builds - the commit before this work (`0a8d6e74`), the commit with
the first three fixes (`dcf3b4a6`), and the working tree with the redenomination
and the machine lots on top - each built by `model/dsk_sfc/build.sh` from its own
copy of the tree, all reading the same `dsk_sfc_inputs.json` (md5 checked equal),
600 periods. One batch ran them in that order and a second in the reverse order,
because this machine drifts several per cent over a session. Single runs: seed 1,
eight and twelve rounds, the first two discarded, median reported. Throughput:
eighty runs, eight processes of ten seeds each, four rounds, median reported.
Nothing else was running.

| | before | first three fixes | redenomination and lots |
|---|---:|---:|---:|
| one run, forward batch | 547 ms | 549 ms | 563 ms |
| one run, reversed batch | 521 ms | 529 ms | 528 ms |
| throughput, forward batch | 394 a minute | 378 | 378 |
| throughput, reversed batch | 390 a minute | 379 | 377 |

Read it this way. The single-run differences change sign between the two batches
and are noise at this resolution. The throughput difference does not: the first
three fixes cost about 3% of saturated throughput in both orders and in every
round, whatever position the build ran in, and the redenomination and the lots
add about 1% more, which is at the edge of what these four rounds resolve. Over
the million-run experiment 4% is about an hour and a half.

Where that 3% sits has not been traced. The fixes add a comparison per carbon
iteration, a comparison per rationing event and one changed constant, none of
which accounts for it on its own; a shift in code layout under link-time
optimisation would. It was measured rather than argued, and it buys runs that
do not stop at period 2,788.

## Two costs that decide whether a long run is practical

### Memory

Four arrays are sized by the full horizon: `g`, `gtemp` and `g_price` hold one
double per period, capital-good firm and consumption-good firm, and `acquired`
one int, allocated in `dsk_sfc_main.cpp` at the start of the run. With
`N1 = 20` and `N2 = 200` that is 112 KB per period of horizon, 44.8 GB at
`T = 400,200`. The measured peak resident memory grows by 0.108 MB per period of
horizon between `T = 600` and `T = 4,800`, which agrees. The machine has 7.5 GB.

### Time

An earlier version of this section said the vintage loops cover every period
since the start. That was wrong. `OVERBOOST()` advances `t0`, the oldest vintage
they cover, past every vintage nobody holds, at the end of every period.

With all three changes above, `cop_0191`, seed 1, a scratch build printing
`t - t0 + 1` and the elapsed time every 250 periods: 21 vintages held at every
sample from period 250 to period 12,750, and 27.8 s from period 250 to period 12,500, about
2.2 ms a period, rising slowly (1.7 s per thousand periods near the start, 2.9 s
near period 12,000). The 41 ms a period measured before the labour change was
that loop. Periods between samples were not measured, so the largest window is
not known.
