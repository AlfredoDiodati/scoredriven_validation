# The simulator this project runs, and how it differs from upstream

`model/dsk_sfc` holds the source code of the DSK stock-flow-consistent
agent-based model.

Upstream is `https://github.com/CoMoS-SA/Reissl_2025.git` at commit
`611ff9cb44348baa55be1bc315eefe2c117ccd44`. Five of its files differ here -
`dsk_sfc_main.cpp`, `dsk_sfc_globalvars.h`, `modules/module_finance_sfc.cpp`,
`modules/module_finance_sfc.h` and `CMakeLists.txt` - and
`model/dsk_sfc/upstream/` holds their original versions so the difference can
be built and compared without going back to the network. Two headers are new,
`dsk_sfc_vintage.h` and `dsk_sfc_reductions.h`; nothing upstream includes
them.

Nothing compiles or links against `model/dsk_sfc/upstream/`. `build.sh
--upstream` copies a scratch tree, drops those five files over their modified
counterparts in it, and builds that, so the reference binary is upstream's code
at upstream's flags and nothing in the working tree is disturbed.

    make model              the simulator this project runs
    make model-upstream     bin/dsk_SFC_upstream, unmodified, upstream's flags
    make test-dsk_long_path            the bug fix has its own test
    make test-dsk_build_equivalence    the speed work changes no number

## The program, and how it is invoked

One executable. It reads a JSON file of parameters, initial values and flags,
simulates, and writes an 83-column whitespace-separated results file and an
error log into an `output/` directory beside the path it was invoked as.

    dsk_SFC INPUTS.json -r RUN_NAME -s SEED -f FULL -c CONSOLE -v VERBOSE

`-r` names the run and is appended, with the seed, to every file it writes.
`-s` is a positive integer seed. `-f 1` writes the model's extended output
instead of the single results file, `-c 1` prints error messages to the
console as well as the log, `-v 1` prints progress. All five are optional and
default to `test`, 1, 0, 0, 0. `applications/abm_system_simulate.c` passes
`-f 0 -c 0 -v 0` and one seed per run.

That the output lands beside the *invocation* path rather than the working
directory is why the simulation driver reaches the model through a symlink in
a per-process scratch directory, and it is what made the filename buffer bug
below reachable.

## What was left out of the copy

`model/dsk_sfc` is upstream's source and nothing else. Four things upstream
ships were dropped because they describe someone else's desk rather than this
program:

- `.vscode/`, a VS Code launch and task configuration.
- `CMakeSettings.json`, Visual Studio's WSL settings, carrying the original
  author's Windows user path.
- upstream's own `.gitignore`, which would sit nested inside this repository's
  and hide files from it for build layouts this project does not use.
- upstream's `README.md`, whose substance was compilation instructions for
  Windows, WSL and MinGW, plus the command-line arguments recorded above.

Nothing removed is read by the compiler or the program, which
`make test-dsk_build_equivalence` confirms by still matching upstream byte for
byte after the removals.

One thing to check before this repository goes anywhere public: upstream ships
no licence file, and neither do the three libraries whose source it carries
inside itself (newmat10, rapidjson, CLI11). That is upstream's omission, not something to
paper over here.

## Why any of this

The experiment is a million runs. At the speed upstream ships, that is 73 days
of this machine; at the speed here it is 5.3. Nothing else about the model is
touched, and nothing that touches a number is allowed: every change below
produces output identical to upstream's, byte for byte, which is what
`tests/dsk_build_equivalence.c` checks rather than assumes.

## What was changed

### The build was unoptimised

`CMakeLists.txt` sets `CMAKE_BUILD_TYPE Debug` and adds no optimisation flag of
its own, so a default build carries none at all. It is now Release at `-O2`
with link-time optimisation.

| build | one 600-step run | speedup |
|---|---:|---:|
| as upstream ships it | 37.20 s | 1.00x |
| `-O2 -flto` | 9.70 s | 3.84x |
| the `MACH()` loop order | 7.29 s | 5.10x |
| column reductions in place | 4.71 s | 7.90x |
| `COSTPROD()`'s vintage costs hoisted | 4.13 s | 9.01x |
| raw storage in `MACH()`, cached second-hand minimum | 3.59 s | 10.36x |
| the vintage arrays as one block each | 3.61 s* | 9.8x* |

*The last row was measured against upstream in the same thermal state rather
than against the earlier absolute numbers: an hour of continuous benchmarking
warms this machine and every absolute time drifts about 10% over a session.
Alternating the two builds and comparing, upstream 40.55 s against 4.15 s, is
the measurement to trust; the flattening itself is 8.9% over the row above it,
measured the same way.

Seed 1, the baseline calibration, three repetitions, median reported, nothing
else running. Every row produces the same output as the row above it and as
upstream, byte for byte.

A second round took it further. Its rows are each their own alternating A/B
between the build before the change and the build after it, three or more runs
each, median, because an hour of benchmarking warms this machine and the
absolute times drift several per cent across a session; the ratio inside a row
is the measurement and the level between rows is not.

| change | before | after |
|---|---:|---:|
| `SCRAPPING()` for every firm at once | 2.883 s | 2.604 s |
| the empty firm-vintage pairs skipped | 2.604 s | 2.227 s |
| `COSTPROD()` searching only what the firm holds | 2.270 s | 2.060 s |
| `PRODMACH()`'s ages cleared for every firm at once | 2.067 s | 1.946 s |
| `CANCMACH()` working through `SCRAPPING()`'s list | 1.946 s | 1.918 s |
| `UPDATE()`'s ages, `LOANRATES()`'s column maximum | 1.914 s | 1.846 s |
| `MACH()`'s four copies as whole blocks | 1.853 s | 1.789 s |
| `g_pb` and `C_pb` never cleared | 1.827 s | 1.614 s |
| the second-hand minimum over the live vintages | 1.614 s | 1.557 s |

End to end, upstream and this build alternated three times each in the same
thermal state: 34.840 s against 1.570 s, 22.2x.

Two flags are deliberately absent. `-O3` measured no faster. `-march=native`
lets the compiler contract a multiply and an add into one instruction, which
changed 356 KB of the 3 MB output; equality with upstream is worth more here
than the nothing it bought.

### Two loops in `MACH()` read memory the wrong way round

`MACH()` was 31.6% of the runtime, the largest single share. Its two hot loops
walk the machine vintage arrays with the firm index outermost, and those arrays
have the firm index innermost, so every read landed on its own cache line: one
loop copies four arrays that way, the other accumulates each firm's unit cost
over every vintage that way, and between them they did roughly 48 million
strided reads a period.

The firm loop now runs inside the vintage loops, and the per-vintage
productivities are read once for all firms instead of once per firm. Each firm
still accumulates exactly the same terms in exactly the same order, and the
expressions are unchanged down to where the divisions sit, which is why the
output does not move. That is the 9.70 to 7.29 second step.

### Ranking borrowers rebuilt a column of the matrix every time

The largest single cost in the run was not a model calculation. Two loops rank
each bank's borrowers by debt service, and both ask for the maximum of a column
of a matrix: `DebtServiceToSales2_bank.Column(i).Maximum()` in `LOANRATES()`,
`DS2_rating.Column(i)` twice in `ALLOCATECREDIT()`. newmat evaluates
`Column(i)` into a temporary before reducing it, and a `Matrix` is stored
row-major, so building that temporary copies the column out one row at a time.
It happened once per firm per bank per period: 264 million row copies a run,
about 30% of the time, to find three numbers.

`dsk_sfc_reductions.h` walks the column where it lies instead, by its stride
through the matrix's own storage. The scan order, the comparisons and the
tie-breaking are newmat's own, copied out of `GeneralMatrix::Maximum` and
`Minimum1`, so the value returned and the index reported are the ones the
library returned. 7.29 seconds to 4.71.

### A unit cost recomputed once per machine instead of once per firm

`COSTPROD()` walks a firm's machine vintages from cheapest to dearest until its
desired output is covered. Each pass scanned all 400 vintages and computed the
unit cost each one implies, twice per vintage - once to compare and once to
assign - and the model does several passes per firm. But that cost does not
change while the firm draws its machines down; only the counts do. It is now
worked out once per call and the passes read it. 4.71 seconds to 4.13.

### Two more of the same kind

`MACH()`'s inner loop runs 48 million times a run and reached eight
bounds-checked newmat accessors each time round for the five running totals and
three per-firm inputs; it now takes their storage once per vintage.
`ENTRYEXIT()` rescanned the whole second-hand price matrix, 12,000 elements, to
find the cheapest machine on offer for every vintage of every exiting firm -
71,629 full scans a run - when that minimum only moves when a machine is
actually taken; it is now kept and refreshed at the point it changes. Together,
4.13 seconds to 3.59.

### Ninety-eight per cent of the machine grid is empty

Every firm-vintage pair the model can hold is a cell of the arrays indexed
`[period][supplier][firm]`, and four functions sweep all of them every period.
A firm holds units of very few vintages, though: instrumenting a run at seed 1
counts 52,076,000 firm-vintage pairs visited in `MACH()`'s weighted-average
loop and 51,143,520 of them, 98.2%, with a machine count of zero.

A pair with no machines contributes a productivity multiplied by zero to each
of five running totals, which leaves each total where it was, so those pairs
are skipped. The skip is only taken while the divisors are non-zero, since with
a zero divisor the term would be a NaN and a NaN does not leave a total where
it was; a check over the firms once per period decides which loop runs.
`SCRAPPING()` gets the same treatment, since a pair with no machines enters
neither of its branches. 2.604 s to 2.227 s.

`COSTPROD()` walks a firm's vintages from cheapest to dearest, and it did that
by scanning all 434 of them on every pass and working out the unit cost of each
one. It now collects the vintages the firm actually holds once per call, in the
order the scan visited them, and the passes walk that. Measured on the same
run: 115,116 calls, 723,656 passes, 49,925,260 unit costs worked out a run
against 7,074,653 shortlist entries scanned. 2.270 s to 2.060 s.

### One firm at a time, in arrays where the firm is innermost

Four loops did their work for one firm at a time, and the firm is the innermost
index, so each of them read one element out of each of 434 rows of firms, 48
million times a run. All four now run the firm loop inside the vintage loops,
which reads along memory instead of across it.

- `SCRAPPING()` was called once per firm from inside `INVEST()`'s firm loop. It
  is now called once, before that loop. Nothing `INVEST()` does between firms
  writes what it reads, and `ORD()` sets its own supplier index, so each firm
  still decides on the same information in the same order. 2.883 s to 2.604 s.
- `PRODMACH()` cleared the age of every vintage a firm no longer holds, inside
  its firm loop. That loop is now three loops: the scrapping cancellations, the
  ages, then the new machines. 2.067 s to 1.946 s.
- `UPDATE()` aged every machine a firm holds, inside a loop that otherwise only
  copies per-firm scalars. The ageing is now its own sweep.
- `MACH()` copies `gtemp` into `g` and then into `g_c`, `g_c2` and `g_c3`. All
  four take the same values and a period's suppliers and firms lie together, so
  it is four copies of one block per period rather than four per supplier.
  1.853 s to 1.789 s.

### Two arrays cleared for everything and filled in for almost nothing

`g_pb` and `C_pb` hold what each firm wants to scrap and what it costs to run.
`SCRAPPING()` cleared both for every firm and vintage, then filled in the two
per cent it marked, and `CANCMACH()` read them back by scanning the whole grid
for the marks. That is 833 MB of writes and a full strided scan a run to carry
a few thousand numbers.

`SCRAPPING()` now records what it marks, per firm, in the order it marks it,
and `CANCMACH()` works through that list. Entries outside the list are never
read, so neither array is cleared at all. The two changes were measured
separately: the list is 1.946 s to 1.918 s, dropping the clearing is 1.827 s to
1.614 s.

### Two more reductions that were being redone

`LOANRATES()` pushes a bank's non-customers to the end of its ranking by
setting each one to the maximum of the column plus one, and it rescanned the
column for each of the roughly 180 of them. Each of those writes is itself the
new maximum, so the next is the one before it plus one; the column is now
scanned once per bank.

`ENTRYEXIT()` finds the cheapest machine on the second-hand market by taking
the minimum of a 600 by 20 matrix, of which only the rows for vintages still in
use hold anything but infinity. The scan covers those 434 entries. The column
maximum was measured together with `UPDATE()`'s ageing, 1.914 s to 1.846 s; the
second-hand minimum on its own, 1.614 s to 1.557 s.

### What was tried and did not work

What follows was measured and thrown away, which is worth writing down so
nobody spends the afternoon again. Everything down to the two build-level ideas
comes from the round that ended at 3.59 s, where `SCRAPPING()` and `PRODMACH()`
were the top of the profile at 15% and 14%; the last two are from the round
after it.

*Hoisting `SCRAPPING()`'s invariants.* Its vintage loop reached seven newmat
accessors and did three divisions that depend on the firm's supplier but not on
the vintage. Hoisting all of it out of the loop: 3.680 s against 3.696 s,
alternating five runs each. No difference. The compiler had already done it -
`-flto` inlines the accessors, and once inlined they are plainly loop
invariant. That is also why the profile overstates what is left: the profiling
build has no link-time optimisation, so it attributes real time to calls that
do not exist in the build that ships.

*Gathering `COSTPROD()`'s machine counts.* The firm index is the innermost
dimension of the vintage arrays, so one firm's 400 machine counts sit 400
different cache lines apart, and `COSTPROD()` walks them once per pass.
Gathering them into one contiguous scratch array per call and walking that
instead: 3.711 s against 3.733 s. No difference either, because the function
averages about one pass per call, so there is nothing for the gather to
amortise - it just moves the same scattered reads earlier.

Both results pointed at the array layout, so the layout was taken apart in
three steps, one kept and two thrown away.

*One block per array, same index order.* The nine vintage arrays were a vector
of vectors of vectors, so every row of firms was its own heap allocation and
consecutive rows landed wherever the allocator put them. The walk that matters
reads one element of each row in turn, which was 400 loads from 400 addresses
the processor cannot guess. `dsk_sfc_vintage.h` lays each array out as one
block, keeping the index order, so the same walk reads addresses a fixed
distance apart - a stride the prefetcher recognises. 8.9% faster, and the 185
element accesses in the model did not have to change: the subscripts return
small proxies, so `X[tt-1][i-1][j-1]` still parses and still means the same
thing. Kept.

*Turning the order round so the firm is outermost.* The obvious next step, and
wrong: measured 24.5% slower, 4.72 s against 3.79 s alternating. It does what
was expected for `COSTPROD()`, `SCRAPPING()` and `PRODMACH()`, which walk
vintages for one firm, and ruins `MACH()`, which copies whole rows of firms and
accumulates across firms for one vintage. The two access patterns want opposite
layouts and `MACH()` wins. This also retracts an estimate made earlier in this
file's history, that transposition was worth about 1.5x: it is worth -24%.

*Keeping both layouts at once.* If one order suits `MACH()` and the other suits
the per-firm sweeps, hold `g_c` twice: as it is, and again with the vintages
adjacent, the second built inside the copy loop `MACH()` already runs over
every element. Measured, and it is 10.3% slower on one run and worse again on
throughput, 67.6 runs a minute against 80.1 at sixteen concurrent, with a run
growing from 161 MB to 180 MB. `out/dsk_dual_layout_experiment.txt` has the
whole run; the output was identical on all five seeds, so the copy was being
kept in step correctly. It simply does not pay, and the reason is arithmetic
rather than implementation: transposing 80,000 elements a period costs the same
scattered memory traffic that the per-firm sweep costs, and each array is swept
exactly once a period, so the copy costs what the copy saves.

That closes the layout question. The order the model already uses is the right
one for it, and the only thing worth doing to it was making each array
contiguous, which is done.

Two build-level ideas were tried and gave nothing. Profile-guided optimisation,
trained on a full run at a different seed, measured 4.04 s against 3.88 s -
inside the noise. `-O3` and `-march=native` were tested earlier with the same
result.

Two more from the round that followed, both of which looked like the changes
above and neither of which measured.

*Collecting every firm's held vintages in one sweep.* `COSTPROD()` builds its
shortlist by reading one element out of each of 434 rows, which is the access
pattern every other change here was about removing. Building all 200 firms'
shortlists in one contiguous sweep at the top of `INVEST()` instead: six
alternating pairs, 1.561 s against 1.535 s median, mean difference 0.8%, and
two of the six pairs went the other way. The strided reads cover 694 KB, which
sits in L2, so there was less to save than the pattern suggests, and the sweep
costs a pass of its own. Not kept.

*Writing the two age loops as selects rather than branches.* Eight alternating
pairs, 1.5965 s against 1.5889 s, 0.5%. Not kept.

### Filename buffers too small for a real path

This one is a bug, not a speed change, and it is upstream's.

Every output filename is built with `strcpy` and `strcat` into a fixed 64-byte
buffer, starting from the directory of the path the executable was invoked as.
Twenty-six of those buffers, plus three stack arrays sized exactly one byte
short of the terminator `strcpy` writes. The longest name is the error file's:
the executable's directory, plus `output`, plus `/errors/Errors`, plus the run
name and seed, plus `.txt`. Past roughly 26 characters of directory the writes
run off the end into whatever global follows.

- Upstream's own unoptimised build: silent past the limit, segmentation fault
  at 101 bytes of error-file name.
- The optimised build: segmentation fault at 66 bytes.

`applications/abm_system_simulate.c` invokes the model through a symlink in a
scratch directory, which on a cluster is `$SLURM_TMPDIR/dsk_<pid>` — squarely
in the range that fails. A million runs would have come back as nonzero exit
statuses, and the failure log would have recorded a crash with no cause.

The buffers are now `PATH_MAX` and the three stack arrays one byte longer.
`tests/dsk_long_path.c` runs the model from a path whose error-file name is 137
bytes and requires both that it completes and that its output matches a run
from a short path, since a name that overflowed into a neighbouring global
could have changed a number as easily as crashed.

## What the tests establish

`tests/dsk_build_equivalence.c` runs this project's build and the unmodified one
over the same seeds and compares every byte of the results file. Equality, not
similarity: two runs that agree byte for byte are indistinguishable under any
test that could be applied to them, so no test statistic is needed while it
holds. When it stops holding, the test reports the worst absolute and relative
gap with the period and column it sits at, which is what separates arithmetic
that moved slightly from a model that moved.

It establishes that at one set of flag settings: the seven shock flags at zero,
which is how `dsk_sfc_inputs.json` ships and how every run of this project's
experiment is configured, since `applications/abm_system_simulate.c` writes
only the `params` block. `docs/ABM_SYSTEM_SIMULATION.md` says why they should
stay there.

Two of the changes have a branch that only runs with a shock flag on, and
neither branch has been executed by any test. `MACH()` skips the firm-vintage
pairs with no machines only while every firm's two shock factors are non-zero,
and reverts to adding every term when one is not; `shocks_labprod2` and
`shocks_eneff2` are written only inside `if(flag_prodshocks2==1)` and `==2` in
`modules/module_climate_sfc.cpp`, so with that flag at zero they hold the zero
they are initialised to and the skip is always taken. `PRODMACH()`'s
restructured loop sits beside the `flag_capshocks` block, which never executes
either. Before any run with a shock flag on, add seeds at those settings to
`tests/dsk_build_equivalence.c` and check that byte equality still holds.

## Parallelising inside a single run: tried, measured, not kept

The one place in this model where firms can be split across threads without
any doubt is `MACH()`'s cost loop: each firm accumulates only its own
totals, nothing there draws a random number, and the loop calls nothing, so
the global indices `i`, `j` and `tt` can be replaced by local ones inside the
region without reaching into any other function. Threads take contiguous
blocks of firms, so every firm still receives its terms in the order it
received them serially and the result cannot depend on the thread count.

It was implemented and it works. `tests/dsk_build_equivalence.c` passed and
the output hash was the same at 1, 4 and 16 threads. It is also slower, on the
build as it stood then:

| threads | one 600-step run |
|---:|---:|
| 1 | 6.18 s |
| 4 | 6.39 s |
| 16 | 6.93 s |

The reason is the size of the model rather than anything about the
parallelisation. There are 200 consumption-good firms, so the arrays being
written are 200 doubles, 1,600 bytes, about 25 cache lines. Sixteen threads
writing into that are writing into each other's cache lines, and the loop
holds only about 1.5 ms of work per period to divide in the first place. The
contention costs more than the division saves. At several thousand firms the
balance would go the other way; at this calibration it does not.

So the code went back to the serial loop. What this establishes is not that
the model cannot be parallelised but that its own size is what stops it
paying, which is worth knowing before anyone spends longer on it.

The other candidates are worse bets. `COSTPROD()` and `SCRAPPING()` are called
once per firm from a loop in `INVEST()`, and they read and write the firm index
`j` and a long list of scalar temporaries as globals. Splitting that loop means
making every one of those globals thread-private, and missing one is a data
race - which is exactly the kind of fault the equivalence test cannot be
trusted to catch, because a race is nondeterministic and can pass three seeds
and fail on the fourth. The test is a sound guard for a deterministic change
and not for that one.

## Why it would not have helped anyway

Even a parallelisation that did speed up one run would not shorten this
experiment. Concurrent independent runs on this machine - a Ryzen 7 4800H, 8
cores with two threads each, 7 GB - in runs per minute, ten runs per concurrent
process, seeds 1 upward, nothing else running:

| concurrent runs | this build | upstream's |
|---:|---:|---:|
| 1 | 38.2 | 1.4 |
| 4 | 105.8 | 3.5 |
| 8 | 130.2 | 3.6 |
| 16 | 120.2 | 9.5 |

The machine saturates near 130 runs a minute at eight concurrent, which is a
million runs in 5.3 days. Sixteen is slower than eight, which is what running
two threads on each of eight cores does when both threads want memory. The
experiment is a million independent runs, so the cores are already full and
making one run use several of them cannot raise that number. It would only
shorten a single run's latency, which nothing here needs.

Note what does not carry over. One run is 22.2 times faster than upstream's,
and saturated throughput is 1.7 times what it was before this work - 130
against the 76 measured at the same point in the history. The gap is memory
bandwidth: sixteen runs at once want their vintage arrays in cache at the same
time, and there is one memory system between them. That also means further work
on a single run buys steadily less of what this experiment actually costs.

The upstream column was measured earlier, at the same concurrencies but before
the second round of changes, and it is noisier than this one - its 8-way row is
slower than its 16-way, which cannot be right and means something else was on
the machine during that batch.

Memory is the other limit, and it is closer than it looks. A run peaks at 165
MB resident, both builds, so sixteen at once is 2.6 GB. On a 7 GB machine
that is part of why the 16-way row gains so little over the 8-way, and on a
cluster it is what decides how many jobs fit on a node.
