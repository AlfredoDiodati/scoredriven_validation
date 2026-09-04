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

    make test-dsk_build_equivalence         the aggregate output matches upstream
    make test-dsk_full_output_equivalence   so do the per-firm files under -f 1
    make test-dsk_design_equivalence        so does the model away from the baseline
    make test-dsk_long_path                 the bug fix has its own test

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
of this machine; at the speed here it is 2.2. Nothing else about the model is
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

A third round, measured the same way. It was directed by the loop-level cycle
counts in the section after next rather than by the per-function profile, which
by then had stopped agreeing with the A/B results.

| change | before | after |
|---|---:|---:|
| `COSTPROD()`'s scan walks a pointer | 1.722 s | 1.670 s |
| `PRODMACH()`'s clearing of unheld ages dropped | 1.655 s | 1.542 s |
| `g_c2` and `g_c3` copied only when read | 1.361 s | 1.257 s |
| `UPDATE()` ages every entry, testing none | 1.259 s | 1.174 s |
| `SCRAPPING()` reads the holders `MACH()` recorded | 1.174 s | 1.109 s |
| an entrant's arrays cleared three ways, not six | 1.100 s | 1.043 s |
| the consumption cap's row sums taken in place, once | 1.056 s | 1.018 s |
| `COMPET2()`'s count of surviving firms hoisted | 1.035 s | 1.006 s |
| `MACH()`'s weights test four firms at a time | 1.029 s | 0.932 s |
| `COSTPROD()` reads the vintages `MACH()` recorded | 0.895 s | 0.830 s |
| the borrower ranking sorts packed pairs | 0.846 s | 0.763 s |
| `ALLOCATECREDIT()` reads the ranking `LOANRATES()` made | 0.778 s | 0.740 s |
| `g_c` replaced by a per-firm scratch | 0.717 s | 0.663 s |
| `g_pb` and `C_pb` too, `g_c2` and `g_c3` not allocated | 0.664 s | 0.642 s |

End to end, upstream and this build alternated three times each in the same
thermal state: 35.28 s against 0.677 s, 52.1x.

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

### Measuring the loops instead of the functions

A third round started by replacing the per-function timing above with cycle
counts taken with `__rdtsc` around individual loops, because the per-function
figures had stopped agreeing with what the A/B measurements showed. Bracketing a
function that is called 115,000 times a run costs more than the thing being
measured, and the timer calls also stop the compiler moving code across them.
The loop counts, one run at seed 1, as a share of the whole run:

| loop | share of the run |
|---|---:|
| `SCRAPPING()`'s scan for the firms holding a vintage | 10.5% |
| `UPDATE()`'s ageing | 9.7% |
| `COSTPROD()`'s scan for the vintages a firm holds | 6.9% |
| the rest of `COSTPROD()` | 5.5% |
| `MACH()`'s block copies | 3.8% |
| `CANCMACH()` | 0.4% |

Two of the three largest were sweeps of all 52 million firm-vintage pairs that
act on the 2% which are not empty, and one turned out not to be needed at all.

### Work done for entries nobody holds

`PRODMACH()` cleared the age of every machine a firm no longer holds, over the
whole grid, every period. Nothing reads those ages. Every read of an age -
`SCRAPPING()`, `CANCMACH()`, `ENTRYEXIT()` - asks for it only where the firm's
count of that machine is positive, and the three places a count rises from zero
all set the age at the same entry: the purchase in `PRODMACH()` and the two
second-hand transfers in `ENTRYEXIT()`. The sweep is gone. 6.6%.

`UPDATE()` ages the machines a firm holds, and tested every count to find them.
By the same argument the test is unnecessary: ageing an entry nobody holds
changes no value that is ever read, and the entry's age is set when it is next
acquired. It now ages every entry in the window, which reads none of the counts
and is a plain increment over contiguous integers. 6.8%.

`MACH()` copies the machine counts into three working arrays, `g_c`, `g_c2` and
`g_c3`. The second and third are read in one function, `ADJUSTEMISSENLAB()`,
which is reached only from the two `flag_capshocks` branches in `PRODMACH()`.
With that flag off nothing reads them, so they are copied only when something
will. 8.5% - the largest single change of the round, for a two-line condition.

### More work whose result nobody reads

`ENTRYEXIT()` clears six arrays over every vintage for each firm that enters,
and three of them do not need it: `g_c2` and `g_c3` are read only under
`flag_capshocks`, and an age is read only where a count is positive, which the
cleared counts make none of them. 3.7% of the run to about half that, and the
run from 1.100 s to 1.043 s.

### Two more reductions taken once instead of many times

`ALLOC()` caps consumption at what households can pay by scaling every firm's
sales down, and asked for the row total twice for each of the 200 firms, each
time through `S2.Row(1).Sum()`, which copies the row before summing it. The row
is already contiguous, and the total is the same in the test and in the
subtraction because nothing changes between them. `dsk_sfc_reductions.h` sums it
where it lies, once per firm. 1.056 s to 1.018 s.

`COMPET2()` divides by the number of surviving firms twice for each firm, and
counted them by summing a 200-element vector each time, inside a loop that
cannot change which firms are exiting. 1.035 s to 1.006 s.

### Finding the same set twice

`SCRAPPING()` reads every firm's count of every vintage to find the few that are
not zero. `MACH()`'s weighted-average loop reads exactly the same counts a few
steps earlier, and nothing writes them in between, so it records which firms
hold each vintage as it goes and `SCRAPPING()` reads that list. 5.5%.

Recording them inside that loop rather than in a pass of its own is not an
incidental choice: the separate pass was measured and is 6.8% *slower* than
recording in place, because it reads all 52 million counts a second time and
that costs more than a store in the middle of the existing loop.

### Arrays that existed to be destroyed

The round above ran into a wall that was not in the code: one run kept getting
faster and the throughput of eight concurrent runs stopped moving. Two changes
worth 9.8% and 4.9% each left saturated throughput where it was. Eight runs at
once were thrashing a shared 8 MB cache with about 4 MB of live machine arrays
apiece, so the constraint had stopped being instructions and become memory.

Three of the nine arrays indexed `[period][supplier][firm]` - 19 MB each, 2.4
million doubles - turned out to exist only so that one function could destroy a
column of them.

- `g_c` was a copy of the machine counts, refreshed from `gtemp` every period.
  `COSTPROD()` is its only reader, and it reads only the vintages one firm
  holds, drawing them down as it assigns production. `ENTRYEXIT()` wrote to it,
  but `MACH()` overwrote those writes before anything could read them. It is now
  a per-firm scratch of about eight numbers.
- `g_pb` and `C_pb` hold what each firm wants to scrap and what each of those
  machines costs to run. `SCRAPPING()` writes them only for the entries it
  marks, and `CANCMACH()` reads them only through that same marked list, which
  is already carried alongside. Both are now columns of that list.

Two more, `g_c2` and `g_c3`, are read only by `ADJUSTEMISSENLAB()`, which only
the two `flag_capshocks` branches reach. They are no longer allocated when that
flag is off, which is every run of this project's experiment: 38 MB that was
being faulted in and zeroed at the start of every run for a branch that never
executes.

| measure | before | after |
|---|---:|---:|
| one 600-step run | 0.717 s | 0.642 s |
| saturated throughput, eight concurrent | 221.1 a minute | 320.8 a minute |
| peak resident memory | 163 MB | 74 MB |

The first two rows are the point. A change worth 7.5% on one run bought 41% of
throughput, because it removed memory traffic rather than instructions, and
throughput is what a million runs is billed in. The memory figure travels
further than this machine: at 74 MB instead of 163 twice as many runs fit on a
cluster node, which is what decides how a job array is packed.

### Ranking borrowers, again

Two functions rank each bank's borrowers by debt service, and between them they
were 18% of the run.

`LOANRATES()` sorts 200 firm indices with a comparison that fetches each firm's
ratio through its index, twice for each of the roughly 1400 comparisons a bank's
ranking takes. Letting the ratio travel with the index - sorting pairs rather
than indices - took that block from 9.2% of the run to 3.4%. Nothing about the
order changes: ties break on the firm index, which is distinct, so the ordering
is total and every correct sort returns the same permutation.

`ALLOCATECREDIT()` then serves each bank's customers in that order, and found
the order again from scratch: scan the whole column of ranks for its smallest,
push that entry above all the others so it is not found twice, repeat. Two scans
of 200 entries for each customer of each bank, about 48 million reads a run, for
an order `LOANRATES()` had just computed and nothing had touched since. It now
reads the ranking directly. 0.778 s to 0.740 s.

### A branch the processor cannot predict

`MACH()`'s weighted-average loop reads every firm's count of every vintage, 52
million a run, and acts on the 2% that are not zero. Which 2% is not a pattern
the branch predictor can learn, so the test costs a mispredict most times it
matters. Testing four firms at a time - if all four hold nothing, skip all four
- replaces four unpredictable branches with one that is almost always not
taken. The loop went from 13.0% of the run to 4.9%, and the run from 1.029 s to
0.932 s. Each firm still meets its own terms in the same order and its five
running totals are its own, so nothing accumulates differently.

That measurement also settles what kind of loop it is. The same treatment
applied to `COSTPROD()`'s scan does nothing at all, because that one reads a
firm's counts a machine row apart, four separate cache lines at a time, and it
is the reads rather than the branches that cost. Two sparse scans, two
different limits, and the fix for one is worthless on the other.

### Recording what the next function needs

`COSTPROD()`'s scan is the read-bound one, and the way to make a read-bound
scan cheaper is not to do it. `MACH()`'s weighted-average loop already reads
those counts, so it records, for each firm, the vintages it holds. `COSTPROD()`
reads that list. Its scan went from 9.0% of the run to 0.7%, and the run from
0.895 s to 0.830 s.

This is the same change that was measured and thrown away twice, and it is
worth saying why it works now. The recording is a store inside `MACH()`'s
weighted-average loop. While that loop was branch-bound, a store in the middle
of it cost about what the scan cost, and the two cancelled. Once the four-at-a-
time test made the loop cheap, the store became cheap with it and the saving
came through. A change that measures nothing can be worth retrying after the
code around it changes character - but only after, and only by measuring again.

### An address worked out 52 million times

`COSTPROD()` scans a firm's vintages by subscripting `g_c`, which works the
address out from the array's base and its two strides on every one of the 52
million reads, and reloads all three each time round because the loop also
writes to arrays the compiler cannot prove sit elsewhere in memory. One vintage
on is a fixed distance, so the scan now walks a pointer. 4.4%.

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
element accesses in the model at the time did not have to change: the subscripts
return small proxies, so `X[tt-1][i-1][j-1]` still parses and still means the
same thing. Kept. Three of those nine arrays have since gone entirely, and the
count of accesses with them; the section above has that.

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

Three more from the third round.

*Taking `SCRAPPING()`'s four vintage rows once instead of subscripting.* The
same change that is worth 4.4% in `COSTPROD()`, applied where the firm index is
already innermost: 0.85% slower over six alternating pairs. The address there
advances by one element, which the compiler strength-reduces on its own, and the
hoisting is paid on every vintage row while 98% of the rows do nothing.

*Recording the holders in a pass of its own* rather than inside `MACH()`'s
existing loop: 6.8% slower over nine pairs, 1.099 s against 1.174 s.

*Building `COSTPROD()`'s per-firm lists in `MACH()` too.* The obvious companion
to the change that worked for `SCRAPPING()`, and at the time it measured
nothing: ten alternating pairs each way, and the sign of the difference changed
with the order the two builds were run in. It is in the model now. What changed
is the loop it records from, which was branch-bound then and is not now; the
section above has the account. The lesson is not that the measurement was wrong,
it is that a measurement is of one version of the code and does not carry over
to another.

*Testing four vintages at a time in `COSTPROD()`'s scan.* The change that is
worth 9% in `MACH()`, applied to the other sparse scan: 8.97% of the run before,
9.30% after, which is nothing. The two scans are limited by different things.

*Working out each supplier's offer once in `BROCHURE()`.* The comparison that
picks a machine supplier rebuilds a four-division cost expression for both
candidates on each of 4000 firm-supplier pairs a period, and every operand is
fixed while the firms choose. Computing the 20 suppliers' offers once instead:
nine alternating pairs said 0.7% slower, ten in the reverse order said 0.1%
faster. Nothing, and it rewrites an expression of the published model for no
return, so it is out.

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

Two more tests widen that in the two directions the aggregate comparison leaves
open, and both are worth running before any further change to the model.

`tests/dsk_full_output_equivalence.c` runs the model with `-f 1`, which writes
twelve files of per-firm quantities instead of the one file of economy-wide
totals: each K-firm's and C-firm's productivity, energy efficiency,
environmental friendliness and net worth, each bank's net worth, and each
C-firm's debt. 13 files and 28 MB a run against 3 MB, and per firm rather than
summed, so a change that moved one firm's machines to another firm while leaving
the totals alone would pass the aggregate comparison and fail this one. That is
the failure this test exists for, because several of the changes above rest on
an argument about which values are ever read.

`tests/dsk_design_equivalence.c` runs the model away from the baseline
calibration. Every change above was measured and checked at the one parameter
vector the model ships with, and the experiment runs a thousand others, where
different numbers of firms enter and exit and the machine arrays are sparse in
different places - which is what several of those changes turn on. It takes four
points of `dataset/abm_system_design.csv`: the design's column-wise minima and
maxima, which are the corners of the box it actually reached, and two rows the
experiment will really simulate. Each point is compared on three things: the
exit status, the results file and the error log. A parameter vector the model
refuses to simulate is a valid outcome and still a comparison, as long as both
builds refuse it the same way.

    make test-dsk_full_output_equivalence
    make test-dsk_design_equivalence

`make asan` builds every test with AddressSanitizer and UndefinedBehaviorSanitizer
and runs it, these three included. That is what catches the class of mistake
where a test itself is wrong rather than the model: writing
`tests/dsk_design_equivalence.c` produced one, freeing a `Mat` that
`df_col_numeric` had returned as a view into the dataframe rather than as
memory of its own, and under the sanitizer it reports as `bad-free` at the line
that did it instead of as a bare abort.

### The one thing these tests do not cover

All three run with the model's seven shock flags at zero. That is how
`dsk_sfc_inputs.json` ships and how every run of this project's experiment is
configured - `applications/abm_system_simulate.c` writes only the `params`
block, never `flags` - and `docs/ABM_SYSTEM_SIMULATION.md` sets out why they
should stay there. The gap is deliberate and is not going to be closed, because
a test at settings the experiment will not use would guard code the experiment
will not run.

What that leaves untested is specific and worth naming, in case someone later
wants the shock scenarios. Two changes have a branch that only executes with a
flag on:

- `MACH()` skips the firm-vintage pairs with no machines only while every
  firm's two shock factors are non-zero, and adds every term as before when one
  is not. `shocks_labprod2` and `shocks_eneff2` are written only inside
  `if(flag_prodshocks2==1)` and `==2` in `modules/module_climate_sfc.cpp`, so
  with that flag at zero they hold the zero they are initialised to and the skip
  is always taken. The fallback has never run.
- `PRODMACH()`'s restructured loop sits beside the `flag_capshocks` block, and
  `ENTRYEXIT()` and `MACH()` both skip work on `g_c2` and `g_c3` under the same
  flag. None of those branches has run either.

Anyone turning a shock flag on should add points at those flag settings to these
tests first, and check that byte equality still holds there.

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
| 1 | 89.8 | 1.4 |
| 4 | 279.4 | 3.5 |
| 8 | 320.8 | 3.6 |
| 16 | 281.5 | 9.5 |

The machine saturates near 321 runs a minute at eight concurrent, which is a
million runs in 2.2 days. Sixteen is slower than eight, which is what running
two threads on each of eight cores does when both threads want memory. The
experiment is a million independent runs, so the cores are already full and
making one run use several of them cannot raise that number. It would only
shorten a single run's latency, which nothing here needs.

Note what does not carry over, because this is now the binding constraint. One
run is 52 times faster than upstream's, and saturated throughput is 4.2 times
what it was before this work - 321 against the 76 measured at the same point in
the history.

The two do not move together, and which one a change moves says what it removed.
Work that took out instructions raised the single-run time steadily and then
stopped raising throughput at all: two changes worth 9.8% and 4.9% on one run
moved saturated throughput from 227 to 221, which is noise. Eight runs at once
were competing for one memory system, so the constraint had stopped being
instructions. Work that took out memory moved it again: dropping three of the
nine machine arrays is worth 7.5% on one run and 41% on throughput, and the
footprint went from 163 MB to 74. That is also what decides how many runs fit on
a cluster node, which no single-run figure shows. The gap is memory bandwidth: eight
runs at once want their vintage arrays in cache at the same time, and there is
one memory system between them. That also means further work on a single run
buys steadily less of what this experiment actually costs.

The upstream column was measured earlier, at the same concurrencies but before
the second round of changes, and it is noisier than this one - its 8-way row is
slower than its 16-way, which cannot be right and means something else was on
the machine during that batch.

Memory is the other limit, and for most of this work it was the binding one. A
run peaked at 165 MB resident in both builds, so sixteen at once was 2.6 GB on a
7 GB machine, and eight at once held far more live machine array than the shared
8 MB cache. Dropping three of the nine machine arrays took this build to 74 MB
against upstream's 166, which is where the jump from 221 to 321 runs a minute
came from. On a cluster the same figure decides how many jobs fit on a node.
