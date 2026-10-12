# Documentation audit, 2026-10-11

Scope: README.md, every file in docs/, and the comment block at the top of the
Makefile. Checked by reading every doc, by searching README.md and docs/ for the
name of every tracked source file and every make target, by checking that every
cited path exists, and by spot-checking claims against the code. Line numbers
refer to the working tree on 2026-10-11, including the uncommitted edit to
docs/MCS_VALIDATION_PROCEDURE.md.

Items are ordered by how much they would mislead a reader, worst first within
each section.


## 1. Undocumented code

### Not mentioned in any doc

- studies/qvarma_stuck_fits.c, studies/qvarma_conditioning.c,
  studies/qvarma_convergence_criteria.c, studies/us_qvarma_nu_sensitivity.c,
  studies/abm_system_winner_nu_profile.c,
  studies/abm_system_winner_nu_likelihood_scan.c,
  studies/abm_system_winner_tail_comparison.c,
  studies/abm_system_winner_diagnostics_plots.py.
  All eight are in `make study` and their outputs are tracked under out/
  (out/qvarma_*.txt, out/abm_system_winner_nu_*, out/abm_system_winner_tail_comparison.txt,
  out/us_qvarma_nu_sensitivity*, out/abm_system_winner_diagnostics_plots/),
  but no doc says what question each answers or what it found.
  studies/abm_system_winner_normality.c is mentioned only as the place Mardia's
  test lives (ABM_SYSTEM_MCS_VALIDATION.md:240); its result is not reported.
- benchmarks/: all seven files (qvarma_fit_cost, qvarma_fit_io,
  qvarma_iteration_budget, qvarma_taped_vs_fused, qvarma_thread_scaling,
  qvarma_process_scaling, fit_contention_source). The only description is the
  Makefile comment above BENCH_STEMS. `make bench` is not in the README.
- montecarlo/fit_metric_pairwise.c and montecarlo/fit_metric_table.h. The
  pairwise program is not in MONTECARLO_STEMS, so no make target builds it;
  its header comment points at a section "Pairwise tests" of
  docs/MONTECARLO_FIT_METRIC_COMPARISON.md that no longer exists, and that doc
  now says it uses one test over all five models "rather than pairwise tests".
  fit_metric_table.h says it serves `bin/fit_metric_comparison`, which does not
  exist. The tracked log montecarlo/out/fit_metric_pairwise.log holds only
  "Failed to find executable /usr/bin/time". Decide whether this program is
  alive; if not, it belongs in old/ or nowhere.
- The throughput study (applications/throughput_*.c, out/throughput_*,
  out/fit_speedup_options.txt): one bullet in README.md:631-633, no write-up of
  what was measured or found.

### Documented incompletely

- README.md "Tests" (659-728) lists two data tests and twelve simulator tests.
  `make test` also runs tests/marginal_laws_correctness.c,
  tests/lp_system_transform.c and tests/fit_metric_correctness.c (TEST_STEMS in
  the Makefile). They are described in other docs but missing from the list
  that claims to be what `make test` runs.
- README.md "Layout" (221-244) has no entry for montecarlo/ (a whole pipeline
  with its own out/) or bin/.
- README.md lists no study by name; "make study" (718-723) is the only mention.
- Tracked outputs no doc mentions: out/refactor_equality_check/,
  out/winner_normality_rerun/, montecarlo/out/lp_growth_run.log.


## 2. Outdated docs

### Wrong about what the pipeline does

- docs/DATA_DOCUMENTATION.md, whole file. It says GDP, CPI and energy are kept
  as log levels because "differencing a co-integrated variable throws the
  relation away" (34-38), and that the growth rate was dropped "because the
  model uses the level, not the growth rate" (119-121, 135-136). The t-QVARMA
  and the local projections are fitted on GDP growth, energy growth, employment
  change, inflation and the interest rate: applications/us_data.h
  `us_system_growth_block`, applications/us_qvarma_spec_choice.c `build_block`,
  ABM_SYSTEM_LP_VALIDATION.md:172-175, MONTECARLO_VALIDATION.md:62-64. The doc
  describes out/us_system.csv correctly but never mentions the differencing
  step applied to it, so it argues for the opposite of what is fitted.
- docs/MCS_VALIDATION_PROCEDURE.md:19 (the uncommitted edit) says the procedure
  fits "the t-QVARMA(1,1,2), the linear projection, and the nonlinear
  projection with combined states" to the benchmark. On the US data the
  combined-states set is not computed (ABM_SYSTEM_LP_VALIDATION.md:232-234), and
  the rest of the same doc still describes a t-QVARMA-only procedure
  (49-51 "with QVARMA standing in for LP", 128-160, column names at 139-141).
- docs/ABM_SYSTEM_MCS_VALIDATION.md:589-595 and 807-814 say the absolute
  validation score (v = 1/(1+s)) "is not implemented here". It is implemented,
  oracle and bootstrap versions, in montecarlo/fit_metric_*.c and reported in
  docs/MONTECARLO_FIT_METRIC.md; what is missing is its value on the US data.
  The formula given there, s = d_i / sigma^2_rw, divides by a variance, while
  MONTECARLO_FIT_METRIC.md:78-95 divides by a standard deviation and explains
  why.
- docs/ABM_SYSTEM_MCS_VALIDATION.md:801-802: "64.8% of the fits did not
  converge, 62.1% of them because the line search could not lower the
  objective". 621,485 is 62.1% of all fits, which is 95.9% of the unconverged
  ones.
- docs/ABM_SYSTEM_SIMULATION.md:231 "Eight are kept", followed by a table of
  five; 247 says the other three are not stored.
- docs/ABM_SYSTEM_SIMULATION.md:394-399 projects 10.6 GB from the old .Rdata
  replications; the run stored 15 GB (same doc, 580).
- docs/ABM_SYSTEM_SIMULATION.md:589-595 recommends a second pilot "before
  committing two days". The full run finished on 2026-09-07.
- docs/MONTECARLO_FIT_METRIC.md:226-229 cites mean losses "published" in
  docs/MONTECARLO_VALIDATION.md for cop_0191 run 706. That single benchmark was
  moved to docs/OLD_MONTECARLO_BENCHMARK_SUBSETS.md and
  docs/OLD_MONTECARLO_LP_VALIDATION.md.
- docs/MONTECARLO_FIT_METRIC.md:335-340 and MONTECARLO_FIT_METRIC_COMPARISON.md:181
  say nothing of the comparison is implemented, while
  montecarlo/fit_metric_pairwise.c implements a pairwise version (see section 1).
- docs/OLD_MONTECARLO_LP_VALIDATION.md:22-28 describes
  docs/MONTECARLO_VALIDATION.md as a one-run experiment; it is now the
  1000 by 1000 grid.
- docs/MONTECARLO_VALIDATION.md:244-248: the tracked figures and
  sweep_grid_identifiability.csv cover two of the five models. Either rerun the
  script or keep saying so in the output table, not only in the prose.
- README.md:490-491 "The figures use Python ...; everything else is C".
  montecarlo/sweep_grid_cop_0191.py, sweep_grid_presence.py and
  compressed_response_learn.py / _report.py produce tables and learned maps.
- README.md:262-269 Requirements list polars, plotly, kaleido only.
  The compressed-response scripts need numpy, scipy and pandas
  (MONTECARLO_COMPRESSED_RESPONSE.md:272-273).

### Internal inconsistencies and stale counts

- README.md:119-123 "Eleven tests ... the first six ... the last five", followed
  by nine bullets; README.md:671 "eleven of the simulator"; README.md:686 "The
  simulator's twelve". The Makefile has eleven dsk stems plus
  dsk_bulk_cancellation_distribution, twelve.
- docs/DSK_MODEL_CHANGES.md:1061 "None of the ten tests".
- Machine memory: 7.5 GB (README.md:291, DSK_MODEL_CHANGES.md:126) against
  7.1 GB (ABM_SYSTEM_MCS_VALIDATION.md:776, MCS_VALIDATION_PROCEDURE.md:184).
  OLD_MONTECARLO_BENCHMARK_SUBSETS.md:226 says "16 cores"; the machine has 8
  cores, 16 threads.
- docs/DSK_LONG_HORIZON.md line references into the source are stale:
  `dsk_sfc_main.cpp:302` (75; the message is at 325), `:5286` (133; the
  assignment is now at 5344), `:2109` (150; now 2163),
  `module_macro_sfc.cpp:46` (162; that loop was replaced by the bulk
  cancellation). Name the function instead of the line.
- Makefile:6-8 says every test is the DSK simulator's; four are not.
- Makefile:30-32 says `make applications` fits the model to the real dataset;
  APPLICATION_STEMS is us_prepare_data and abm_system_design, neither fits
  anything.
- Makefile:473 "with the three changes docs/DSK_MODEL_CHANGES.md records"; that
  doc records more than twenty.
- docs/OLD_MONTECARLO_LP_VALIDATION.md:199 cites docs/LP_DOCUMENTATION.md, which
  does not exist.
- docs/MONTECARLO_FIT_METRIC_COMPARISON.md announces two nulls (193-196) and
  has a section only for null 2. The file ends mid-thought (282).


## 3. Unprofessional docs

### Citing files a reader cannot get

- `_temp/` is gitignored, yet these are cited as the source of a method:
  `_temp/Note on non-lin LP.pdf` (README.md:566, MONTECARLO_VALIDATION.md:55,80,
  MONTECARLO_COP_0191.md:269, MONTECARLO_COMPRESSED_RESPONSE.md:204),
  `_temp/main_code.R` (README.md:565, OLD_MONTECARLO_LP_VALIDATION.md:27),
  `_temp/functions.R` (OLD_MONTECARLO_LP_VALIDATION.md:27), `_temp/estimate.R`
  (ABM_SYSTEM_LP_VALIDATION.md:180). Also `rw_data_prep.R`
  (OLD_LP_R_LEVELS.md:150).
- "The thesis" / "Fabiano" (an unpublished thesis, as DATA_DOCUMENTATION.md:23
  admits): MCS_VALIDATION_PROCEDURE.md:52-73, 84, 96;
  DATA_DOCUMENTATION.md:20-24, 34-42, 60-66; ABM_SYSTEM_LP_VALIDATION.md:232-239;
  OLD_LP_R_LEVELS.md:150; ABM_SYSTEM_MCS_VALIDATION.md:592, 807.
  MCS_VALIDATION_PROCEDURE.md "The procedure this replicates" (34-63) is a
  summary of that thesis's protocol, including steps this project does not
  follow (MSE, local projections). Either cite it with full reference and
  a link, or state the method in this project's own terms.
- "the collaborator", "the collaborator's R pipeline", "R's transform_data()",
  "the scripts received" (README.md:565, ABM_SYSTEM_LP_VALIDATION.md:166,
  OLD_EXPERIMENTS.md:82, OLD_LP_R_LEVELS.md:108, OLD_MONTECARLO_LP_VALIDATION.md:26,
  194). No name, no reference, no code in the repository. A reader cannot tell
  what "R's settings" are except by the numbers repeated next to them.
- `qvarma_data.txt` and the transformation-search script, neither in the
  repository: DATA_DOCUMENTATION.md:10-16 and the whole background section
  79-160, which says itself it "is not something that can be re-run here".
- The old .Rdata dataset ("Dati simulati.zip", gitignored), named
  `EstimationSeriesSample1_N`: ABM_SYSTEM_MCS_VALIDATION.md:753-769,
  ABM_SYSTEM_SIMULATION.md:394-399, MCS_VALIDATION_PROCEDURE.md:111-126, 144-152.
- Gitignored logs cited as the evidence for timings: out/sweep_grid.log and
  out/lp_nl_run.log (MONTECARLO_VALIDATION.md:220-224).
- Things outside the repository: `../et_al./KNOWN_ISSUES.md` (Makefile:18),
  "a scratch copy of the source outside this repository"
  (DSK_LONG_HORIZON.md:52-55), "A second pipeline wants one run ... for 400,200
  periods" (DSK_LONG_HORIZON.md:4-5), which no doc or code here describes.
- MCS_VALIDATION_PROCEDURE.md:128-135: "Its code is not part of this
  repository" for t-QVARMAd.

### Terms and symbols used without definition

- The nine design parameters (Gamma, chi, psi1, psi3, alfa, taylor1, taylor2,
  taylor, kappa) are never explained in any doc. ABM_SYSTEM_SIMULATION.md:48-54
  declines to, and says the two guesses it offers are not "verified against the
  model's source here", although the source is in model/dsk_sfc.
- MCS_TR and MCS_TMAX are defined only in MCS_VALIDATION_PROCEDURE.md:164-168,
  then used bare in ABM_SYSTEM_LP_VALIDATION.md, MONTECARLO_VALIDATION.md,
  MONTECARLO_COP_0191.md, MONTECARLO_COMPRESSED_RESPONSE.md and every OLD_ doc.
- `LM` as a label for the weighted score loss is defined only in
  OLD_SCORE_LOSS.md:111; OLD_MONTECARLO_BENCHMARK_SUBSETS.md:131-192 uses it in
  tables without definition.
- `L_bar_i`, `sigma_i^2` (MCS_VALIDATION_PROCEDURE.md:58), `sigma^2_rw`
  (ABM_SYSTEM_MCS_VALIDATION.md:808), never defined.
- Model code identifiers in prose: `GDP_r(1) = Q1tot * dim_mach + Q2tot`
  (ABM_SYSTEM_MCS_VALIDATION.md:245), K-firms and C-firms, `flatten_estimated`,
  `stats_mse` (MCS_VALIDATION_PROCEDURE.md:77-79).
- "HP cycle" (ABM_SYSTEM_LP_VALIDATION.md:182), spelled out only in
  MONTECARLO_VALIDATION.md:73.
- "the older dataset" (ABM_SYSTEM_SIMULATION.md:267, 275, 278, 328), never said
  what it is in that doc.
- "Sims et al. 1990; Cochrane 1997" (DATA_DOCUMENTATION.md:37), and other
  author-year citations, with no reference list anywhere.
- "single-benchmark experiment" (MONTECARLO_COP_0191.md:243) without saying it
  is discarded and where it is.
- `qvarma.h` (ABM_SYSTEM_MCS_VALIDATION.md:17), `us_real.csv`
  (DATA_DOCUMENTATION.md:5) without their paths.

### History and argument instead of description

These paragraphs narrate how the file or the code got here, or answer an
objection, rather than describe what is. They belong in git history or nowhere.

- MCS_VALIDATION_PROCEDURE.md:75-79 (an earlier version compared parameters),
  109-126 (the MSE episode, "verified directly rather than assumed"),
  128-135 (t-QVARMAd), 144-160 (p1q1r4, "it answered"), 181-186 (an old et_al
  allocation). The decisions can stay as one sentence each; the stories cannot.
- DSK_MODEL_CHANGES.md:19-26 (a bug in build.sh's cleanup trap),
  1025-1098 "The tests were reviewed as well" (a changelog of test fixes),
  1201 "Stated plainly, because a reader should not have to work it out."
- README.md:250-256 and 653-657 ("each of them was violated somewhere before",
  "were removed once the runs that replaced them were on disk").
- ABM_SYSTEM_SIMULATION.md:213-216 "Any claim that it does should be read
  against this measurement."
- MONTECARLO_VALIDATION.md:145 "which is cheaper than arguing it away".
- OLD_MONTECARLO_BENCHMARK_SUBSETS.md:108-111 says cop_0191 being the US winner
  is "a consequence of the convergence criterion rather than a choice", but
  requirement one (60-63) restricts the pool to the US set, which is cop_0191
  alone. The sentence defends the choice and gets the logic backwards.

### Filler and rhetorical phrasing

- "Measured rather than argued" (OLD_SCORE_LOSS.md:89, DSK_LONG_HORIZON.md:423),
  "Measured rather than assumed" (README.md:141), "rather than assumed away"
  (OLD_SCORE_LOSS.md:166).
- OLD_SCORE_LOSS.md:80-81 "it is not a small one", 125-126 "This is the whole
  reason for the weighting; the null distribution below is a bonus", 144 "with
  no sign to get wrong", 150 "the entire job".
- OLD_MONTECARLO_BENCHMARK_SUBSETS.md:150 "This is not a marginal miss", 190
  "for what it is worth", 192 "It is simply not enough to win".
- ABM_SYSTEM_MCS_VALIDATION.md:811 "the largest of the open questions rather
  than one of them".
- README.md:145-146 and DSK_MODEL_CHANGES.md:1000 "no room left ... to hide".
- DSK_LONG_HORIZON.md:300 "honest".
- ABM_SYSTEM_SIMULATION.md:315-318 "Simulation is the whole pipeline's only
  step that is not C, and it is not one either", which contradicts itself.
- MONTECARLO_FIT_METRIC_COMPARISON.md:282: one unwrapped line, "Montecarlo",
  "normal friedman test", ends without a conclusion.

### Structure

- ABM_SYSTEM_MCS_VALIDATION.md:589-595 ("None of this makes the confidence set
  wrong") sits at the end of the correctness-test section and refers to
  "What the confidence set does not say", 490 lines earlier.
- ABM_SYSTEM_MCS_VALIDATION.md is titled as the result on the US data but holds
  four tail and distribution studies (166-595) that are about the simulator,
  not the validation. A reader looking for the result has to skip 500 lines.
- docs/Model_Simulation.Rmd is a third party's file, tracked, with the
  author's name in it (line 3). README.md:43-45 says it is not this project's
  writing. Confirm it may be published before the repository goes public.
- DSK_MODEL_CHANGES.md:94-97 leaves an open to-do in the doc ("One thing to
  check before this repository goes anywhere public: upstream ships no licence
  file").
