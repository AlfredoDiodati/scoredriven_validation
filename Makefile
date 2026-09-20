# t-QVARMA(p,q,r), Blazsek, Escribano and Licht (2023).
#
# Header-only against et_al.; nothing here builds a library of its own.
#
#   make                          build every test binary
#   make test                     every test script. All of them are the DSK
#                                 simulator's: the auxiliary model is et_al.'s
#                                 and is tested there, not here.
#   make test-stress              every test script, slow checks included
#   make test-<stem>              one of them on its own, one stem per entry in
#                                 TEST_STEMS
#   make test                     float64 by default. The analytic gradient is
#                                 identical in both builds; float64 is the
#                                 default because et_al.'s _syevd (used by
#                                 standard_errors) has been observed to fail
#                                 to converge under float32 on a fitted
#                                 t-QVARMA Hessian - see
#                                 ../et_al./KNOWN_ISSUES.md.
#   make MAT_DOUBLE=0 test        float32, not recommended; see above.
#   make study                    the robustness studies on the pipeline's own
#   make study-robustness         result, the same target under the name that
#                                 says what they are. Every one of them reads
#                                 what the pipeline wrote to out/ and dataset/ -
#                                 15 GB of simulation and a million cached fits -
#                                 and none rebuilds it, so none runs on a fresh
#                                 clone.
#   make study-<stem>             one study on its own
#   make bench                    run the benchmarks (never part of test)
#   make bench-performance ETAL_DEV=1
#                                 the same, built against the development et_al.
#   make applications             fit the model to the real dataset, results to
#                                 out/. Always built in float64. One app-<stem>
#                                 target per stem in APPLICATION_STEMS.
#   make app-<stem>                one application on its own, EXPERIMENT_STEMS
#                                 included - one-off structural searches
#                                 already answered (see that variable's own
#                                 comment), not part of `make applications`
#                                 so a routine build does not redo minutes of
#                                 already-settled work every time.
#   make asan                     under AddressSanitizer and UndefinedBehaviorSanitizer
#
# Three kinds of script, three directories, one stem list each. A stem maps to
# <directory>/<stem>.c and bin/<stem>, and the targets below are generated from
# the lists rather than written by hand.
#
#   tests/       pass or fail. A failure stops `make test` and the build.
#   benchmarks/  how long something takes. Never part of `make test`: a function
#                that returns the wrong answer quickly is not fast.
#   studies/     a question about behaviour, answered into out/. No verdict, so
#                nothing here can gate anything.
#
# The directory is what says which of the three a script is, so adding one means
# putting it in the right tree and adding its stem to the matching list.

CC ?= gcc
CFLAGS ?= -O2 -march=native -Wall -Wextra -std=c11
ETAL_CFLAGS := $(shell pkg-config --cflags et_al.-core)
ETAL_LIBS := $(shell pkg-config --libs et_al.-core)
# et_al.'s own -I covers <et_al./...> and the bare <frame/...> spellings. The
# repository root is added for applications/abm_system.h and for the tests' own
# headers, which are included by their directory-qualified path so a reader can
# see which tree a header comes from.
INCLUDES := -I.

# Every header under the installed et_al., so that reinstalling et_al. (its own
# `sudo make install`) after a change is what makes the next `make` here
# actually recompile, rather than seeing every tracked prerequisite older than
# an existing binary and silently doing nothing. et_al. is header-only, so this
# is the only place its own changes take effect at all.
ETAL_INSTALLED_HEADERS := $(shell find $(shell pkg-config --variable=includedir et_al.-core) \
                                    -name '*.h' 2>/dev/null)

# float64 by default; make MAT_DOUBLE=0 for float32.
MAT_DOUBLE ?= 1
ifneq ($(MAT_DOUBLE),0)
CFLAGS += -DMAT_DOUBLE
endif

# The model, the solver, the unit root tests, Ljung-Box and the QLR test all
# live in the installed et_al. now, so this project carries no header of its
# own beyond the two under applications/ that describe its own data. A rebuild
# on an et_al. change is already covered by ETAL_INSTALLED_HEADERS above.
HEADERS :=
TEST_HEADERS :=
# Every stem here is the DSK simulator's. The auxiliary model lives in the
# installed et_al. and its correctness suite lives there with it; a copy kept
# here drifted out of date against the library it was testing and failed on a
# contract et_al. had since changed, which is worse than not having one.
# abm_system_layout is this project's own and tests the one header both dataset
# writers and every reader agree about the stored layout through.
TEST_STEMS := abm_system_layout dsk_long_path dsk_build_equivalence \
               dsk_full_output_equivalence dsk_design_equivalence \
               dsk_memory_safety dsk_ulp_sensitivity dsk_redenomination_invariance \
               dsk_machine_lot_rebase dsk_good_unit_invariance \
               dsk_dataset_reproduction
# Where the wall time of a t-QVARMA fit goes, and what each way of speeding it
# up is worth. Measured 2026-08-29 against a 500,000-fit run of
# abm_system_fit_qvarma; out/fit_speedup_options.txt collects the numbers and
# the setup behind them. fit_contention_source is the only one that touches
# neither et_al. nor this project: it times cblas and malloc on their own,
# because that is what isolates which of the two the taped filter contends on.
BENCH_STEMS := qvarma_fit_cost qvarma_fit_io qvarma_iteration_budget \
                qvarma_taped_vs_fused qvarma_thread_scaling qvarma_process_scaling \
                fit_contention_source
# Robustness checks on the result of the main pipeline, not steps of it. Every
# one of them reads what the pipeline already wrote to out/ and dataset/ and
# never rebuilds it, so none can run on a fresh clone: `study` checks for the
# dataset and the fit cache first and says so rather than letting the first
# binary fail on an assert.
#
# The order is load-bearing in two places, which is why this is a list and not a
# set: abm_system_winner_diagnostics writes the theta table
# abm_system_winner_normality reads, and us_qvarma_nu_sensitivity writes the
# held-nu US fits abm_system_winner_tail_comparison reads. Each producer comes
# before its consumer here and `study` runs them in this order.
STUDY_STEMS := qvarma_stuck_fits qvarma_conditioning qvarma_convergence_criteria \
               us_qvarma_nu_sensitivity abm_system_winner_diagnostics \
               abm_system_winner_nu_profile abm_system_winner_nu_likelihood_scan \
               abm_system_winner_tail_comparison abm_system_tail_origin \
               abm_system_seed_correlation abm_system_winner_normality
# us_prepare_data has to come first: it is the only script that reads
# us_real.csv's Unemployment and Des_Energy_demand columns, and it writes
# out/us_system.csv, which us_data.h reads back. The order here is what makes
# `applications`' shell loop run it first; the app-<name> targets that need it
# also depend on it directly below, for a standalone `make app-<name>` to
# regenerate the CSV rather than read a stale one.
#
# Only data preparation is routine. A model fit already takes real time and
# does not change on a routine basis, so the fitting scripts and the whole
# ABM chain sit in EXPERIMENT_STEMS and are run one at a time via
# `make app-<name>`.
APPLICATION_STEMS := us_prepare_data abm_system_design
#
# throughput_fit is deliberately not a stem here. Its solver budget
# is a compile-time constant that names every file it writes, so one binary per
# budget is built below from THROUGHPUT_ITERATION_CAPS instead, and the generic
# app-<stem> rule - which would build one unnamed binary at the default budget -
# would give a second way to write the same tree.
EXPERIMENT_STEMS := us_qvarma_spec_choice \
                     abm_system_convert_rdata abm_system_fit_qvarma \
                     abm_system_irf_loss abm_system_mcs \
                     abm_system_mcs_statistic_comparison \
                     abm_system_winner_irf \
                     throughput_dataset
# Whatever the application scripts share, so editing it rebuilds them.
APPLICATION_HEADERS := applications/us_data.h applications/abm_system.h
BIN := bin
OUT := out

TEST_BINARIES := $(addprefix $(BIN)/,$(TEST_STEMS))
# The one C++ test, built by its own rule further down.
TEST_BINARIES += $(BIN)/dsk_bulk_cancellation_distribution
BENCH_BINARIES := $(addprefix $(BIN)/,$(BENCH_STEMS))
STUDY_BINARIES := $(addprefix $(BIN)/,$(STUDY_STEMS))
APPLICATION_BINARIES := $(addprefix $(BIN)/,$(APPLICATION_STEMS))
EXPERIMENT_BINARIES := $(addprefix $(BIN)/,$(EXPERIMENT_STEMS))

# Benchmarks need clock_gettime, which -std=c11 hides, and they are never part
# of `test`: a function that returns the wrong answer quickly is not fast.
# -fopenmp because the scaling benchmarks time one thread against four.
BENCH_CFLAGS := $(CFLAGS) -D_POSIX_C_SOURCE=199309L -fopenmp

# ETAL_DEV=1 builds against the development et_al. at ETAL_DEV_PATH instead of
# the installed headers, which is how a change to et_al. is measured before and
# after without installing it. Correctness targets always use the installed
# copy, since that is what the project actually ships against.
ETAL_DEV_PATH ?= /home/dioda/Documents/_py
ifdef ETAL_DEV
BENCH_ETAL_CFLAGS := -I$(ETAL_DEV_PATH) $(shell pkg-config --cflags openblas 2>/dev/null) -fopenmp
BENCH_ETAL_HEADERS := $(shell find $(ETAL_DEV_PATH)/et_al. -name '*.h' 2>/dev/null)
else
BENCH_ETAL_CFLAGS := $(ETAL_CFLAGS)
BENCH_ETAL_HEADERS := $(ETAL_INSTALLED_HEADERS)
endif

.PHONY: all test test-stress bench bench-performance study applications asan clean

all: $(TEST_BINARIES)

$(BIN) $(OUT):
	mkdir -p $@

$(BIN)/%: tests/%.c $(HEADERS) $(TEST_HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) $(ETAL_CFLAGS) $(INCLUDES) $< -o $@ $(ETAL_LIBS)

# One pair of targets per stem. These are generated rather than pattern rules
# because make skips the implicit rule search for phony targets, so a
# `test-%:` pattern would silently report nothing to be done.
define test_targets_for_stem
.PHONY: test-$(1) test-$(1)-stress
test-$(1): $(BIN)/$(1) | $(OUT)
	./$(BIN)/$(1)
test-$(1)-stress: $(BIN)/$(1) | $(OUT)
	STRESS=1 ./$(BIN)/$(1)
endef
$(foreach stem,$(TEST_STEMS),$(eval $(call test_targets_for_stem,$(stem))))

# Benchmarks are built separately because they need different flags and must
# not be swept into the aggregate test targets.
define bench_target_for_stem
.PHONY: bench-$(1)
$(BIN)/$(1): benchmarks/$(1).c $(HEADERS) $(BENCH_ETAL_HEADERS) | $(BIN)
	$(CC) $(BENCH_CFLAGS) $(BENCH_ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
bench-$(1): $(BIN)/$(1) | $(OUT)
	./$(BIN)/$(1)
endef
$(foreach stem,$(BENCH_STEMS),$(eval $(call bench_target_for_stem,$(stem))))

bench-performance: $(BENCH_BINARIES) | $(OUT)
	@for binary in $(BENCH_BINARIES); do ./$$binary || exit 1; done

# Studies are neither correctness gates nor speed measurements: they answer a
# question about the model's behaviour and write their answer to out/. Kept out
# of `test` because they take minutes and have no pass or fail.
define study_target_for_stem
.PHONY: study-$(1)
$(BIN)/$(1): studies/$(1).c $(HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -fopenmp $(ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
study-$(1): $(BIN)/$(1) | $(OUT)
	./$(BIN)/$(1)
endef
$(foreach stem,$(STUDY_STEMS),$(eval $(call study_target_for_stem,$(stem))))

# What every study reads and no study rebuilds. Checked once here because the
# alternative is the first binary aborting on an assert several minutes in.
study-inputs:
	@test -d dataset/abm_system || { \
	  echo "dataset/abm_system/ is missing. The studies read the simulated dataset"; \
	  echo "and never rebuild it; applications/abm_system_simulate_all.sh produces it."; \
	  exit 1; }
	@test -d out/abm_system_fit_qvarma || { \
	  echo "out/abm_system_fit_qvarma/ is missing. The studies read the fit cache and"; \
	  echo "never rebuild it; make app-abm_system_fit_qvarma produces it."; \
	  exit 1; }
	@test -f out/abm_system_mcs.csv || { \
	  echo "out/abm_system_mcs.csv is missing. The studies read the confidence"; \
	  echo "set and never rebuild it; make app-abm_system_mcs produces it."; \
	  exit 1; }

.PHONY: study-inputs study-robustness
# One name per thing: study-robustness says what these are, study is the name
# the README and the habit already use.
study-robustness: study
study: study-inputs $(STUDY_BINARIES) | $(OUT)
	@for binary in $(STUDY_BINARIES); do ./$$binary || exit 1; done

# Applications fit the model to a real dataset and write their results to out/.
# Neither a test nor a study on simulated data, so they are their own category
# and are never swept into `test`. -DMAT_DOUBLE is unconditional here rather
# than a switch: the float32 filter aborts once mu_dag's random walk leaves
# float32's range, which is not a risk worth leaving to a command line.
define application_target_for_stem
.PHONY: app-$(1)
$(BIN)/$(1): applications/$(1).c $(HEADERS) $(APPLICATION_HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -DMAT_DOUBLE -fopenmp $(ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
app-$(1): $(BIN)/$(1) | $(OUT)
	./$(BIN)/$(1)
endef
$(foreach stem,$(APPLICATION_STEMS),$(eval $(call application_target_for_stem,$(stem))))
$(foreach stem,$(EXPERIMENT_STEMS),$(eval $(call application_target_for_stem,$(stem))))

# abm_system_simulate is built but never run by a target of its own: with no
# arguments it simulates the whole design, a million runs and about a fortnight
# of this machine, which is not something a make target should start.
# app-abm_system_design draws the design it reads, and keeps an existing one
# rather than renumbering configurations under results already on disk.
$(BIN)/abm_system_simulate: applications/abm_system_simulate.c $(HEADERS) $(APPLICATION_HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -DMAT_DOUBLE -fopenmp $(ETAL_CFLAGS) $(INCLUDES) $< -o $@ $(ETAL_LIBS)

# us_data.h reads out/us_system.csv rather than the raw file, so a standalone
# run has to regenerate it first rather than trust whatever an earlier run
# left behind.
app-us_qvarma_spec_choice: app-us_prepare_data

# app-abm_system_fit_qvarma has no prerequisite that builds its dataset. It fits
# every subdirectory of dataset/abm_system whatever wrote it, and the two writers
# cannot share that directory: applications/abm_system_convert_rdata.c produces
# EstimationSeries* from the .Rdata files under dataset/simulated, and
# applications/abm_system_simulate.c produces cop_* from the parameter design.
# Naming either one here would add its dataset on top of whichever is already
# there, and the fits would run over both without saying so. The design route
# also takes about forty hours, which is not something make should start;
# applications/abm_system_simulate_all.sh runs it.

# throughput_dataset deliberately does not appear as a prerequisite of
# app-throughput_fit: it writes 21 GB and takes tens of minutes, so
# rerunning the throughput test must not rebuild the dataset it reads.
app-abm_system_irf_loss: app-us_qvarma_spec_choice
app-abm_system_mcs: app-abm_system_irf_loss
app-abm_system_winner_irf: app-abm_system_mcs

# One throughput run per solver budget. 86.85% of the 500,000 fits at a cap of
# 2000 stopped at the cap, so what a larger budget costs and what it moves is
# measured rather than assumed, which needs both budgets on disk at once. The
# cap is in the binary's own name as well as in every file it writes, so make
# cannot hand back a binary built at one budget for a request at another.
THROUGHPUT_ITERATION_CAPS := 2000 4000 8000

# One comparison binary per pair of budgets, named for the pair, because the
# budgets it reads are compile-time constants and its report is named for them
# too. A pair is added here rather than passed on the command line so that make
# cannot serve a binary built for one pair against a request for another.
THROUGHPUT_BUDGET_PAIRS := 2000_4000 4000_8000

define scale_comparison_for_pair
.PHONY: app-throughput_iteration_budget-i$(1)
$(BIN)/throughput_iteration_budget_i$(1): applications/throughput_iteration_budget.c \
                                                    $(APPLICATION_HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -DMAT_DOUBLE -fopenmp \
	      -DBASE_ITERATIONS=$(word 1,$(subst _, ,$(1))) -DHIGH_ITERATIONS=$(word 2,$(subst _, ,$(1))) \
	      $(ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
app-throughput_iteration_budget-i$(1): $(BIN)/throughput_iteration_budget_i$(1) | $(OUT)
	./$(BIN)/throughput_iteration_budget_i$(1)
endef
$(foreach pair,$(THROUGHPUT_BUDGET_PAIRS),$(eval $(call scale_comparison_for_pair,$(pair))))

define scale_fit_target_for_cap
.PHONY: app-throughput_fit-i$(1)
$(BIN)/throughput_fit_i$(1): applications/throughput_fit.c \
                                          $(APPLICATION_HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -DMAT_DOUBLE -DMAX_ITERATIONS=$(1) -fopenmp $(ETAL_CFLAGS) $(INCLUDES) \
	      $$< -o $$@ $(ETAL_LIBS)
app-throughput_fit-i$(1): $(BIN)/throughput_fit_i$(1) | $(OUT)
	./$(BIN)/throughput_fit_i$(1)
endef
$(foreach cap,$(THROUGHPUT_ITERATION_CAPS),$(eval $(call scale_fit_target_for_cap,$(cap))))

# The DSK simulator itself, vendored under model/dsk_sfc with the three
# changes docs/DSK_MODEL_CHANGES.md records. Its own script builds it rather
# than a rule here: it is 35 translation units of someone else's C++ with its
# own vendored dependencies, and restating that here would be a second place
# for it to go wrong. model-upstream builds the unmodified reference the
# equivalence test compares against, from the copies in model/dsk_sfc/upstream.
.PHONY: model model-upstream model-sanitized
model:
	./model/dsk_sfc/build.sh

model-upstream: | $(BIN)
	./model/dsk_sfc/build.sh --upstream $(BIN)/dsk_SFC_upstream

# This project's own code under AddressSanitizer and UndefinedBehaviorSanitizer,
# for the memory-safety test: comparing outputs cannot find an index that runs
# off the end of an array and reads a plausible number.
model-sanitized: | $(BIN)
	./model/dsk_sfc/build.sh --sanitize $(BIN)/dsk_SFC_sanitized

# dsk_bulk_cancellation_distribution does not run the simulator. It calls the
# bulk cancellation header and the model's own random number generators directly,
# which are C++, so it is compiled here rather than by the rule for tests/%.c.
DSK_RANDOM_SOURCES := model/dsk_sfc/auxiliary/ran1.cpp model/dsk_sfc/auxiliary/bnldev.cpp \
                      model/dsk_sfc/auxiliary/gammln.cpp
$(BIN)/dsk_bulk_cancellation_distribution: tests/dsk_bulk_cancellation_distribution.cpp \
                                           model/dsk_sfc/dsk_sfc_bulk_cancellation.h $(DSK_RANDOM_SOURCES) | $(BIN)
	g++ -std=c++11 -O2 -I. $< $(DSK_RANDOM_SOURCES) -o $@
.PHONY: test-dsk_bulk_cancellation_distribution
test-dsk_bulk_cancellation_distribution: $(BIN)/dsk_bulk_cancellation_distribution | $(OUT)
	./$(BIN)/dsk_bulk_cancellation_distribution

# dsk_dataset_reproduction runs both builds over a sample of the design and
# compares against the archives the experiment produced. The upstream build is
# unoptimised, about 35 seconds a run, so the sample is spread over the cores
# with OpenMP, which the rule for tests/%.c does not enable.
$(BIN)/dsk_dataset_reproduction: tests/dsk_dataset_reproduction.c $(APPLICATION_HEADERS) \
                                 $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -fopenmp $(ETAL_CFLAGS) $(INCLUDES) $< -o $@ $(ETAL_LIBS)

# Every DSK test runs the simulator, so the binaries have to exist first.
test-dsk_long_path: model
test-dsk_build_equivalence: model model-upstream
test-dsk_full_output_equivalence: model model-upstream
test-dsk_design_equivalence: model model-upstream app-abm_system_design
test-dsk_memory_safety: model model-sanitized app-abm_system_design
test-dsk_ulp_sensitivity: model
test-dsk_redenomination_invariance: model
test-dsk_machine_lot_rebase: model
test-dsk_good_unit_invariance: model
test-dsk_dataset_reproduction: model model-upstream

applications: $(APPLICATION_BINARIES) | $(OUT)
	@for binary in $(APPLICATION_BINARIES); do ./$$binary || exit 1; done

bench: $(BENCH_BINARIES) | $(OUT)
	@for binary in $(BENCH_BINARIES); do ./$$binary || exit 1; done

# Ten of the eleven stems run the simulator, so every build they reach for has
# to exist before the loop starts. The per-stem targets below declare the same
# prerequisites, but an aggregate that only depended on the test binaries would
# pass on a tree where a simulator build happened to be left over and fail on a
# clean one.
TEST_PREREQUISITES := model model-upstream model-sanitized app-abm_system_design

test: $(TEST_PREREQUISITES) $(TEST_BINARIES) | $(OUT)
	@for binary in $(TEST_BINARIES); do ./$$binary || exit 1; done

test-stress: $(TEST_PREREQUISITES) $(TEST_BINARIES) | $(OUT)
	@for binary in $(TEST_BINARIES); do STRESS=1 ./$$binary || exit 1; done

# Sanitizers, per et_al.'s testing policy for allocation-heavy code. CFLAGS has
# to be a make argument rather than a shell prefix, since the assignment above
# is unconditional and would override an inherited environment variable.
# Most of the stems run the simulator, so every build they reach for has to
# exist first, the same list `test` depends on.
asan: $(TEST_PREREQUISITES) | $(BIN) $(OUT)
	@for stem in $(TEST_STEMS); do \
	  $(CC) -fsanitize=address,undefined -g -O1 -std=c11 -DMAT_DOUBLE \
	        $(ETAL_CFLAGS) $(INCLUDES) tests/$$stem.c \
	        -o $(BIN)/$${stem}_asan $(ETAL_LIBS) || exit 1; \
	  STRESS=1 ./$(BIN)/$${stem}_asan || exit 1; \
	done

clean:
	rm -rf $(BIN)
