# t-QVARMA(p,q,r), Blazsek, Escribano and Licht (2023).
#
# Header-only against et_al.; nothing here builds a library of its own.
#
#   make                          build every test binary
#   make test                     every test script
#   make test-stress              every test script, slow checks included
#   make test-qvarma_correctness  does the implementation compute what it claims
#   make test-qvarma_correctness-stress
#                                 the same plus the slow checks
#   make test-qvarma_identification
#                                 which parameters the data can pin down
#   make test-lbfgs_correctness   does the solver find the minimum it claims to
#   make test                     float64 by default. The analytic gradient is
#                                 identical in both builds; float64 is the
#                                 default because et_al.'s _syevd (used by
#                                 standard_errors) has been observed to fail
#                                 to converge under float32 on a fitted
#                                 t-QVARMA Hessian - see
#                                 ../et_al./KNOWN_ISSUES.md.
#   make MAT_DOUBLE=0 test        float32, not recommended; see above.
#   make study-qvarma_recovery_study
#                                 Monte Carlo recovery study, writes out/, prints
#                                 nothing. REPLICATIONS sets draws per cell,
#                                 default 12; MAX_ITERATIONS the solver budget.
#   make bench                    run the benchmarks (never part of test)
#   make bench-performance ETAL_DEV=1
#                                 the same, built against the development et_al.
#   make examples                 run the example scripts, output to the
#                                 terminal only, nothing written to out/
#   make applications             fit the model to the real dataset, results to
#                                 out/. Always built in float64. One app-<stem>
#                                 target per stem in APPLICATION_STEMS, which is
#                                 empty while the specification is rebuilt.
#   make app-<stem>                one application on its own, EXPERIMENT_STEMS
#                                 included - one-off structural searches
#                                 already answered (see that variable's own
#                                 comment), not part of `make applications`
#                                 so a routine build does not redo minutes of
#                                 already-settled work every time.
#   make asan                     under AddressSanitizer and UndefinedBehaviorSanitizer
#
# Adding a test script means adding its stem to TEST_STEMS. A stem maps to
# tests/qvarma_<stem>.c, bin/qvarma_<stem>, and the targets
# test-<stem> and test-<stem>-stress, which are generated below rather than
# written by hand.

CC ?= gcc
CFLAGS ?= -O2 -march=native -Wall -Wextra -std=c11
ETAL_CFLAGS := $(shell pkg-config --cflags et_al.-core)
ETAL_LIBS := $(shell pkg-config --libs et_al.-core)
# Nothing to add: et_al.'s own -I covers <et_al./...> and the bare <frame/...>
# spellings, and a "..." include next to the source resolves on its own.
INCLUDES :=

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
TEST_STEMS := qvarma_correctness
BENCH_STEMS :=
STUDY_STEMS := qvarma_recovery_study
EXAMPLE_STEMS :=
# _old/ holds the previous attempt and is deliberately not built; _old/README.md
# says why it was left.
#
# us_prepare_data has to come first: it is the only script that reads
# qvarma_data.txt and us_real.csv's Unemployment and Des_Energy_demand columns,
# and it writes out/us_system.csv and out/us_system_author_full.csv, which
# us_data.h reads back for the specification-grid scripts. The order here is
# what makes `applications`' shell loop run it first; the app-<name> targets
# that need it also depend on it directly below, for a standalone
# `make app-<name>` to regenerate the CSVs rather than read stale ones.
#
# Only data preparation is routine. A model fit already takes real time and
# does not change on a routine basis, so the fitting scripts and the whole
# ABM chain sit in EXPERIMENT_STEMS and are run one at a time via
# `make app-<name>`.
#
# not_used/ holds everything outside the chain docs/call_24082026.md
# describes, in the same directory layout, and is deliberately not built.
APPLICATION_STEMS := us_prepare_data
EXPERIMENT_STEMS := us_qvarma_employment_change \
                     abm_system_extract abm_system_fit_qvarma \
                     abm_system_mse_qvarma abm_system_mcs \
                     abm_system_winner_irf
# Whatever the application scripts share, so editing it rebuilds them.
APPLICATION_HEADERS := applications/us_data.h applications/abm_system.h
BIN := bin
OUT := out

TEST_BINARIES := $(addprefix $(BIN)/,$(TEST_STEMS))
BENCH_BINARIES := $(addprefix $(BIN)/,$(BENCH_STEMS))
STUDY_BINARIES := $(addprefix $(BIN)/,$(STUDY_STEMS))
EXAMPLE_BINARIES := $(addprefix $(BIN)/,$(EXAMPLE_STEMS))
APPLICATION_BINARIES := $(addprefix $(BIN)/,$(APPLICATION_STEMS))
EXPERIMENT_BINARIES := $(addprefix $(BIN)/,$(EXPERIMENT_STEMS))

# Benchmarks need clock_gettime, which -std=c11 hides, and they are never part
# of `test`: a function that returns the wrong answer quickly is not fast.
BENCH_CFLAGS := $(CFLAGS) -D_POSIX_C_SOURCE=199309L

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

.PHONY: all test test-stress bench bench-performance study examples applications asan clean

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
$(BIN)/$(1): tests/$(1).c $(HEADERS) $(BENCH_ETAL_HEADERS) | $(BIN)
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
$(BIN)/$(1): tests/$(1).c $(HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) -fopenmp $(ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
study-$(1): $(BIN)/$(1) | $(OUT)
	./$(BIN)/$(1)
endef
$(foreach stem,$(STUDY_STEMS),$(eval $(call study_target_for_stem,$(stem))))

study: $(STUDY_BINARIES) | $(OUT)
	@for binary in $(STUDY_BINARIES); do ./$$binary || exit 1; done

# Examples exist to be read and to print to the terminal, not to gate
# anything or to write to out/, so unlike test/study they take no -stress
# variant and no OUT dependency.
define example_target_for_stem
.PHONY: example-$(1)
$(BIN)/$(1): examples/$(1).c $(HEADERS) $(ETAL_INSTALLED_HEADERS) | $(BIN)
	$(CC) $(CFLAGS) $(ETAL_CFLAGS) $(INCLUDES) $$< -o $$@ $(ETAL_LIBS)
example-$(1): $(BIN)/$(1)
	./$(BIN)/$(1)
endef
$(foreach stem,$(EXAMPLE_STEMS),$(eval $(call example_target_for_stem,$(stem))))

examples: $(EXAMPLE_BINARIES)
	@for binary in $(EXAMPLE_BINARIES); do ./$$binary || exit 1; done

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

# us_data.h reads out/us_system.csv rather than the raw file, so a standalone
# run has to regenerate it first rather than trust whatever an earlier run
# left behind.
app-us_qvarma_employment_change: app-us_prepare_data
app-abm_system_fit_qvarma: app-abm_system_extract
app-abm_system_mse_qvarma: app-us_qvarma_employment_change
app-abm_system_mcs: app-abm_system_mse_qvarma
app-abm_system_winner_irf: app-abm_system_mcs

applications: $(APPLICATION_BINARIES) | $(OUT)
	@for binary in $(APPLICATION_BINARIES); do ./$$binary || exit 1; done

bench: $(BENCH_BINARIES) | $(OUT)
	@for binary in $(BENCH_BINARIES); do ./$$binary || exit 1; done

test: $(TEST_BINARIES) | $(OUT)
	@for binary in $(TEST_BINARIES); do ./$$binary || exit 1; done

test-stress: $(TEST_BINARIES) | $(OUT)
	@for binary in $(TEST_BINARIES); do STRESS=1 ./$$binary || exit 1; done

# Sanitizers, per et_al.'s testing policy for allocation-heavy code. CFLAGS has
# to be a make argument rather than a shell prefix, since the assignment above
# is unconditional and would override an inherited environment variable.
asan: | $(BIN) $(OUT)
	@for stem in $(TEST_STEMS); do \
	  $(CC) -fsanitize=address,undefined -g -O1 -std=c11 -DMAT_DOUBLE \
	        $(ETAL_CFLAGS) $(INCLUDES) tests/$$stem.c \
	        -o $(BIN)/$${stem}_asan $(ETAL_LIBS) || exit 1; \
	  STRESS=1 ./$(BIN)/$${stem}_asan || exit 1; \
	done

clean:
	rm -rf $(BIN)
