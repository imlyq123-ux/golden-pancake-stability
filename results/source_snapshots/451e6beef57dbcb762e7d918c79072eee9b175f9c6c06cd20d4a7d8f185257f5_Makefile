PYTHON ?= .venv/bin/python
CXX ?= g++
MPLCONFIGDIR ?= /tmp/threebody-mpl
# Official MPFR/GMP headers: system packages (libmpfr-dev, libgmp-dev) or the local copy in .deps/
# (fetched without root by `make deps`). The Taylor integrator also builds with tools/mpfr_abi.hpp.
DEPS_INC = $(if $(wildcard /usr/include/mpfr.h),,-I.deps/usr/include -I.deps/usr/include/x86_64-linux-gnu)
MPFR_LIBS = -Wl,-l:libmpfr.so.6 -Wl,-l:libgmp.so.10
export OPENBLAS_NUM_THREADS = 1

all: tools/taylor_mpfr tools/gbs_mpfr

tools/taylor_mpfr: tools/taylor_mpfr.cpp tools/mpfr_abi.hpp
	$(CXX) -std=c++17 -O3 tools/taylor_mpfr.cpp -Wl,-l:libmpfr.so.6 -o $@

tools/gbs_mpfr: tools/gbs_mpfr.cpp
	$(CXX) -std=c++17 -O3 -Wall -Wextra $(DEPS_INC) tools/gbs_mpfr.cpp $(MPFR_LIBS) -o $@

# Download the Debian/Ubuntu development headers matching the installed MPFR/GMP runtime (no root needed).
deps:
	mkdir -p .deps/debs && cd .deps/debs && apt-get download libmpfr-dev libgmp-dev
	for d in .deps/debs/*.deb; do dpkg -x $$d .deps; done

test: tools/taylor_mpfr tools/gbs_mpfr
	$(PYTHON) -m pytest tests -q --junitxml=results/tests_all.xml

# ---- Phase A (binary64) -------------------------------------------------------------------
phase-a:
	$(PYTHON) threebody_expand.py --mass-list 0.613,0.612,0.611,0.610,0.608,0.605 --stability --out results/branch_recomputed_double.json
	$(PYTHON) pancake_convergence_check.py
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) make_pancake_results.py

dense-sweep:
	$(PYTHON) dense_sweep_double.py

SEED613 = .151481541079005834,.140303019570684087,8.846608333332379,.613
SEED056 = 0.2654808237951244,0.13176080185416522,11.449696087028153,0.5603560388149639
SCAN = $(PYTHON) isola_scan_double.py --max-steps 1500 --newton-tol 2e-10 --rtol 3e-13 --atol 1e-15
# Output locations (override to run elsewhere, e.g. in tests). Parallel jobs go through
# tools/run_parallel.sh, which waits on every PID and fails the target if any job fails.
SCAN_DIR ?= results/isola_scan
REFINE_DIR ?= results/boundary_refinement
SUMMARY_DIR ?= results
FIGURE_DIR ?= figures
RUN_PARALLEL = bash tools/run_parallel.sh
comma := ,
scan_cmd = '$(SCAN) --seed $(2) --direction $(3) --label $(1) --output $(SCAN_DIR)/$(1).json > $(SCAN_DIR)/$(1).log 2>&1'

isola-scan:
	mkdir -p $(SCAN_DIR)
	$(RUN_PARALLEL) $(call scan_cmd,W_down,$(SEED613),-1) $(call scan_cmd,W_up,$(SEED613),1) \
	  $(call scan_cmd,C_down,$(SEED056),-1) $(call scan_cmd,C_up,$(SEED056),1)
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) isola_scan_summary.py --scan-dir $(SCAN_DIR) \
	  --output $(SUMMARY_DIR)/isola_scan_summary.json --figure-dir $(FIGURE_DIR)

phase-e:
	$(PYTHON) asymmetric_double.py --label main_window --mu-center 0.613 --uvT 0.151481541079005834,0.140303019570684087,8.846608333332379 --deltas 0:0.012:0.0005 --mu-range 0.594,0.624 --mu-step 0.0005 --workers 8 > results/asymmetric_main_window.log 2>&1
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) asymmetric_summary.py

figures:
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) summarize_transitions.py --refinement-dir $(REFINE_DIR) --output-dir $(SUMMARY_DIR)
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) make_figures.py

# ---- Phases B/C, first session (MPFR Taylor + Decimal quotient) -----------------------------
phase-bc: tools/taylor_mpfr
	$(PYTHON) high_precision.py
	$(PYTHON) high_precision.py --digits 85 --order 50 --correction-tolerance 1e-54 --verification-tolerance 1e-60 --output results/high_precision_orbits_85.json
	$(PYTHON) high_precision.py --masses 0.59904378,0.59904379,0.61784426,0.61784428 --seed-file results/high_precision_orbits.json --output results/boundary_orbits_65.json
	$(PYTHON) high_precision.py --masses 0.59904378,0.59904379,0.61784426,0.61784428 --digits 85 --order 50 --target 1e-48 --correction-tolerance 1e-54 --verification-tolerance 1e-60 --output results/boundary_orbits_85.json
	$(PYTHON) verify_time_reversal.py
	MPLCONFIGDIR=$(MPLCONFIGDIR) $(PYTHON) finalize_results.py

# ---- Second session: independent GBS + mpmath pipeline, refinements, classification -------
MASSES15 = 0.613,0.5990,0.5991,0.6178,0.6179,0.597,0.605,0.618,0.59904378,0.59904379,0.61784426,0.61784428,0.5876435519022455,0.612356088827994835162485140041506525569967,0.630314929124448

independent: tools/gbs_mpfr
	$(PYTHON) independent_mp.py --masses $(MASSES15) --digits 80 --stages 20 --verify-stages 24 --newton-tol 1e-62 --verify-tol 1e-66 --target 1e-55 --workers 8 --output results/independent_gbs_orbits.json --raw-dir results/gbs_raw
	$(PYTHON) crosscheck_report.py
	$(PYTHON) half_period_symmetry.py

refine_cmd = '$(PYTHON) refine_boundaries.py --boundary $(1) --pipeline $(2) --output-dir $(REFINE_DIR) $(3) > $(REFINE_DIR)/$(1)_$(2).log 2>&1'
REFINE_MAIN = $(foreach b,lower upper,$(foreach p,taylor50 taylor70 taylor90 gbs80,$(call refine_cmd,$(b),$(p),)))
REFINE_NEW = $(foreach p,gbs80 taylor50,\
  $(call refine_cmd,K0_lowwindow_top,$(p),--bracket 0.5800$(comma)0.5805 --kind krein) \
  $(call refine_cmd,K3_nearfold,$(p),--bracket 0.63417975$(comma)0.63435086 --kind krein) \
  $(call refine_cmd,T1_abab_postfold,$(p),--bracket 0.4756077$(comma)0.4763384 --kind tangent --seed-scan $(SCAN_DIR)/C_down.json --seed-arc-min 0.2))

refine: tools/taylor_mpfr tools/gbs_mpfr
	mkdir -p $(REFINE_DIR)
	$(RUN_PARALLEL) $(REFINE_MAIN)
	$(RUN_PARALLEL) $(REFINE_NEW)
	$(PYTHON) summarize_transitions.py --refinement-dir $(REFINE_DIR) --output-dir $(SUMMARY_DIR)

transitions: tools/gbs_mpfr
	$(PYTHON) transition_analysis.py --workers 8

.PHONY: all deps test phase-a dense-sweep isola-scan phase-e figures phase-bc independent refine transitions
