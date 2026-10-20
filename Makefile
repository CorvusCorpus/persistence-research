# Top-level entry points. The scripts themselves run from work/ (base 10, small bases) or
# work/gsieve/ (chain sieve) because they use relative paths; see CLAUDE.md.
#
#   make                build everything (CPU helpers, zf7, chain sieve)
#   make test           chain-sieve validation suite + base-10 smoke test + record checks
#   make gsieve-test    chain-sieve validation suite only (work/gsieve/Makefile `test`)
#   make smoke          base-10 smoke test only (scripts/smoke_base10.py: prove_cum.py 98 147 from scratch)
#   make records-check  SHA-256 manifest of the record files + coverage of the production logs
#   make manifest       regenerate MANIFEST.sha256 (only after adding new records; review the diff)
#   make logs-tarball   release tarball of all run logs in dist/ (for a GitHub release)
#
# Chain-sieve build variant. The production main path (gsieve.c default) uses AVX-512 IFMA; on CPUs
# without avx512ifma the only buildable variant is the scalar reference path (-DSCALAR). Override with
#   make GSIEVE_VARIANT=-DSCALAR   or   make GSIEVE_VARIANT=
HAVE_IFMA      := $(shell grep -qw avx512ifma /proc/cpuinfo 2>/dev/null && echo 1)
GSIEVE_VARIANT ?= $(if $(HAVE_IFMA),,-DSCALAR)
GSIEVE_CFLAGS  ?= -O3 -march=native -Wall -Wno-unused-function $(GSIEVE_VARIANT)
PYTHON         ?= python3

all: helpers gsieve

helpers:
	$(MAKE) -C work helpers

gsieve:
	@echo "gsieve build: CFLAGS='$(GSIEVE_CFLAGS)'"
	$(MAKE) -C work/gsieve CFLAGS="$(GSIEVE_CFLAGS)"

gsieve-test: gsieve
	$(MAKE) -C work/gsieve test CFLAGS="$(GSIEVE_CFLAGS)"

smoke: helpers
	$(PYTHON) scripts/smoke_base10.py 98 147

records-check:
	sha256sum --quiet -c MANIFEST.sha256 && echo "MANIFEST.sha256: all record files intact"
	$(PYTHON) scripts/check_coverage.py

test: gsieve-test smoke records-check
	@echo "ALL TESTS PASSED"

manifest:
	scripts/make_manifest.sh > MANIFEST.sha256
	@git diff --stat -- MANIFEST.sha256 || true

logs-tarball:
	scripts/logs_tarball.sh

clean:
	$(MAKE) -C work clean

.PHONY: all helpers gsieve gsieve-test smoke records-check test manifest logs-tarball clean
