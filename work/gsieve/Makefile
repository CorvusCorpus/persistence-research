# Chain sieve builds; see CHAIN_SIEVE.md. Stage-2 variants of the main path (gsieve.c):
#   ifma    AVX-512 F/DQ/BW/VL + IFMA (the production build of 19-20 Oct 2026)
#   avx512  AVX-512 F/DQ/BW/VL without IFMA
#   avx2    AVX2
#   scalar  portable reference
# `make` builds the default names (gsieve_k13 ...) with the first variant the target CPU supports.
# `make variants` builds every variant this CPU can run, with suffixes (gsieve_k13_avx2 ...).
# `make test` validates the default build; `make test V=avx2` validates the suffixed build of a variant.
# `make test-prod` re-runs a logged production invocation (main N=36000, c in [7540,7898), 8.6e11 points) and
# requires every counter to equal the production DONE line (about 2 min on 4 cores).
CC     ?= gcc
CFLAGS ?= -O3 -march=native -Wall -Wno-unused-function
LIBS    = -lgmp -lm -lpthread
BINS    = gsieve_k13 gsieve_k12 gsieve_k13_nl7 gsieve_k12_nl7
VARIANTS = ifma avx512 avx2 scalar
VFLAG_ifma = -DS2_IFMA
VFLAG_avx512 = -DS2_AVX512
VFLAG_avx2 = -DS2_AVX2
VFLAG_scalar = -DSCALAR
CPUFLAGS := $(shell grep -m1 '^flags' /proc/cpuinfo 2>/dev/null)
RUNNABLE := $(if $(filter avx512ifma,$(CPUFLAGS)),ifma) $(if $(filter avx512f,$(CPUFLAGS)),avx512) $(if $(filter avx2,$(CPUFLAGS)),avx2) scalar
V      ?=
SFX     = $(if $(V),_$(V))

all: $(BINS)

gsieve_k13:     gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=13 -o $@ gsieve.c $(LIBS)
gsieve_k12:     gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=12 -o $@ gsieve.c $(LIBS)
gsieve_k13_nl7: gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=13 -DNL=7 -o $@ gsieve.c $(LIBS)
gsieve_k12_nl7: gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=12 -DNL=7 -o $@ gsieve.c $(LIBS)

define VRULES
gsieve_k13_$(1):     gsieve.c gconst.h ; $$(CC) $$(CFLAGS) $$(VFLAG_$(1)) -DKK=13 -o $$@ gsieve.c $$(LIBS)
gsieve_k12_$(1):     gsieve.c gconst.h ; $$(CC) $$(CFLAGS) $$(VFLAG_$(1)) -DKK=12 -o $$@ gsieve.c $$(LIBS)
gsieve_k13_nl7_$(1): gsieve.c gconst.h ; $$(CC) $$(CFLAGS) $$(VFLAG_$(1)) -DKK=13 -DNL=7 -o $$@ gsieve.c $$(LIBS)
gsieve_k12_nl7_$(1): gsieve.c gconst.h ; $$(CC) $$(CFLAGS) $$(VFLAG_$(1)) -DKK=12 -DNL=7 -o $$@ gsieve.c $$(LIBS)
endef
$(foreach v,$(VARIANTS),$(eval $(call VRULES,$(v))))

variants: $(foreach v,$(RUNNABLE),$(BINS:%=%_$(v)))

gconst.h: gconst.py
	python3 gconst.py

val_kernel: val_kernel.c ; $(CC) -O2 -Wall -o $@ val_kernel.c

# validation suite (each prints OK / exact match; any assertion failure is a bug)
test: $(BINS:%=%$(SFX)) val_kernel
	./val_kernel
	python3 tiny.py 131
	python3 val_dump.py 13 G 150 ./gsieve_k13$(SFX)
	python3 val_dump.py 13 E 150 ./gsieve_k13$(SFX)
	python3 val_dump.py 12 E 150 ./gsieve_k12$(SFX)
	python3 val_dump.py 12 G 150 ./gsieve_k12$(SFX)
	python3 val_e2e.py ./gsieve_k13_nl7$(SFX) 160 7
	python3 val_e2e.py ./gsieve_k12_nl7$(SFX) 160 7
	python3 val_direct_odd.py ./gsieve_k13_nl7$(SFX) 13 160 7
	python3 val_prodrow.py 36000 35500 ./gsieve_k13$(SFX)

NPROC := $(shell nproc 2>/dev/null || echo 1)
test-prod: gsieve_k13$(SFX)
	python3 val_proddone.py ./gsieve_k13$(SFX) $(NPROC) --only main:7540:7898 logs/prod_N36000_n0-1_main_0_7898.log

clean:
	rm -f $(BINS) $(foreach v,$(VARIANTS),$(BINS:%=%_$(v))) val_kernel

.PHONY: all variants test test-prod clean
