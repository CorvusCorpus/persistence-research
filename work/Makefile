# Build all compiled helpers. Run from work/.
CC      ?= gcc
CFLAGS  ?= -O2 -march=native
LIBS     = -lmpfr -lgmp -lm

HELPERS = csrc/topc csrc/topc3 csrc/bottomc csrc/bottomc2 csrc/bottomc3 csrc/topb csrc/bottomb csrc/top csrc/bottom

all: $(HELPERS) zf7 gsieve

# everything except the chain sieve (which needs AVX-512 IFMA unless built with -DSCALAR)
helpers: $(HELPERS) zf7

csrc/topc:     csrc/top.c     ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/topc3:    csrc/top3.c    ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/bottomc:  csrc/bottom.c  ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/bottomc2: csrc/bottom2.c ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/bottomc3: csrc/bottom3.c ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/topb:     csrc/topb.c    ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/bottomb:  csrc/bottomb.c ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
# legacy gap-chain helpers (prove_iter3.py); the original binaries were built from these sources
csrc/top:      csrc/top_v1.c    ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)
csrc/bottom:   csrc/bottom_v1.c ; $(CC) $(CFLAGS) -o $@ $< $(LIBS)

zf7: zf7.c ; $(CC) -O3 -march=native -fopenmp -o $@ $< -lm

gsieve:
	$(MAKE) -C gsieve

clean:
	rm -f $(HELPERS) zf7
	$(MAKE) -C gsieve clean

.PHONY: all helpers gsieve clean
