# Chain sieve builds. Requires AVX-512 (F, DQ, BW, VL, IFMA, VBMI2) and BMI2; see CHAIN_SIEVE.md.
CC     ?= gcc
CFLAGS ?= -O3 -march=native -Wall -Wno-unused-function
LIBS    = -lgmp -lm
BINS    = gsieve_k13 gsieve_k12 gsieve_k13_nl7 gsieve_k12_nl7

all: $(BINS)

gsieve_k13:     gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=13 -o $@ gsieve.c $(LIBS)
gsieve_k12:     gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=12 -o $@ gsieve.c $(LIBS)
gsieve_k13_nl7: gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=13 -DNL=7 -o $@ gsieve.c $(LIBS)
gsieve_k12_nl7: gsieve.c gconst.h ; $(CC) $(CFLAGS) -DKK=12 -DNL=7 -o $@ gsieve.c $(LIBS)

gconst.h: gconst.py
	python3 gconst.py

# validation suite (each prints OK / exact match; any assertion failure is a bug)
test: all
	python3 tiny.py 131
	python3 val_dump.py 13 G 150
	python3 val_dump.py 13 E 150
	python3 val_dump.py 12 E 150
	python3 val_dump.py 12 G 150
	python3 val_e2e.py ./gsieve_k13_nl7 160 7
	python3 val_e2e.py ./gsieve_k12_nl7 160 7
	python3 val_direct_odd.py ./gsieve_k13_nl7 13 160 7
	python3 val_prodrow.py 36000 35500

clean:
	rm -f $(BINS)

.PHONY: all test clean
