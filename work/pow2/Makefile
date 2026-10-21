# z2sieve: decimal zeros of powers of two (OEIS A007377), see docs/ZEROLESS_POW2.md
CC ?= gcc
CFLAGS = -O3 -march=x86-64-v2 -std=gnu11 -Wall -Wextra
LIBS = -lmpfr -lgmp -lpthread -lm

all: z2sieve z2sieve_k4

# production: K = 13 (stages 1-2: the last 26 digits), NL = 18 (stage 4: the last 162 digits), TB = 24, TD = 7
z2sieve: z2sieve.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

# positive control: K = 4, NL = 2 (the last 18 digits), TD = 3, NMIN = 60, so that the known terms 67..86 are in range
z2sieve_k4: z2sieve.c
	$(CC) $(CFLAGS) -DKK=4 -DNL=2 -DTB=18 -DTD=3 -DNMIN=60 -o $@ $< $(LIBS)

test: all
	./val_z2.sh

clean:
	rm -f z2sieve z2sieve_k4

.PHONY: all test clean
