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

# CUDA variant: stages 2-4 on an NVIDIA GPU (z2sieve_cuda.cu), the same decisions as the CPU build. Not part of
# `make all` (needs nvcc). CUDA_ARCH=sm_86 is the RTX 3060 Ti.
NVCC      ?= nvcc
CUDA_ARCH ?= sm_86
NVFLAGS   ?= -O3 -arch=$(CUDA_ARCH)
CUDA_HOME ?= $(patsubst %/bin/nvcc,%,$(shell command -v $(NVCC) 2>/dev/null))
CUDALIBS   = -L$(CUDA_HOME)/lib64 -lcudart_static -ldl -lrt -lstdc++
K4 = -DKK=4 -DNL=2 -DTB=18 -DTD=3
z2sieve_cuda.o: z2sieve_cuda.cu z2sieve_cuda.h
	$(NVCC) $(NVFLAGS) -c -o $@ z2sieve_cuda.cu
z2sieve_k4_cuda.o: z2sieve_cuda.cu z2sieve_cuda.h
	$(NVCC) $(NVFLAGS) $(K4) -c -o $@ z2sieve_cuda.cu
z2sieve_cuda: z2sieve.c z2sieve_cuda.h z2sieve_cuda.o
	$(CC) $(CFLAGS) -DZ2_CUDA -o $@ z2sieve.c z2sieve_cuda.o $(LIBS) $(CUDALIBS)
z2sieve_k4_cuda: z2sieve.c z2sieve_cuda.h z2sieve_k4_cuda.o
	$(CC) $(CFLAGS) -DZ2_CUDA $(K4) -DNMIN=60 -o $@ z2sieve.c z2sieve_k4_cuda.o $(LIBS) $(CUDALIBS)
cuda: z2sieve_cuda z2sieve_k4_cuda

test: all
	./val_z2.sh

clean:
	rm -f z2sieve z2sieve_k4 z2sieve_cuda z2sieve_k4_cuda z2sieve_cuda.o z2sieve_k4_cuda.o

.PHONY: all cuda test clean
