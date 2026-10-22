# Build all compiled helpers. Run from work/.
CC      ?= gcc
CFLAGS  ?= -O2 -march=native
LIBS     = -lmpfr -lgmp -lm

HELPERS = csrc/topc csrc/topc3 csrc/bottomc csrc/bottomc2 csrc/bottomc3 csrc/topb csrc/bottomb csrc/top csrc/bottom

all: $(HELPERS) zf7 gsieve

# everything except the chain sieve
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

# GPU bottom solver (drop-in for csrc/bottomc3 in 'cum' mode; see csrc/bottomgpu.cu). `make cuda` needs nvcc.
NVCC      ?= nvcc
CUDA_ARCH ?= sm_86
BGDEPS = csrc/bottomgpu_host.c csrc/bottomgpu.cu csrc/bottomgpu.h csrc/bigfix.h csrc/bottom3.c
csrc/bottomgpu_host.o: $(BGDEPS) ; $(CC) $(CFLAGS) -Wno-unused-function -c -o $@ csrc/bottomgpu_host.c
csrc/bottomgpu_cuda.o: $(BGDEPS) ; $(NVCC) -O3 -arch=$(CUDA_ARCH) -c -o $@ csrc/bottomgpu.cu
csrc/bottomgpu_emu.o: $(BGDEPS) ; g++ -O2 -march=native -DBG_EMU -x c++ -c -o $@ csrc/bottomgpu.cu
csrc/bottomgpu: csrc/bottomgpu_host.o csrc/bottomgpu_cuda.o ; $(NVCC) -o $@ $^ $(LIBS)
csrc/bottomgpu_emu: csrc/bottomgpu_host.o csrc/bottomgpu_emu.o ; g++ -o $@ $^ $(LIBS)
# GPU top filter (drop-in for csrc/topc3; see csrc/topgpu.cu)
TPDEPS = csrc/topgpu_host.c csrc/topgpu.cu csrc/topgpu.h csrc/bigfix.h
csrc/topgpu_host.o: $(TPDEPS) ; $(CC) $(CFLAGS) -c -o $@ csrc/topgpu_host.c
csrc/topgpu_cuda.o: $(TPDEPS) ; $(NVCC) -O3 -arch=$(CUDA_ARCH) -c -o $@ csrc/topgpu.cu
csrc/topgpu_emu.o: $(TPDEPS) ; g++ -O2 -march=native -DTP_EMU -x c++ -c -o $@ csrc/topgpu.cu
csrc/topgpu: csrc/topgpu_host.o csrc/topgpu_cuda.o ; $(NVCC) -o $@ $^ $(LIBS)
csrc/topgpu_emu: csrc/topgpu_host.o csrc/topgpu_emu.o ; g++ -o $@ $^ $(LIBS)
cuda: csrc/bottomgpu csrc/bottomgpu_emu csrc/topgpu csrc/topgpu_emu

gsieve:
	$(MAKE) -C gsieve

clean:
	rm -f $(HELPERS) zf7 csrc/bottomgpu csrc/bottomgpu_emu csrc/topgpu csrc/topgpu_emu csrc/*.o
	$(MAKE) -C gsieve clean

.PHONY: all helpers gsieve clean cuda
