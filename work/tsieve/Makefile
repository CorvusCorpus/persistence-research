# Tree sieve for the chain sieve's main path (see tsieve6.c, docs/NOTEBOOK.md section T).
#   make            tsieve6 (CPU), tsieve6_nl7 (positive control), tsieve_gpu_emu (the GPU kernel compiled for the CPU)
#   make cuda       tsieve_gpu, tsieve_gpu_nl7 (needs nvcc; CUDA_ARCH defaults to sm_86, the RTX 3060 Ti)
#   make test       validation that runs anywhere (CPU and emulated GPU kernel against gsieve, 314/314)
#   make test-gpu   ./val_tsieve_gpu.sh: the GPU acceptance test (desktop)
CC     ?= cc
CFLAGS ?= -O3 -march=native -Wall -Wno-unused-function -Wno-unused-variable -Werror=implicit-function-declaration
NVCC   ?= nvcc
CUDA_ARCH ?= sm_86
GS     = ../gsieve
LIBS   = -lgmp -lm -lpthread
DEPS   = tsieve6.c $(GS)/gsieve.c $(GS)/gconst.h

all: tp13.txt tsieve6 tsieve6_nl7 tsieve_gpu_emu tsieve_gpu_emu_nl7
cuda: tp13.txt tsieve_gpu tsieve_gpu_nl7

tp13.txt: tparams.py ; python3 tparams.py 13 $@
tsieve6: $(DEPS) ; $(CC) $(CFLAGS) -DKK=13 -I$(GS) -o $@ tsieve6.c $(LIBS)
tsieve6_nl7: $(DEPS) ; $(CC) $(CFLAGS) -DKK=13 -DNL=7 -I$(GS) -o $@ tsieve6.c $(LIBS)

tg_host.o: tsieve_gpu.c tsieve_gpu.h $(DEPS) ; $(CC) $(CFLAGS) -DKK=13 -I$(GS) -c -o $@ tsieve_gpu.c
tg_host_nl7.o: tsieve_gpu.c tsieve_gpu.h $(DEPS) ; $(CC) $(CFLAGS) -DKK=13 -DNL=7 -I$(GS) -c -o $@ tsieve_gpu.c
tg_emu.o: tsieve_gpu.cu tsieve_gpu.h ; g++ -O2 -march=native -DTG_EMU -x c++ -c -o $@ tsieve_gpu.cu
tg_cuda.o: tsieve_gpu.cu tsieve_gpu.h ; $(NVCC) -O3 -arch=$(CUDA_ARCH) -c -o $@ tsieve_gpu.cu
tsieve_gpu_emu: tg_host.o tg_emu.o ; g++ -o $@ tg_host.o tg_emu.o $(LIBS)
tsieve_gpu_emu_nl7: tg_host_nl7.o tg_emu.o ; g++ -o $@ tg_host_nl7.o tg_emu.o $(LIBS)
tsieve_gpu: tg_host.o tg_cuda.o ; $(NVCC) -o $@ tg_host.o tg_cuda.o $(LIBS)
tsieve_gpu_nl7: tg_host_nl7.o tg_cuda.o ; $(NVCC) -o $@ tg_host_nl7.o tg_cuda.o $(LIBS)

test: all
	./val_tsieve6.sh tp13.txt G 4 300 -1 4
	./val_tsieve6.sh tp13.txt G 3 900 600 4
	./val_tsieve6.sh tp13.txt E 5 600 -1 4
	python3 val_e2e_tree.py ./tsieve6_nl7 $(GS)/gsieve_k13_nl7 tp13.txt 4 160 7 4
	./val_tsieve6.sh tp13.txt G 4 300 -1 4 $(GS)/gsieve_k13 ./tsieve_gpu_emu
	python3 val_e2e_tree.py ./tsieve_gpu_emu_nl7 $(GS)/gsieve_k13_nl7 tp13.txt 4 160 7 4
	@echo "TSIEVE TESTS PASSED"
test-gpu: cuda tsieve6
	./val_tsieve_gpu.sh

clean:
	rm -f tsieve tsieve_check tsieve4 tsieve6 tsieve6_nl7 tsieve_gpu tsieve_gpu_nl7 tsieve_gpu_emu tsieve_gpu_emu_nl7 *.o t4_* t6_*

.PHONY: all cuda test test-gpu clean
