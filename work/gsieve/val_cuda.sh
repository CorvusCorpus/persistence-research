#!/bin/bash
# Validation of the CUDA variant (stages 1-2 of the main path on the GPU) on this machine: run it after building, on
# the exact binaries used for production. One line per check; logs go to LOGDIR. About 10 min (RTX 3060 Ti, 20 CPUs).
#   K: kernel tests: the GPU's stage-1/2 decisions on the whole table for many lift indices j against the SCALAR
#      formula, and on 2^20 random (e, j) against 2^E mod 10^2K (GMP); K = 13 and 12, modes G and E
#   A: the validation suite (make test V=cuda) at 1 and 8 threads (E-mode end to end must find exactly 314)
#   C: production slices against exact arithmetic (val_prodrow.py), N = 36000 and the shell 36000 < m <= 40000
#   D: cross-variant dumps (val_xvariant.py): the CUDA binary at 1 and 8 threads against the CPU variants
#   P: a logged production invocation (val_proddone.py, main N = 36000, c in [7540, 7898), 8.6e11 points)
#   Z: the zero-free paths: stages 3-4 on the GPU (GSIEVE_GPU34=1: val_dump, dumps byte-identical to avx2, val_e2e,
#      val_prodrow, val_proddone) and odd2 (GSIEVE_ODD2=1: val_odd2.py against exact arithmetic)
# usage: val_cuda.sh LOGDIR      (from work/gsieve; exit status 0 only if every check passes)
set -u
L=$(mkdir -p "$1" && cd "$1" && pwd) || exit 1
cd "$(dirname "$0")" || exit 1
st=0
step() {          # step LABEL LOGFILE CMD...: run CMD with output to LOGFILE, print LABEL, exit code and time
    local label=$1 log=$2; shift 2
    local t0=$(date +%s); "$@" > "$log" 2>&1; local rc=$?
    [ $rc -eq 0 ] || st=1
    echo "$label exit=$rc $(( $(date +%s) - t0 ))s"
}
echo "== start $(date -u) commit $(git rev-parse --short HEAD) gsieve.c $(sha256sum gsieve.c | cut -c1-16)" \
     "gsieve_cuda.cu $(sha256sum gsieve_cuda.cu | cut -c1-16)"
make cuda all variants val_kernel > "$L/build.log" 2>&1 || { echo "build failed, see $L/build.log"; exit 1; }
sha256sum gsieve_k1* > "$L/binaries.sha256"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader > "$L/gpu.txt" 2>&1
F=" $(grep -m1 '^flags' /proc/cpuinfo) "; VARS=""          # CPU variants this machine runs (RUNNABLE in the Makefile)
for v in ifma:avx512ifma avx512:avx512f avx2:avx2; do case "$F" in *" ${v#*:} "*) VARS="$VARS ${v%%:*}";; esac; done
VARS="$VARS scalar"; echo "gpu: $(cat "$L/gpu.txt"); CPU variants:$VARS"
# K
for spec in "13 G 1024 1" "13 E 16 2" "12 G 256 3" "12 E 16 4"; do
    set -- $spec
    step "K kerneltest K=$1 $2 NJ=$3" "$L/K_kerneltest_k$1_$2.log" ./gsieve_k$1_cuda kerneltest $2 $3 $4 threads=$(nproc)
    echo "  $(grep -h 'mismatches' "$L/K_kerneltest_k$1_$2.log" | sed 's/.*mismatches/mismatches/' | tr '\n' ' ')"
done
# A
for T in 1 8; do step "A suite V=cuda T=$T" "$L/A_suite_cuda_T$T.log" env GSIEVE_THREADS=$T make test V=cuda; done
# C
for spec in "36000 34200" "36000 33000" "36000 30000" "36000 25000" "40000 39500 36000" "40000 38200 36000" \
            "40000 35990 36000" "40000 30000 36000"; do
    set -- $spec
    step "C val_prodrow N=$1 c=$2 n0=${3:--1}" "$L/C_prodrow_$1_$2.log" env GSIEVE_THREADS=8 python3 val_prodrow.py $1 $2 ./gsieve_k13_cuda ${3:+$3}
    echo "  $(grep -E '^c=' "$L/C_prodrow_$1_$2.log")"
done
# D
X="./gsieve_k13_cuda:1 ./gsieve_k13_cuda:8"; for V in $VARS; do X="$X ./gsieve_k13_$V:8"; done
for r in "main G 1000 0 1001" "main E 400 0 401" "main G 36000 20000 20001" "main G 40000 35990 36000 n0=36000" \
         "main G 36000 7540 7541" "main G 40000 11824 11826 n0=36000" "direct G 3000 0 3001" "direct E 400 0 401" \
         "odd G 3000 0 3001" "odd E 400 0 401"; do
    step "D xvariant [$r]" "$L/D_xvariant_$(echo $r | tr ' =' '__').log" python3 val_xvariant.py "$r" $X
done
# P
step "P val_proddone main:7540:7898" "$L/P_proddone.log" python3 val_proddone.py ./gsieve_k13_cuda $(nproc) \
     --only main:7540:7898 logs/prod_N36000_n0-1_main_0_7898.log
echo "  $(grep -o 'variant=.*' "$L/P_proddone.log")"
# Z: the zero-free extension's paths: stages 3-4 on the GPU (GSIEVE_GPU34=1) and odd2 (GSIEVE_ODD2=1)
same() {          # same "ARGS": avx2 and cuda34 dumps (all non-DONE lines) must be byte-identical
    local t; t=$(mktemp -d)
    ./gsieve_k13 $1 dump prog threads=8 > "$t/a" 2>/dev/null; GSIEVE_GPU34=1 ./gsieve_k13_cuda $1 dump prog threads=8 > "$t/b" 2>/dev/null
    local n; n=$(grep -vc '^DONE' "$t/a")
    if [ "$n" -gt 0 ] && cmp -s <(grep -v '^DONE' "$t/a") <(grep -v '^DONE' "$t/b"); then echo "identical: $n lines"; rm -rf "$t"; return 0; fi
    echo "DIFFERENT"; rm -rf "$t"; return 1
}
for spec in "13 G" "13 E" "12 E" "12 G"; do
    set -- $spec
    step "Z val_dump K=$1 $2 N=150 (GPU34)" "$L/Z_dump_k$1_$2.log" env GSIEVE_GPU34=1 python3 val_dump.py $1 $2 150 ./gsieve_k$1_cuda
done
for r in "main G 1000 0 1001" "main E 1000 0 1001" "main G 40000 35990 36000 n0=36000"; do
    step "Z same dump avx2 = cuda34 [$r]" "$L/Z_same_$(echo $r | tr ' =' '__').log" same "$r"
    echo "  $(cat "$L/Z_same_$(echo $r | tr ' =' '__').log")"
done
for k in 13 12; do
    step "Z val_e2e K=$k NL=7 (GPU34 + ODD2)" "$L/Z_e2e_k$k.log" env GSIEVE_GPU34=1 GSIEVE_ODD2=1 python3 val_e2e.py ./gsieve_k${k}_nl7_cuda 160 7
    echo "  $(grep -E '^expected' "$L/Z_e2e_k$k.log")"
done
step "Z val_prodrow N=36000 c=33000 (GPU34)" "$L/Z_prodrow_36000_33000.log" env GSIEVE_GPU34=1 GSIEVE_THREADS=8 python3 val_prodrow.py 36000 33000 ./gsieve_k13_cuda
step "Z val_proddone main:7540:7898 (GPU34)" "$L/Z_proddone.log" env GSIEVE_GPU34=1 python3 val_proddone.py ./gsieve_k13_cuda $(nproc) \
     --only main:7540:7898 logs/prod_N36000_n0-1_main_0_7898.log
for spec in "400 ./gsieve_k13_cuda" "500 ./gsieve_k13_cuda" "400 ./gsieve_k13_cuda 250" "300 ./gsieve_k13_nl7_cuda -1 7"; do
    step "Z val_odd2 $spec" "$L/Z_odd2_$(echo $spec | tr ' ./' '___').log" python3 val_odd2.py $spec
done
echo "== end $(date -u): $([ $st -eq 0 ] && echo 'ALL CHECKS PASSED' || echo 'SOME CHECKS FAILED')"
exit $st
