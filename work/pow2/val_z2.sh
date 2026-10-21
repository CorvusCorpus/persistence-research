#!/bin/bash
# Validation of z2sieve on this machine (docs/ZEROLESS_POW2.md, "Validation"). One line per check; logs go to LOGDIR.
#   K: kernel tests: every stage against GMP/MPFR at 10^6 random n < 2^56 (both builds)
#   E: exact region checks (val_z2.py): every decision of the dump against Python integers and ball arithmetic, the
#      S2 set in both directions, CDONE coverage and the DONE counters; regions at the start of the range, across the
#      first j boundary, at chunk boundaries (10^12 = 1024 P, 10^15 = 1024000 P), at 3.7 10^13 and just below 2^61
#   C: positive control: the K = 4 build over [60, 2 10^6) must report exactly the known terms 67 72 76 77 81 86
#   T: the dumps do not depend on the number of threads
#   G: the CUDA builds (if nvcc is available): the exact region checks and the positive control on the GPU
#   X: CPU against GPU: sorted dumps identical on regions, and identical CDONE counters on a production-scale slice
# usage: val_z2.sh [LOGDIR]      (from work/pow2; exit status 0 only if every check passes)
set -u
L=$(mkdir -p "${1:-val_logs}" && cd "${1:-val_logs}" && pwd) || exit 1
cd "$(dirname "$0")" || exit 1
st=0
step() {          # step LABEL LOGFILE CMD...: run CMD with output to LOGFILE, print LABEL, exit code and time
    local label=$1 log=$2; shift 2
    local t0=$(date +%s); "$@" > "$log" 2>&1; local rc=$?
    [ $rc -eq 0 ] || st=1
    echo "$label exit=$rc $(( $(date +%s) - t0 ))s"
}
echo "== start $(date -u) commit $(git rev-parse --short HEAD) z2sieve.c $(sha256sum z2sieve.c | cut -c1-16)" \
     "z2sieve_cuda.cu $(sha256sum z2sieve_cuda.cu | cut -c1-16)"
make all > "$L/build.log" 2>&1 || { echo "build failed, see $L/build.log"; exit 1; }
CUDA=0; command -v nvcc > /dev/null && make cuda >> "$L/build.log" 2>&1 && CUDA=1
sha256sum z2sieve z2sieve_k4 $([ $CUDA = 1 ] && echo z2sieve_cuda z2sieve_k4_cuda) > "$L/binaries.sha256"
[ $CUDA = 1 ] && nvidia-smi --query-gpu=name,driver_version --format=csv,noheader > "$L/gpu.txt" 2>&1
T=${THREADS:-8}
# K
step "K kerneltest K=13" "$L/K_k13.log" ./z2sieve kerneltest 1000000 11
echo "  $(grep -h '^kerneltest' "$L/K_k13.log")"
step "K kerneltest K=4" "$L/K_k4.log" ./z2sieve_k4 kerneltest 1000000 12
echo "  $(grep -h '^kerneltest' "$L/K_k4.log")"
# E
for r in "1000 2001000" "975562500 977562500" "1000000000000 1000002000000" "999999999000007 1000000001000007" \
         "37000000000000 37000002000000" "2305843009211693952 2305843009213693952"; do
    step "E val_z2 K=13 [$r]" "$L/E_k13_$(echo $r | tr ' ' '_').log" python3 val_z2.py ./z2sieve $r $T
    echo "  $(grep -h '^val_z2' "$L/E_k13_$(echo $r | tr ' ' '_').log" | sed 's/.*K=13 NL=18 TB=24 TD=7: //')"
done
for r in "1000000000000 1000001000000" "2305843009212693952 2305843009213693952"; do
    step "E val_z2 K=4 [$r]" "$L/E_k4_$(echo $r | tr ' ' '_').log" python3 val_z2.py ./z2sieve_k4 $r $T
    echo "  $(grep -h '^val_z2' "$L/E_k4_$(echo $r | tr ' ' '_').log" | sed 's/.*TD=3: //')"
done
# C
step "C positive control K=4 [60, 2000000)" "$L/C_k4_60_2000000.log" python3 val_z2.py ./z2sieve_k4 60 2000000 $T
got=$(grep -h '^val_z2' "$L/C_k4_60_2000000.log" | grep -o 'members=\[[^]]*\]')
echo "  $got"
[ "$got" = "members=[67, 72, 76, 77, 81, 86]" ] || { echo "  POSITIVE CONTROL FAILED"; st=1; }
# T
same() {          # same EXE N_LO N_HI: dumps at 1 and $T threads identical (all lines but DONE)
    local a b; a=$("$1" "$2" "$3" dump threads=1 2>/dev/null | grep -v '^DONE' | sha256sum); b=$("$1" "$2" "$3" dump threads=$T 2>/dev/null | grep -v '^DONE' | sha256sum)
    [ -n "$a" ] && [ "$a" = "$b" ] && { echo "identical: $a"; return 0; }; echo "DIFFERENT"; return 1
}
step "T threads 1 = $T [K=13 999999999000007 1000000001000007]" "$L/T_k13.log" same ./z2sieve 999999999000007 1000000001000007
step "T threads 1 = $T [K=4 60 2000000]" "$L/T_k4.log" same ./z2sieve_k4 60 2000000
if [ $CUDA = 1 ]; then
    # G
    for r in "1000 2001000" "975562500 977562500" "1000000000000 1000002000000" "999999999000007 1000000001000007" \
             "37000000000000 37000002000000" "2305843009211693952 2305843009213693952"; do
        step "G val_z2 cuda K=13 [$r]" "$L/G_k13_$(echo $r | tr ' ' '_').log" python3 val_z2.py ./z2sieve_cuda $r $T
        echo "  $(grep -h '^val_z2' "$L/G_k13_$(echo $r | tr ' ' '_').log" | sed 's/.*K=13 NL=18 TB=24 TD=7: //')"
    done
    step "G positive control cuda K=4 [60, 2000000)" "$L/G_k4_60_2000000.log" python3 val_z2.py ./z2sieve_k4_cuda 60 2000000 $T
    got=$(grep -h '^val_z2' "$L/G_k4_60_2000000.log" | grep -o 'members=\[[^]]*\]')
    echo "  $got"
    [ "$got" = "members=[67, 72, 76, 77, 81, 86]" ] || { echo "  POSITIVE CONTROL FAILED"; st=1; }
    step "G threads 1 = $T cuda [K=13 999999999000007 1000000001000007]" "$L/G_T_k13.log" same ./z2sieve_cuda 999999999000007 1000000001000007
    # X
    xdump() {     # xdump EXE_A EXE_B N_LO N_HI: sorted dump lines (all but DONE) identical
        local a b; a=$("$1" "$3" "$4" dump threads=$T 2>/dev/null | grep -v '^DONE' | sort | sha256sum); b=$("$2" "$3" "$4" dump threads=$T 2>/dev/null | grep -v '^DONE' | sort | sha256sum)
        [ -n "$a" ] && [ "$a" = "$b" ] && { echo "identical: $a"; return 0; }; echo "DIFFERENT"; return 1
    }
    step "X dumps cpu = cuda [K=13 999999999000007 1000000001000007]" "$L/X_k13.log" xdump ./z2sieve ./z2sieve_cuda 999999999000007 1000000001000007
    step "X dumps cpu = cuda [K=4 60 2000000]" "$L/X_k4.log" xdump ./z2sieve_k4 ./z2sieve_k4_cuda 60 2000000
    xprod() {     # CDONE counters of the CPU and GPU builds on a production-scale slice (2.5 10^11 exponents)
        local a b; a=$(./z2sieve 4000000000000 4250000000000 threads=$T 2>/dev/null | grep '^CDONE'); b=$(./z2sieve_cuda 4000000000000 4250000000000 threads=$T 2>/dev/null | grep '^CDONE')
        echo "cpu:  $a"; echo "cuda: $b"; [ -n "$a" ] && [ "$a" = "$b" ]
    }
    step "X counters cpu = cuda [4 10^12, 4.25 10^12)" "$L/X_prod.log" xprod
    echo "  $(grep -h '^cuda' "$L/X_prod.log")"
fi
echo "== end $(date -u): $([ $st -eq 0 ] && echo 'ALL CHECKS PASSED' || echo 'SOME CHECKS FAILED')"
exit $st
