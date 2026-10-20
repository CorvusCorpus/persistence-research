#!/bin/bash
# Full validation of the chain-sieve builds on this machine (about 15 min on 4 cores). Run it on the machine that
# will do production runs, after building, before trusting new binaries. One line per check; logs go to LOGDIR.
#   A: validation suite (make test) for every variant this CPU runs (4 threads) and for the default build (1 thread)
#   B: strict-ISA builds (-march=haswell -> avx2 variant, -march=x86-64 -> scalar) in scratch copies
#   C: production slices against exact arithmetic (val_prodrow.py), including the shell boundary with n0=36000
#   D: cross-variant dumps on regions too large for the Python checkers (val_xvariant.py)
# usage: val_all.sh LOGDIR      (from work/gsieve; exit status 0 only if every check passes)
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
echo "== start $(date -u) commit $(git rev-parse --short HEAD) gsieve.c $(sha256sum gsieve.c | cut -c1-16)"
make all variants val_kernel > "$L/build.log" 2>&1 || { echo "build failed, see $L/build.log"; exit 1; }
sha256sum gsieve_k1* > "$L/binaries.sha256"
F=" $(grep -m1 '^flags' /proc/cpuinfo) "; VARS=""          # variants this CPU can run (RUNNABLE in the Makefile)
for v in ifma:avx512ifma avx512:avx512f avx2:avx2; do case "$F" in *" ${v#*:} "*) VARS="$VARS ${v%%:*}";; esac; done
VARS="$VARS scalar"; echo "variants:$VARS; default build: $(./gsieve_k13 direct G 1 0 1 2>/dev/null | grep -o 'variant=[a-z0-9]*')"
# A
for V in $VARS; do step "A suite V=$V T=4" "$L/A_suite_${V}_T4.log" env GSIEVE_THREADS=4 make test V=$V; done
step "A suite default build T=1" "$L/A_suite_default_T1.log" make test
# B
for arch in haswell x86-64; do
    T=$(mktemp -d); mkdir -p "$T/work/gsieve"; cp ../Z.json "$T/work/"
    cp gsieve.c gconst.h gconst.py Makefile val_kernel.c tiny.py val_dump.py val_e2e.py val_direct_odd.py val_prodrow.py \
       G_small.json "$T/work/gsieve/"
    step "B suite -march=$arch" "$L/B_suite_$arch.log" env GSIEVE_THREADS=4 make -C "$T/work/gsieve" test \
         CFLAGS="-O3 -march=$arch -Wall -Wno-unused-function"
    echo "  variant $(grep -o 'variant=[a-z0-9]*' "$L/B_suite_$arch.log" | sort -u | tr '\n' ' ')AVX-512 instructions:" \
         "$(objdump -d "$T/work/gsieve/gsieve_k13" | grep -cE 'zmm|%k[0-7]')"
    rm -rf "$T"
done
# C
for spec in "36000 34200" "36000 33000" "40000 39500 36000" "40000 38200 36000" "40000 35990 36000"; do
    set -- $spec
    step "C val_prodrow N=$1 c=$2 n0=${3:--1}" "$L/C_prodrow_$1_$2.log" env GSIEVE_THREADS=4 python3 val_prodrow.py $1 $2 ./gsieve_k13 ${3:+$3}
    echo "  $(grep -E '^c=' "$L/C_prodrow_$1_$2.log")"
done
# D
X="./gsieve_k13:1"; for V in $VARS; do X="$X ./gsieve_k13_$V:4"; done
for r in "main G 1000 0 1001" "main E 400 0 401" "main G 36000 20000 20001" "main G 40000 35990 36000 n0=36000" \
         "direct G 3000 0 3001" "direct E 400 0 401" "odd G 3000 0 3001" "odd E 400 0 401"; do
    step "D xvariant [$r]" "$L/D_xvariant_$(echo $r | tr ' =' '__').log" python3 val_xvariant.py "$r" $X
done
echo "== end $(date -u): $([ $st -eq 0 ] && echo 'ALL CHECKS PASSED' || echo 'SOME CHECKS FAILED')"
exit $st
