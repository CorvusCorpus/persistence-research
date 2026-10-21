#!/bin/bash
# Validation of z2sieve on this machine (docs/ZEROLESS_POW2.md, "Validation"). One line per check; logs go to LOGDIR.
#   K: kernel tests: every stage against GMP/MPFR at 10^6 random n < 2^56 (both builds)
#   E: exact region checks (val_z2.py): every decision of the dump against Python integers and ball arithmetic, the
#      S2 set in both directions, CDONE coverage and the DONE counters; regions at the start of the range, across the
#      first j boundary, at chunk boundaries (10^12 = 1024 P, 10^15 = 1024000 P), at 3.7 10^13 and just below 2^61
#   C: positive control: the K = 4 build over [60, 2 10^6) must report exactly the known terms 67 72 76 77 81 86
#   T: the dumps do not depend on the number of threads
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
echo "== start $(date -u) commit $(git rev-parse --short HEAD) z2sieve.c $(sha256sum z2sieve.c | cut -c1-16)"
make all > "$L/build.log" 2>&1 || { echo "build failed, see $L/build.log"; exit 1; }
sha256sum z2sieve z2sieve_k4 > "$L/binaries.sha256"
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
echo "== end $(date -u): $([ $st -eq 0 ] && echo 'ALL CHECKS PASSED' || echo 'SOME CHECKS FAILED')"
exit $st
