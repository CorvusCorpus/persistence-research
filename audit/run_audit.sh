#!/bin/bash
# Standard audit: independent reference (audit/ref_sieve.py) against the fast binary on full small regions and
# on production slices that contain survivors. Prints one line per check; exit status 0 only if all agree.
# usage: audit/run_audit.sh [BIN] [JOBS]      (default BIN work/gsieve/gsieve_k13, JOBS = number of CPUs)
cd "$(dirname "$0")/.." || exit 1
BIN=${1:-work/gsieve/gsieve_k13}; J=${2:-$(nproc)}
R="python3 audit/ref_sieve.py"; st=0
run() {           # run one region check; print its verdict line and the elapsed time
    local t0=$(date +%s) out
    out=$($R "$@" --exe "$BIN" --threads "$J" 2>&1); local rc=$?
    printf '%-62s %s (%ss)\n' "$*" "$(echo "$out" | head -1 | grep -o 'points=.*')" "$(( $(date +%s) - t0 ))"
    echo "$out" | grep -E "MISMATCH|DISAGREES|Error|Traceback" ; [ $rc -eq 0 ] || st=1
}
lem=$($R lemmas) || st=1; echo "$lem" | tail -1
# full small regions, G and E modes (E mode has survivors and members at this size)
for m in main direct odd; do for dm in G E; do
    n=$([ $m = main ] && echo 150 || echo 160); run $m $dm $n 0 $((n + 1)) --full
done; done
# production slices of the logged runs (N = 36000, and the shell 36000 < m <= 40000 with n0)
# --full (all digits of every survivor, ~14 ms per 30,000-digit number) only where the survivor count is small
run main G 36000 35500 35501 --full
run main G 36000 34200 34201
run main G 40000 39500 39501 n0=36000 --full
run main G 40000 35990 35991 n0=36000
# direct path: every slice has survivors (mostly from a = 0, x = 3^b 7^c)
for c in 35000 20000 1000; do run direct G 36000 $c $((c + 1)) --full; done
# odd path: slices holding thread candidates found by the production runs (30 for N = 36000, 11 for the shell)
for c in 153 1058 2164 29958; do run odd G 36000 $c $((c + 1)) --full; done
for c in 1787 36343; do run odd G 40000 $c $((c + 1)) n0=36000 --full; done
run thread G 36000 0 400 --jobs "$J" --full                  # holds (2490, 20092, 153)
run thread G 36000 1000 1400 --jobs "$J" --full              # holds c = 1058, 1301
run thread G 40000 1700 1800 n0=36000 --jobs "$J" --full     # holds c = 1787
[ $st -eq 0 ] && echo "AUDIT: reference agrees with $BIN on every region" || echo "AUDIT FAILED"
exit $st
