#!/bin/bash
# Acceptance test of the GPU tree sieve (desktop, after `make cuda tsieve6` in work/tsieve and `make` in work/gsieve).
#   A. Exact regions: tsieve_gpu against gsieve, every DONE counter and every S2/R3/R4/S4/MEMBER dump line
#      (G and E modes, full regions and shells).
#   B. Positive control: tsieve_gpu_nl7 on the main path + gsieve's direct and odd paths = the 314 known zero-free numbers.
#   C. Production scale: tsieve_gpu against tsieve6 (the CPU tree) on a slice of the N = 10^5 cosets, every counter.
#   D. Full region N = 16000 against gsieve's DONE line, and timings.
# Exit status 0 only if A-D pass. Usage: ./val_tsieve_gpu.sh [LOGDIR]
set -u
cd "$(dirname "$0")" || exit 1
L=${1:-gpu_val}; mkdir -p "$L"
P=tp13.txt; G=../gsieve/gsieve_k13; G7=../gsieve/gsieve_k13_nl7; X=./tsieve_gpu; X7=./tsieve_gpu_nl7; T=$(nproc)
st=0
echo "== val_tsieve_gpu $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD) gpu $(sha256sum $X | cut -c1-16) threads $T" | tee "$L/summary.txt"
echo "== A. exact regions" | tee -a "$L/summary.txt"
for a in "G 4 300 -1" "G 5 500 -1" "G 3 900 600" "G 6 1200 1100" "G 7 2000 -1" "E 4 300 -1" "E 5 600 -1" "E 6 900 700"; do
    ./val_tsieve6.sh $P $a $T $G $X 2>/dev/null | tee -a "$L/summary.txt" | grep -q IDENTICAL || st=1
done
echo "== B. positive control (NL = 7, N = 160)" | tee -a "$L/summary.txt"
python3 val_e2e_tree.py $X7 $G7 $P 4 160 7 $T 2>&1 | tail -2 | tee -a "$L/summary.txt" | grep -q "exact match" || st=1
echo "== C. production scale: GPU tree vs CPU tree, N = 100000, roots 0:16 of 16384" | tee -a "$L/summary.txt"
for m in 6 7; do
    c=$(./tsieve6 $P G $m 100000 roots=0:16 threads=$T | grep ^DONE); g=$($X $P G $m 100000 roots=0:16 threads=$T 2>/dev/null | grep ^DONE)
    keys="s1 s2 s3 s4 gmp full members leaves rows hits"; same=1
    for k in $keys; do [ "$(echo "$c" | grep -o " $k=[0-9]*")" = "$(echo "$g" | grep -o " $k=[0-9]*")" ] || same=0; done
    echo "M=$m CPU: $(echo "$c" | grep -o 's1=.*hits=[0-9]*') wall=$(echo "$c" | grep -o 'wall=.*')" | tee -a "$L/summary.txt"
    echo "M=$m GPU: $(echo "$g" | grep -o 's1=.*hits=[0-9]*') wall=$(echo "$g" | grep -o 'wall=.*') -> $([ $same = 1 ] && echo MATCH || echo MISMATCH)" | tee -a "$L/summary.txt"
    [ $same = 1 ] || st=1
done
echo "== D. N = 16000 against gsieve (main G)" | tee -a "$L/summary.txt"
gd=$($G main G 16000 0 16001 threads=$T 2>/dev/null | grep ^DONE); td=$($X $P G 6 16000 threads=$T 2>/dev/null | grep ^DONE)
same=1; for k in points s1 s2 s3 s4 gmp full members; do [ "$(echo "$gd" | grep -o " $k=[0-9]*")" = "$(echo "$td" | grep -o " $k=[0-9]*")" ] || same=0; done
echo "gsieve: $(echo "$gd" | grep -o 'points=.*members=[0-9]*') $(echo "$gd" | grep -o 'wall=.*')" | tee -a "$L/summary.txt"
echo "tree:   $(echo "$td" | grep -o 'points=.*members=[0-9]*') $(echo "$td" | grep -o 'wall=.*') -> $([ $same = 1 ] && echo MATCH || echo MISMATCH)" | tee -a "$L/summary.txt"
[ $same = 1 ] || st=1
echo "== $([ $st = 0 ] && echo 'ALL PASSED' || echo 'SOMETHING FAILED') $(date -u +%FT%TZ)" | tee -a "$L/summary.txt"
exit $st
