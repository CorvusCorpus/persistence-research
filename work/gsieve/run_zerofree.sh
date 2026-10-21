#!/bin/bash
# Zero-free extension: every 7-smooth x < 10^D without a digit 0 is known (OEIS A238985; Chaffin: none in (10^140, 10^20000)).
# A zero-free x is not divisible by 10, so x = 2^a 3^b 7^c (form i) or 3^a 5^b 7^c, b >= 1 (form iii). With the digit
# length m: m <= log10(x) / log10(7) in form i and m <= log10(x) / log10(5) in form iii, so the regions m <= N1 with
# N1 = ceil(D / log10 7) + 1 (main + direct paths) and m <= N3 with N3 = ceil(D / log10 5) + 1 (odd path) contain every
# x < 10^D of the two forms that has >= 130 digits; smaller x are in work/Z.json (complete below 10^5000).
# Runs, in parallel, with the CUDA binary in mode E (zero-free) and checkpointed prod.py:
#   form i:   main (stages 1-4 on the GPU, GSIEVE_GPU34=1) and direct, N = N1, in RUNDIR/zf_formi_N<N1>/
#   form iii: odd (2-adic bitmap, GSIEVE_ODD2=1), N = N3, in RUNDIR/zf_formiii_N<N3>/
# Expected members (the only known zero-free 7-smooth numbers with >= 130 digits): 2^77 3^64 7^98 (main path, 137
# digits) and 2^25 3^227 7^28 (direct path, 140 digits). Tripwire: any other MEMBER line stops everything at once,
# is recorded in RUNDIR/ALERT_MEMBER.txt and checked by check_member.py; exit status 42. Afterwards: coverage of every
# path, members = exactly the expected two, logs archived (scripts/archive_logs.sh) to LOGDEST (logs/).
# Resumable: run it again with the same arguments.
# usage: run_zerofree.sh RUNDIR D [THREADS_FORM_I THREADS_FORM_III]      (from work/gsieve)
set -u
R=$(mkdir -p "$1" && cd "$1" && pwd) || exit 1
D=$2; T1=${3:-4}; T3=${4:-16}      # the GPU pipeline needs few CPU threads; the odd path takes the rest
cd "$(dirname "$0")" || exit 1
GS=$(pwd); REPO=$(cd ../.. && pwd); EXE=$GS/gsieve_k13_cuda
LOGDEST=${LOGDEST:-$GS/logs}; mkdir -p "$LOGDEST"
read N1 N3 < <(python3 -c "import math; D=$D; print(math.ceil(D / math.log10(7)) + 1, math.ceil(D / math.log10(5)) + 1)")
DI=$R/zf_formi_N$N1; DO=$R/zf_formiii_N$N3; mkdir -p "$DI" "$DO"
EXPECT='MEMBER form=0 a=25 b=227 c=28 digits=140
MEMBER form=0 a=77 b=64 c=98 digits=137'
echo "== run_zerofree D=$D (N1=$N1, N3=$N3) start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD)" \
     "binary $(sha256sum "$EXE" | cut -c1-16) threads $T1+$T3"
unexpected() { grep -hs '^MEMBER' "$DI"/*.log "$DO"/*.log | sort -u | grep -vxF "$EXPECT"; }
alarm() {
    echo "!!!!!! UNEXPECTED MEMBER: A NEW ZERO-FREE 7-SMOOTH NUMBER? ALL RUNS STOPPED. $(date -u +%FT%TZ) !!!!!!"
    for p in $P1 $P3; do pkill -TERM -s "$p"; done; sleep 2; for p in $P1 $P3; do pkill -KILL -s "$p"; done
    unexpected | tee "$R/ALERT_MEMBER.txt"
    python3 check_member.py "$R/ALERT_MEMBER.txt" 2>&1 | tee -a "$R/ALERT_MEMBER.txt"
    exit 42
}
setsid bash -c "cd '$DI' && DM=E EXE='$EXE' GSIEVE_GPU34=1 GSIEVE_PIECE=67108864 THREADS=$T1 JOBS=main:0:$((N1 + 1)),direct:0:$((N1 + 1)) \
    python3 '$GS/prod.py' $N1 ALL 172800" > "$DI/prod.out" 2>&1 &
P1=$!
setsid bash -c "cd '$DO' && DM=E EXE='$EXE' GSIEVE_ODD2=1 THREADS=$T3 JOBS=odd:0:$((N3 + 1)) \
    python3 '$GS/prod.py' $N3 ALL 172800" > "$DO/prod.out" 2>&1 &
P3=$!
while kill -0 $P1 2>/dev/null || kill -0 $P3 2>/dev/null; do
    [ -n "$(unexpected)" ] && alarm
    sleep 5
done
wait $P1; r1=$?; wait $P3; r3=$?
[ -n "$(unexpected)" ] && alarm
st=0; [ $r1 -eq 0 ] && [ $r3 -eq 0 ] || st=1
echo "== form i (prod.py exit $r1):"; cat "$DI/prod.out"; grep -h '^DONE' "$DI"/*.log | cut -d' ' -f2-
echo "== form iii (prod.py exit $r3):"; cat "$DO/prod.out"; grep -h '^DONE' "$DO"/*.log | cut -d' ' -f2-
o1=$(python3 aggregate.py $N1 -1 "$DI"); o3=$(python3 aggregate.py $N3 -1 "$DO"); echo "$o1"; echo "$o3"
echo "$o1" | grep -qx "main: c covered $((N1 + 1))/$((N1 + 1)), missing \[\]" || st=1
echo "$o1" | grep -qx "direct: c covered $((N1 + 1))/$((N1 + 1)), missing \[\]" || st=1
echo "$o3" | grep -qx "odd: c covered $((N3 + 1))/$((N3 + 1)), missing \[\]" || st=1
got=$(grep -hs '^MEMBER' "$DI"/*.log "$DO"/*.log | sort -u)
if [ "$got" = "$(echo "$EXPECT" | sort -u)" ]; then echo "members: exactly the two expected known numbers"; else echo "MEMBERS DIFFER FROM EXPECTED:"; echo "$got"; st=1; fi
for d in "$DI" "$DO"; do sh "$REPO/scripts/archive_logs.sh" "$d" && cp "$d.tar.xz" "$d.sha256" "$LOGDEST/" || st=1; done
echo "== run_zerofree end $(date -u +%FT%TZ): $([ $st -eq 0 ] && echo 'ALL PASSED' || echo 'SOMETHING FAILED')"
exit $st
