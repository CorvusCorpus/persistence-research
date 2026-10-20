#!/bin/bash
# Task 4 driver: extend the chain sieve in shells N0 < m <= N, then give every shell an independent second pass.
#   gpu:  checkpointed production with prod.py (worker ALL: main, direct, odd) and the CUDA binary, in
#         RUNDIR/prod_N<N>_n0<N0>/; then aggregate.py (full coverage, no members), the main-path survival rates against
#         0.8^12 and 0.8^13, and a reproducible archive of the logs (scripts/archive_logs.sh) copied to LOGDEST (logs/).
#   cpu:  second pass with the CPU binary (an independent implementation of stages 1-2): every complete production
#         invocation is re-run and must reproduce all eight DONE counters (val_proddone.py); then aggregate.py and
#         compare_cdone.py (per-c counts) against the production logs; logs in RUNDIR/check_N<N>_n0<N0>/, archived to
#         CHECKDEST (../../validation_pc/task4_check/).
# Tripwire: a MEMBER line is a new chain number (more than 27,950 digits here). The moment one appears in any log under
# RUNDIR, everything is stopped, the lines are written to RUNDIR/ALERT_MEMBER.txt, check_member.py recomputes each number
# and its persistence orbit exactly, and the driver exits with status 42.
# Resumable: run it again with the same arguments and it continues (prod.py checkpoints; complete logs are reused).
# usage: run_shells.sh RUNDIR "N0:N N0:N ..." [gpu|cpu|both]      (from work/gsieve; exit status 0 only if all pass)
set -u
R=$(mkdir -p "$1" && cd "$1" && pwd) || exit 1
SHELLS=$2; WHAT=${3:-both}
cd "$(dirname "$0")" || exit 1
GS=$(pwd); REPO=$(cd ../.. && pwd); GPU=$GS/gsieve_k13_cuda; CPU=$GS/gsieve_k13; T=$(nproc)
LOGDEST=${LOGDEST:-$GS/logs}; CHECKDEST=${CHECKDEST:-$REPO/validation_pc/task4_check}     # where the archives go
mkdir -p "$LOGDEST" "$CHECKDEST"
echo "== run_shells $WHAT [$SHELLS] start $(date -u +%FT%TZ) commit $(git rev-parse --short HEAD)" \
     "gpu $(sha256sum "$GPU" | cut -c1-16) cpu $(sha256sum "$CPU" | cut -c1-16) threads $T"
st=0
alarm() {
    echo "!!!!!! MEMBER LINE FOUND: A NEW CHAIN NUMBER. ALL RUNS STOPPED. $(date -u +%FT%TZ) !!!!!!"
    grep -h '^MEMBER' "$R"/*/*.log | sort -u | tee "$R/ALERT_MEMBER.txt"
    python3 check_member.py "$R/ALERT_MEMBER.txt" 2>&1 | tee -a "$R/ALERT_MEMBER.txt"
    exit 42
}
run_watched() {   # run a command in its own session; stop all of it at once if a MEMBER line appears under RUNDIR
    setsid "$@" &
    local pid=$!     # the session id: prod.py's `timeout` moves the binary to another process group, not another session
    while kill -0 "$pid" 2>/dev/null; do
        if grep -qs '^MEMBER' "$R"/*/*.log; then
            pkill -TERM -s "$pid"; sleep 2; pkill -KILL -s "$pid"; alarm
        fi
        sleep 5
    done
    wait "$pid"; local rc=$?
    grep -qs '^MEMBER' "$R"/*/*.log && alarm
    return $rc
}
survival() {      # main-path stage survival over all main DONE lines of a log directory, against 0.8^12 and 0.8^13
    python3 - "$1" <<'EOF'
import glob, math, sys
pts = s1 = s2 = 0
for p in glob.glob(sys.argv[1] + "/*_main_*.log"):
    for line in open(p):
        if line.startswith("DONE "):
            d = dict(kv.split("=", 1) for kv in line.split()[1:])
            pts += int(d["points"]); s1 += int(d["s1"]); s2 += int(d["s2"])
p1, p2 = 0.8 ** 12, 0.8 ** 13
r1, r2 = s1 / pts, s2 / s1
z1, z2 = (r1 - p1) / math.sqrt(p1 * (1 - p1) / pts), (r2 - p2) / math.sqrt(p2 * (1 - p2) / s1)   # binomial scale
ok = abs(z1) < 6 and abs(z2) < 6
print(f"survival: s1/points = {r1:.8f} (0.8^12 = {p1:.8f}, z = {z1:+.2f}), s2/s1 = {r2:.8f} (0.8^13 = {p2:.8f},"
      f" z = {z2:+.2f}): {'match theory' if ok else 'DEVIATE FROM THEORY'}")
sys.exit(0 if ok else 1)
EOF
}
covered() {       # aggregate.py N N0 DIR must show every c covered on all three paths and no members
    local out; out=$(python3 aggregate.py "$@"); echo "$out"
    [ "$(echo "$out" | grep -c 'c covered .*, missing \[\]$')" -eq 3 ] && echo "$out" | grep -q '^members: none$'
}
if [ "$WHAT" != cpu ]; then
    for s in $SHELLS; do
        N0=${s%%:*}; N=${s##*:}; D=$R/prod_N${N}_n0$N0; mkdir -p "$D"
        echo "== production $N0 < m <= $N start $(date -u +%FT%TZ)"
        run_watched bash -c "cd '$D' && EXE='$GPU' THREADS=$T python3 '$GS/prod.py' $N ALL 86400 $N0" || st=1
        covered "$N" "$N0" "$D" || st=1
        grep -h '^DONE' "$D"/*.log | cut -d' ' -f2-6,9-16,18-
        survival "$D" || st=1
        sh "$REPO/scripts/archive_logs.sh" "$D" && cp "$D.tar.xz" "$D.sha256" "$LOGDEST/" || st=1
        echo "== production $N0 < m <= $N end $(date -u +%FT%TZ)"
    done
fi
if [ "$WHAT" != gpu ]; then
    for s in $SHELLS; do
        N0=${s%%:*}; N=${s##*:}; D=$R/prod_N${N}_n0$N0; C=$R/check_N${N}_n0$N0
        echo "== second pass (CPU) $N0 < m <= $N start $(date -u +%FT%TZ)"
        run_watched python3 val_proddone.py "$CPU" "$T" --logdir "$C" "$D"/prod_N${N}_n0${N0}_*.log || st=1
        covered "$N" "$N0" "$C" || st=1
        python3 compare_cdone.py "$N" "$N0" "$D" "$C" || st=1
        sh "$REPO/scripts/archive_logs.sh" "$C" && cp "$C.tar.xz" "$C.sha256" "$CHECKDEST/" || st=1
        echo "== second pass (CPU) $N0 < m <= $N end $(date -u +%FT%TZ)"
    done
fi
echo "== run_shells end $(date -u +%FT%TZ): $([ $st -eq 0 ] && echo 'ALL PASSED' || echo 'SOMETHING FAILED')"
exit $st
