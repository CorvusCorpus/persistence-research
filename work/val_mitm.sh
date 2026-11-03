#!/bin/bash
# Validation of the MITM top steps (docs/MITM_STEPS.md): writes LOGDIR/val_mitm.out and the e2e logs.
#   1. mitm_ref self-test (hash join against brute force, planted genuine numbers);
#   2. the CPU programs against the reference (mitm_fast.py test, test_top);
#   3. the GPU probe against mitm_top and mpmath (mitm_fast.py test_gpu: csrc/mitmgpu_emu and csrc/mitmgpu);
#   4. NODE re-proved from scratch with MITM steps (val_twochain_e2e.sh with PC_MITM=1), once with the CPU probe (MITM_TOP =
#      csrc/mitm_top) and once with the GPU probe; every bag record of the two runs must be identical.
# usage (from work/): ./val_mitm.sh LOGDIR [NODE] [JOBS]      (default node 24192, 8 workers)
set -u
cd "$(dirname "$0")" || exit 1
LOG=$1; NODE=${2:-24192}; J=${3:-8}; W=$(pwd)
mkdir -p "$LOG"; OUT="$LOG/val_mitm.out"; bad=0
{
  echo "val_mitm.sh $NODE $J  commit $(git rev-parse HEAD)$(git diff --quiet HEAD -- . || echo ' (working tree differs)')  $(date -u +%FT%TZ)"
  make -s csrc/mitm_bottom csrc/mitm_top csrc/mitmgpu csrc/mitmgpu_emu >/dev/null || echo "BUILD FAILED"
  sha256sum csrc/mitm_bottom csrc/mitm_top csrc/mitmgpu csrc/mitmgpu_emu
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null
} > "$OUT" 2>&1
run() { echo "== $*" >> "$OUT"; "$@" >> "$OUT" 2>&1 || { echo "FAILED: $*" >> "$OUT"; bad=$((bad + 1)); }; }
run python3 mitm_ref.py
run python3 mitm_fast.py test
run python3 mitm_fast.py test_top
run python3 mitm_fast.py test_gpu csrc/mitmgpu_emu
run python3 mitm_fast.py test_gpu csrc/mitmgpu
for probe in cpu gpu; do
  exe=$W/csrc/mitm_top; [ $probe = gpu ] && exe=$W/csrc/mitmgpu
  S=/root/scratch/val_mitm_${probe}_$NODE
  echo "== e2e $NODE, $probe probe ($exe)" >> "$OUT"
  PC_MITM=1 PC_GB_FLOOR=1 MITM_TOP=$exe OMP_NUM_THREADS=2 ./val_twochain_e2e.sh "$NODE" "$J" "$S" >> "$OUT" 2>&1 || bad=$((bad + 1))
  cp "$S/run.log" "$LOG/e2e_${NODE}_${probe}_run.log"
  echo "MITM steps $(grep -c 'MITM top' "$S/run.log"), failed $(grep -c FAILED "$S/run.log")" >> "$OUT"
done
echo "== bag records, CPU probe against GPU probe" >> "$OUT"
python3 - /root/scratch/val_mitm_cpu_$NODE/bagresC /root/scratch/val_mitm_gpu_$NODE/bagresC >> "$OUT" 2>&1 <<'EOF' || bad=$((bad + 1))
import json, glob, os, sys
A, B = sys.argv[1], sys.argv[2]
fa = sorted(os.path.basename(f) for f in glob.glob(A + "/*.json")); fb = sorted(os.path.basename(f) for f in glob.glob(B + "/*.json"))
bad = 0 if fa == fb and fa else 1
for f in sorted(set(fa) & set(fb)):
    a, b = json.load(open(os.path.join(A, f))), json.load(open(os.path.join(B, f)))
    if a != b: bad += 1; print("DIFF", f, {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)})
print(f"{len(fa)} / {len(fb)} bag records; {bad} differences")
sys.exit(1 if bad else 0)
EOF
echo "VAL_MITM $([ $bad = 0 ] && echo PASSED || echo "FAILED ($bad)")" >> "$OUT"
tail -1 "$OUT"
