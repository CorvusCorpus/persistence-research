#!/bin/bash
# End-to-end check of the GPU paths inside the prover: prove NODE from scratch twice, in empty scratch directories (the
# repository's ckptC/, bagresC/, bagres/ and proof files are not read or written), with the case loops on the CPU
# (cbottom4 -> bottomc3, ctop3 -> topc3) and on the GPU (cbottom5 -> bottomgpu, ctop4 -> topgpu). Every bag's facts,
# length bound and step counters (top cases, bottom solves, survivor runs, resolved candidates) and the node's proofC file
# must be identical.
# usage (from work/): ./val_prove_gpu.sh NODE [JOBS] [SCRATCH]
set -u
cd "$(dirname "$0")" || exit 1
NODE=$1; J=${2:-8}; S=${3:-/root/scratch/e2e_$NODE}; W=$(pwd)
rm -rf "$S"; mkdir -p "$S"
for side in cpu gpu; do
    mkdir -p "$S/$side/ckptC" "$S/$side/bagresC" "$S/$side/bagres"; ln -s "$W/Z.json" "$S/$side/Z.json"
    flag=$([ $side = cpu ] && echo "--cpu-bottom --cpu-top" || echo)
    ( cd "$S/$side" && /usr/bin/time -f "$side wall %e s" python3 "$W/prove_cum_gpu.py" --jobs "$J" $flag "$NODE" > run.log 2> time.log ) || { echo "$side run FAILED"; cat "$S/$side/run.log" | tail -5; exit 1; }
    tail -1 "$S/$side/run.log"; cat "$S/$side/time.log"
done
python3 - "$S" "$NODE" <<'EOF'
import json, sys, os, glob
S, node = sys.argv[1], sys.argv[2]
keys = ('top', 'bottom', 'enum', 'proofs', 'cost', 'resolved', 'resolved_bottom', 'l1fail')
bad = 0; n = 0
for f in sorted(glob.glob(f"{S}/cpu/bagresC/*.json")):
    a = json.load(open(f)); b = json.load(open(f.replace('/cpu/', '/gpu/')))
    n += 1
    for k in ('bag', 'Lmax', 'top', 'bot'):
        if a[k] != b[k]: bad += 1; print("DIFF", os.path.basename(f), k, a[k], b[k])
    for k in keys:
        if a['stats'].get(k) != b['stats'].get(k): bad += 1; print("DIFF", os.path.basename(f), 'stats', k, a['stats'].get(k), b['stats'].get(k))
pa = json.load(open(f"{S}/cpu/proofC_{node}.json")); pb = json.load(open(f"{S}/gpu/proofC_{node}.json"))
if pa != pb: bad += 1; print("DIFF proofC")
print(f"{n} bags compared: {'IDENTICAL' if not bad else f'{bad} DIFFERENCES'}")
sys.exit(1 if bad else 0)
EOF
