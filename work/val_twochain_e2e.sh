#!/bin/bash
# Re-prove NODE from scratch with the two-chain bootstrap (PC_BOOTSTRAP=twochain) on the GPU case loops, in an empty scratch
# directory (the repository's ckptC/, bagresC/, bagres/ and proof files are not read or written). Checks:
#   - every bag record carries bootstrap 'twochain' and Lmax <= 5000, and so does every bag entry of the new proofC file;
#   - the preimage set equals that of the committed proofC_NODE.json;
#   - every genuine solution of every bag (Z.json) satisfies every proved fact D_j < top[j], a_i < bot[i];
#   - the committed records' cost (units of the cost model) next to the new run's, for the record.
# usage (from work/): ./val_twochain_e2e.sh NODE [JOBS] [SCRATCH]
set -u
cd "$(dirname "$0")" || exit 1
NODE=$1; J=${2:-8}; S=${3:-/root/scratch/twochain_$NODE}; W=$(pwd)
rm -rf "$S"; mkdir -p "$S/ckptC" "$S/bagresC" "$S/bagres"; ln -s "$W/Z.json" "$S/Z.json"
( cd "$S" && PC_BOOTSTRAP=twochain /usr/bin/time -f "wall %e s" python3 "$W/prove_cum_gpu.py" --jobs "$J" "$NODE" > run.log 2> time.log ) \
    || { echo "run FAILED"; tail -5 "$S/run.log"; cat "$S/time.log"; exit 1; }
tail -1 "$S/run.log"; cat "$S/time.log"
python3 - "$S" "$NODE" "$W" <<'EOF'
import json, sys, os, glob
S, node, W = sys.argv[1], sys.argv[2], sys.argv[3]
bad = 0
def fail(m):
    global bad; bad += 1; print("FAIL", m)
new = json.load(open(f"{S}/proofC_{node}.json"))[node]; old = json.load(open(f"{W}/proofC_{node}.json"))[node]
if new['preimages'] != old['preimages']: fail(f"preimages {new['preimages']} != committed {old['preimages']}")
if not new['Lbound'] <= 5000: fail(f"Lbound {new['Lbound']}")
for rec in new['bags']:
    if len(rec) != 5 or rec[4] != 'twochain': fail(f"proofC bag entry {rec[0]} without bootstrap tag")
Z = json.load(open(f"{W}/Z.json"))
cost_new = cost_old = 0.0; nb = nviol = nsol = 0; missing = 0
for f in sorted(glob.glob(f"{S}/bagresC/*.json")):
    r = json.load(open(f)); nb += 1
    if r.get('bootstrap') != 'twochain' or not r['Lmax'] <= 5000: fail(f"{os.path.basename(f)}: {r.get('bootstrap')} {r['Lmax']}")
    bag = tuple(r['bag']); cost_new += r['stats']['cost']
    o = os.path.join(W, 'bagresC', os.path.basename(f))
    if os.path.exists(o): cost_old += json.load(open(o))['stats']['cost']
    else: missing += 1
    for y in Z:
        s = str(y)
        if y < 10 or '0' in s or tuple(sorted(int(c) for c in s if c != '1')) != bag: continue
        nsol += 1; L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; a = [L] + pos
        for j, T in r['top'].items():
            if not L - pos[int(j)] < T: nviol += 1
        for i, H in r['bot'].items():
            if not a[int(i)] < H: nviol += 1
if nviol: fail(f"{nviol} fact violations by genuine solutions")
print(f"{nb} bags: twochain cost {cost_new:.3e} units; committed records {cost_old:.3e}" + (f" ({missing} bags without a record)" if missing else f" ({cost_old / max(cost_new, 1):.1f}x)")
      + f"; {nsol} genuine solutions checked against the facts; preimages {new['preimages'] or 'none'}")
print("E2E OK" if not bad else f"{bad} FAILURES"); sys.exit(1 if bad else 0)
EOF
