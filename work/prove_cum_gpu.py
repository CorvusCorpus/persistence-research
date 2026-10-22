"""prove_cum.py with its case loops on the GPU (bottom: cbottom5 -> csrc/bottomgpu; top: ctop4 -> csrc/topgpu), bags in
parallel processes.

Every statement, fact, checkpoint and output file is prove_cum's: prove_bag runs unchanged; only its step modules are
cbottom5 / ctop4 (the same decisions and counters as bottomc3 / topc3: val_bottomgpu.py, val_topgpu.py) instead of
cbottom4 / ctop3, and the bags of all the given nodes are proved by one pool of worker processes, the costliest first
(estimates from logs/node_costs.json; each bag writes its own bagresC/ckptC file, as runner4/prove_cum do, so a stopped
run resumes). A node's proofC_<node>.json is written exactly as prove_cum.py writes it, once every bag of it is proved.
usage (from work/): python3 prove_cum_gpu.py [--jobs J] [--cpu-bottom] [--cpu-top] NODE ...
  --cpu-bottom / --cpu-top keep cbottom4 / ctop3 (the CPU paths), for comparisons."""
import sys
_here = sys.path.pop(0)           # work/struct.py (the Z.json script) would shadow the standard struct that pickle and
import struct, pickle             # multiprocessing need: import them with work/ off the path, then put it back
from multiprocessing import Pool
sys.path.insert(0, _here)
import os, json, time

def _save_retry(path, obj):
    """prove_cum._save (temp file + atomic rename), retrying the rename: on the Windows file system (WSL drvfs) a scanner or
    indexer briefly holding the target makes os.replace fail with EACCES (seen once with 12 workers)"""
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh: json.dump(obj, fh)
    for i in range(60):
        try:
            os.replace(tmp, path); return
        except PermissionError:
            if i == 59: raise
            time.sleep(0.5)

def _prove(args):
    bag, cpu_b, cpu_t = args
    import prove_cum, cbottom5, ctop4
    prove_cum._save = _save_retry
    if not cpu_b: prove_cum.cbottom4 = cbottom5
    if not cpu_t: prove_cum.ctop3 = ctop4
    out = []
    Lmax, res = prove_cum.prove_bag(bag, log=lambda x: out.append(x))
    return bag, Lmax, res, out

def _est_costs():
    """estimated cost per bag (logs/node_costs.json, node_costs.py's dry runs), for the order of execution only"""
    try:
        d = json.load(open(os.path.join('logs', 'node_costs.json')))
        return {pb['bag']: pb['cost'] for v in d.values() for pb in v.get('per_bag', [])}
    except (OSError, ValueError, KeyError):
        return {}

if __name__ == "__main__":
    args = sys.argv[1:]; jobs = 4; cpu = False; cpu_t = False
    if '--jobs' in args: i = args.index('--jobs'); jobs = int(args[i + 1]); del args[i:i + 2]
    if '--cpu-bottom' in args: args.remove('--cpu-bottom'); cpu = True
    if '--cpu-top' in args: args.remove('--cpu-top'); cpu_t = True
    import prove_cum
    from bags import bags_of
    Z = json.load(open('Z.json'))
    nodes = list(map(int, args)); t0 = time.time()
    bags = {s: sorted(bags_of(s), key=len) for s in nodes}   # prove_cum.py's order (records of proofC)
    # one pool for the bags of all nodes, the costliest first (estimates, else k), so that no node ends with one bag alone
    est = _est_costs()
    todo = sorted(((s, tuple(b)) for s in nodes for b in bags[s]),
                  key=lambda sb: (-est.get('_'.join(map(str, sorted(sb[1]))), 0), -len(sb[1])))
    for s in nodes:
        print(f"node {s}: {len(bags[s])} digit multisets", flush=True)
    print(f"{len(todo)} bags, {jobs} workers, bottom on {'CPU' if cpu else 'GPU'}, top on {'CPU' if cpu_t else 'GPU'}", flush=True)
    owner = {tuple(sorted(b)): s for s, b in todo}; left = {s: len(bags[s]) for s in nodes}; results = {}
    with Pool(jobs) as pool:
        for bag, Lmax, res, out in pool.imap_unordered(_prove, [(b, cpu, cpu_t) for s, b in todo]):
            for x in out: print(x, flush=True)
            key = tuple(sorted(bag)); results[key] = (Lmax, res); s = owner[key]; left[s] -= 1
            if left[s]: continue
            Lb = 0; recs = []
            for b in bags[s]:                             # prove_cum.py's order and record format
                Lmax, res = results[tuple(sorted(b))]
                Lb = max(Lb, Lmax); recs.append([list(b), Lmax, res['top'], res['bot']])
            assert Lb <= 5000
            bset = {tuple(sorted(b)) for b in bags[s]}
            pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                         and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
            print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
            _save_retry(f'proofC_{s}.json', {s: {'preimages': pre, 'Lbound': Lb, 'bags': recs, 'method': 'cumulative'}})
