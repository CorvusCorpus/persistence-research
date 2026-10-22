"""prove_cum.py with its case loops on the GPU (bottom: cbottom5 -> csrc/bottomgpu; top: ctop4 -> csrc/topgpu), bags in
parallel processes.

Every statement, fact, checkpoint and output file is prove_cum's: prove_bag runs unchanged; only its step modules are
cbottom5 / ctop4 (the same decisions and counters as bottomc3 / topc3: val_bottomgpu.py, val_topgpu.py) instead of
cbottom4 / ctop3, and the bags of a node are proved by a pool of worker processes (each bag writes its own bagresC/ckptC
file, as runner4/prove_cum do). The node's proofC_<node>.json is written exactly as prove_cum.py writes it, once every bag
is proved.
usage (from work/): python3 prove_cum_gpu.py [--jobs J] [--cpu-bottom] [--cpu-top] NODE ...
  --cpu-bottom / --cpu-top keep cbottom4 / ctop3 (the CPU paths), for comparisons."""
import sys, os, json, time
from multiprocessing import Pool

def _prove(args):
    bag, cpu_b, cpu_t = args
    import prove_cum, cbottom5, ctop4
    if not cpu_b: prove_cum.cbottom4 = cbottom5
    if not cpu_t: prove_cum.ctop3 = ctop4
    out = []
    Lmax, res = prove_cum.prove_bag(bag, log=lambda x: out.append(x))
    return bag, Lmax, res, out

if __name__ == "__main__":
    args = sys.argv[1:]; jobs = 4; cpu = False; cpu_t = False
    if '--jobs' in args: i = args.index('--jobs'); jobs = int(args[i + 1]); del args[i:i + 2]
    if '--cpu-bottom' in args: args.remove('--cpu-bottom'); cpu = True
    if '--cpu-top' in args: args.remove('--cpu-top'); cpu_t = True
    import prove_cum
    from bags import bags_of
    Z = json.load(open('Z.json'))
    for s in map(int, args):
        t0 = time.time()
        bs = sorted(bags_of(s), key=len)
        print(f"node {s}: {len(bs)} digit multisets, {jobs} workers, bottom on {'CPU' if cpu else 'GPU'}, top on {'CPU' if cpu_t else 'GPU'}", flush=True)
        results = {}
        with Pool(jobs) as pool:
            for bag, Lmax, res, out in pool.imap_unordered(_prove, [(tuple(b), cpu, cpu_t) for b in bs]):
                for x in out: print(x, flush=True)
                results[tuple(sorted(bag))] = (Lmax, res)
        Lb = 0; recs = []
        for b in bs:                                      # prove_cum.py's order and record format
            Lmax, res = results[tuple(sorted(b))]
            Lb = max(Lb, Lmax); recs.append([list(b), Lmax, res['top'], res['bot']])
        assert Lb <= 5000
        bset = {tuple(sorted(b)) for b in bs}
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                     and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
        prove_cum._save(f'proofC_{s}.json', {s: {'preimages': pre, 'Lbound': Lb, 'bags': recs, 'method': 'cumulative'}})
