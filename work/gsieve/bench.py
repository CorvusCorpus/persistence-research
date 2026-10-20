# Throughput benchmark of the chain-sieve main path: candidates (region points) per second by stage-2 variant
# and thread count. Default workload: N = 36000, c in [7540, 7560) (4.86e10 points, 5 c values per thread at
# T = 4). Table construction is excluded (DONE `wall`/`time` start after it). Best of REPEAT runs.
# usage: bench.py [--variants avx512,avx2,scalar] [--threads 1,2,4] [--repeat 2] [--range N:c0:c1]
import subprocess, sys
opt = {"variants": "avx512,avx2,scalar", "threads": "1,2,4", "repeat": "2", "range": "36000:7540:7560"}
a = sys.argv[1:]
while a: opt[a[0].lstrip('-')] = a[1]; a = a[2:]
N, c0, c1 = opt["range"].split(":")
print(f"workload: main G N={N} c=[{c0},{c1})")
print(f"{'variant':8} {'T':>2} {'points':>15} {'wall s':>8} {'cpu s':>8} {'Gpts/s':>8} {'ns/pt/thread':>13} {'s1':>13} {'s2':>11}")
ref = None
for v in opt["variants"].split(","):
    for t in map(int, opt["threads"].split(",")):
        best = None
        for _ in range(int(opt["repeat"])):
            out = subprocess.run([f"./gsieve_k13_{v}", "main", "G", N, c0, c1, f"threads={t}"],
                                 capture_output=True, text=True, check=True).stdout
            d = dict(kv.split("=", 1) for kv in [l for l in out.splitlines() if l.startswith("DONE")][0].split()[1:])
            if best is None or float(d["wall"]) < float(best["wall"]): best = d
        d = best; pts = int(d["points"]); wall = float(d["wall"]); cpu = float(d["time"])
        key = (d["points"], d["s1"], d["s2"], d["s3"], d["s4"])
        if ref is None: ref = key
        assert key == ref, ("counters differ between runs", key, ref)
        print(f"{v:8} {t:>2} {pts:>15,} {wall:>8.2f} {cpu:>8.2f} {pts / wall / 1e9:>8.3f} {1e9 * wall * t / pts:>13.3f} "
              f"{int(d['s1']):>13,} {int(d['s2']):>11,}", flush=True)
print("counters identical across all variants and thread counts")
