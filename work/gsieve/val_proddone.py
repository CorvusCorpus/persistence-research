# Re-run complete production invocations and compare every counter with the production DONE line.
# The production logs record, per completed invocation, a DONE line with the region point count and the
# survivor counts of every stage (s1 bitmap, s2 affine lift, s3 leading digits, s4 limb engine, gmp, full,
# members). Equal counters on the same c-range are a production-scale check that a binary (another stage-2
# variant, another thread count, another machine) makes the same decisions as the production binary.
# usage: val_proddone.py EXE THREADS [--only mode:c0:c1,...] [--max-points P] [--logdir DIR] LOG [LOG ...]
#   --logdir DIR: also run with `prog` and keep each re-run's output as DIR/prod_N<N>_n0<N0>_<mode>_<c0>_<c1>.log
#   (aggregate.py-compatible; an existing complete log there is reused instead of re-running).
import os, re, subprocess, sys, time
args = sys.argv[1:]
exe, threads = args[0], int(args[1]); args = args[2:]
only, maxp, logdir = None, None, None
while args and args[0].startswith('--'):
    if args[0] == '--only': only = set(args[1].split(',')); args = args[2:]
    elif args[0] == '--max-points': maxp = int(float(args[1])); args = args[2:]
    elif args[0] == '--logdir': logdir = args[1]; os.makedirs(logdir, exist_ok=True); args = args[2:]
    else: sys.exit(f"unknown option {args[0]}")
KEYS = ("points", "s1", "s2", "s3", "s4", "gmp", "full", "members")
def parse(line):
    d = dict(kv.split("=", 1) for kv in line.split()[1:])
    d["c0"], d["c1"] = map(int, re.match(r"\[(\d+),(\d+)\)", d["c"]).groups())
    return d
jobs = []
for log in args:
    for line in open(log):
        if line.startswith("DONE "):
            d = parse(line)
            key = f"{d['mode']}:{d['c0']}:{d['c1']}"
            if only is not None and key not in only: continue
            if maxp is not None and int(d["points"]) > maxp: continue
            jobs.append((log, d))
if not jobs: sys.exit("no matching DONE lines")
ok = True
for log, d in jobs:
    cmd = [exe, d["mode"], d["dm"], d["N"], str(d["c0"]), str(d["c1"]), f"threads={threads}"]
    if int(d["N0"]) >= 0: cmd.append(f"n0={d['N0']}")
    t0 = time.time()
    if logdir:
        path = os.path.join(logdir, f"prod_N{d['N']}_n0{d['N0']}_{d['mode']}_{d['c0']}_{d['c1']}.log")
        if not (os.path.exists(path) and any(l.startswith("DONE ") for l in open(path))):
            with open(path, "w") as f:
                f.write(f"# invocation start c={d['c0']} at {time.ctime()} exe={exe} threads={threads} (re-run of {log})\n"); f.flush()
                subprocess.run(cmd + ["prog"], stdout=f, stderr=subprocess.STDOUT, check=True)
        out = open(path).read()
    else:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    r = parse([l for l in out.splitlines() if l.startswith("DONE ")][-1])
    same = all(r[k] == d[k] for k in KEYS) and r["K"] == d["K"] and r["NL"] == d["NL"]
    ok &= same
    print(f"{d['mode']} N={d['N']} N0={d['N0']} c=[{d['c0']},{d['c1']}): " +
          " ".join(f"{k}={r[k]}" + ("" if r[k] == d[k] else f"(prod {d[k]})") for k in KEYS) +
          f" | variant={r.get('variant', '?')} threads={r.get('threads', '?')} wall={time.time() - t0:.0f}s"
          f" | prod time={d['time']}s -> {'MATCH' if same else 'MISMATCH'}", flush=True)
print("all counters match the production DONE lines" if ok else "MISMATCH FOUND")
sys.exit(0 if ok else 1)
