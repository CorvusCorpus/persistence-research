# Aggregate production logs: verify every c in [0, N] was completed for each path, sum counts, list members.
import re, sys, glob, collections
N = int(sys.argv[1]); N0 = int(sys.argv[2]) if len(sys.argv) > 2 else -1
LOGDIR = sys.argv[3] if len(sys.argv) > 3 else ('logs' if glob.glob(f'logs/prod_N{N}_*') else '.')
tot = collections.Counter(); members = []; covered = {"main": set(), "direct": set(), "odd": set()}
for path in sorted(glob.glob(f"{LOGDIR}/prod_N{N}_n0{N0}_*.log")):
    mode = path.split("/")[-1].split("_")[3]
    start = None; last = None; lastline = None
    def flush():
        if lastline is None: return
        d = dict(kv.split("=") for kv in lastline.split()[2:])
        tot[mode + "_points"] += int(d["points"]); tot[mode + "_s4"] += int(d["s4"]); tot[mode + "_members"] += int(d["members"])
        covered[mode].update(range(start, last + 1))
    for line in open(path):
        if line.startswith("# invocation start c="):
            flush(); start = int(line.split("c=")[1].split()[0]); last = None; lastline = None
        elif line.startswith("CDONE"):
            last = int(line.split()[1]); lastline = line
        elif line.startswith("MEMBER"): members.append(line.strip())
    flush()
for mode, s in covered.items():
    missing = sorted(set(range(0, N + 1)) - s)
    print(f"{mode}: c covered {len(s)}/{N+1}, missing {missing[:5]}{'...' if len(missing) > 5 else ''}")
for k in sorted(tot): print(f"  {k} = {tot[k]}")
print("points total:", sum(v for k, v in tot.items() if k.endswith("_points")))
print("members:", members if members else "none")
