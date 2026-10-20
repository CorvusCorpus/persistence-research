# Production driver for the chain sieve: runs a worker's job list with per-c checkpoints and resume.
import subprocess, sys, os, time, re
N = int(sys.argv[1]); worker = sys.argv[2]; budget = float(sys.argv[3])   # seconds this driver may run
N0 = int(sys.argv[4]) if len(sys.argv) > 4 else -1
def vol(c):                    # main-path points for a given c (shell N0 < m <= N)
    v = 3 * (N - c) ** 2
    if N0 >= 0 and c < N0: v -= 3 * (N0 - c) ** 2
    return max(v, 0)
_tot = sum(vol(c) for c in range(N + 1)); _extra = 0.05 * _tot   # direct+odd on worker B (approx.)
_acc = 0; CM = 0
while _acc < (_tot + _extra) / 2: _acc += vol(CM); CM += 1
JOBS = {"A": [("main", 0, CM)], "B": [("main", CM, N + 1), ("direct", 0, N + 1), ("odd", 0, N + 1)]}.get(worker)
if os.environ.get("JOBS"):                      # explicit job list: mode:c0:c1,mode:c0:c1
    JOBS = [(j.split(":")[0], int(j.split(":")[1]), int(j.split(":")[2])) for j in os.environ["JOBS"].split(",")]
EXE = "./gsieve_k13"
t_start = time.time()
def logname(mode, c0, c1): return f"prod_N{N}_n0{N0}_{mode}_{c0}_{c1}.log"
def resume_point(path, c0):
    last = c0 - 1
    if os.path.exists(path):
        for line in open(path):
            m = re.match(r"CDONE (\d+) ", line)
            if m: last = max(last, int(m.group(1)))
    return last + 1
for (mode, c0, c1) in JOBS:
    path = logname(mode, c0, c1)
    while True:
        start = resume_point(path, c0)
        if start >= c1: break
        left = budget - (time.time() - t_start)
        if left < 120: print("budget exhausted; rerun to resume", flush=True); sys.exit(0)
        with open(path, "a") as f:
            f.write(f"# invocation start c={start} at {time.ctime()}\n"); f.flush()
            args = [EXE, mode, "G", str(N), str(start), str(c1), "prog"] + ([f"n0={N0}"] if N0 >= 0 else [])
            subprocess.run(["timeout", str(int(left - 60)), *args], stdout=f, stderr=subprocess.STDOUT)
    print(f"job {mode} [{c0},{c1}) complete", flush=True)
print(f"worker {worker} finished all jobs", flush=True)
