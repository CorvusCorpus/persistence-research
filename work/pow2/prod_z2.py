#!/usr/bin/env python3
"""Production driver for z2sieve: covers [N_LO, N_HI) with checkpointed invocations, resuming after its last CDONE line.

Every CDONE line is written only after the chunk's integrity checks passed, so any CDONE line can be trusted, also
from an invocation that was killed. Tripwire: a CAND or MEMBER line (no 0 found in 2^n) stops the run at once; the
lines go to ALERT_z2.txt and the exit status is 42.
usage: prod_z2.py N_LO N_HI BUDGET_SECONDS      (env THREADS (default: all CPUs), EXE (default ./z2sieve))
log: prod_z2_<N_LO>_<N_HI>.log in the current directory
"""
import os, re, subprocess, sys, time

nlo, nhi, budget = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
THREADS = int(os.environ.get("THREADS", os.cpu_count() or 1))
EXE = os.environ.get("EXE", "./z2sieve")
LOG = f"prod_z2_{nlo}_{nhi}.log"
t_start = time.time()


def covered_to():
    """End of the contiguous coverage of [nlo, ...) by the CDONE lines of the log."""
    iv = []
    if os.path.exists(LOG):
        for line in open(LOG):
            m = re.match(r"CDONE j=\[\d+,\d+\) n=\[(\d+),(\d+)\)", line)
            if m:
                iv.append((int(m.group(1)), int(m.group(2))))
    end = nlo
    for a, b in sorted(iv):
        if a <= end:
            end = max(end, b)
    return end


def alerts():
    if not os.path.exists(LOG):
        return []
    return [l for l in open(LOG) if l.startswith(("CAND", "MEMBER"))]


while True:
    start = covered_to()
    if start >= nhi:
        break
    left = budget - (time.time() - t_start)
    if left < 120:
        print("budget exhausted; rerun to resume", flush=True)
        sys.exit(0)
    with open(LOG, "a") as f:
        f.write(f"# invocation start n={start} at {time.ctime()} exe={EXE} threads={THREADS}\n")
        f.flush()
        p = subprocess.Popen(["timeout", str(int(left - 60)), EXE, str(start), str(nhi), f"threads={THREADS}"],
                             stdout=f, stderr=subprocess.STDOUT)
        while p.poll() is None:
            time.sleep(5)
            if alerts():
                p.terminate()
                time.sleep(2)
                p.kill()
                break
    if alerts():
        with open("ALERT_z2.txt", "w") as a:
            a.writelines(alerts())
        print("!!!!!! z2sieve found no 0 in some 2^n: run stopped, see ALERT_z2.txt !!!!!!", flush=True)
        sys.exit(42)
    if p.returncode not in (0, 124):
        print(f"z2sieve exited with status {p.returncode}; see {LOG}", flush=True)
        if p.returncode == 3 or covered_to() == start:
            sys.exit(p.returncode or 1)       # integrity check failed, or no progress: do not retry blindly
print(f"[{nlo},{nhi}) complete", flush=True)
