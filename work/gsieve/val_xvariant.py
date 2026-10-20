# Cross-variant check on regions beyond the reach of the Python checkers: every binary (stage-2 variant,
# thread count) must produce the same set of dump lines (S2/R3/R4/S4, D4/T4, MEMBER) and the same DONE counters.
# Lines are compared as sorted multisets because the scalar path interleaves S2 and R3/R4/S4 lines per candidate.
# usage: val_xvariant.py "MODE DM N C0 C1 [n0=M]" EXE[:THREADS] [EXE[:THREADS] ...]
import hashlib, subprocess, sys, time
spec = sys.argv[1].split()
ref = None
for e in sys.argv[2:]:
    exe, _, t = e.partition(":")
    t0 = time.time()
    p = subprocess.Popen([exe, *spec, "dump", f"threads={t or 1}"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    lines = []; done = None
    for line in p.stdout:
        if line.startswith("DONE"): done = line
        else: lines.append(line)
    assert p.wait() == 0, f"{exe} failed"
    lines.sort()
    h = hashlib.sha256("".join(lines).encode()).hexdigest()[:16]
    d = dict(kv.split("=", 1) for kv in done.split()[1:])
    counters = tuple(d[k] for k in ("points", "s1", "s2", "s3", "s4", "gmp", "full", "members"))
    print(f"{e}: variant={d.get('variant')} lines={len(lines)} sha256={h} counters={counters} ({time.time() - t0:.0f}s)", flush=True)
    if ref is None: ref = (h, counters)
    assert (h, counters) == ref, "MISMATCH between binaries"
print(f"identical dumps and counters across {len(sys.argv) - 2} binaries for: {sys.argv[1]}")
