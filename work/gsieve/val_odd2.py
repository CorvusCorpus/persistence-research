# Validate the odd2 path (GSIEVE_ODD2=1: zero-free mode, x = 3^a 5^b 7^c, b >= 26 by the 2-adic bitmap) against exact
# big-integer arithmetic on a full region: the point count; the set of points whose last 26 digits have no 0 (S2 lines,
# b >= 26) and of points whose last 54 digits have no 0 (D4 lines, b <= 25, the unchanged point-by-point code); every
# leading-digit rejection (R3: a 0 among the first 5 digits); every limb-engine decision (R4 iff a 0 among digits
# 26..9NL-1); and the members against the zero-free numbers found by exact arithmetic among the S4 points.
# usage: val_odd2.py N EXE [n0] [NL]      (EXE = a CUDA build, e.g. ./gsieve_k13_cuda; NL = 14 unless an _nl7 build)
import math, os, subprocess, sys
sys.set_int_max_str_digits(0)
N, exe = int(sys.argv[1]), sys.argv[2]
N0 = int(sys.argv[3]) if len(sys.argv) > 3 else -1
NL = int(sys.argv[4]) if len(sys.argv) > 4 else 14
XMIN = 9 * NL + 4
cmd = [exe, "odd", "E", str(N), "0", str(N + 1), "dump", "threads=8"] + ([f"n0={N0}"] if N0 >= 0 else [])
out = subprocess.run(cmd, capture_output=True, text=True, check=True, env=dict(os.environ, GSIEVE_ODD2="1"))
S = {k: set() for k in ("S2", "R3", "R4", "S4", "D4")}; members = set(); done = None
for line in out.stdout.splitlines():
    p = line.split()
    if p[0] in ("S2", "R3", "R4", "S4"): S[p[0]].add(tuple(map(int, p[1:4])))
    elif p[0] == "D4": S["D4"].add(tuple(map(int, p[2:5])))
    elif p[0] == "MEMBER":
        d = dict(kv.split("=") for kv in p[1:]); members.add((int(d["a"]), int(d["b"]), int(d["c"])))
    elif p[0] == "DONE": done = dict(kv.split("=", 1) for kv in p[1:])
assert "odd2:" in out.stderr, "odd2 was not active"
LG3, LG5, LG7 = 0.47712125471966243730, 0.69897000433601880479, 0.84509804001425683071
m26, m54, mL = 10 ** 26, 10 ** 54, 10 ** (9 * NL)
exp_S2, exp_D4, npts = set(), set(), 0
for c in range(N + 1):
    for b in range(1, N + 1):                      # the point set of odd_c in mode E: every b, a in [alo, amax]
        amax = 2 * N - 2 * b - 2 * c
        if amax < 0: break
        need = (XMIN - b * LG5 - c * LG7) / LG3
        alo = math.ceil(need - 1e-6) if need > 0 else 0
        if N0 >= 0 and 2 * N0 - 2 * b - 2 * c >= 0: alo = max(alo, 2 * N0 - 2 * b - 2 * c + 1)
        if alo > amax: continue
        m = m26 if b >= 26 else m54
        x = pow(3, alo, m) * pow(5, b, m) * pow(7, c, m) % m
        for a in range(alo, amax + 1):
            npts += 1
            s = str(x).zfill(26 if b >= 26 else 54)
            if "0" not in s: (exp_S2 if b >= 26 else exp_D4).add((a, b, c))
            x = 3 * x % m
print(f"DONE {' '.join(f'{k}={done[k]}' for k in ('points', 's1', 's2', 's3', 's4', 'gmp', 'full', 'members'))}")
assert int(done["points"]) == npts, ("point count", done["points"], npts)
print(f"points {npts}: equal")
assert S["S2"] == exp_S2, ("S2 set", len(S["S2"]), len(exp_S2), sorted(S["S2"] ^ exp_S2)[:5])
assert S["D4"] == exp_D4, ("D4 set", len(S["D4"]), len(exp_D4))
print(f"S2 (b >= 26, last 26 digits zero-free): {len(exp_S2)} = exact; D4 (b <= 25, last 54 digits): {len(exp_D4)} = exact")
assert S["R3"] | S["R4"] | S["S4"] == S["S2"] and not (S["R3"] & S["R4"]) and not (S["R3"] & S["S4"]) and not (S["R4"] & S["S4"])
for (a, b, c) in S["R3"]:
    assert "0" in str(3 ** a * 5 ** b * 7 ** c)[:5], ("leading-digit table rejected a zero-free prefix", a, b, c)
exp_mem = set()
for (a, b, c) in S["R4"] | S["S4"]:
    x = pow(3, a, mL) * pow(5, b, mL) * pow(7, c, mL) % mL
    bad = "0" in str(x).zfill(9 * NL)[: 9 * NL - 26]       # digits 26 .. 9NL-1 from the bottom
    assert bad == ((a, b, c) in S["R4"]), ("limb engine", a, b, c)
    if not bad and "0" not in str(3 ** a * 5 ** b * 7 ** c): exp_mem.add((a, b, c))
for (a, b, c) in S["D4"]:
    if "0" not in str(3 ** a * 5 ** b * 7 ** c): exp_mem.add((a, b, c))
print(f"R3 {len(S['R3'])} (every one has a 0 among its first 5 digits), R4 {len(S['R4'])}, S4 {len(S['S4'])}: verified")
assert members == exp_mem, ("members", sorted(members ^ exp_mem)[:5])
print(f"members: {len(members)} = exact")
print("odd2 decisions verified")
