# Validate main-path stage decisions of gsieve against exact big-integer arithmetic.
# usage: val_dump.py K G|E|P N [exe]   (exe defaults to ./gsieve_k<K>)
import subprocess, sys, math
K = int(sys.argv[1]); dm = sys.argv[2]; N = int(sys.argv[3]); NL = 14
XMIN = 9 * NL + 4; AMIN = 2 * K
exe = sys.argv[4] if len(sys.argv) > 4 else f"./gsieve_k{K}"
out = subprocess.run([exe, "main", dm, str(N), "0", str(N + 1), "dump"], capture_output=True, text=True, check=True).stdout
S = {k: set() for k in ("S2", "R3", "R4", "S4")}; members = []; done = None
for line in out.splitlines():
    p = line.split()
    if p[0] in S: S[p[0]].add(tuple(map(int, p[1:4])))
    elif p[0] == "MEMBER": members.append(line)
    elif p[0] == "DONE": done = line
bad = set('05') if dm == 'G' else set('0') if dm == 'E' else set()
LG2, LG3, LG7 = 0.30102999566398119521, 0.47712125471966243730, 0.84509804001425683071
m26 = 10 ** (2 * K); M = 10 ** (9 * NL)
exp_S2 = set(); npts = 0
for c in range(N + 1):
    bmax = (6 * N - 6 * c - 2 * AMIN) // 3
    for b in range(bmax + 1):
        ahi = (6 * N - 3 * b - 6 * c) // 2
        need = (XMIN - b * LG3 - c * LG7) / LG2
        alo = math.ceil(need - 1e-6) if need > AMIN else AMIN
        if alo > ahi: continue
        base = (pow(3, b, m26) * pow(7, c, m26)) % m26
        x = (pow(2, alo, m26) * base) % m26
        for a in range(alo, ahi + 1):
            npts += 1
            if not (set(str(x).zfill(2 * K)) & bad): exp_S2.add((a, b, c))
            x = (2 * x) % m26
print(done)
print(f"region points (python) = {npts}")
assert f"points={npts} " in done, "point count mismatch"
print(f"S2: C={len(S['S2'])} python={len(exp_S2)}", "OK" if S['S2'] == exp_S2 else "MISMATCH")
assert S['S2'] == exp_S2
assert S['R3'] | S['R4'] | S['S4'] == S['S2'] and not (S['R3'] & S['R4']) and not (S['R3'] & S['S4']) and not (S['R4'] & S['S4'])
for (a, b, c) in S['R3']:
    s = str(2 ** a * 3 ** b * 7 ** c)
    assert set(s[:7]) & bad, ("top-table rejected a good prefix", a, b, c, s[:10])
nR4 = 0
for (a, b, c) in S['R4'] | S['S4']:
    x = (pow(2, a, M) * pow(3, b, M) * pow(7, c, M)) % M
    s = str(x).zfill(9 * NL)          # digit i from the bottom = s[-1-i]
    window = s[: 9 * NL - 2 * K]      # digits 2K .. 9NL-1 from the bottom
    isbad = bool(set(window) & bad)
    assert isbad == ((a, b, c) in S['R4']), ("limb engine wrong", a, b, c)
print(f"R3={len(S['R3'])} R4={len(S['R4'])} S4={len(S['S4'])}: top-table and limb-engine decisions verified")
for m in members: print(m)
