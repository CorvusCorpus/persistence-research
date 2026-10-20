# Spot-check the production binary at production scale: one full c-slice (e.g. N=36000, or a shell N0 < m <= N).
# usage: val_prodrow.py N c [exe] [n0]   (exe defaults to ./gsieve_k13; n0 restricts to the shell as gsieve n0=)
import subprocess, sys, math, random
from mpmath import mp, mpf, log10 as mlog10, floor as mfloor, power
mp.dps = 60
N, c = int(sys.argv[1]), int(sys.argv[2]); K = 13; NL = 14; XMIN = 9 * NL + 4; AMIN = 26
exe = sys.argv[3] if len(sys.argv) > 3 else "./gsieve_k13"
N0 = int(sys.argv[4]) if len(sys.argv) > 4 else -1
out = subprocess.run([exe, "main", "G", str(N), str(c), str(c + 1), "dump"] + ([f"n0={N0}"] if N0 >= 0 else []),
                     capture_output=True, text=True, check=True).stdout
S = {k: set() for k in ("S2", "R3", "R4", "S4")}
for line in out.splitlines():
    p = line.split()
    if p[0] in S: S[p[0]].add(tuple(map(int, p[1:4])))
    if p[0] == "DONE": print(line)
LG2, LG3, LG7 = 0.30102999566398119521, 0.47712125471966243730, 0.84509804001425683071
m26 = 10 ** 26; exp = set(); n = 0
p7 = pow(7, c, m26)
for b in range((6 * N - 6 * c - 2 * AMIN) // 3 + 1):
    ahi = (6 * N - 3 * b - 6 * c) // 2
    need = (XMIN - b * LG3 - c * LG7) / LG2
    alo = math.ceil(need - 1e-6) if need > AMIN else AMIN
    if N0 >= 0 and 6 * N0 - 3 * b - 6 * c >= 0: alo = max(alo, (6 * N0 - 3 * b - 6 * c) // 2 + 1)
    if alo > ahi: continue
    x = pow(2, alo, m26) * pow(3, b, m26) * p7 % m26
    for a in range(alo, ahi + 1):
        n += 1
        s = str(x).zfill(26)
        if '0' not in s and '5' not in s: exp.add((a, b, c))
        x = 2 * x % m26
done = [l for l in out.splitlines() if l.startswith("DONE")][0]
assert f" points={n} " in done, "point count mismatch"
print(f"c={c}: region points {n}, S2 C={len(S['S2'])} python={len(exp)}", "OK" if exp == S['S2'] else "MISMATCH")
assert exp == S['S2']
assert S['R3'] | S['R4'] | S['S4'] == S['S2']
L2, L3, L7 = mlog10(2), mlog10(3), mlog10(7)
for (a, b, cc) in random.sample(sorted(S['R3']), min(300, len(S['R3']))):
    f = a * L2 + b * L3 + cc * L7; f -= mfloor(f)
    lead = str(int(mfloor(power(10, f + 9))))[:5]
    assert set(lead) & set('05'), ("bad top rejection", a, b, cc, lead)
M = 10 ** 126
for (a, b, cc) in S['R4'] | S['S4']:
    s = str(pow(2, a, M) * pow(3, b, M) * pow(7, cc, M) % M).zfill(126)[: 126 - 26]
    assert bool(set(s) & set('05')) == ((a, b, cc) in S['R4']), ("limb engine", a, b, cc)
print(f"R3 sample verified ({min(300, len(S['R3']))} of {len(S['R3'])}), R4={len(S['R4'])} S4={len(S['S4'])} verified")
