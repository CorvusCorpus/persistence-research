# Validate G-mode direct path (a < 2K) and odd path (b <= 39 direct, b >= 40 thread) decisions.
import subprocess, sys, math
exe, K, N, NL = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
XMIN = 9 * NL + 4; AMIN = 2 * K
LG2, LG3, LG5, LG7 = 0.30102999566398119521, 0.47712125471966243730, 0.69897000433601880479, 0.84509804001425683071
m54 = 10 ** 54
def gbad_even(s): return bool(set(s) & set('05'))
def gbad_oddn5(s): return '0' in s or ('5' in s and bool(set(s) & set('2468')))
def gbad_odd5(s): return bool(set(s) & set('02468'))
D = {}
for mode in ("direct", "odd"):
    r = subprocess.run([exe, mode, "G", str(N), "0", str(N + 1), "dump"], capture_output=True, text=True, check=True)
    D[mode] = {tuple(map(int, l.split()[1:5])) for l in r.stdout.splitlines() if l[:2] in ("D4",)}
    D[mode + "T"] = {tuple(map(int, l.split()[1:5])) for l in r.stdout.splitlines() if l[:2] == "T4"}
    print([l for l in r.stdout.splitlines() if l.startswith("DONE")][0])
# direct path brute force
exp = set(); n = 0
for c in range(N + 1):
    for a in range(AMIN):
        bmax = (6 * N - 2 * a - 6 * c) // 3
        need = (XMIN - a * LG2 - c * LG7) / LG3
        blo = math.ceil(need - 1e-6) if need > 0 else 0
        for b in range(blo, bmax + 1):
            n += 1
            s = str((pow(2, a, m54) * pow(3, b, m54) * pow(7, c, m54)) % m54).zfill(54)
            if not (gbad_oddn5(s) if a == 0 else gbad_even(s)): exp.add((0, a, b, c))
print(f"direct: points {n}, C survivors {len(D['direct'])}, python {len(exp)}", "OK" if exp == D['direct'] else "MISMATCH")
assert exp == D['direct']
# odd path, b <= 39 brute force
exp = set(); n = 0
for c in range(N + 1):
    for b in range(1, 40):
        amax = 2 * N - 2 * b - 2 * c
        if amax < 0: break
        need = (XMIN - b * LG5 - c * LG7) / LG3
        alo = math.ceil(need - 1e-6) if need > 0 else 0
        for a in range(alo, amax + 1):
            n += 1
            s = str((pow(3, a, m54) * pow(5, b, m54) * pow(7, c, m54)) % m54).zfill(54)
            if not gbad_odd5(s): exp.add((1, a, b, c))
print(f"odd b<=39: points {n}, C survivors {len(D['odd'])}, python {len(exp)}", "OK" if exp == D['odd'] else "MISMATCH")
assert exp == D['odd']
# thread candidates: each must have all-odd bottom 40 digits; and completeness on b>=40 by brute force (bottom 40 digits)
m40 = 10 ** 40; full = set(); n = 0
for c in range(N + 1):
    for b in range(40, N + 1):
        amax = 2 * N - 2 * b - 2 * c
        if amax < 0: break
        for a in range(0, amax + 1):
            if a * LG3 + b * LG5 + c * LG7 < XMIN: continue
            n += 1
            s = str((pow(3, a, m40) * pow(5, b, m40) * pow(7, c, m40)) % m40).zfill(40)
            if not gbad_odd5(s): full.add((1, a, b, c))
print(f"thread: brute-force points {n}, all-odd-40 points {len(full)}, thread candidates {len(D['oddT'])}")
assert full == D['oddT']
print("G-mode direct/odd/thread decisions verified")
