# Brute-force reference for tsieve: every point of the region, exact integers. Checks the point count with K
# admissible last digits (tsieve "points") and the exact set of points with 2K admissible last digits ("S" lines).
# usage: val_tsieve.py PARAMFILE N N0 AMIN TSIEVE_EXE
import subprocess, sys
par, N, N0, AMIN, exe = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
K = int(open(par).readline().split()[1])
m1, m2 = 10 ** K, 10 ** (2 * K)
def ok(v, n): return all(d not in "05" for d in str(v).zfill(n)[-n:])
okK = {}                                    # memo on the residue mod 10^K
npts = nK = 0; S = set()
for c in range(N + 1):
    for b in range(0, (6 * N - 6 * c) // 3 + 1):
        amax = (6 * N - 3 * b - 6 * c) // 2
        a = AMIN
        if amax < a: continue
        x = pow(2, a, m2) * pow(3, b, m2) * pow(7, c, m2) % m2
        while a <= amax:
            if N0 < 0 or 2 * a + 3 * b + 6 * c > 6 * N0:
                npts += 1
                lo = x % m1
                k = okK.get(lo)
                if k is None: k = okK[lo] = ok(lo, K)
                if k:
                    nK += 1
                    if ok(x, 2 * K): S.add((a, b, c))
            x = x * 2 % m2; a += 1
out = subprocess.run([exe, par, str(N), str(N0), str(AMIN), "dump"], capture_output=True, text=True, check=True).stdout
T = set(); done = None
for line in out.splitlines():
    if line.startswith("S "):
        t = tuple(map(int, line.split()[1:]))
        if t in T: sys.exit(f"duplicate point {t}")
        T.add(t)
    elif line.startswith("DONE"): done = dict(kv.split("=") for kv in line.split()[1:])
print(f"region points {npts}, K-digit survivors {nK} (tsieve points={done['points']}), "
      f"2K-digit survivors {len(S)} (tsieve {len(T)}), leaves {done['leaves']} (4^K = {4 ** K})")
ok_all = int(done["points"]) == nK and S == T and int(done["leaves"]) == 4 ** K and int(done["surv"]) == len(S)
print("EXACT MATCH" if ok_all else f"MISMATCH: missing {sorted(S - T)[:5]} extra {sorted(T - S)[:5]}")
sys.exit(0 if ok_all else 1)
