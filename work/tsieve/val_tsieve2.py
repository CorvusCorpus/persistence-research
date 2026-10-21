# Brute-force reference for tsieve4 (and the hit-list versions): exact integers over every point of the region. Checks the number of points
# whose last K + M digits are admissible ("hits") and the exact set with 2K admissible digits ("S" lines).
# usage: val_tsieve2.py PARAMFILE M N N0 AMIN EXE      (EXE: tsieve4)
import subprocess, sys
par, M, N, N0, AMIN, exe = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
K = int(open(par).readline().split()[1])
m2 = 10 ** (2 * K)
def ok(v, n): return all(d not in "05" for d in str(v % 10 ** n).zfill(n))
npts = nKM = 0; S = set()
for c in range(N + 1):
    for b in range(0, (6 * N - 6 * c) // 3 + 1):
        amax = (6 * N - 3 * b - 6 * c) // 2
        a = AMIN
        if amax < a: continue
        x = pow(2, a, m2) * pow(3, b, m2) * pow(7, c, m2) % m2
        while a <= amax:
            if N0 < 0 or 2 * a + 3 * b + 6 * c > 6 * N0:
                npts += 1
                if ok(x, K + M):
                    nKM += 1
                    if ok(x, 2 * K): S.add((a, b, c))
            x = x * 2 % m2; a += 1
out = subprocess.run([exe, par, str(M), str(N), str(N0), str(AMIN), "dump"], capture_output=True, text=True, check=True).stdout
T = set(); done = None
for line in out.splitlines():
    if line.startswith("S "):
        t = tuple(map(int, line.split()[1:]))
        if t in T: sys.exit(f"duplicate point {t}")
        T.add(t)
    elif line.startswith("DONE"): done = dict(kv.split("=") for kv in line.split()[1:])
print(f"K={K} M={M} N={N} N0={N0}: region points {npts}, (K+M)-digit survivors {nKM} (hits={done['hits']}), "
      f"2K-digit survivors {len(S)} (sieve {len(T)}), leaves {done['leaves']}, rows {done['rows']}")
good = int(done["hits"]) == nKM and S == T and int(done["leaves"]) == 4 ** K and int(done["surv"]) == len(S)
print("EXACT MATCH" if good else f"MISMATCH: missing {sorted(S - T)[:5]} extra {sorted(T - S)[:5]}")
sys.exit(0 if good else 1)
