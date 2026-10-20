# End-to-end E-mode validation: the three paths must reproduce exactly the zero-free 7-smooth numbers
# with >= XMIN digits in the region (Z.json is complete below 10^5000).
import subprocess, sys, json
exe, N = sys.argv[1], int(sys.argv[2]); NL = int(sys.argv[3]); XMIN = 9 * NL + 4
Z = json.load(open('../Z.json'))
def fac(x):
    e = {}
    for p in (2, 3, 5, 7):
        while x % p == 0: x //= p; e[p] = e.get(p, 0) + 1
    return e
exp = set()
for x in Z:
    if x < 10 or '0' in str(x) or len(str(x)) < XMIN + 1: continue
    f = fac(x)
    if 5 in f: key = (1, f.get(3, 0), f[5], f.get(7, 0)); ok = f.get(3, 0) + 2 * f[5] + 2 * f.get(7, 0) <= 2 * N
    else: key = (0, f.get(2, 0), f.get(3, 0), f.get(7, 0)); ok = 2 * f.get(2, 0) + 3 * f.get(3, 0) + 6 * f.get(7, 0) <= 6 * N
    if ok: exp.add(key)
got = set(); tot = {}
for mode in ("main", "direct", "odd"):
    r = subprocess.run([exe, mode, "E", str(N), "0", str(N + 1)], capture_output=True, text=True, check=True)
    for line in r.stdout.splitlines():
        if line.startswith("MEMBER"):
            d = dict(kv.split('=') for kv in line.split()[1:])
            key = (int(d['form']), int(d['a']), int(d['b']), int(d['c']))
            assert key not in got, ("duplicate", key)
            got.add(key)
        elif line.startswith("DONE"): print(line)
print(f"expected {len(exp)}  found {len(got)}  missing {len(exp - got)}  extra {len(got - exp)}")
assert got == exp
print("E-mode end-to-end: exact match")
