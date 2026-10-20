"""Exact companion to zf7.c (pure Python big integers, no floating point).

  python3 verify.py small K            -> exact list of zero-free 7-smooth x < 10^(K+1)
  python3 verify.py validate K D FILE  -> brute force: compare zf7 output with exact set
  python3 verify.py final K FILE       -> exact check of zf7 survivors, then persistence graph
"""
import sys
from math import prod

sys.set_int_max_str_digits(0)


def family(bound):
    """Yield (tag, a, b, c, x) for every x = 2^a3^b7^c or 5^a3^b7^c (a>=1) with x < bound."""
    for tag, base, amin in (("A", 2, 0), ("B", 5, 1)):
        p7, c = 1, 0
        while p7 < bound:
            p37, b = p7, 0
            while p37 < bound:
                x, a = p37 * base ** amin, amin
                while x < bound:
                    yield tag, a, b, c, x
                    x *= base; a += 1
                p37 *= 3; b += 1
            p7 *= 7; c += 1


def value(tag, a, b, c):
    return (2 if tag == "A" else 5) ** a * 3 ** b * 7 ** c


def zero_free(x):
    return "0" not in str(x)


def P(x):
    return prod(int(d) for d in str(x))


def persistence(x):
    s = 0
    while x > 9:
        x = P(x); s += 1
    return s


def small(K):
    bound = 10 ** (K + 1)
    out, n = set(), 0
    for _, _, _, _, x in family(bound):
        n += 1
        if zero_free(x):
            out.add(x)
    return out, n


def validate(K, D, fname):
    lo, hi, mod = 10 ** (K + 1), 10 ** D, 10 ** K
    exact = set()
    for tag, a, b, c, x in family(hi):
        if x >= lo and "0" not in str(x % mod).zfill(K):
            exact.add((tag, a, b, c))
    got = set()
    for line in open(fname):
        tag, a, b, c = line.split()
        t = (tag, int(a), int(b), int(c))
        if lo <= value(*t) < hi:
            got.add(t)
    print(f"validate K={K} D={D}: brute-force {len(exact)}  zf7 {len(got)}  "
          f"missing {len(exact - got)}  spurious {len(got - exact)}")
    return exact == got


def final(K, fname):
    zf, n = small(K)
    print(f"exact pass: {n} family members below 10^{K+1}, {len(zf)} zero-free")
    surv = 0
    for line in open(fname):
        tag, a, b, c = line.split()
        x = value(tag, int(a), int(b), int(c))
        surv += 1
        if zero_free(x):
            zf.add(x)
    zf |= set(range(1, 10))
    print(f"zf7 survivors checked exactly: {surv}")
    zl = sorted(zf)
    print(f"zero-free 7-smooth numbers found: {len(zl)}; largest = {zl[-1]} ({len(str(zl[-1]))} digits)")
    pers = {x: persistence(x) for x in zl}
    mx = max(pers.values())
    print("persistence histogram:", {k: sum(1 for v in pers.values() if v == k) for k in range(mx + 1)})
    for x in zl:
        if pers[x] >= 9:
            print(f"  persistence {pers[x]}: {x}")
    return zl, pers


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "small":
        zf, n = small(int(sys.argv[2]))
        print(n, len(zf), max(zf))
    elif mode == "validate":
        ok = validate(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
        sys.exit(0 if ok else 1)
    elif mode == "final":
        final(int(sys.argv[2]), sys.argv[3])
