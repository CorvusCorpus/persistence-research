# Constants for the G-sieve (chain sieve). Exact integer arithmetic only.
from mpmath import mp, mpf, log10, nint
mp.dps = 60

def dlog2_mod5(target, m):
    """x mod 4*5^(m-1) with 2^x == target (mod 5^m); 2 is a primitive root mod 5^m."""
    x = next(i for i in range(4) if pow(2, i, 5) == target % 5)
    per = 4
    for i in range(1, m):
        mod = 5 ** (i + 1)
        x = next(x + j * per for j in range(5) if pow(2, x + j * per, mod) == target % mod)
        per *= 5
    assert pow(2, x, 5 ** m) == target % 5 ** m
    return x

def dlog5_mod2(target, K):
    """f mod 2^(K-2) with 5^f == target (mod 2^K); target == 1 (mod 4)."""
    assert target % 4 == 1
    f = 0
    for i in range(0, K - 2):          # determine bit i using modulus 2^(i+3)
        mod = 2 ** (i + 3)
        if pow(5, f, mod) != target % mod:
            f += 2 ** i
        assert pow(5, f, mod) == target % mod
    assert pow(5, f, 2 ** K) == target % 2 ** K
    return f

out = []
def D(name, val): out.append(f"#define {name} {val}ULL")

for k in (12, 13):
    P1 = 4 * 5 ** (k - 1); Q = 5 ** k; m2 = 2 * k; P2 = 4 * 5 ** (m2 - 1); Q2 = 5 ** m2
    mu3 = dlog2_mod5(3, m2); mu7 = dlog2_mod5(7, m2)
    g = pow(2, P1, Q2); assert (g - 1) % Q == 0
    W = ((g - 1) // Q) % Q
    y0 = pow(pow(2, m2, Q2), -1, Q2)          # 2^(0 - 2k) mod 5^(2k)
    x1 = pow(2, P1, 10 ** k)                  # bottom-k window of 2^(P1+0)
    p = f"K{k}_"
    D(p + "P1", P1); D(p + "Q", Q); D(p + "P2", P2); D(p + "Q2", Q2)
    D(p + "MU3", mu3); D(p + "MU7", mu7); D(p + "W", W); D(p + "Y0", y0); D(p + "X1", x1)
    D(p + "TEN", 10 ** k); D(p + "TWO", 2 ** k)
    # sanity: affine lift formula on random exponents
    import random
    for _ in range(2000):
        e13 = random.randrange(P1); j = random.randrange(Q)
        e26 = e13 + P1 * j
        x = pow(2, e26 + P2, 10 ** m2)        # last 2k digits of 2^(e26 + P2), exponent >= 2k
        yy = pow(2, e13 - m2, Q2) if e13 >= m2 else (pow(2, e13, Q2) * y0) % Q2
        Alo, Ahi = yy % Q, yy // Q
        h0 = (2 ** k * Alo) // Q
        beta = (Alo * W) % Q
        Y = (Ahi + beta * j) % Q
        assert x % 10 ** k == pow(2, e13 + P1, 10 ** k)
        assert x // 10 ** k == h0 + 2 ** k * Y, (e13, j)
    # sanity: suffix collapse for 2^a 3^b 7^c, a >= 2k
    for _ in range(500):
        a = random.randrange(m2, 300); b = random.randrange(300); c = random.randrange(300)
        e = (a + b * mu3 + c * mu7) % P2
        assert (2 ** a * 3 ** b * 7 ** c) % 10 ** m2 == pow(2, e + P2, 10 ** m2)
    print(f"K1={k}: P1={P1} Q={Q} W={W} affine-lift and suffix-collapse checks passed")

# fixed-point log10 constants: round(2^64 * log10 p)
for pr in (2, 3, 5, 7):
    v = int(nint(log10(mpf(pr)) * mpf(2) ** 64))
    D(f"L{pr}", v)

# odd branch: x = 3^a 5^b 7^c, x mod 2^K = (-1)^(a+c) 5^(b + a f3 + c f7)
K = 40
f3 = dlog5_mod2((-3) % 2 ** K, K); f7 = dlog5_mod2((-7) % 2 ** K, K)
D("OF3", f3); D("OF7", f7)
# thread: residues e mod 2^(K-2) such that (eps * 5^e') mod 10^K has all odd digits, e' == e, e' >= K
def odd_digits(v, n):
    s = str(v).zfill(n)
    return all(ch in '13579' for ch in s)
threads = {}
for eps in (1, -1):
    T = [e for e in range(2) if odd_digits((eps * pow(5, e + 2 * 4, 1000)) % 1000, 3)]   # K=3, mod 2^(1)
    for KK in range(3, K):
        per = 2 ** (KK - 2); new = []
        for e in T:
            for cand in (e, e + per):
                ee = cand + 2 ** (KK - 1) * (KK + 2)        # exponent >= KK+1, same residue mod 2^(KK-1)
                if odd_digits((eps * pow(5, ee, 10 ** (KK + 1))) % 10 ** (KK + 1), KK + 1):
                    new.append(cand)
        T = new
    threads[eps] = T
    print(f"eps={eps:+d}: {len(T)} thread(s) mod 2^{K-2}: {T}")
assert len(threads[1]) == 0 and len(threads[-1]) == 1
t = threads[-1][0]
D("OTHREAD", t)
w = (-pow(5, t + 2 ** (K - 2) * 2, 10 ** K)) % 10 ** K
print("odd thread last 40 digits:", str(w).zfill(K))
# verify thread formula on planted cases
import random
for _ in range(300):
    a = random.randrange(0, 400); c = random.randrange(0, 400)
    if (a + c) % 2 == 0: c += 1
    b0 = (t - a * f3 - c * f7) % 2 ** (K - 2)
    b = b0 + 2 ** (K - 2) * 0
    if b < K: b += 2 ** (K - 2)
    x = (3 ** a * 7 ** c * pow(5, b, 10 ** K)) % 10 ** K
    assert odd_digits(x, K), (a, b, c)
print("odd-thread planted check passed")
open("gconst.h", "w").write("\n".join(out) + "\n")
print(open("gconst.h").read())
