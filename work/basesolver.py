"""Rigorous solver for bases with ONE free prime (b = 3, 4, 6):
find all y = B*q^c (B a product of bounded powers of the primes dividing b, q the only prime < b not
dividing b) whose base-b digits are a fixed multiset D of non-1 digits plus any number of 1s.

Notation: L digits; non-1 digits d_1..d_k at positions a_1 > ... > a_k; gaps g_0 = L - a_1,
g_i = a_i - a_{i+1} (1 <= i < k) and the trailing run g_k = a_k >= 0.
   (b-1) y = b^L - 1 + sum c_i b^{a_i},  c_i = (b-1)(d_i - 1).
Linear forms: for 0 <= j < k,
   Lam_j = log((b-1)B) + c log q - L log b - log(1+z_j),  z_j = sum_{i<=j} c_i b^{-(g_0+..+g_{i-1})},
   0 < Lam_j < eps0 * b^{-g_j},  eps0 = (b-2) b.
Matveev (n = 3 logarithms: q, b, (b-1)B/(1+z_j)) bounds every gap; proved bounds are substituted.
Top test of gap j: 2-dim lattice; bottom test of gap j (j >= 1, incl. trailing run j = k): y mod b^M
known, which pins c modulo N >> X via p-adic logarithms for p | b."""
import math, itertools, json, sys, time
from fractions import Fraction
import flint
from solver import plog

LN = math.log
C3 = 1.4 * 30**6 * 3**4.5

def primes_below(b): return [p for p in range(2, b) if all(p % r for r in range(2, p))]
def vp(n, p):
    n = abs(n)
    if n == 0: return 10**9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k

class Base:
    def __init__(self, b):
        self.b = b
        S = primes_below(b)
        self.div = [p for p in S if b % p == 0]                  # primes of y that divide b (bounded)
        self.bdiv = [p for p in range(2, b + 1) if b % p == 0 and all(p % r for r in range(2, p))]
        free = [p for p in S if b % p]
        assert len(free) == 1, "solver needs exactly one free prime"
        self.q = free[0]
        self.eps0 = (b - 2) * b

def pbound(bag, b, p):
    """Exact max v_p(y) over all y with non-1 digit multiset `bag` (Fonga Thm 4.1, base b, p | b)."""
    e = vp(b, p); bag = tuple(sorted(bag))
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def V(rem, a, T):
        best = vp(T, p); lim = vp(T, p) // e
        if rem:
            for i, d in enumerate(rem):
                if i and rem[i] == rem[i-1]: continue
                nrem = rem[:i] + rem[i+1:]
                for t in range(a + 1, lim + 1):
                    best = max(best, V(nrem, t, T + (b-1)*(d-1)*b**t))
        else:
            for L in range(a + 1, lim + 1):
                best = max(best, vp(T + b**L, p))
        return best
    best = 0
    for i, d in enumerate(bag):
        if (i == 0 or bag[i] != bag[i-1]) and vp((b-1)*(d-1) - 1, p) > 0:
            best = max(best, V(bag[:i] + bag[i+1:], 0, (b-1)*(d-1) - 1))
    return best

def cofactors(bag, base):
    """All allowed B = prod p^{e_p} (p | b) with e_p <= pbound and b not dividing B*q^c."""
    b = base.b
    bounds = {p: pbound(bag, b, p) for p in base.div}
    out = []
    for es in itertools.product(*[range(bounds[p] + 1) for p in base.div]):
        B = 1
        for p, e in zip(base.div, es): B *= p**e
        if B % b == 0: continue                      # y would end in 0
        out.append(B)
    return out, bounds

def matveev_partial(k, base, Bmax, known):
    """Bound X_L on L given proved bounds known[j] (g_j < known[j]) for j in 0..k (k = trailing run)."""
    b, q = base.b, base.q
    lq, lb = LN(q), LN(b)
    cst = LN((b - 1) * Bmax) + LN(b)
    def f(Bc):
        beta = 1 + LN(Bc); s = 0.0
        for j in range(k):
            if j in known: g = known[j] - 1
            else:
                A3 = s * lb + cst
                g = LN(base.eps0) / lb + C3 * beta * lq * A3
            s += g
        if k in known: L = s + known[k] - 1
        else:
            A3 = s * lb + cst
            L = 1 + LN(2) / lb + C3 * beta * lq * A3
            L = max(L, s)
        coef = L * lb / lq + 2
        return max(coef, L, 10.0), L
    Bc = 1e300
    for _ in range(300):
        nb, _ = f(Bc)
        if nb >= Bc * (1 - 1e-9): break
        Bc = nb
    Bc = max(Bc * 1.01, 10.0)
    nb, L = f(Bc)
    assert nb <= Bc
    return max(4, math.ceil(L * 1.01)), Bc

class Top2:
    """2-dim lattice test for one gap threshold G and length bound XL."""
    def __init__(self, base, XL, G):
        b, q = base.b, base.q
        self.base, self.XL, self.G = base, XL, G
        self.Xc = int(XL * LN(b) / LN(q)) + 2
        S = self.Xc
        self.C = S * b**G // base.eps0 + 1
        self.prec = max(256, self.C.bit_length() + 128)
        flint.ctx.prec = self.prec
        def rnd(x):
            v = flint.arb(self.C) * x
            m = int(v.mid().floor().unique_fmpz())
            if v - m > 0.5: m += 1
            assert abs(v - m) < 0.51
            return m
        self.rq = rnd(flint.arb(q).log()); self.rb = rnd(flint.arb(b).log())
        self.rnd = rnd
        B = flint.fmpz_mat([[1, self.rq], [0, -self.rb]]).lll()
        self.Bm = [[int(B[i, j]) for j in range(2)] for i in range(2)]
        b1, b2 = [Fraction(x) for x in self.Bm[0]], [Fraction(x) for x in self.Bm[1]]
        mu = (b1[0]*b2[0] + b1[1]*b2[1]) / (b1[0]**2 + b1[1]**2)
        self.mu = mu
        self.g1 = b1[0]**2 + b1[1]**2
        b2s = [b2[0] - mu*b1[0], b2[1] - mu*b1[1]]
        self.g2 = b2s[0]**2 + b2s[1]**2
        det = self.Bm[0][0]*self.Bm[1][1] - self.Bm[0][1]*self.Bm[1][0]
        self.det = det
        E = Fraction(51, 100) * (self.Xc + XL + 2)
        self.R2 = Fraction(self.Xc)**2 + (S + E + 1)**2

    def survivors(self, kappa_arb, check):
        """Lattice points within R of target (0,-m), m = round(C*kappa); returns genuine (c,L)."""
        m = self.rnd(kappa_arb)
        y = (0, -m)
        # coordinates z of y in basis rows: y = z1*B[0] + z2*B[1]
        B = self.Bm; det = self.det
        z1 = Fraction(y[0]*B[1][1] - y[1]*B[1][0], det)
        z2 = Fraction(B[0][0]*y[1] - B[0][1]*y[0], det)
        # quick test on the second Gram-Schmidt coordinate
        f = z2 - math.floor(z2); d = min(f, 1 - f)
        if d * d * self.g2 > self.R2: return []
        out = []
        s2 = math.isqrt(int(self.R2 / self.g2)) + 2
        for x2 in range(math.floor(z2) - s2, math.floor(z2) + s2 + 2):
            t2 = (x2 - z2)**2 * self.g2
            if t2 > self.R2: continue
            c1 = z1 - self.mu * (x2 - z2)
            s1 = math.isqrt(int((self.R2 - t2) / self.g1)) + 2
            for x1 in range(math.floor(c1) - s1, math.floor(c1) + s1 + 2):
                if (x1 - c1)**2 * self.g1 + t2 > self.R2: continue
                v = (x1*B[0][0] + x2*B[1][0], x1*B[0][1] + x2*B[1][1])
                c = v[0]
                num = c * self.rq - v[1]
                if num % self.rb: continue
                L = num // self.rb
                if 0 <= c <= self.Xc and 1 <= L <= self.XL and check(c, L):
                    out.append((c, L))
        return out

def lam_ok(base, B, c, L, one_plus_z, G):
    """Decide 0 < Lam < eps0 b^-G exactly enough (interval arithmetic); undecided counts as True."""
    b, q = base.b, base.q
    flint.ctx.prec = max(flint.ctx.prec, 64 + int(G * LN(b) / LN(2)) + int(LN(max(c, L, 2)) / LN(2)) + 64)
    lam = (flint.arb((b - 1) * B).log() + c * flint.arb(q).log() - L * flint.arb(b).log()
           - (flint.arb(one_plus_z.numerator) / one_plus_z.denominator).log())
    bound = base.eps0 * flint.arb(b)**(-G)
    if lam <= 0 or lam >= bound: return False
    return True

# ------------------------------------------------------------- bottom (1-dim congruence)
def solve_c(base, B, r, M, Xc):
    """All c in [0, Xc] with B*q^c = r (mod b^M). Raises if more than one could fit."""
    b, q = base.b, base.q
    congr = []                                     # list of (residue, modulus)
    for p in base.bdiv:
        e = M * vp(b, p); beta = vp(B, p)
        if vp(r, p) != beta and not (vp(r, p) >= e and beta >= e): return []
        n = e - beta
        if n <= 0: continue
        rho = (r // p**beta) * pow(B // p**beta, -1, p**n) % p**n
        if p == 2:
            if n < 4: return None                  # too little information (never happens for large M)
            sq = 1 if q % 4 == 1 else -1; sr = 1 if rho % 4 == 1 else -1
            if sq == 1 and sr == -1: return []
            if sq == -1: congr.append((0 if sr == 1 else 1, 2))
            lq = plog((sq * q) % 2**(n + 2), 2, n); lr = plog((sr * rho) % 2**(n + 2), 2, n)
            v = vp(lq, 2)
            if vp(lr, 2) < v and lr % 2**n: return []
            mod = 2**(n - v)
            congr.append(((lr >> v) * pow(lq >> v, -1, mod) % mod, mod))
        else:
            # residue mod p: brute force over c mod ord_p(q)
            o = 1
            while pow(q, o, p) != 1: o += 1
            res = [c for c in range(o) if pow(q, c, p) == rho % p]
            if not res: return []
            assert len(res) == 1
            congr.append((res[0], o))
            lq = plog(pow(q, p - 1, p**(n + 1)), p, n); lr = plog(pow(rho, p - 1, p**(n + 1)), p, n)
            v = vp(lq, p)
            if lr % p**n and vp(lr, p) < v: return []
            mod = p**(n - v)
            congr.append(((lr // p**v) * pow(lq // p**v, -1, mod) % mod, mod))
    # generalized CRT
    c0, N = 0, 1
    for r_, m_ in congr:
        g = math.gcd(N, m_)
        if (r_ - c0) % g: return []
        l = N // g * m_
        k = ((r_ - c0) // g) * pow(N // g, -1, m_ // g) % (m_ // g)
        c0 = (c0 + N * k) % l; N = l
    if N <= Xc: raise RuntimeError("congruence too weak")
    return [c0] if c0 <= Xc else []

# ------------------------------------------------------------- bootstrapped driver
def prove_bag(bag, base, limit, log=print):
    bag = tuple(sorted(bag)); k = len(bag); b = base.b
    Bs, pb = cofactors(bag, base)
    Bmax = max(Bs)
    known = {}
    t0 = time.time(); st = {'top': 0, 'bottom': 0}
    while True:
        XL, Bc = matveev_partial(k, base, Bmax, known)
        Xc = XL * LN(b) / LN(base.q) + 2
        G = math.ceil((LN(Xc) + 7 * LN(10)) / LN(b)) + 2               # b^G >> Xc * 10^7
        if len(known) == k + 1 and sum(v - 1 for v in known.values()) <= limit:
            Lmax = sum(v - 1 for v in known.values())
            log(f"  base {b} bag {bag}: cofactors {Bs} p-bounds {pb} gap bounds {[known[j] for j in range(k+1)]} "
                f"-> every solution has L <= {Lmax} (top {st['top']}, bottom {st['bottom']}, {time.time()-t0:.1f}s)")
            return Lmax
        # choose an action: top for gaps 0..k-1 (needs all smaller known), bottom for gaps 1..k (all larger known)
        best = None
        for j in range(k + 1):
            if j in known and known[j] <= G + 3: continue
            if j < k and all(i in known for i in range(j)):
                cost = 1
                for i in range(j): cost *= known[i] - 1
                cost *= max(1, len(set(itertools.permutations(bag, j))))
                if best is None or cost < best[0]: best = (cost, j, 'top')
            if all(i in known for i in range(j + 1, k + 1)):
                cost = 5
                for i in range(j + 1, k + 1): cost *= max(known[i] - (0 if i == k else 1), 1)
                cost *= max(1, len(set(itertools.permutations(bag, k - j)))) * len(Bs)
                if best is None or cost < best[0]: best = (cost, j, 'bottom')
        assert best is not None, (known, G)
        _, j, side = best
        if side == 'top':
            ok = top_test(j, bag, base, Bs, XL, G, known, st)
        else:
            ok = bottom_test(j, bag, base, Bs, XL, G, known, st)
        if not ok: raise RuntimeError(f"base {b} bag {bag}: gap {j} ({side}) survives at G={G}")
        known[j] = min(known.get(j, 10**18), G)

def top_test(j, bag, base, Bs, XL, G, known, st):
    b = base.b
    T = Top2(base, XL, G)
    ranges = [range(1, known[i]) for i in range(j)]
    for digs in sorted(set(itertools.permutations(bag, j))):
        cs = [(b - 1) * (d - 1) for d in digs]
        for gaps in itertools.product(*ranges):
            s = sum(gaps); N = b**s; acc = 0
            for c, g in zip(cs, gaps):
                acc += g; N += c * b**(s - acc)
            opz = Fraction(N, b**s)
            for B in Bs:
                st['top'] += 1
                flint.ctx.prec = T.prec
                kappa = flint.arb((b - 1) * B).log() - (flint.arb(N) / flint.arb(b)**s).log()
                surv = T.survivors(kappa, lambda c, L: lam_ok(base, B, c, L, opz, G))
                if surv: return False
    return True

def bottom_test(j, bag, base, Bs, XL, G, known, st):
    b = base.b; k = len(bag)
    Xc = int(XL * LN(b) / LN(base.q)) + 2
    nb = k - j                                       # digits below gap j
    for digs in sorted(set(itertools.permutations(bag, nb))):    # top-to-bottom d_{j+1}..d_k
        # gaps below: g_{j+1}..g_{k-1} (>=1) and trailing g_k (>=0)
        rngs = [range(1, known[i]) for i in range(j + 1, k)] + ([range(0, known[k])] if j < k else [])
        for gaps in itertools.product(*rngs):
            # positions from the bottom: a_k = g_k, a_{k-1} = a_k + g_{k-1}, ...
            pos = []
            if nb:
                a = gaps[-1]; pos = [a]
                for g in reversed(gaps[:-1]):
                    a += g; pos.append(a)
                pos = pos[::-1]                              # aligned with digs (top to bottom)
                M = pos[0] + G
            else:
                M = G
            r = (b**M - 1) // (b - 1) + sum((d - 1) * b**p for d, p in zip(digs, pos))
            for B in Bs:
                st['bottom'] += 1
                sol = solve_c(base, B, r % b**M, M, Xc)
                if sol is None: return False
                if sol:
                    print("    bottom candidate", digs, pos, B, sol); return False
    return True

# ------------------------------------------------------------- whole-base driver
def prove_base(b, D, log=print):
    from basetrees import trees, bags_of, fprod
    base = Base(b)
    n, nz, nodes = trees(b, D)
    zf_all = sorted({y for d in nodes for y in nodes[d]})          # zero-free nodes below b^D
    # every zero-free smooth number below b^D (not only nonzero-target ones) for preimage lookup
    from basetrees import smooth_not_div
    zf_list = [y for y in smooth_not_div(b, D) if '0' not in flint.fmpz(y).str(b)]
    result = {}
    for d in sorted(nodes):
        for s in sorted(nodes[d]):
            bs = bags_of(s, b)
            Lb = 0
            for bag in bs:
                Lb = max(Lb, prove_bag(bag, base, D - 1, log=log))
            pre = sorted(y for y in zf_list if y >= b and fprod(y, b) == s)
            result[s] = {'target': d, 'depth': nodes[d][s], 'preimages': pre, 'Lbound': Lb}
            log(f"base {b}: node {s} (target {d}): all smooth y with digit product {s}: {pre}")
    # closure: every preimage of a node is itself a node of the same tree
    for s, v in result.items():
        for y in v['preimages']:
            assert y in result and result[y]['target'] == v['target'], (s, y)
    return result

if __name__ == "__main__":
    b, D = int(sys.argv[1]), int(sys.argv[2])
    res = prove_base(b, D)
    json.dump({str(k): v for k, v in res.items()}, open(f'base{b}_proof.json', 'w'))
    for d in sorted({v['target'] for v in res.values()}):
        H = max(v['depth'] for v in res.values() if v['target'] == d)
        print(f"THEOREM base {b}, final digit {d}: tree complete ({sum(1 for v in res.values() if v['target']==d)} nodes); "
              f"max persistence {H+1}")
