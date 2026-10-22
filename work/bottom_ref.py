"""Reference (Python, exact) of the cumulative bottom solver csrc/bottom3.c ('cum' mode), in the restructured form that the
GPU kernel uses, plus a direct transcription of bottom3.c's formulas to check the restructuring against.

The restructuring (no exact divisions, every value a residue modulo the final modulus):
  log1p(x) with x = p^v y:  sum_{n >= 1} (-1)^(n+1) p^(n v - e) y^n / n'  mod p^M,   n = p^e n'  (p ∤ n'),
  over the n with n v - e < M (the other terms vanish mod p^M). bottom3.c sums x^n / p^e mod p^(M+12) instead; both are the
  exact value of the convergent series mod p^M. plog(x) = log(x^(p^s)) / p^s is treated the same way with z = x^(p^s) - 1.
usage (from work/):
  python3 bottom_ref.py selfcheck [N]            random checks: restructured == bottom3.c formulas; planted (t, u', w)
  python3 bottom_ref.py job JOBFILE              run a bottomc3 'cum' job in Python; prints CAND/DONE lines like bottomc3
"""
import sys, math, random

DL5 = (-1, 0, 1, 3, 2)                       # dlog base 2 mod 5: 1 -> 0, 2 -> 1, 4 -> 2, 3 -> 3

def v_p(x, p):
    if x == 0: return 10**9
    e = 0
    while x % p == 0: x //= p; e += 1
    return e

def floorlog(n, p):
    lg, q = 0, p
    while q <= n: q *= p; lg += 1
    return lg

# ---------------------------------------------------------------- transcription of bottom3.c
def log1p_series_c(x, p, prec, vx):
    assert vx >= (1 if p > 2 else 2) and x % p**vx == 0
    out = 0; xp = 1; modx = p**(prec + 12); n = 1
    while True:
        lg = floorlog(n, p)
        if n * vx - lg >= prec + 2 and n > 1: break
        e = v_p(n, p); nn = n // p**e
        assert e <= 12
        xp = xp * x % modx
        term = xp // p**e * pow(nn, -1, p**prec) % p**prec
        out = (out + term) % p**prec if n & 1 else (out - term) % p**prec
        n += 1
    return out

def plog_c(x, p, M):
    s = math.isqrt(M) + 2; N = M + s; modN = p**N
    X = pow(x % modN, p**s, modN); z = (X - 1) % modN
    vz = s + (1 if p > 2 else 2)
    assert z % p**vz == 0
    total = 0; zn = 1; modext = p**(N + 12); n = 1
    while True:
        lg = floorlog(n, p)
        if n * vz - lg >= N + 5 and n > 2: break
        e = v_p(n, p); nn = n // p**e
        zn = zn * z % modext
        term = zn // p**e * pow(nn, -1, modN) % modN
        total = (total + term) % modN if n & 1 else (total - term) % modN
        n += 1
    assert total % p**s == 0
    return total // p**s % p**M

# ---------------------------------------------------------------- restructured (the GPU form)
def log1p_r(y, p, v, M):
    """log1p(p^v y) mod p^M, as sum (-1)^(n+1) p^(n v - e) y^n n'^-1 over n v - e < M (no exact division)."""
    m = p**M; out = 0; yn = 1; n = 1
    while True:
        e = v_p(n, p)
        if n * v - floorlog(n, p) >= M: break            # every later term has valuation n v - e >= n v - lg >= M
        yn = yn * y % m
        if n * v - e < M:
            term = pow(p, n * v - e, m) * yn % m * pow(n // p**e, -1, m) % m
            out = (out + term) % m if n & 1 else (out - term) % m
        n += 1
    return out

def plog_r(x, p, M):
    """plog(x) mod p^M (x = 1 mod p, or mod 4 for p = 2) = log(x^(p^s)) / p^s; z = x^(p^s) - 1 = p^vz z'."""
    s = math.isqrt(M) + 2; N = M + s
    X = pow(x % p**N, p**s, p**N)
    vz = s + (1 if p > 2 else 2)
    z = (X - 1) % p**N
    assert z % p**vz == 0
    zp = z // p**vz                                      # known mod p^(N - vz) = p^(M + s - vz): enough, see below
    m = p**M; out = 0; zn = 1; n = 1
    while True:
        e = v_p(n, p)
        if n * vz - floorlog(n, p) - s >= M: break
        zn = zn * zp % m
        k = n * vz - e - s                               # >= 1: valuation of the term after the division by p^s
        if k < M:
            term = pow(p, k, m) * zn % m * pow(n // p**e, -1, m) % m
            out = (out + term) % m if n & 1 else (out - term) % m
        n += 1
    return out

# ---------------------------------------------------------------- the solver
class Ctx:
    def __init__(self, G, tm, U, W, Mmax, logs):
        self.G, self.tm, self.U, self.W, self.Mmax = G, tm, U, W, Mmax
        self.L5_3, self.L5_7, self.L5_16, self.L2_3, self.L2_7 = logs
        self.mtc = {}
    def mt(self, M, t):
        key = (M, t)
        if key in self.mtc: return self.mtc[key]
        n = M - t
        assert n >= 6
        q2, q5 = 2**(n - 2), 5**(M - 1); Q = q2 * q5
        A = (self.L2_3 % 2**n) >> 2; Bc = (self.L2_7 % 2**n) >> 2
        A5 = (self.L5_3 % 5**M) // 5; B5 = (self.L5_7 % 5**M) // 5
        Ainv = pow(A, -1, q2); A5inv = pow(A5, -1, q5)
        alpha2 = (-Bc * Ainv) % q2; alpha5 = (-B5 * A5inv) % q5
        e2 = q5 * pow(q5, -1, q2); e5 = q2 * pow(q2, -1, q5)
        alpha = (alpha2 * e2 + alpha5 * e5) % Q
        b1 = [4 * alpha % Q, 4]; b2 = [Q, 0]
        while True:
            n1 = b1[0]**2 + b1[1]**2; n2 = b2[0]**2 + b2[1]**2
            if n1 > n2: b1, b2 = b2, b1; n1, n2 = n2, n1
            num = b1[0] * b2[0] + b1[1] * b2[1]
            k = (2 * num + n1) // (2 * n1)
            if k == 0: break
            b2 = [b2[0] - k * b1[0], b2[1] - k * b1[1]]
        assert n1 > self.U**2 + self.W**2, "lattice too fine"
        r = dict(q2=q2, q5=q5, Q=Q, Ainv=Ainv, A5inv=A5inv, alpha2=alpha2, alpha=alpha, e2=e2, e5=e5,
                 b1=b1, b2=b2, det=b1[0] * b2[1] - b1[1] * b2[0], q2inv5=pow(q2, -1, q5))
        self.mtc[key] = r
        return r

    def low(self, Flow, restructured=True):
        """per-configuration data (bottom3.c low()); restructured: only the precision the solves need (mod 5^M, 2^M)"""
        tl = v_p(Flow, 2)
        odd = Flow >> tl if Flow > 0 else -((-Flow) >> tl)
        sig = 1 if odd % 4 == 1 else -1
        P = self.G if restructured else self.Mmax
        l2low = (plog_r if restructured else plog_c)((sig * odd) % 2**(P + 2), 2, P)
        l5low = (plog_r if restructured else plog_c)(pow(Flow % 5**(P + 1), 4, 5**(P + 1)), 5, P)
        inv5 = pow(Flow % 5**P, -1, 5**P)
        inv2 = pow(odd % 2**(P + 2), -1, 2**(P + 2))
        return dict(Flow=Flow, tl=tl, sig=sig, l2low=l2low, l5low=l5low, inv5=inv5, inv2=inv2)

    def beta(self, lw, c, p, M, restructured=True):
        """returns None (no solution possible before the lattice step) or (t, sigma, F, beta2, beta, r)"""
        Flow = lw['Flow']; F = c * 10**p + Flow
        t = v_p(F, 2)
        if t >= M or t > self.tm or F % 5 == 0: return None
        n = M - t
        if p >= 1:
            if restructured:
                y = c * pow(2, p, 5**M) * lw['inv5'] % 5**M
                l5F = (4 * log1p_r(y, 5, p, M) + lw['l5low']) % 5**M
            else:
                x5 = c * 10**p * lw['inv5'] % 5**(M + 10)
                l5F = (4 * log1p_series_c(x5, 5, M, p) + lw['l5low']) % 5**M
        else:
            l5F = (plog_r if restructured else plog_c)(pow(F % 5**(M + 1), 4, 5**(M + 1)), 5, M)
        vx2 = p + v_p(c, 2) - lw['tl']
        if lw['tl'] == t and vx2 >= 2:
            if restructured:
                codd = c >> v_p(c, 2)
                y2 = codd * pow(5, p, 2**n) * lw['inv2'] % 2**n
                l2F = (log1p_r(y2, 2, vx2, n) + lw['l2low']) % 2**n
            else:
                x2 = (c * 10**p >> lw['tl']) * lw['inv2'] % 2**(n + 40)
                l2F = (log1p_series_c(x2, 2, n, vx2) + lw['l2low']) % 2**n
            sigma = lw['sig']
        else:
            T2 = (F >> t) % 2**(n + 2)
            sigma = 1 if T2 % 4 == 1 else -1
            l2F = (plog_r if restructured else plog_c)(sigma * T2 % 2**(n + 2), 2, n)
        r = self.mt(M, t)
        T = (l5F - t * self.L5_16) % 5**M
        assert T % 5 == 0
        D5 = T // 5; D2 = (l2F % 2**n) >> 2
        beta5 = D5 * r['A5inv'] % r['q5']; beta2 = D2 * r['Ainv'] % r['q2']
        if restructured:     # Garner: beta = beta2 + q2 ((beta5 - beta2) q2^-1 mod q5)
            beta = beta2 + r['q2'] * ((beta5 - beta2) * r['q2inv5'] % r['q5'])
        else:
            beta = (beta2 * r['e2'] + beta5 * r['e5']) % r['Q']
        return t, sigma, F, beta2, beta, r

    def candidates(self, t, sigma, F, beta2, beta, r):
        rho5 = F * pow(3, t, 5) % 5
        a2m4, b2m4 = r['alpha2'] % 4, beta2 % 4
        b1x, b1y = r['b1']; b2x, b2y = r['b2']; tq = 2 * r['det']; out = []
        for w0 in range(4):
            u0 = (a2m4 * w0 + b2m4) % 4
            if (u0 + w0) % 2 != (0 if sigma == 1 else 1): continue
            if (3 * u0 + w0) % 4 != DL5[rho5]: continue
            off0 = (r['alpha'] * w0 + beta) % r['Q']
            cx2 = self.U - 2 * off0; cy2 = self.W - 2 * w0
            ci = (cx2 * b2y - cy2 * b2x) // tq; cj = (b1x * cy2 - b1y * cx2) // tq
            for di in (0, 1):
                for dj in (0, 1):
                    up = off0 + (ci + di) * b1x + (cj + dj) * b2x; ww = w0 + (ci + di) * b1y + (cj + dj) * b2y
                    if 2 <= up <= self.U and 0 <= ww <= self.W: out.append((t, up, ww))
        return out

def read_job(path):
    tok = open(path).read().split(); i = 0
    def nxt():
        nonlocal i; i += 1; return tok[i - 1]
    mode = nxt(); assert mode == 'cum'
    kv = {}
    for name in ('G', 'tm', 'U', 'W', 'Mmax', 'L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7', 'k', 'j'):
        assert nxt() == name; kv[name] = int(nxt())
    assert nxt() == 'heff'; b = kv['k'] - kv['j']; heff = [int(nxt()) for _ in range(b)]
    assert nxt() == 'nord'; nord = int(nxt())
    orders = [[int(nxt()) for _ in range(b)] for _ in range(nord)]
    for name in ('skip', 'aklo', 'akhi'):
        assert nxt() == name; kv[name] = int(nxt())
    return kv, heff, orders

def run_job(path):
    kv, heff, orders = read_job(path)
    ctx = Ctx(kv['G'], kv['tm'], kv['U'], kv['W'], kv['Mmax'], (kv['L5_3'], kv['L5_7'], kv['L5_16'], kv['L2_3'], kv['L2_7']))
    G = kv['G']; b = len(heff); nsolve = nd5 = nd2 = 0
    for o in range(kv['skip'], len(orders)):
        cc = orders[o]; pos = [0] * b
        def rec(r):
            nonlocal nsolve, nd5, nd2
            if r == 0:
                lo = 0 if b == 1 else pos[1] + 1; phi = heff[0]
                if b == 1: lo = max(lo, kv['aklo']); phi = min(phi, kv['akhi'] + 1)
                if lo >= phi: return
                Flow = -1 + sum(cc[q] * 10**pos[q] for q in range(1, b))
                lw = ctx.low(Flow)
                for p in range(lo, phi):
                    M = max(p + 1, G); nsolve += 1
                    res = ctx.beta(lw, cc[0], p, M)
                    if res is None: continue                   # bottom3.c counts direct paths only past the early exit
                    t, sigma, F = res[0], res[1], res[2]
                    if p == 0: nd5 += 1
                    if not (lw['tl'] == t and p + v_p(cc[0], 2) - lw['tl'] >= 2): nd2 += 1
                    for (tt, up, ww) in ctx.candidates(*res):
                        print(f"CAND {o} pos {' '.join(str(pos[q]) for q in range(1, b))} p {p} M {M} t {tt} up {up} w {ww}")
                return
            lo = 0 if r == b - 1 else pos[r + 1] + 1; hi = heff[r]
            if r == b - 1: lo = max(lo, kv['aklo']); hi = min(hi, kv['akhi'] + 1)
            for a in range(lo, hi):
                pos[r] = a; rec(r - 1)
        rec(b - 1)
        print(f"PROG {o + 1} {nsolve}")
    print(f"DONE {nsolve} {nd5} {nd2}")

def selfcheck(N=300, seed=1):
    """restructured == transcription on random inputs; planted (t, u', w) are recovered exactly"""
    rnd = random.Random(seed)
    for G in (20, 37, 61):
        tm = 12; XL = 10**6; U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2); Mmax = 4 * G + 8
        from cbottom import consts
        cs = consts(Mmax)
        ctx = Ctx(G, tm, U, W, Mmax, (cs['L5_3'], cs['L5_7'], cs['L5_16'], cs['L2_3'], cs['L2_7']))
        for p_ in (2, 5):
            for M in (G, G + 3):
                for _ in range(N // 6):
                    x = rnd.randrange(5**M) * (p_ if p_ == 5 else 4) + 1
                    assert plog_r(x, p_, M) == plog_c(x, p_, M), ("plog", p_, M)
                    v = rnd.randrange(1 if p_ == 5 else 2, 9); y = rnd.randrange(1, 5**M)
                    assert log1p_r(y, p_, v, M) == log1p_series_c(p_**v * y % p_**(M + 10), p_, M, v), ("log1p", p_, v, M)
        nb = nc = 0
        for _ in range(N):
            k = rnd.randrange(1, 4); cc = [9 * (rnd.choice((2, 3, 4, 6, 7, 8, 9)) - 1) for _ in range(k + 1)]
            pos = sorted(rnd.sample(range(0, G - 1), k))
            Flow = -1 + sum(cc[q + 1] * 10**pos[q] for q in range(k))
            if Flow % 5 == 0: continue
            lr, lc = ctx.low(Flow, True), ctx.low(Flow, False)
            for key in ('tl', 'sig'): assert lr[key] == lc[key]
            for p in range(max(pos) + 1, G):
                a = ctx.beta(lr, cc[0], p, G, True); b_ = ctx.beta(lc, cc[0], p, G, False)
                assert (a is None) == (b_ is None)
                if a is None: continue
                assert a[:4] == b_[:4] and a[4] == b_[4], "beta"
                assert ctx.candidates(*a) == ctx.candidates(*b_)
                nb += 1
        # planted: F = 2^t 3^u' 7^w mod 10^G with (u', w) in the box must give exactly that candidate
        for _ in range(N):
            t = rnd.randrange(0, tm + 1); up = rnd.randrange(2, U + 1); w = rnd.randrange(0, W + 1)
            F = pow(2, t, 10**G) * pow(3, up, 10**G) * pow(7, w, 10**G) % 10**G
            # write F = c 10^p + Flow with a digit-like split: p = G - 2, c = F div 10^p, Flow = F mod 10^p (any split works)
            p = G - 2; c = F // 10**p; Flow = F % 10**p
            if c == 0 or Flow == 0 or Flow % 5 == 0: continue
            lw = ctx.low(Flow)
            res = ctx.beta(lw, c, p, G)
            assert res is not None and res[0] == t
            cands = ctx.candidates(*res)
            assert (t, up, w) in cands, ("planted not found", t, up, w, cands)
            nc += 1
        print(f"G={G}: plog/log1p restructured == bottom3.c; {nb} (Flow, p) betas identical; planted {nc}/{nc} found")
    print("SELFCHECK OK")

if __name__ == '__main__':
    if sys.argv[1] == 'selfcheck': selfcheck(int(sys.argv[2]) if len(sys.argv) > 2 else 300)
    elif sys.argv[1] == 'job': run_job(sys.argv[2])
