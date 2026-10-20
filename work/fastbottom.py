"""Fast exact 10-adic test for bottom configurations.

For a bottom configuration (digits d_i at positions p_i, i = top..bottom of the configuration),
with a run of >= G ones above the highest one, y mod 10^M is known (M = p_top + G), and
    9y = F (mod 10^M),  F = -1 + sum_i 9(d_i - 1) 10^{p_i}     (since 10^M = 0 mod 10^M).
So 2^t 3^{u'} 7^w = F (mod 10^M) with u' = u + 2.  We decide whether any (t,u',w) with
t <= tm, 2 <= u' <= U, 0 <= w <= W satisfies this, using 2-adic/5-adic logarithms.
Speed-up: F = F_low + c 10^p (p = topmost position), log F = log F_low + log(1 + c 10^p / F_low).
"""
import math
from fractions import Fraction
from solver import plog, v2

DL5 = {1: 0, 2: 1, 4: 2, 3: 3}

class BottomContext:
    def __init__(self, G, tm, U, W):
        self.G, self.tm, self.U, self.W = G, tm, U, W
        self.Mmax = 4 * G + 8
        Mm = self.Mmax
        self.p5 = 5**Mm; self.p2 = 2**(Mm + 2)
        self.L5_3 = plog(pow(3, 4, 5**(Mm + 1)), 5, Mm)
        self.L5_7 = plog(pow(7, 4, 5**(Mm + 1)), 5, Mm)
        self.L5_16 = plog(16, 5, Mm)
        self.L2_3 = plog(-3 % 2**(Mm + 2), 2, Mm)
        self.L2_7 = plog(-7 % 2**(Mm + 2), 2, Mm)
        self._mt = {}
        self._low = {}

    def mt(self, M, t):
        key = (M, t)
        if key in self._mt: return self._mt[key]
        n = M - t; q2 = 2**(n - 2); q5 = 5**(M - 1); Q = q2 * q5
        A = (self.L2_3 % 2**n) >> 2; Bc = (self.L2_7 % 2**n) >> 2
        A5 = (self.L5_3 % 5**M) // 5; B5 = (self.L5_7 % 5**M) // 5
        Ainv = pow(A, -1, q2); A5inv = pow(A5, -1, q5)
        alpha2 = (-Bc * Ainv) % q2; alpha5 = (-B5 * A5inv) % q5
        e2 = q5 * pow(q5, -1, q2); e5 = q2 * pow(q2, -1, q5)
        alpha = (alpha2 * e2 + alpha5 * e5) % Q
        b1, b2 = [(4 * alpha) % Q, 4], [Q, 0]
        n2 = lambda v: v[0] * v[0] + v[1] * v[1]
        while True:
            if n2(b1) > n2(b2): b1, b2 = b2, b1
            num = b1[0] * b2[0] + b1[1] * b2[1]; den = n2(b1)
            k = (2 * num + den) // (2 * den)          # nearest integer to num/den
            if k == 0: break
            b2 = [b2[0] - k * b1[0], b2[1] - k * b1[1]]
        diam2 = self.U * self.U + self.W * self.W
        assert n2(b1) > diam2, "lattice too fine"
        det = b1[0] * b2[1] - b1[1] * b2[0]
        rec = dict(n=n, q2=q2, q5=q5, Q=Q, Ainv=Ainv, A5inv=A5inv, alpha2=alpha2, alpha=alpha,
                   e2=e2, e5=e5, b1=b1, b2=b2, det=det)
        self._mt[key] = rec
        return rec

    def low(self, Flow):
        """Cache logs of a fixed lower part F_low (a unit at 5; 2-part separated)."""
        if Flow in self._low: return self._low[Flow]
        Mm = self.Mmax
        tl = v2(Flow)
        odd = (Flow >> tl) if Flow > 0 else -((-Flow) >> tl)
        sig = 1 if odd % 4 == 1 else -1
        l2 = plog((sig * odd) % 2**(Mm + 2), 2, Mm)
        l5 = plog(pow(Flow % 5**(Mm + 1), 4, 5**(Mm + 1)), 5, Mm)
        inv5 = pow(Flow % self.p5, -1, self.p5)
        inv2 = pow(odd % self.p2, -1, self.p2)
        rec = (tl, sig, l2, l5, inv5, inv2)
        self._low[Flow] = rec
        return rec

    @staticmethod
    def _log1p_series(x, p, prec, vx):
        """log(1+x) mod p^prec for integer x with v_p(x) >= vx >= 1 (>= 2 if p == 2)."""
        assert vx >= (1 if p > 2 else 2) and x % p**vx == 0
        mod = p**prec; total = 0; n = 1; xp = 1
        while True:
            e = 0; nn = n
            while nn % p == 0: nn //= p; e += 1
            lg = 0
            while p**(lg + 1) <= n: lg += 1          # floor(log_p n) >= v_p(m) for all m <= n... used as safe bound
            if n * vx - lg >= prec + 2 and n > 1: break
            xp_full = pow(x, n, p**(prec + e))
            term = (xp_full // p**e) * pow(nn, -1, mod) % mod
            total = (total + term) % mod if n % 2 else (total - term) % mod
            n += 1
        return total

    def solutions(self, Flow, c, p, M):
        """All (t,u',w) in the box with 2^t 3^u' 7^w = Flow + c*10^p (mod 10^M)."""
        assert M <= self.Mmax, "p-adic precision exceeded (M > Mmax)"
        F = Flow + c * 10**p
        t = v2(F)
        if t >= M or t > self.tm or F % 5 == 0: return []
        tl, sig, l2low, l5low, inv5, inv2 = self.low(Flow)
        n = M - t
        # 5-adic log of F^4
        if p >= 1:
            x5 = c * 10**p * inv5 % 5**(M + 10)
            l5F = (l5low + 4 * self._log1p_series(x5, 5, M, p)) % 5**M
        else:                                          # no 5-adic smallness: direct
            l5F = plog(pow(F % 5**(M + 1), 4, 5**(M + 1)), 5, M)
        # 2-adic: need v2(F) = v2(Flow) and v2(x2) >= 2
        vx2 = p + v2(c) - tl
        if tl == t and vx2 >= 2:
            x2 = (c * 10**p >> tl) * inv2 % 2**(n + 40)
            l2F = (l2low + self._log1p_series(x2, 2, n, vx2)) % 2**n
            sigma = sig
        else:                                         # rare: direct computation
            odd = (F >> t) % 2**(n + 2)
            sigma = 1 if odd % 4 == 1 else -1
            l2F = plog((sigma * odd) % 2**(n + 2), 2, n)
        rec = self.mt(M, t)
        D5 = ((l5F - t * self.L5_16) % 5**M) // 5
        D2 = (l2F % 2**n) >> 2
        beta5 = D5 * rec['A5inv'] % rec['q5']; beta2 = D2 * rec['Ainv'] % rec['q2']
        Q = rec['Q']; beta = (beta2 * rec['e2'] + beta5 * rec['e5']) % Q
        rho5 = F * pow(2, -t, 5) % 5
        sols = []
        b1, b2, det = rec['b1'], rec['b2'], rec['det']
        for w0 in range(4):
            u0 = (rec['alpha2'] * w0 + beta2) % 4
            if (u0 + w0) % 2 != (0 if sigma == 1 else 1): continue
            if (3 * u0 + w0) % 4 != DL5[rho5]: continue
            off = ((rec['alpha'] * w0 + beta) % Q, w0)
            cx = Fraction(self.U, 2) - off[0]; cy = Fraction(self.W, 2) - off[1]
            ci = (cx * b2[1] - cy * b2[0]) / det; cj = (b1[0] * cy - b1[1] * cx) / det
            fi, fj = math.floor(ci), math.floor(cj)
            for i in range(fi - 2, fi + 4):
                for j in range(fj - 2, fj + 4):
                    up = off[0] + i * b1[0] + j * b2[0]; w = off[1] + i * b1[1] + j * b2[1]
                    if 2 <= up <= self.U and 0 <= w <= self.W:
                        sols.append((t, up, w))
        return sols

if __name__ == "__main__":
    import random, time
    random.seed(3)
    G, tm = 60, 15
    U, W = 10**25, 6 * 10**24
    ctx = BottomContext(G, tm, U, W)
    # planted solutions: pick (t,u',w), compute F = 2^t 3^u' 7^w mod 10^M, split as Flow + c*10^p
    ok = 0; tot = 0; t0 = time.time()
    for trial in range(300):
        t = random.randint(0, tm); up = random.randint(2, U); w = random.randint(0, W)
        p = random.randint(2, 3 * G); M = p + G
        Fm = pow(2, t, 10**M) * pow(3, up, 10**M) * pow(7, w, 10**M) % 10**M
        c = random.choice([9, 18, 27, 45, 54, 63, 72])
        Flow = (Fm - c * 10**p) % 10**M          # so that Flow + c*10^p = Fm (mod 10^M)
        if Flow % 5 == 0: continue
        tot += 1
        if (t, up, w) in ctx.solutions(Flow, c, p, M): ok += 1
    print(f"planted recovered {ok}/{tot} in {time.time()-t0:.1f}s")
