"""Rigorous solver for: find all y = 2^t 3^u 7^w whose decimal digits are a fixed multiset
`bag` of non-1 digits (from 2,3,4,6,7,8,9) plus any number of 1s.

Notation: y has L digits; non-1 digits d_1..d_k at positions a_1 > ... > a_k;
gaps g_0 = L - a_1, g_i = a_i - a_{i+1}; c_i = 9(d_i - 1).
    9y = 10^L - 1 + sum c_i 10^{a_i}.
Facts used (each proved in the write-up):
  (F1) a_k <= 1   (3^u 7^w mod 100 has even tens digit; even y has even last digit)
  (F2) t <= tmax(bag)        (twoadic.py)
  (F3) for 0 <= j < k:  0 < Lam_j < 80*10^(-g_j), where
        Lam_j = t log2 + (u+2) log3 + w log7 - L log10 - log(1+z_j),
        z_j = sum_{i<=j} c_i 10^{-(g_0+...+g_{i-1})}   (z_0 = 0)
  (F4) Matveev's theorem bounds L by X_L (matveev_bound).
Case split with threshold G: if some gap >= G, let j_top/j_bot be the first/last such gap.
  top case  (j_top <= JT): lattice test on Lam_{j_top} with known z_{j_top}.
  bottom case (j_bot > JT): y mod 10^M known, M = a_{j_bot+1} + G; 2-adic & 5-adic logs.
If every case is impossible, all gaps < G, so L <= k(G-1)+1.
"""
import math, itertools
from fractions import Fraction
import flint
from twoadic import tmax, v2

# ---------------------------------------------------------------- Matveev (real, D=1)
def matveev_bound(k, tm):
    """Return (X_L, B): rigorous upper bounds for L and for max coefficient B."""
    ln = math.log
    C4 = 1.4 * 30**7 * 4**4.5
    C5 = 1.4 * 30**8 * 5**4.5
    Pi = ln(2) * ln(3) * ln(7) * ln(10) * (1 + 1e-12)
    def f(B):
        beta = 1 + ln(B)
        s = (ln(80) + C4 * beta * Pi) / ln(10)                 # g_0 bound
        for j in range(1, k):
            g = (ln(80) + C5 * beta * Pi * (s * ln(10) + ln(9))) / ln(10)
            s += g
        Lmax = s + 1
        return max(2.0960 * Lmax + 2.2, tm + 1), Lmax
    B = 1e300
    for _ in range(200):
        nb, _ = f(B)
        if nb >= B * (1 - 1e-9): break
        B = nb
    B *= 1.01                                                  # safety margin
    nb, Lmax = f(B)
    assert nb <= B
    return math.ceil(Lmax * 1.01), B

# ---------------------------------------------------------------- top-side lattice
class TopLattice:
    def __init__(self, XL, T, G):
        flint.ctx.prec = 2400
        self.XL, self.T, self.G = XL, T, G
        U = int(2.0960 * XL + 3)              # bound for u' = u+2
        W = int(1.1834 * XL + 2)              # bound for w
        S = U
        lw = max(1, round(S / W)); lt = max(1, S // max(T, 1))
        self.S_t, self.S_u, self.S_w = lt * T, U, lw * W
        C = S * 10**G // 80 + 1
        self.C = C
        flint.ctx.prec = max(2400, C.bit_length() + 200)      # enough for C*log p to well within 0.01
        self.rC = {}
        for p in (2, 3, 7, 10):
            v = flint.arb(C) * flint.arb(p).log()
            m = int(v.mid().floor().unique_fmpz())
            if v - m > 0.5: m += 1
            assert abs(v - m) < 0.51
            self.rC[p] = m
        rows = [[lt, 0, 0, self.rC[2]],
                [0, 1, 0, self.rC[3]],
                [0, 0, lw, self.rC[7]],
                [0, 0, 0, -self.rC[10]]]
        if T == 0:                             # t fixed to 0: 3-dimensional lattice
            rows = [r[1:] for r in rows[1:]]
        self.dim = len(rows)
        Bm = flint.fmpz_mat(rows).lll()
        self.B = [[int(Bm[i, j]) for j in range(self.dim)] for i in range(self.dim)]
        # exact Gram-Schmidt norms and inverse
        Bq = [[Fraction(x) for x in r] for r in self.B]
        star = []
        for i in range(self.dim):
            v = Bq[i][:]
            for s in star:
                mu = sum(a*b for a, b in zip(Bq[i], s)) / sum(a*a for a in s)
                v = [a - mu*b for a, b in zip(v, s)]
            star.append(v)
        self.gs2 = [sum(a*a for a in s) for s in star]       # squared GS norms
        # error budget
        E = Fraction(51, 100) * (T + U + W + XL + 2)
        self.R2 = self.S_t**2 + self.S_u**2 + self.S_w**2 + (S + E + 1)**2   # squared radius
        self.min_gs2 = min(self.gs2)

    def kills_homogeneous(self):
        # shortest nonzero vector >= min GS norm
        return self.min_gs2 > self.R2

    def kills(self, m):
        """Target y = (0,..,0,-m). Return True if no lattice vector lies within R of y."""
        # z = y B^{-1}; z_last = -m * (B^{-1})[last, last] ... computed via last column of B^{-1}
        q = self.Binv_last_col_last()
        z = -m * q
        frac = z - math.floor(z)
        d = min(frac, 1 - frac)
        return d * d * self.gs2[-1] > self.R2

    def Binv_last_col_last(self):
        if not hasattr(self, "_q"):
            Binv = flint.fmpq_mat(flint.fmpz_mat(self.B)).inv()
            n = self.dim
            self._q = Fraction(int(Binv[n-1, n-1].p), int(Binv[n-1, n-1].q))
        return self._q

    def target_m(self, one_plus_z: Fraction):
        """m = nearest integer to C*kappa, kappa = -log(1+z); rigorous within 0.51."""
        v = flint.arb(self.C) * (-(flint.arb(one_plus_z.numerator) / one_plus_z.denominator).log())
        m = int(v.mid().floor().unique_fmpz())
        if v - m > 0.5: m += 1
        assert abs(v - m) < 0.51
        return m

# ---------------------------------------------------------------- p-adic logarithms
def plog(x, p, M):
    """log_p(x) mod p^M for x = 1 mod p (p odd) or x = 1 mod 4 (p = 2). Exact integer result."""
    s = math.isqrt(M) + 2
    N = M + s
    modN = p**N
    X = pow(x, p**s, modN)
    z = (X - 1) % modN
    assert z % p**(s + (1 if p > 2 else 2)) == 0
    total = 0
    n = 1
    vz = s + (1 if p > 2 else 2)               # lower bound for v_p(z)
    while True:
        e = 0; nn = n
        while nn % p == 0: nn //= p; e += 1
        lg = 0
        while p**(lg + 1) <= n: lg += 1
        if n * vz - lg >= N + 5 and n > 2: break
        term = pow(z, n, p**(N + e)) // p**e    # exact: v(z^n) >= n*vz >= e
        term = term * pow(nn, -1, modN) % modN
        total = (total + term) % modN if n % 2 == 1 else (total - term) % modN
        n += 1
    assert total % p**s == 0
    return (total // p**s) % p**M

_logcache = {}
def unit_logs(M):
    """Precompute logs of 3 and 7 for modulus parameters M (5-adic) and 2-adic up to M."""
    if M not in _logcache:
        L5_3 = plog(pow(3, 4, 5**(M+1)), 5, M); L5_7 = plog(pow(7, 4, 5**(M+1)), 5, M)
        L2_3 = plog(-3 % 2**(M+2), 2, M); L2_7 = plog(-7 % 2**(M+2), 2, M)
        _logcache[M] = (L5_3, L5_7, L2_3, L2_7)
    return _logcache[M]

DL5 = {1: 0, 2: 1, 4: 2, 3: 3}                 # discrete log base 2 mod 5

def solve_bottom(r, M, tm, Umax, Wmax):
    """All (t,u,w) with t<=tm, 0<=u<=Umax, 0<=w<=Wmax and 2^t 3^u 7^w = r (mod 10^M).
    Returns list of (t,u,w); raises if the box could contain more than one point per coset
    without being enumerated."""
    t = v2(r)
    if t >= M or t > tm or r % 5 == 0: return []
    n = M - t
    assert n >= 6
    rho2 = (r >> t) % 2**n
    rho5 = r * pow(2, -t, 5**M) % 5**M
    L5_3, L5_7, L2_3, L2_7 = unit_logs(M)
    # 2-adic
    sigma = 1 if rho2 % 4 == 1 else -1
    l_rho2 = plog((sigma * rho2) % 2**(n+2), 2, n)
    A, Bc, Dd = (L2_3 % 2**n) >> 2, (L2_7 % 2**n) >> 2, l_rho2 >> 2      # mod 2^(n-2)
    assert A % 2 == 1 and (L2_3 % 4 == 0) and (l_rho2 % 4 == 0)
    q2 = 2**(n-2)
    alpha2 = (-Bc * pow(A, -1, q2)) % q2; beta2 = (Dd * pow(A, -1, q2)) % q2
    # 5-adic
    l_rho5 = plog(pow(rho5, 4, 5**(M+1)), 5, M)
    A5, B5, D5 = L5_3 // 5, L5_7 // 5, l_rho5 // 5                          # mod 5^(M-1)
    assert L5_3 % 5 == 0 and L5_7 % 5 == 0 and l_rho5 % 5 == 0 and A5 % 5 != 0
    q5 = 5**(M-1)
    alpha5 = (-B5 * pow(A5, -1, q5)) % q5; beta5 = (D5 * pow(A5, -1, q5)) % q5
    # CRT: u = alpha*w + beta (mod Q)
    Q = q2 * q5
    def crt(a2, a5):
        return (a2 * q5 * pow(q5, -1, q2) + a5 * q2 * pow(q2, -1, q5)) % Q
    alpha, beta = crt(alpha2, alpha5), crt(beta2, beta5)
    sols = []
    for w0 in range(4):
        u0 = (alpha2 * w0 + beta2) % 4
        if (u0 + w0) % 2 != (0 if sigma == 1 else 1): continue
        if (3 * u0 + w0) % 4 != DL5[rho5 % 5]: continue
        # coset: (u,w) = (alpha*w0+beta, w0) + i*(4*alpha mod Q, 4) + j*(Q, 0)
        off = ((alpha * w0 + beta) % Q, w0)
        b1, b2 = [(4 * alpha) % Q, 4], [Q, 0]
        # Gauss reduction
        def n2(v): return v[0]*v[0] + v[1]*v[1]
        while True:
            if n2(b1) > n2(b2): b1, b2 = b2, b1
            mu = Fraction(b1[0]*b2[0] + b1[1]*b2[1], n2(b1))
            k = round(mu)
            if k == 0: break
            b2 = [b2[0] - k*b1[0], b2[1] - k*b1[1]]
        diam2 = Umax*Umax + Wmax*Wmax
        if n2(b1) <= diam2:
            raise RuntimeError("lattice too fine for unique-point argument")
        # at most one lattice point in the box; find it by exact enumeration near the box centre
        cx, cy = Fraction(Umax, 2) - off[0], Fraction(Wmax, 2) - off[1]
        det = b1[0]*b2[1] - b1[1]*b2[0]
        ci = (cx*b2[1] - cy*b2[0]) / det; cj = (b1[0]*cy - b1[1]*cx) / det
        for i in range(math.floor(ci) - 2, math.floor(ci) + 4):
            for j in range(math.floor(cj) - 2, math.floor(cj) + 4):
                u = off[0] + i*b1[0] + j*b2[0]; w = off[1] + i*b1[1] + j*b2[1]
                if 0 <= u <= Umax and 0 <= w <= Wmax:
                    sols.append((t, u, w))
    return sols

# ---------------------------------------------------------------- exact enumeration fallback
def _gs(B):
    n = len(B); Bq = [[Fraction(x) for x in r] for r in B]
    star, mu = [], [[Fraction(0)]*n for _ in range(n)]
    for i in range(n):
        v = Bq[i][:]
        for j in range(i):
            mu[i][j] = sum(a*b for a, b in zip(Bq[i], star[j])) / sum(a*a for a in star[j])
            v = [a - mu[i][j]*b for a, b in zip(v, star[j])]
        star.append(v)
    return [sum(a*a for a in s) for s in star], mu

def _isqrt_up(q):           # integer >= sqrt(q) for a nonnegative Fraction q
    return math.isqrt(q.numerator // q.denominator) + 1

def close_points(B, gs2, mu, Binv, y, R2):
    """All lattice vectors v (as integer coordinate vectors in the ambient space) with
    ||v - y||^2 <= R2. Exact (Fincke-Pohst)."""
    n = len(B)
    z = [sum(Fraction(y[j]) * Binv[j][i] for j in range(n)) for i in range(n)]
    out = []
    x = [0]*n
    def rec(i, rem):
        # c_i = z_i - sum_{j>i} mu[j][i] (x_j - z_j)
        c = z[i] - sum(mu[j][i] * (x[j] - z[j]) for j in range(i+1, n))
        s = _isqrt_up(rem / gs2[i])
        for xi in range(math.ceil(c - s), math.floor(c + s) + 1):
            term = (xi - c)**2 * gs2[i]
            if term <= rem:
                x[i] = xi
                if i == 0:
                    out.append([sum(x[k]*B[k][j] for k in range(n)) for j in range(n)])
                else:
                    rec(i - 1, rem - term)
    rec(n - 1, Fraction(R2))
    return out

def _attach(lat):
    if hasattr(lat, '_mu'): return
    gs2, mu = _gs(lat.B)
    lat._mu = mu
    Bi = flint.fmpq_mat(flint.fmpz_mat(lat.B)).inv()
    lat._Binv = [[Fraction(int(Bi[i, j].p), int(Bi[i, j].q)) for j in range(lat.dim)] for i in range(lat.dim)]
    U = int(2.0960 * lat.XL + 3); W = int(1.1834 * lat.XL + 2)
    lat._U, lat._W = U, W
    lat._lw = max(1, round(U / W)); lat._lt = max(1, U // max(lat.T, 1))

def valid_points(lat, m, one_plus_z):
    """Exact check of a top branch: return list of genuine (t,u',w,L) in the box with
    0 < Lam < 80*10^-G, where Lam = t log2 + u' log3 + w log7 - L log10 - log(1+z)."""
    _attach(lat)
    n = lat.dim
    y = [0]*(n-1) + [-m]
    pts = close_points(lat.B, lat.gs2, lat._mu, lat._Binv, y, lat.R2)
    good = []
    for v in pts:
        if n == 4:
            if v[0] % lat._lt or v[2] % lat._lw: continue
            t, up, w, last = v[0] // lat._lt, v[1], v[2] // lat._lw, v[3]
        else:
            if v[1] % lat._lw: continue
            t, up, w, last = 0, v[0], v[1] // lat._lw, v[2]
        num = t*lat.rC[2] + up*lat.rC[3] + w*lat.rC[7] - last
        if num % lat.rC[10]: continue
        L = num // lat.rC[10]
        if not (0 <= t <= lat.T and 2 <= up <= lat._U and 0 <= w <= lat._W and 1 <= L <= lat.XL):
            continue
        lam = (t*flint.arb(2).log() + up*flint.arb(3).log() + w*flint.arb(7).log() - L*flint.arb(10).log()
               - (flint.arb(one_plus_z.numerator) / one_plus_z.denominator).log())
        if lam > 0 and lam < 80 * flint.arb(10)**(-lat.G):
            good.append((t, up, w, L))
        elif not (lam <= 0 or lam >= 80 * flint.arb(10)**(-lat.G)):
            good.append((t, up, w, L, 'undecided'))
    return good
