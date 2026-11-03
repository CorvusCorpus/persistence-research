"""Reference implementation (exact, slow) of the meet-in-the-middle top step (docs/MITM_STEPS.md section 2).

Step at top level j with thresholds (Gs, Gp).  Claim: every solution with D_j >= Gs and L >= max(top[j-1], Gs) + Gp is
impossible; together with the trivial case this gives  D_j < Gs  or  L < max(top[j-1], Gs) + Gp.
For such a solution, with T its top j digits (ranges from the top facts) and B the digits at heights < Gp (the lowest b
digits, b >= 1 since a_k <= 1):
  top     0 < t log2 + u' log3 + w log7 - L log10 - log(1 + z_T) < eps = 80 10^-Gs          (prove_cum's (F3'))
  bottom  2^t 3^u' 7^w = F_B (mod 10^Gp), F_B = -1 + sum_B c 10^h; hence t = v2(F_B) when v2(F_B) < Gp (else t >= Gp),
          and 3^u' 7^w = F_B 2^-t (mod 5^Gp), i.e.  u' + e7 w = dl_B (mod phi), phi = 4 5^(Gp-1), e7 = dlog_3 7.
With (u', w) = v_B + m1 b1 + m2 b2 (b1, b2 a reduced basis of {u' + e7 w = 0 mod phi}, v_B a reduced representative) the
pair condition is: some integer m = (m1, m2, L) has m M - tau in K = [-hu, hu] x [-1/2, 1/2] x [-1/2, 1/2], where
  M rows  (b1u/U, b1w/W, l(b1)/eps), (b2u/U, b2w/W, l(b2)/eps), (0, 0, -log10/eps),   l(v) = vu log3 + vw log7,
  tau     = ((c_u - vu)/U, (W/2 - vw)/W, (y - l(v_B))/eps + 1/2),   y = log(1 + z_T) - t log2,   c_u = (U + 2)/2.
With Mr = Tr M LLL-reduced (Tr unimodular) this reads  P_T - Q_B in E = K Mr^-1 (mod Z^3), where
  P_T = (log(1 + z_T)/eps) e3 Mr^-1  and  Q_B = -tau_B Mr^-1,  tau_B = tau with log(1 + z_T) removed.
Bottom entries are hashed by cells at least as wide as E's bounding box (plus the numerical margin); each top configuration
probes the 27 neighbouring cells; every pair found is decided exactly (all integer points m in the box, then the top
inequality with interval arithmetic, the full congruence mod 10^Gp, the side condition and the candidate's digits).
usage (from work/): python3 mitm_ref.py   (self-test: planted genuine numbers, and the join against all pairs)"""
import math, itertools, json, random, sys
import mpmath as mp
import flint
import prove_cum
from twoadic import v2

LN = math.log

def dlog3(x, n):
    """e in [0, 4*5^(n-1)) with 3^e = x (mod 5^n); 3 is a primitive root mod 5^n (Pohlig-Hellman)"""
    M = 5 ** n; ph = 4 * 5 ** (n - 1); x %= M
    assert x % 5
    e4 = next(e for e in range(4) if pow(3, e * 5 ** (n - 1), M) == pow(x, 5 ** (n - 1), M))
    g = pow(3, 4, M); h = pow(x, 4, M); m = n - 1; e5 = 0
    if m:
        gm = pow(g, 5 ** (m - 1), M); tab = {pow(gm, d, M): d for d in range(5)}
        for i in range(m):
            e5 += tab[pow(h * pow(g, -e5, M) % M, 5 ** (m - 1 - i), M)] * 5 ** i
    e = next(e for e in range(e5, ph, 5 ** (n - 1)) if e % 4 == e4)
    assert pow(3, e, M) == x
    return e

class MitmLattice:
    def __init__(self, XL, tm, Gs, Gp):
        self.XL, self.tm, self.Gs, self.Gp = XL, tm, Gs, Gp
        self.U = int(2.0960 * XL + 3); self.W = int(1.1834 * XL + 2)
        self.phi = 4 * 5 ** (Gp - 1); self.e7 = dlog3(7, Gp); self.M5 = 5 ** Gp
        a, b = [self.phi, 0], [(-self.e7) % self.phi, 1]               # Gauss reduction of {u + e7 w = 0 mod phi}
        while True:
            if a[0] ** 2 + a[1] ** 2 < b[0] ** 2 + b[1] ** 2: a, b = b, a
            mu = round(Fraction_div(a[0] * b[0] + a[1] * b[1], b[0] ** 2 + b[1] ** 2))
            if mu == 0: break
            a = [a[0] - mu * b[0], a[1] - mu * b[1]]
        self.b1, self.b2 = a, b
        self.det2 = a[0] * b[1] - a[1] * b[0]
        assert abs(self.det2) == self.phi
        mp.mp.dps = Gs + Gp + 80
        self.eps = 80 * mp.mpf(10) ** (-Gs)
        self.l3, self.l7, self.l10, self.l2 = mp.log(3), mp.log(7), mp.log(10), mp.log(2)
        U, W, eps = self.U, self.W, self.eps
        l = lambda v: v[0] * self.l3 + v[1] * self.l7
        M = mp.matrix([[mp.mpf(a[0]) / U, mp.mpf(a[1]) / W, l(a) / eps],
                       [mp.mpf(b[0]) / U, mp.mpf(b[1]) / W, l(b) / eps],
                       [0, 0, -self.l10 / eps]])
        scale = mp.mpf(2) ** 80 * U                                     # integer image for LLL (rounding is harmless: the
        Mi = flint.fmpz_mat([[int(mp.nint(M[i, k] * scale)) for k in range(3)] for i in range(3)])   # unimodular Tr is exact)
        R = Mi.lll()
        Tq = flint.fmpq_mat(R) * flint.fmpq_mat(Mi).inv()
        self.Tr = [[int(Tq[i, k]) for k in range(3)] for i in range(3)]
        assert all(Tq[i, k] == self.Tr[i][k] for i in range(3) for k in range(3))
        assert abs(flint.fmpz_mat(self.Tr).det()) == 1
        self.Mr = mp.matrix(3, 3)
        for i in range(3):
            for k in range(3):
                self.Mr[i, k] = sum(self.Tr[i][r] * M[r, k] for r in range(3))
        self.Minv = self.Mr ** -1
        self.K = [mp.mpf(U - 2) / (2 * U), mp.mpf(1) / 2, mp.mpf(1) / 2]
        self.hw = [sum(self.K[r] * abs(self.Minv[r, i]) for r in range(3)) for i in range(3)]
        self.margin = mp.mpf(10) ** (-30)                               # numerical error of P_T, Q_B (dps >> 30)
        self.ncell = [max(1, int(1 / (2 * (h + self.margin)))) for h in self.hw]
        self.det3 = abs(mp.det(self.Mr))

    def v_rep(self, dl):
        """representative of {(u', w): u' + e7 w = dl mod phi} near the origin"""
        v = [dl, 0]; a, b = self.b1, self.b2
        c1 = round(Fraction_div(v[0] * b[1] - v[1] * b[0], self.det2)); c2 = round(Fraction_div(a[0] * v[1] - a[1] * v[0], self.det2))
        r = [v[0] - c1 * a[0] - c2 * b[0], v[1] - c1 * a[1] - c2 * b[1]]
        assert (r[0] + self.e7 * r[1] - dl) % self.phi == 0
        return r

    def tau_B(self, vB, t):
        U, W, eps = self.U, self.W, self.eps
        return [mp.mpf(U + 2 - 2 * vB[0]) / (2 * U), mp.mpf(W - 2 * vB[1]) / (2 * W),
                (-t * self.l2 - vB[0] * self.l3 - vB[1] * self.l7) / eps + mp.mpf(1) / 2]

    def coords(self, vec):
        return [sum(vec[r] * self.Minv[r, i] for r in range(3)) for i in range(3)]

    def cell(self, p):
        return tuple(int(mp.floor((p[i] - mp.floor(p[i])) * self.ncell[i])) % self.ncell[i] for i in range(3))

def Fraction_div(a, b):
    from fractions import Fraction
    return Fraction(a, b)

def bottom_entries(lat, bag, bot, bmax):
    """every bottom window: the lowest b <= bmax digits at heights < Gp (a_k <= 1, the bottom facts), with their values;
    one entry per admissible t.  Returns list of (Q_B cell, entry)."""
    k = len(bag); Gp = lat.Gp; out = []
    for b in range(1, bmax + 1):
        if b < k and bot.get(k - b, 10 ** 9) <= Gp: continue              # digit k-b must lie inside the window
        hi = []
        for i in range(k - b + 1, k + 1):                                   # bounds a_i < h for the window's digits
            h = Gp - (i - (k - b + 1))
            for ip in range(k - b + 1, i + 1):
                if ip in bot: h = min(h, bot[ip] - (i - ip))
            hi.append(h)
        if min(hi) <= 0: continue
        orders = prove_cum.ordered(bag, b)
        def heights(r, top):
            if r == b: yield []; return
            for a in range(top - 1, -1, -1):
                if a < hi[r] and (r < b - 1 or a <= 1):
                    for rest in heights(r + 1, a): yield [a] + rest
        for hs in heights(0, Gp):
            for o in orders:
                F = -1 + sum(9 * (d - 1) * 10 ** h for d, h in zip(o, hs))
                tv = v2(F)
                ts = [tv] if tv < Gp else list(range(Gp, lat.tm + 1))
                for t in ts:
                    if t > lat.tm: continue
                    dl = dlog3(F * pow(2, -t, lat.M5) % lat.M5, Gp)
                    vB = lat.v_rep(dl)
                    q = [-c for c in lat.coords(lat.tau_B(vB, t))]
                    out.append((lat.cell(q), {'heights': hs, 'digits': o, 'F': F, 't': t, 'vB': vB, 'Q': q}))
    return out

def top_configs(bag, j, ranges):
    """(digits, depths D_0 < ... < D_(j-1)) with D_i < ranges[i]"""
    def deps(r, last):
        if r == j: yield []; return
        for D in range(last + 1, ranges[r]):
            for rest in deps(r + 1, D): yield [D] + rest
    for o in prove_cum.ordered(bag, j):
        for ds in deps(0, 0): yield o, ds

def lam_T(o, ds):
    return mp.log(1 + sum(9 * (d - 1) * mp.mpf(10) ** (-D) for d, D in zip(o, ds)))

def exact_pairs(lat, lam, e):
    """every (t, u', w, L) in the box with both relations for top value lam = log(1+z_T) and bottom entry e"""
    tau = lat.tau_B(e['vB'], e['t']); tau[2] += lam / lat.eps
    c = lat.coords(tau); out = []
    rngs = [range(int(mp.floor(c[i] - lat.hw[i] - lat.margin)), int(mp.ceil(c[i] + lat.hw[i] + lat.margin)) + 1) for i in range(3)]
    for mr in itertools.product(*rngs):
        m = [sum(mr[r] * lat.Tr[r][i] for r in range(3)) for i in range(3)]
        u = e['vB'][0] + m[0] * lat.b1[0] + m[1] * lat.b2[0]; w = e['vB'][1] + m[0] * lat.b1[1] + m[1] * lat.b2[1]; L = m[2]
        t = e['t']
        if not (2 <= u <= lat.U and 0 <= w <= lat.W): continue
        val = t * lat.l2 + u * lat.l3 + w * lat.l7 - L * lat.l10 - lam
        if not (0 < val < lat.eps): continue
        if (pow(2, t, 10 ** lat.Gp) * pow(3, u, 10 ** lat.Gp) * pow(7, w, 10 ** lat.Gp) - e['F']) % 10 ** lat.Gp: continue
        out.append((t, u, w, L))
    return out

def mitm_join(lat, tops, entries, bag):
    """candidate pairs (top index, entry) whose torus points are within E (cell probing); multiset-compatible only"""
    from collections import Counter
    table = {}
    for key, e in entries: table.setdefault(key, []).append(e)
    full = Counter(bag); found = []; probes = 0
    for ti, (o, ds) in enumerate(tops):
        lam = lam_T(o, ds)
        p = [lam / lat.eps * lat.Minv[2, i] for i in range(3)]
        k0 = lat.cell(p)
        seen = set()
        for d in itertools.product((-1, 0, 1), repeat=3):
            key = tuple((k0[i] + d[i]) % lat.ncell[i] for i in range(3))
            if key in seen: continue
            seen.add(key)
            for e in table.get(key, []):
                probes += 1
                if Counter(o) + Counter(e['digits']) <= full: found.append((ti, e, lam))
    return found, probes

def selftest():
    """1) join against all pairs on a small step; 2) planted genuine numbers from Z.json"""
    bag = (2, 2, 3, 3, 7, 7)
    k = len(bag); ok = True
    from collections import Counter
    full = Counter(bag)
    for XL, tm, Gs, Gp, j, ranges in ((3000, 12, 9, 7, 2, [6, 7]), (3000, 12, 6, 4, 2, [6, 7]), (3000, 12, 8, 6, 3, [7, 8, 9])):
        lat = MitmLattice(XL, tm, Gs, Gp)
        print(f"lattice Gs {Gs} Gp {Gp}: U {lat.U} W {lat.W} phi {lat.phi} b1 {lat.b1} b2 {lat.b2} det3 {mp.nstr(lat.det3, 5)} "
              f"half-widths {[mp.nstr(h, 3) for h in lat.hw]} cells {lat.ncell}")
        tops = list(top_configs(bag, j, ranges))
        ents = bottom_entries(lat, bag, {k: 2}, k - j)
        found, probes = mitm_join(lat, tops, ents, bag)
        fset = {(ti, id(e)) for ti, e, _ in found}
        # every pair with an exact solution must be among the probed pairs
        missed = 0; allsol = 0; nsol = 0
        for ti, (o, ds) in enumerate(tops):
            lam = lam_T(o, ds)
            for _, e in ents:
                if not Counter(o) + Counter(e['digits']) <= full: continue
                s = exact_pairs(lat, lam, e)
                if s:
                    allsol += 1; nsol += len(s)
                    if (ti, id(e)) not in fset: missed += 1
        print(f"1) {len(tops)} top configurations x {len(ents)} bottom entries: {probes} probes, {len(found)} candidate pairs; "
              f"pairs with an exact solution {allsol} ({nsol} lattice points), missed by the join {missed}")
        ok &= missed == 0
    # 2) planted genuine numbers: their own (T, B) pair must be found with their own (t, u', w, L)
    Z = json.load(open('Z.json'))
    tested = hit = 0
    for y in Z:
        s = str(y)
        if y < 10 or '0' in s or '5' in s: continue
        L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch) for ch in s if ch != '1']
        kk = len(pos)
        if kk < 2 or L < 12: continue
        for jj in range(1, kk):
            Dprev = L - pos[jj - 1]; Dj = L - pos[jj]
            Gp_ = min(6, pos[jj] + 1)                               # the window holds digit j+1 and everything below it
            low = [p for p in pos if p < Gp_]
            if Dj < 3 or not low or len(low) > kk - jj or L < max(Dprev, Dj) + Gp_: continue
            Gs_ = min(Dj, 14)
            t = v2(y); r = y >> t; u = 0
            while r % 3 == 0: r //= 3; u += 1
            w = 0
            while r % 7 == 0: r //= 7; w += 1
            lat2 = MitmLattice(max(200, L + 50), t + 2, Gs_, Gp_)
            o = tuple(digs[:jj]); ds = [L - p for p in pos[:jj]]
            lam = lam_T(o, ds)
            hs = [p for p in pos if p < Gp_]; od = tuple(digs[len(pos) - len(hs):])
            F = -1 + sum(9 * (d - 1) * 10 ** h for d, h in zip(od, hs)); tt = v2(F)
            if tt != t: continue
            dl = dlog3(F * pow(2, -tt, lat2.M5) % lat2.M5, Gp_); vB = lat2.v_rep(dl)
            q = [-c for c in lat2.coords(lat2.tau_B(vB, tt))]
            e = {'heights': hs, 'digits': od, 'F': F, 't': tt, 'vB': vB, 'Q': q}
            p = [lam / lat2.eps * lat2.Minv[2, i] for i in range(3)]
            d = [(p[i] - q[i]) - mp.nint(p[i] - q[i]) for i in range(3)]
            inside = all(abs(d[i]) <= lat2.hw[i] + lat2.margin for i in range(3))
            sols = exact_pairs(lat2, lam, e)
            tested += 1; good = inside and (t, u + 2, w, L) in sols
            hit += good
            if not good: print("   planted MISS", y, jj, inside, sols[:3], (t, u + 2, w, L))
            break
        if tested >= 40: break
    print(f"2) planted genuine numbers: {hit}/{tested} found (torus condition and exact pair)")
    ok &= hit == tested and tested > 0
    print("MITM REFERENCE SELF-TEST", "PASSED" if ok else "FAILED")
    return ok

if __name__ == '__main__':
    sys.exit(0 if selftest() else 1)
