"""Bootstrapped case-split prover for bags with up to ~9 non-1 digits.

Facts F1-F3 as in solver.py.  We prove 'g_j < G' gap by gap:
  * top proof of gap j: enumerate the top j digits and gaps g_0..g_{j-1} < G (already proved)
    and kill each branch with the lattice test on Lam_j (requires g_j >= G otherwise);
  * bottom proof of gap j: enumerate the k-j digits below gap j with gaps < G (already proved)
    and a_k in {0,1}; a run of >= G ones above them fixes y mod 10^(a_{j+1}+G); exact 10-adic test.
After each proved gap, Matveev's chain is recomputed with the proved bounds substituted, which
shrinks the length bound X and hence the threshold G; when G drops we restart at the new G.
When every gap is < G, L <= k(G-1)+1 and the exhaustive list (L <= 5000) finishes the job."""
import sys, math, json, time, itertools
from fractions import Fraction
import flint
from solver import TopLattice, valid_points
from fastbottom import BottomContext
from twoadic import tmax
from bags import bags_of

LN = math.log
C4 = 1.4 * 30**7 * 4**4.5
C5 = 1.4 * 30**8 * 5**4.5
PI = LN(2) * LN(3) * LN(7) * LN(10) * (1 + 1e-12)

def matveev_partial(k, tm, known):
    """Upper bound X_L on L given proved bounds known[j] (meaning g_j < known[j])."""
    def f(B):
        beta = 1 + LN(B); s = 0.0
        for j in range(k):
            if j in known:
                g = known[j] - 1
            elif j == 0:
                g = (LN(80) + C4 * beta * PI) / LN(10)
            else:
                g = (LN(80) + C5 * beta * PI * (s * LN(10) + LN(9))) / LN(10)
            s += g
        Lmax = s + 1
        return max(2.0960 * Lmax + 2.2, tm + 1), Lmax
    B = 1e300
    for _ in range(300):
        nb, _ = f(B)
        if nb >= B * (1 - 1e-9): break
        B = nb
    B = max(B * 1.01, 10.0)
    nb, Lmax = f(B)
    assert nb <= B
    return max(10, math.ceil(Lmax * 1.01))

def ordered(bag, r):
    return sorted(set(itertools.permutations(bag, r)))

class FastTop:
    def __init__(self, XL, tm, G):
        self.lat = TopLattice(XL, tm, G)
        lat = self.lat
        q = lat.Binv_last_col_last(); self.qa, self.qb = q.numerator, q.denominator
        g = lat.gs2[-1]; R2 = Fraction(lat.R2)
        self.K1 = g.numerator * R2.denominator
        self.K2 = R2.numerator * g.denominator * self.qb * self.qb
        self.C = lat.C
        self.prec = max(256, self.C.bit_length() + 96)   # enough for C*log to +-0.01 and for valid_points
        flint.ctx.prec = self.prec
        self.log10 = flint.arb(10).log()
    def target(self, N, s):
        # m = round(C * (s log10 - log N)) = round(-C log(1+z)), with 1+z = N/10^s
        if flint.ctx.prec != self.prec: flint.ctx.prec = self.prec
        v = flint.arb(self.C) * (s * self.log10 - flint.arb(N).log())
        m = int(v.mid().floor().unique_fmpz())
        if v - m > 0.5: m += 1
        assert abs(v - m) < 0.51
        return m
    def kills(self, m):
        r = (-m * self.qa) % self.qb
        d = min(r, self.qb - r)
        return d * d * self.K1 > self.K2

def prove_bag(bag, log=print):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    known = {}                                     # j -> bound (g_j < bound), any G level
    XL = matveev_partial(k, tm, known)
    G = math.ceil(2 * math.log10(XL)) + 12
    stats = {'top': 0, 'bottom': 0, 'enum': 0, 'restarts': 0}
    t0 = time.time()
    while True:
        top = FastTop(XL, tm, G)
        U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
        ctx = BottomContext(G, tm, U, W)
        proven = set()                             # gaps proved < G at this G
        restart = False
        while len(proven) < k:
            # candidate gaps: next from the top (all smaller indices proven), next from the bottom
            jt = next(j for j in range(k) if j not in proven)
            jb = max(j for j in range(k) if j not in proven)
            ct = (G ** jt) * max(1, len(ordered(bag, jt)))
            cb = 2 * (G ** max(k - jb - 1, 0)) * len(ordered(bag, k - jb)) * 30   # bottom cases ~30x costlier
            if ct <= cb:
                ok = prove_top(jt, bag, G, top, stats); j = jt; side = 'top'
            else:
                ok = prove_bottom(jb, bag, G, ctx, stats); j = jb; side = 'bottom'
            if not ok:
                raise RuntimeError(f"gap {j} ({side}) not provable at G={G} for bag {bag}")
            proven.add(j); known[j] = min(known.get(j, 10**18), G)
            XLn = matveev_partial(k, tm, known)
            Gn = math.ceil(2 * math.log10(XLn)) + 12
            if Gn <= G - 4 and len(proven) < k:
                stats['restarts'] += 1
                XL, G = XLn, Gn; restart = True
                break
        if not restart:
            Lmax = k * (G - 1) + 1
            log(f"  bag {bag}: k={k} tmax={tm} final G={G} XL={XL:.2e} -> every solution has L <= {Lmax} "
                f"(top {stats['top']}, bottom {stats['bottom']}, enum {stats['enum']}, restarts {stats['restarts']}, "
                f"{time.time()-t0:.0f}s)")
            return Lmax

def prove_top(j, bag, G, top, stats):
    lat = top.lat
    if j == 0:
        if lat.kills_homogeneous(): return True
        stats['enum'] += 1
        return not valid_points(lat, 0, Fraction(1))
    for digs in ordered(bag, j):
        cs = [9 * (d - 1) for d in digs]
        for gaps in itertools.product(range(1, G), repeat=j):
            s = sum(gaps)
            # N = 10^s (1 + sum c_i 10^{-s_i})
            N = 10**s; acc = 0
            for c, g in zip(cs, gaps):
                acc += g; N += c * 10**(s - acc)
            stats['top'] += 1
            m = top.target(N, s)
            if not top.kills(m):
                stats['enum'] += 1
                if valid_points(lat, m, Fraction(N, 10**s)): return False
    return True

def prove_bottom(j, bag, G, ctx, stats):
    k = len(bag); b = k - j
    for digs in ordered(bag, b):                  # top-to-bottom: d_{j+1} .. d_k
        cs = [9 * (d - 1) for d in digs]
        lower = list(reversed(cs[1:]))            # bottom-up coefficients of d_k, d_{k-1}, ...
        for ak in (0, 1):
            for lowgaps in itertools.product(range(1, G), repeat=max(b - 2, 0)):
                posl = [ak]
                for g in lowgaps: posl.append(posl[-1] + g)
                if b == 1:
                    Flow, pts = -1, [ak]
                else:
                    Flow = -1 + sum(c * 10**p for c, p in zip(lower, posl))
                    pts = [posl[-1] + g for g in range(1, G)]
                for p in pts:
                    stats['bottom'] += 1
                    if ctx.solutions(Flow, cs[0], p, p + G): return False
    return True

if __name__ == "__main__":
    Z = json.load(open('Z.json'))
    for s in map(int, sys.argv[1:]):
        t0 = time.time()
        bs = bags_of(s)
        print(f"node {s}: {len(bs)} digit multisets", flush=True)
        Lb = 0; recs = []
        for b in sorted(bs, key=len):
            Lmax = prove_bag(b, log=lambda x: print(x, flush=True))
            Lb = max(Lb, Lmax); recs.append([list(b), Lmax])
        assert Lb <= 5000
        bset = {tuple(sorted(b)) for b in bs}
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                     and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
        json.dump({s: {'preimages': pre, 'Lbound': Lb, 'bags': recs}}, open(f'proofI_{s}.json', 'w'))
