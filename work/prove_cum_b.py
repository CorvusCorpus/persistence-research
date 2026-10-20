"""Cumulative-depth prover for bases with one free prime (b = 3, 4, 6).
y = B q^c (B = product of bounded powers of the primes dividing b; q the free prime), base-b digits = bag + 1s.
(b-1) y = b^L (1 + z_j) + R_j,  0 < R_j < (b-2) b^(a_{j+1}+1)  (j < k)  =>
   0 < Lam_j = log((b-1)B) + c log q - L log b - log(1+z_j) < eps0 * b^(-D_j),  eps0 = (b-2) b,  D_j = L - a_{j+1}.
Trailing run: (b-1) y = b^L (1+z_k) - 1  =>  |Lam_k| < 2 b^-L.
Facts: top[j] = T: D_j < T (j < k);  bot[i] = H: a_i < H (i = 0..k, a_0 = L).
Top test of level j (threshold G): every prefix allowed by top[0..j-1], every cofactor B: 2-dim lattice test
  (no (c, L) in the box with 0 < Lam_j < eps0 b^-G)  =>  D_j < G.
Bottom test of level j (G_b): every configuration of digits j+1..k at heights allowed by bot[j+1..k]
  (j = k: no digits), M = max(G_b, a_{j+1}+1): y mod b^M is fixed if a_j >= M; no c <= Xc with
  B q^c = r (mod b^M) for any B  =>  a_j < M, hence a_j < G_b.
Matveev (n = 3): D_j < ln(eps0)/lb + C3 beta lq (D_{j-1} lb + cst);  L < 1 + ln2/lb + C3 beta lq (D_{k-1} lb + cst)."""
import sys, math, json, time, os, itertools, subprocess, tempfile
from fractions import Fraction
from functools import lru_cache
import flint
from basesolver import Base, cofactors, Top2, lam_ok, solve_c, C3
from solver import _gs

LN = math.log
HERE = os.path.dirname(os.path.abspath(__file__))
BINB = os.path.join(HERE, 'csrc', 'topb')

@lru_cache(maxsize=None)
def ordered(bag, r):
    vals = sorted(set(bag)); cnt = [bag.count(v) for v in vals]
    out = []; cur = []
    def rec(left):
        if left == 0: out.append(tuple(cur)); return
        for i, v in enumerate(vals):
            if cnt[i]:
                cnt[i] -= 1; cur.append(v); rec(left - 1); cur.pop(); cnt[i] += 1
    rec(r)
    return out

def matveev_cum_b(k, base, Bmax, top, bot):
    b, q = base.b, base.q
    lq, lb = LN(q), LN(b)
    cst = LN((b - 1) * Bmax) + LN(b)
    INF = float('inf')
    def f(Bc):
        beta = 1 + LN(Bc)
        Lb = bot[0] - 1 if 0 in bot else INF
        s = 0.0
        for j in range(k):
            m = LN(base.eps0) / lb + C3 * beta * lq * (s * lb + cst)
            if j in top: m = min(m, top[j] - 1)
            s = m
            if (j + 1) in bot: Lb = min(Lb, s + bot[j + 1] - 1)
        Lk = 1 + LN(2) / lb + C3 * beta * lq * (s * lb + cst)
        Lb = min(Lb, Lk)
        coef = Lb * lb / lq + 2
        return max(coef, Lb, 10.0), Lb
    Bc = 1e300
    for _ in range(300):
        nb, _ = f(Bc)
        if nb >= Bc * (1 - 1e-9): break
        Bc = nb
    Bc = max(Bc * 1.01, 10.0)
    nb, Lb = f(Bc)
    assert nb <= Bc
    return max(4, math.ceil(Lb * 1.01))

def deff_top(j, top):
    return [min(top[ip] - (ip - i) for ip in range(i, j)) for i in range(j)]

def heff_of(j, k, bot):
    return [min(bot[ip] - (i - ip) for ip in range(j + 1, i + 1)) for i in range(j + 1, k + 1)]

def ncomp(bounds):
    cnt = {0: 1}
    for bd in bounds:
        nc = {}
        for s, c in cnt.items():
            for s2 in range(s + 1, bd): nc[s2] = nc.get(s2, 0) + c
        cnt = nc
    return sum(cnt.values())

def nbottom(heff):
    cnt = {-1: 1}
    for r in range(len(heff) - 1, -1, -1):
        nc = {}
        for a0, c in cnt.items():
            for a in range(a0 + 1, heff[r]): nc[a] = nc.get(a, 0) + c
        cnt = nc
    return sum(cnt.values())

class TopB:
    """Exact 2-dim lattice data (Top2) packaged for the C loop (level-1/level-2 kill tests)."""
    def __init__(self, base, XL, G):
        T = Top2(base, XL, G); self.T = T
        self.B = T.Bm; self.dim = 2; self.C = T.C; self.R2 = T.R2
        gs2, mu = _gs(self.B)
        assert Fraction(gs2[0]) == Fraction(T.g1) and Fraction(gs2[1]) == Fraction(T.g2)
        det = self.B[0][0] * self.B[1][1] - self.B[0][1] * self.B[1][0]
        row = [Fraction(-self.B[1][0], det), Fraction(self.B[0][0], det)]       # last row of B^-1
        D = 1
        for x in row: D = D * x.denominator // math.gcd(D, x.denominator)
        a = [int(x * D) for x in row]
        q = row[1]; qa, qb = q.numerator, q.denominator
        R2 = Fraction(T.R2); g = Fraction(gs2[1])
        self.data = dict(C=T.C, qa=qa, qb=qb, K1=g.numerator * R2.denominator, K2=R2.numerator * g.denominator * qb * qb,
                         D=D, a2=a[0], a3=a[1], P=mu[1][0].numerator, Q=mu[1][0].denominator,
                         g2n=Fraction(gs2[0]).numerator, g2d=Fraction(gs2[0]).denominator,
                         g3n=g.numerator, g3d=g.denominator, Rn=R2.numerator, Rd=R2.denominator)
        self.prec = max(T.prec, T.C.bit_length() + 64)

def run_topb(topb, base, Ws, j, bounds, orders, skip=0, g0lo=1, g0hi=10**6):
    d = topb.data
    lines = [f"prec {topb.prec}"] + [f"{kk} {d[kk]}" for kk in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
    lines += [f"base {base.b}", f"ntg {len(Ws)}"] + [f"W{i} {w}" for i, w in enumerate(Ws)]
    lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, bounds)), f"nord {len(orders)}"]
    lines += [" ".join(str((base.b - 1) * (x - 1)) for x in o) for o in orders]
    lines += [f"skip {skip}", "maxsurv 0", "debug 0", f"g0lo {g0lo}", f"g0hi {g0hi}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BINB, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("C topb failed")
    finally:
        os.unlink(path)

def digits_ok(base, bag, B, c, L, K=40):
    """Could y = B q^c (L digits) be a solution with digit multiset `bag`?  Exact check of all digits when y is
    short, else of its last K base-b digits (zero-free, only 1s and bag digits within multiplicities).
    A lattice survivor (c, L) whose y fails this is not a solution, so it can be discarded."""
    b, q = base.b, base.q
    from collections import Counter
    if L <= K + 8:
        y = B * q**c
        ds = flint.fmpz(y).str(b)
        return '0' not in ds and sorted(int(ch, b) for ch in ds if ch != '1') == sorted(bag)
    y = B * pow(q, c, b**K) % b**K
    cnt = Counter(bag)
    for _ in range(K):
        d = y % b; y //= b
        if d == 0: return False
        if d != 1:
            cnt[d] -= 1
            if cnt[d] < 0: return False
    return True

def _resolve(st, base, bag, B, c, L):
    """True if the lattice survivor (c, L) might be a solution (keeps the case open)."""
    if digits_ok(base, bag, B, c, L):
        return True
    st['resolved'] = st.get('resolved', 0) + 1
    return False

def top_test(j, bag, base, Bs, XL, G, top, st, skip=0, progress=lambda d: None):
    b = base.b
    tb = TopB(base, XL, G)
    T = tb.T
    if j == 0:
        for B in Bs:
            flint.ctx.prec = T.prec
            kappa = flint.arb((b - 1) * B).log()
            if T.survivors(kappa, lambda c, L, B=B: lam_ok(base, B, c, L, Fraction(1), G) and _resolve(st, base, bag, B, c, L)): return False
        return True
    orders = ordered(bag, j)
    bounds = deff_top(j, top)
    ok = True
    o0, i0 = (skip, 1) if isinstance(skip, int) else skip
    g0max = bounds[0] - 1
    for oo in range(o0, len(orders)):
      for g0 in range(i0 if oo == o0 else 1, g0max + 1):
        progress((oo, g0))
        for f in run_topb(tb, base, [(b - 1) * B for B in Bs], j, bounds, [orders[oo]], 0, g0, g0):
          if f[0] == 'SURV':
            o = oo; gaps = list(map(int, f[2:2 + j - 1])) + [int(f[-4])]; B = Bs[int(f[-2])]
            st['enum'] += 1
            s = sum(gaps); N = b**s; acc = 0
            for x, g in zip(orders[o], gaps):
                acc += g; N += (b - 1) * (x - 1) * b**(s - acc)
            opz = Fraction(N, b**s)
            flint.ctx.prec = T.prec
            kappa = flint.arb((b - 1) * B).log() - (flint.arb(N) / flint.arb(b)**s).log()
            if T.survivors(kappa, lambda c, L: lam_ok(base, B, c, L, opz, G) and _resolve(st, base, bag, B, c, L)): ok = False
          elif f[0] == 'DONE': st['top'] += int(f[1])
    return ok

def bottom_test(j, bag, base, Bs, XL, Gb, bot, st, skip=0, progress=lambda d: None, pb=None):
    """C implementation (csrc/bottomb); identical statement to bottom_test_py."""
    import cbottomb
    k = len(bag)
    Xc = int(XL * LN(base.b) / LN(base.q)) + 2
    heff = heff_of(j, k, bot)
    orders = ordered(bag, k - j) if k - j else [()]
    nb = k - j
    o0, i0 = (skip, 0) if isinstance(skip, int) else skip
    ok = True
    for oo in range(o0, len(orders)):
        akmax = heff[nb - 1] - 1 if nb else 0
        for ak in range(i0 if oo == o0 else 0, akmax + 1):
            progress((oo, ak))
            if not cbottomb.bottom_test_c(j, bag, base, pb, Xc, Gb, bot, st, [orders[oo]], heff, 0, lambda d: None, ak, ak):
                ok = False
                return ok
    return ok

def bottom_test_py(j, bag, base, Bs, XL, Gb, bot, st, skip=0, progress=lambda d: None):
    b = base.b; k = len(bag)
    Xc = int(XL * LN(b) / LN(base.q)) + 2
    nb = k - j
    heff = heff_of(j, k, bot)
    orders = ordered(bag, nb) if nb else [()]
    for idx, digs in enumerate(orders):
        if idx < skip: continue
        progress(idx)
        # heights: digs[0] at a_{j+1} (top) ... digs[-1] at a_k (bottom); enumerate increasing from the bottom
        def rec(r, lo, pos):
            if r < 0:
                yield pos; return
            for a in range(lo, heff[r]):
                pos[r] = a
                yield from rec(r - 1, a + 1, pos)
        gen = rec(nb - 1, 0, [0] * nb) if nb else iter([[]])
        for pos in gen:
            M = max(Gb, pos[0] + 1) if nb else Gb
            r = (b**M - 1) // (b - 1) + sum((d - 1) * b**a for d, a in zip(digs, pos))
            rm = r % b**M
            for B in Bs:
                st['bottom'] += 1
                sol = solve_c(base, B, rm, M, Xc)
                if sol is None or sol: return False
    return True

def _save(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh: json.dump(obj, fh)
    os.replace(tmp, path)

BOTTOM_COST = 8.0      # C bottom: one solve (~3us) per configuration vs ~0.3us per (top configuration, cofactor)
def prove_bag(bag, base, limit, log=print, dry=False, ckdir='ckptB'):
    os.makedirs(ckdir, exist_ok=True)
    bag = tuple(sorted(bag)); k = len(bag); b = base.b
    Bs, pb = cofactors(bag, base); Bmax = max(Bs)
    Tp = max([pb[p] for p in pb] + [0])
    key = f"b{b}_" + '_'.join(map(str, bag)); cpath = os.path.join(ckdir, key + '.json'); rpath = os.path.join(ckdir, key + '.res.json')
    if not dry and os.path.exists(rpath):
        r = json.load(open(rpath)); log(f"  base {b} bag {bag}: proved earlier -> L <= {r['Lmax']}"); return r['Lmax'], r
    top, bot = {}, {}
    st = {'top': 0, 'bottom': 0, 'enum': 0, 'proofs': 0, 'cost': 0}
    partial = None
    if not dry and os.path.exists(cpath):
        c = json.load(open(cpath)); top = {int(a): v for a, v in c['top'].items()}; bot = {int(a): v for a, v in c['bot'].items()}
        st = c['stats']; partial = c.get('partial')
    t0 = time.time()
    while True:
        XL = matveev_cum_b(k, base, Bmax, top, bot)
        Xc = XL * LN(b) / LN(base.q) + 2
        G = math.ceil((LN(Xc) + 7 * LN(10)) / LN(b)) + 2
        Gb = max(G, Tp + 6)
        if XL <= limit:
            res = {'bag': list(bag), 'Lmax': XL, 'top': top, 'bot': bot, 'stats': st, 'cofactors': Bs, 'pbounds': pb}
            if not dry: _save(rpath, res)
            log(f"  base {b} bag {bag}: cofactors {len(Bs)} p-bounds {pb} top {dict(sorted(top.items()))} bot {dict(sorted(bot.items()))} "
                f"-> every solution has L <= {XL} (top {st['top']}, bottom {st['bottom']}, survivors {st['enum']}, proofs {st['proofs']}, {time.time()-t0:.1f}s)")
            return XL, res
        best = None
        for j in range(k + 1):
            if j < k and all(i in top for i in range(j)) and not (j in top and top[j] <= G + 3):
                c = (ncomp(deff_top(j, top)) if j else 1) * max(1, len(ordered(bag, j))) * len(Bs)
                if best is None or c < best[0]: best = (c, j, 'top')
            if all(i in bot for i in range(j + 1, k + 1)) and not (j in bot and bot[j] <= Gb + 3):
                c = (nbottom(heff_of(j, k, bot)) if j < k else 1) * max(1, len(ordered(bag, k - j))) * BOTTOM_COST
                if best is None or c < best[0]: best = (c, j, 'bottom')
        assert best is not None, (top, bot, G)
        cost, j, side = best
        st['cost'] += cost
        if not dry:
            skip = 0
            if partial and partial[:4] == [j, side, G, XL]:
                skip = (partial[4], partial[5]) if len(partial) > 5 else partial[4]
            partial = None
            last = [time.time()]
            def progress(done, j=j, side=side, G=G, XL=XL):
                if time.time() - last[0] > 30:
                    last[0] = time.time()
                    dd = list(done) if isinstance(done, tuple) else [done]
                    _save(cpath, {'top': top, 'bot': bot, 'stats': st, 'partial': [j, side, G, XL] + dd})
            ok = (top_test(j, bag, base, Bs, XL, G, top, st, skip, progress) if side == 'top'
                  else bottom_test(j, bag, base, Bs, XL, Gb, bot, st, skip, progress, pb=pb))
            if not ok: raise RuntimeError(f"base {b} bag {bag}: level {j} ({side}) survives at G={G}")
        st['proofs'] += 1
        if side == 'top': top[j] = min(top.get(j, 10**18), G)
        else: bot[j] = min(bot.get(j, 10**18), Gb)
        if not dry: _save(cpath, {'top': top, 'bot': bot, 'stats': st, 'partial': None})

if __name__ == "__main__":
    b, D = int(sys.argv[1]), int(sys.argv[2])
    base = Base(b)
    from basetrees import trees, bags_of, fprod, smooth_not_div
    cache = f'trees_b{b}_D{D}.json'
    if os.path.exists(cache):
        cj = json.load(open(cache))
        nodes = {int(d): {int(s): dp for s, dp in v.items()} for d, v in cj['nodes'].items()}
        zf_list = [int(y) for y in cj['zerofree']]
    else:
        n, nz, nodes = trees(b, D)
        zf_list = [y for y in smooth_not_div(b, D) if '0' not in flint.fmpz(y).str(b)]
    only = set(map(int, sys.argv[3:])) if len(sys.argv) > 3 else None
    out = {}
    if os.path.exists(f'base{b}_cum.json'): out = {int(s): v for s, v in json.load(open(f'base{b}_cum.json')).items()}
    order = sorted(((max((len(t) for t in bags_of(s, b)), default=0), s, d) for d in nodes for s in nodes[d]))
    for km, s, d in order:
        if only and s not in only: continue
        if s in out: continue
        t0 = time.time()
        Lb = max([prove_bag(t, base, D - 1, log=lambda x: print(x, flush=True))[0] for t in bags_of(s, b)] + [0])
        pre = sorted(y for y in zf_list if y >= b and fprod(y, b) == s)
        out[s] = {'target': d, 'depth': nodes[d][s], 'preimages': pre, 'Lbound': Lb, 'kmax': km}
        _save(f'base{b}_cum.json', {str(x): v for x, v in out.items()})
        print(f"base {b}: node {s} (target {d}, kmax {km}): all smooth y with digit product {s}: {pre} (L <= {Lb}, {time.time()-t0:.0f}s)", flush=True)
    for d in sorted(nodes):
        ns = nodes[d]
        if all(s in out for s in ns):
            for s in ns:
                for y in out[s]['preimages']: assert y in ns, (s, y)
            print(f"THEOREM base {b}, final digit {d}: tree complete ({len(ns)} nodes), max persistence {max(ns.values()) + 1}")
        else:
            print(f"base {b}, final digit {d}: open nodes {[s for s in ns if s not in out]}")
