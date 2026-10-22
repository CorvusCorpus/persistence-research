"""Cumulative-depth bootstrapped prover for base 10 (bags of non-1 digits from {2,3,4,6,7,8,9}).

Facts (y = 2^t 3^u 7^w with digits = bag + 1s, L digits, non-1 digits d_1..d_k at a_1 > ... > a_k):
  (F1) a_k <= 1.                (F2) t <= tmax(bag).
  (F3') for 0 <= j < k:  0 < Lam_j < 80 * 10^(-D_j),  D_j = L - a_{j+1} = g_0 + ... + g_j,
        Lam_j = t log2 + (u+2) log3 + w log7 - L log10 - log(1 + z_j).
        [9y = 10^L (1+z_j) + R_j with 0 < R_j < 80*10^(a_{j+1})]
  (F4) Matveev (real, D=1): D_0 < (ln80 + C4 beta Pi)/ln10,  D_j < (ln80 + C5 beta Pi (D_{j-1} ln10 + ln9))/ln10.
Proved facts: top[j] = T means D_j < T for every solution;  bot[i] = H means a_i < H (a_0 := L); bot[k] = 2 (F1).
Tests at threshold G = ceil(2 log10 XL) + 12 (XL = current bound on L):
  top level j: every prefix (d_1..d_j; g_0..g_{j-1}) allowed by top[0..j-1] is killed by the lattice test
     (no (t,u',w,L) in the box with 0 < Lam_j < 80*10^-G)  =>  D_j < G.
  bottom level j (G_b = max(G, tm+8)): every bottom configuration (d_{j+1}..d_k at heights allowed by
     bot[j+1..k]) has no (t,u',w) in the box with 2^t 3^u' 7^w = F (mod 10^M), M = max(G_b, a_{j+1}+1)
     => a_j < M, so a_j < G_b (configurations with a_{j+1} >= G_b are impossible).
Stop when XL <= LIMIT: every solution is then in the exhaustive list (all zero-free 7-smooth < 10^5000)."""
import sys, math, json, time, os
from fractions import Fraction
from prove_iter import FastTop
from functools import lru_cache

@lru_cache(maxsize=None)
def ordered(bag, r):
    """Distinct length-r sequences from the multiset `bag`, in lexicographic order
    (identical to sorted(set(itertools.permutations(bag, r))) but without enumerating permutations)."""
    vals = sorted(set(bag)); cnt = [bag.count(v) for v in vals]
    out = []; cur = []
    def rec(left):
        if left == 0: out.append(tuple(cur)); return
        for i, v in enumerate(vals):
            if cnt[i]:
                cnt[i] -= 1; cur.append(v); rec(left - 1); cur.pop(); cnt[i] += 1
    rec(r)
    return out
from solver import valid_points
from twoadic import tmax
import ctop2, cbottom2, cbottom3, ctop3, cbottom4

LN = math.log
C4 = 1.4 * 30**7 * 4**4.5
C5 = 1.4 * 30**8 * 5**4.5
PI = LN(2) * LN(3) * LN(7) * LN(10) * (1 + 1e-12)
LIMIT = 5000
BOTTOM_COST = 4.0
CK = 'ckptC'; RES = 'bagresC'
USE_GAP_RESULTS = True        # accept digit sets already proved by the gap-chain prover
os.makedirs(CK, exist_ok=True); os.makedirs(RES, exist_ok=True)

# Bootstrapping bound: 'matveev' (Matveev's top chain only; every record so far) or 'twochain' (top chain with four
# logarithms and a 5-adic bottom chain, docs/BOOTSTRAP_TWO_CHAINS.md). Records written with 'twochain' carry a 'bootstrap' field.
BOOTSTRAP = os.environ.get('PC_BOOTSTRAP', 'matveev')
assert BOOTSTRAP in ('matveev', 'twochain'), BOOTSTRAP
PI3 = LN(3) * LN(7) * LN(10) * (1 + 1e-12)
# Yu, Acta Math. 211 (2013), Theorem 2 with K = Q, p = 5 (case III.2: c2 = 557 (p-1)/(p-2), a2 = 7, kappa2 = 25), n = 3,
# q = 2, u = 1, e_p = f_p = 1, kappa = 0 (1.2), delta >= 1:  ord_5(Xi - 1) < (n+1) C2 Omega max{log B, log 5} (1.10), (1.26)
YU_C2 = 557 * 4 / 3 * 7**3 * 4**4 / (2 * 2 * LN(5)**3) * max(5, (math.e / 3)**3 * LN(5)**4) * max(4 + LN(4), 1, LN(5))
YU_K = 4 * YU_C2 * LN(3) * LN(7) * (1 + 1e-9)       # (n+1) C2 h(3) h(7)
YU_HMIN = max(3 / LN(5), 1) / 25                     # floor of h^(n) in (1.25)

def twochain_cum(k, tm, top, bot):
    """Bound on L from the facts with two chains.  Top (D_j = L - a_(j+1)): Matveev with n = 4 for j = 0 (2, 3, 7, 10) and
    for j >= 1 (3, 7, 10, N_j/2^t with N_j = 10^(D_(j-1)) (1 + z_j) < 9*10^(D_(j-1)), the same Lam_j).  Bottom: for i < k,
    F = -1 + sum_{m>i} c_m 10^(a_m) is a 5-adic unit and 5^(a_i) | 9y - F = F (Xi - 1), Xi = 3^u' 7^w (2^t/F), so
    a_i < YU_K max(h(2^t/F), YU_HMIN) max(log B, log 5) with h(2^t/F) <= max(log(1 + 80*10^(a_(i+1))), t log 2); a_0 = L."""
    INF = float('inf')
    def f(B):
        beta = 1 + LN(B); lB = max(LN(B), LN(5))
        s = [0.0] * k
        for j in range(k):
            if j == 0: m = (LN(80) + C4 * beta * PI) / LN(10)
            else: m = (LN(80) + C4 * beta * PI3 * max(tm * LN(2), s[j - 1] * LN(10) + LN(9))) / LN(10)
            s[j] = min(m, top[j] - 1) if j in top else m
        b = [INF] * (k + 1); b[k] = bot[k] - 1
        for i in range(k - 1, -1, -1):
            yb = YU_K * max(max(LN(80) + b[i + 1] * LN(10), tm * LN(2)) + LN(2), YU_HMIN) * lB
            b[i] = min(bot[i] - 1, yb) if i in bot else yb
        Lb = min([b[0]] + [s[j] + b[j + 1] for j in range(k)])
        return max(2.0960 * Lb + 2.2, tm + 1), Lb
    B = 1e300
    for _ in range(300):
        nb, _ = f(B)
        if nb >= B * (1 - 1e-9): break
        B = nb
    B = max(B * 1.01, 10.0)
    nb, Lb = f(B)
    assert nb <= B
    return max(10, math.ceil(Lb * 1.01))

def matveev_cum(k, tm, top, bot):
    if BOOTSTRAP == 'twochain': return twochain_cum(k, tm, top, bot)
    INF = float('inf')
    def f(B):
        beta = 1 + LN(B)
        Lb = bot[0] - 1 if 0 in bot else INF
        s = 0.0
        for j in range(k):
            m = (LN(80) + (C4 * beta * PI if j == 0 else C5 * beta * PI * (s * LN(10) + LN(9)))) / LN(10)
            if j in top: m = min(m, top[j] - 1)
            s = m
            if (j + 1) in bot: Lb = min(Lb, s + bot[j + 1] - 1)
        return max(2.0960 * Lb + 2.2, tm + 1), Lb
    B = 1e300
    for _ in range(300):
        nb, _ = f(B)
        if nb >= B * (1 - 1e-9): break
        B = nb
    B = max(B * 1.01, 10.0)
    nb, Lb = f(B)
    assert nb <= B
    return max(10, math.ceil(Lb * 1.01))

def deff_top(j, top):
    return [min(top[ip] - (ip - i) for ip in range(i, j)) for i in range(j)]

def ncomp(bounds):
    """number of (g_0..g_{j-1}), g_i >= 1, with g_0+..+g_i <= bounds[i]-1"""
    cnt = {0: 1}
    for bd in bounds:
        nc = {}
        for s, c in cnt.items():
            for s2 in range(s + 1, bd):
                nc[s2] = nc.get(s2, 0) + c
        cnt = nc
    return sum(cnt.values())

def nbottom(heff):
    """number of heights a_{j+1} > ... > a_k >= 0 with a < heff[r] (r = 0 top .. b-1 bottom)"""
    b = len(heff)
    cnt = {-1: 1}                                  # last placed height (start below 0)
    for r in range(b - 1, -1, -1):
        nc = {}
        for a0, c in cnt.items():
            for a in range(a0 + 1, heff[r]):
                nc[a] = nc.get(a, 0) + c
        cnt = nc
    return sum(cnt.values())

def _save(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh: json.dump(obj, fh)
    os.replace(tmp, path)

def prove_bag(bag, log=print, limit=LIMIT, dry=False):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    key = '_'.join(map(str, bag))
    rpath = os.path.join(RES, key + '.json'); cpath = os.path.join(CK, key + '.json')
    if not dry and os.path.exists(rpath):
        r = json.load(open(rpath)); log(f"  bag {bag}: proved earlier (cumulative) -> L <= {r['Lmax']}"); return r['Lmax'], r
    gpath = os.path.join('bagres', key + '.json')          # proved by the gap-chain prover (prove_iter2/3)
    if not dry and USE_GAP_RESULTS and os.path.exists(gpath):
        r = json.load(open(gpath))
        log(f"  bag {bag}: proved earlier by gap chain ({r['source']}) -> L <= {r['Lmax']}")
        return r['Lmax'], {'bag': list(bag), 'Lmax': r['Lmax'], 'top': {}, 'bot': {}, 'gap_bounds': r['bounds'], 'source': r['source']}
    top, bot = {}, {k: 2}
    stats = {'top': 0, 'bottom': 0, 'enum': 0, 'proofs': 0, 'cost': 0}
    partial = None
    if not dry and os.path.exists(cpath):
        c = json.load(open(cpath))
        top = {int(a): v for a, v in c['top'].items()}; bot = {int(a): v for a, v in c['bot'].items()}
        stats = c['stats']; partial = c.get('partial')
        log(f"  bag {bag}: resuming top={top} bot={bot} partial={partial}")
    t0 = time.time(); cache = {}
    while True:
        XL = matveev_cum(k, tm, top, bot)
        G = math.ceil(2 * math.log10(XL)) + 12
        Gb = max(G, tm + 8)
        if XL <= limit:
            res = {'bag': list(bag), 'Lmax': XL, 'top': top, 'bot': bot, 'stats': stats}
            if BOOTSTRAP != 'matveev': res['bootstrap'] = BOOTSTRAP
            if not dry: _save(rpath, res)
            log(f"  bag {bag}: k={k} tmax={tm} top {dict(sorted(top.items()))} bot {dict(sorted(bot.items()))} -> every solution has "
                f"L <= {XL} (top cases {stats['top']}, bottom solves {stats['bottom']}, survivor runs {stats['enum']}, "
                f"proofs {stats['proofs']}, {time.time()-t0:.0f}s)")
            return XL, res
        best = None
        for j in range(k):
            if all(i in top for i in range(j)) and not (j in top and top[j] <= G + 4):
                c = (ncomp(deff_top(j, top)) if j else 1) * max(1, len(ordered(bag, j)))
                if best is None or c < best[0]: best = (c, j, 'top')
            if all(i in bot for i in range(j + 1, k + 1)) and not (j in bot and bot[j] <= Gb + 4):
                c = nbottom(cbottom2.heff_of(j, k, bot)) * len(ordered(bag, k - j)) * BOTTOM_COST
                if best is None or c < best[0]: best = (c, j, 'bottom')
        assert best is not None
        cost, j, side = best
        stats['cost'] += cost
        if not dry:
            skey = (XL, G)
            if skey not in cache:
                cache.clear(); cache[skey] = FastTop(XL, tm, G)
            topl = cache[skey]
            skip = 0
            if partial and partial[:4] == [j, side, G, XL]:
                skip = (partial[4], partial[5]) if len(partial) > 5 else partial[4]
            partial = None
            last = [time.time()]
            def progress(done, j=j, side=side, G=G, XL=XL):
                if time.time() - last[0] > 30:
                    last[0] = time.time()
                    dd = list(done) if isinstance(done, tuple) else [done]
                    _save(cpath, {'top': top, 'bot': bot, 'stats': stats, 'partial': [j, side, G, XL] + dd})
            if side == 'top':
                if j == 0:
                    ok = topl.lat.kills_homogeneous() or not [pt for pt in valid_points(topl.lat, 0, Fraction(1))
                                                             if ctop2.could_be_solution(bag, pt)]
                else:
                    D = deff_top(j, top)
                    known = {i: D[i] for i in range(j)}
                    ok = ctop3.prove_top_c(j, bag, G, topl, known, stats, ordered(bag, j), skip, progress)
            else:
                U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
                ok = cbottom4.prove_bottom_cum(j, bag, Gb, tm, U, W, bot, stats, ordered(bag, k - j), skip, progress)
            if not ok:
                raise RuntimeError(f"bag {bag}: level {j} ({side}) not provable at G={G}")
        stats['proofs'] += 1
        if side == 'top': top[j] = min(top.get(j, 10**18), G)
        else: bot[j] = min(bot.get(j, 10**18), Gb)
        if not dry: _save(cpath, {'top': top, 'bot': bot, 'stats': stats, 'partial': None})

if __name__ == "__main__":
    from bags import bags_of
    Z = json.load(open('Z.json'))
    for s in map(int, sys.argv[1:]):
        t0 = time.time()
        bs = bags_of(s)
        print(f"node {s}: {len(bs)} digit multisets", flush=True)
        Lb = 0; recs = []
        for b in sorted(bs, key=len):
            Lmax, res = prove_bag(b, log=lambda x: print(x, flush=True))
            Lb = max(Lb, Lmax); recs.append([list(b), Lmax, res['top'], res['bot']] + ([res['bootstrap']] if res.get('bootstrap') else []))
        assert Lb <= 5000
        bset = {tuple(sorted(b)) for b in bs}
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                     and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
        _save(f'proofC_{s}.json', {s: {'preimages': pre, 'Lbound': Lb, 'bags': recs, 'method': 'cumulative'}})
