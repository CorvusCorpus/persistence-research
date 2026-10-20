"""Bootstrapped prover v2.  Same facts and tests as prove_iter.py; differences:
  * every proved fact 'g_j < known[j]' is kept and used as the enumeration range for that gap
    in all later top/bottom case enumerations (valid: it is a proved inequality);
  * the next action is the cheapest available proof among unknown gaps and gaps whose bound can
    be tightened to the current threshold G (G = ceil(2 log10 X) + 12, X from Matveev's chain with
    all proved bounds substituted);
  * stop when sum_j (known[j]-1) + 1 <= LIMIT (5000): then every solution lies in the exhaustive list."""
import sys, math, json, time, itertools
from fractions import Fraction
from prove_iter import matveev_partial, ordered, FastTop
from solver import valid_points
from fastbottom import BottomContext
from twoadic import tmax
from bags import bags_of

LIMIT = 5000
BOTTOM_COST = 25.0                       # one bottom solve ~ 25 top cases

def prod(it):
    r = 1
    for x in it: r *= x
    return r

def prove_bag(bag, log=print, limit=LIMIT):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    known = {}
    stats = {'top': 0, 'bottom': 0, 'enum': 0, 'proofs': 0}
    t0 = time.time(); cache = {}
    while True:
        XL = matveev_partial(k, tm, known)
        G = math.ceil(2 * math.log10(XL)) + 12
        if len(known) == k and sum(v - 1 for v in known.values()) + 1 <= limit:
            Lmax = sum(v - 1 for v in known.values()) + 1
            log(f"  bag {bag}: k={k} tmax={tm} bounds {[known[j] for j in range(k)]} -> every solution has "
                f"L <= {Lmax} (top {stats['top']}, bottom {stats['bottom']}, enum {stats['enum']}, "
                f"proofs {stats['proofs']}, {time.time()-t0:.0f}s)")
            return Lmax, [known[j] for j in range(k)]
        # candidate actions
        best = None
        for j in range(k):
            if j in known and known[j] <= G + 4: continue           # nothing to gain
            if all(i in known for i in range(j)):
                c = prod(known[i] - 1 for i in range(j)) * max(1, len(ordered(bag, j)))
                if best is None or c < best[0]: best = (c, j, 'top')
            if all(i in known for i in range(j + 1, k)):
                b = k - j
                c = 2 * prod(known[i] - 1 for i in range(j + 1, k)) * len(ordered(bag, b)) * BOTTOM_COST
                if best is None or c < best[0]: best = (c, j, 'bottom')
        assert best is not None
        _, j, side = best
        key = (XL, G)
        if key not in cache:
            U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
            cache.clear(); cache[key] = (FastTop(XL, tm, G), BottomContext(G, tm, U, W))
        top, ctx = cache[key]
        ok = prove_top(j, bag, G, top, known, stats) if side == 'top' else prove_bottom(j, bag, G, ctx, known, stats)
        if not ok:
            raise RuntimeError(f"gap {j} ({side}) not provable at G={G} for bag {bag}")
        stats['proofs'] += 1
        known[j] = min(known.get(j, 10**18), G)

def prove_top(j, bag, G, top, known, stats):
    lat = top.lat
    if j == 0:
        if lat.kills_homogeneous(): return True
        stats['enum'] += 1
        return not valid_points(lat, 0, Fraction(1))
    ranges = [range(1, known[i]) for i in range(j)]
    for digs in ordered(bag, j):
        cs = [9 * (d - 1) for d in digs]
        for gaps in itertools.product(*ranges):
            s = sum(gaps)
            N = 10**s; acc = 0
            for c, g in zip(cs, gaps):
                acc += g; N += c * 10**(s - acc)
            stats['top'] += 1
            m = top.target(N, s)
            if not top.kills(m):
                stats['enum'] += 1
                if valid_points(lat, m, Fraction(N, 10**s)): return False
    return True

def prove_bottom(j, bag, G, ctx, known, stats):
    k = len(bag); b = k - j
    # gaps below gap j: g_{j+1} .. g_{k-1}; the topmost of them (g_{j+1}) is the inner loop
    for digs in ordered(bag, b):                       # top-to-bottom: d_{j+1} .. d_k
        cs = [9 * (d - 1) for d in digs]
        lower = list(reversed(cs[1:]))                 # coefficients bottom-up: d_k, d_{k-1}, ..., d_{j+2}
        # lower gaps bottom-up: g_{k-1}, g_{k-2}, ..., g_{j+2}
        lranges = [range(1, known[i]) for i in range(k - 1, j + 1, -1)]
        for ak in (0, 1):
            for lowgaps in itertools.product(*lranges):
                posl = [ak]
                for g in lowgaps: posl.append(posl[-1] + g)
                if b == 1:
                    Flow, pts = -1, [ak]
                else:
                    Flow = -1 + sum(c * 10**p for c, p in zip(lower, posl))
                    pts = [posl[-1] + g for g in range(1, known[j + 1])]
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
            Lmax, bounds = prove_bag(b, log=lambda x: print(x, flush=True))
            Lb = max(Lb, Lmax); recs.append([list(b), Lmax, bounds])
        assert Lb <= 5000
        bset = {tuple(sorted(b)) for b in bs}
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                     and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
        json.dump({s: {'preimages': pre, 'Lbound': Lb, 'bags': recs}}, open(f'proofJ_{s}.json', 'w'))
