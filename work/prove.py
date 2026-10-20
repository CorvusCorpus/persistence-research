"""Prove the complete solution set of 'y = 2^t3^u7^w with non-1 digit multiset = bag' for bags
with k <= 5 digits, via the case split described in solver.py."""
import sys, math, json, time, itertools
from fractions import Fraction
from solver import TopLattice, matveev_bound, solve_bottom
from twoadic import tmax
from bags import bags_of

def ordered(bag, r):
    return sorted(set(itertools.permutations(bag, r)))

def prove_bag(bag, JT=2, verbose=True):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    XL, B = matveev_bound(k, tm)
    G = math.ceil(2 * math.log10(XL)) + 12
    while True:
        lat = TopLattice(XL, tm, G)
        ok = True; stats = {'top': 0, 'bottom': 0}
        # ---- top cases: first large gap at index j <= JT
        if not lat.kills_homogeneous():
            from solver import valid_points
            vp = valid_points(lat, 0, Fraction(1))
            if vp: ok = False; print("  homogeneous case has genuine points:", vp)
            else: stats['enum'] = stats.get('enum', 0) + 1
        for j in range(1, min(JT, k - 1) + 1):
            if not ok: break
            for digs in ordered(bag, j):
                for gaps in itertools.product(range(1, G), repeat=j):
                    z = Fraction(0); s = 0
                    for d, g in zip(digs, gaps):
                        s += g; z += Fraction(9 * (d - 1), 10**s)
                    stats['top'] += 1
                    m = lat.target_m(1 + z)
                    if not lat.kills(m):
                        from solver import valid_points
                        vp = valid_points(lat, m, 1 + z)
                        stats['enum'] = stats.get('enum', 0) + 1
                        if vp:
                            ok = False; print("  top branch has genuine points:", digs, gaps, vp); break
                if not ok: break
        # ---- bottom cases: last large gap at index j_bot > JT (needs k-1 > JT)
        Umax = int(2.0960 * XL + 3); Wmax = int(1.1834 * XL + 2)
        for jb in range(JT + 1, k):
            if not ok: break
            b = k - jb                                   # number of digits below the large gap
            for digs in ordered(bag, b):                 # top-to-bottom order d_{jb+1..k}
                for ak in (0, 1):
                    for gaps in itertools.product(range(1, G), repeat=b - 1):
                        # positions from the bottom: a_k = ak, then add gaps upward
                        posns = [ak]
                        for g in reversed(gaps): posns.append(posns[-1] + g)
                        posns = posns[::-1]              # aligned with digs (top to bottom)
                        M = posns[0] + G
                        r = (10**M - 1) // 9 + sum((d - 1) * 10**p for d, p in zip(digs, posns))
                        stats['bottom'] += 1
                        sols = solve_bottom(r, M, tm, Umax, Wmax)
                        if sols:
                            ok = False; print("  bottom branch has candidates:", digs, posns, sols); break
                    if not ok: break
                if not ok: break
        if ok:
            Lmax = k * (G - 1) + 1
            if verbose:
                print(f"  bag {bag}: k={k} tmax={tm} XL={XL:.2e} G={G} -> every solution has "
                      f"L <= {Lmax}  (top cases {stats['top']}, bottom cases {stats['bottom']}, exact enumerations {stats.get('enum',0)})")
            return Lmax, G, tm, XL, stats
        G += 6
        print(f"  retry with G={G}")
        if G > 400: raise RuntimeError('giving up on bag ' + str(bag))

if __name__ == "__main__":
    Z = json.load(open('Z.json'))
    nodes = [int(x) for x in sys.argv[1:]]
    out = {}
    for s in nodes:
        t0 = time.time()
        bs = bags_of(s)
        print(f"node {s}: {len(bs)} digit multisets")
        Lbound = 0; recs = []
        for b in bs:
            Lmax, G, tm, XL, st = prove_bag(b)
            Lbound = max(Lbound, Lmax); recs.append((b, Lmax, G, tm, st))
        # all solutions are among exhaustively known numbers (< 10^5000)
        assert Lbound <= 5000
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y) and
                     tuple(sorted(int(c) for c in str(y) if c != '1')) in {tuple(sorted(b)) for b in bs})
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)")
        out[s] = {'preimages': pre, 'Lbound': Lbound, 'bags': [[list(b), L, G, tm] for b, L, G, tm, st in recs]}
    json.dump(out, open('proof_' + '_'.join(map(str, nodes)) + '.json', 'w'))
