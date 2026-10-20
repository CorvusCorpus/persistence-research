"""Case-split proof for bags with k <= 6 non-1 digits: top cases j_top <= 2 (lattice),
bottom cases j_bot >= 3 (fast exact 10-adic test, up to 3 digits below the large gap)."""
import sys, math, json, time, itertools
from fractions import Fraction
from solver import TopLattice, matveev_bound, valid_points
from fastbottom import BottomContext
from twoadic import tmax
from bags import bags_of

def ordered(bag, r):
    return sorted(set(itertools.permutations(bag, r)))

def prove_bag(bag, JT=2, verbose=True):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    assert k - JT - 1 <= 3
    XL, B = matveev_bound(k, tm)
    G = math.ceil(2 * math.log10(XL)) + 12
    while True:
        t0 = time.time()
        lat = TopLattice(XL, tm, G)
        ok = True; st = {'top': 0, 'bottom': 0, 'enum': 0}
        if not lat.kills_homogeneous():
            st['enum'] += 1
            vp = valid_points(lat, 0, Fraction(1))
            if vp: ok = False; print("  homogeneous case has genuine points:", vp)
        for j in range(1, min(JT, k - 1) + 1):
            if not ok: break
            for digs in ordered(bag, j):
                cs = [9 * (d - 1) for d in digs]
                for gaps in itertools.product(range(1, G), repeat=j):
                    z = Fraction(0); s = 0
                    for c, g in zip(cs, gaps):
                        s += g; z += Fraction(c, 10**s)
                    st['top'] += 1
                    m = lat.target_m(1 + z)
                    if not lat.kills(m):
                        st['enum'] += 1
                        vp = valid_points(lat, m, 1 + z)
                        if vp: ok = False; print("  top branch has genuine points:", digs, gaps, vp); break
                if not ok: break
        U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
        ctx = BottomContext(G, tm, U, W)
        for jb in range(JT + 1, k):
            if not ok: break
            b = k - jb
            for digs in ordered(bag, b):                    # top-to-bottom
                cs = [9 * (d - 1) for d in digs]
                for ak in (0, 1):
                    # lower positions (all but the topmost bottom digit), built upward from a_k
                    for lowgaps in itertools.product(range(1, G), repeat=max(b - 2, 0)):
                        posl = [ak]
                        for g in lowgaps: posl.append(posl[-1] + g)  # positions of digs[-1], digs[-2], ...
                        if b == 1:
                            Flow, ptop_choices = -1, [ak]
                        else:
                            lower = list(reversed(cs[1:]))            # bottom-up coefficients
                            Flow = -1 + sum(c * 10**p for c, p in zip(lower, posl))
                            ptop_choices = [posl[-1] + g for g in range(1, G)]
                        for p in ptop_choices:
                            st['bottom'] += 1
                            sols = ctx.solutions(Flow, cs[0], p, p + G)
                            if sols:
                                ok = False; print("  bottom branch has candidates:", digs, ak, lowgaps, p, sols); break
                        if not ok: break
                    if not ok: break
                if not ok: break
        if ok:
            Lmax = k * (G - 1) + 1
            if verbose:
                print(f"  bag {bag}: k={k} tmax={tm} XL={XL:.2e} G={G} -> every solution has L <= {Lmax} "
                      f"(top {st['top']}, bottom {st['bottom']}, exact enumerations {st['enum']}, {time.time()-t0:.0f}s)",
                      flush=True)
            return Lmax, G, tm, XL, st
        G += 6
        print(f"  retry with G={G}", flush=True)
        if G > 500: raise RuntimeError('giving up on ' + str(bag))

if __name__ == "__main__":
    Z = json.load(open('Z.json'))
    out = {}
    for s in map(int, sys.argv[1:]):
        t0 = time.time()
        bs = bags_of(s)
        print(f"node {s}: {len(bs)} digit multisets", flush=True)
        Lb = 0; recs = []
        for b in bs:
            Lmax, G, tm, XL, st = prove_bag(b)
            Lb = max(Lb, Lmax); recs.append([list(b), Lmax, G, tm])
        assert Lb <= 5000
        bset = {tuple(sorted(b)) for b in bs}
        pre = sorted(y for y in Z if y >= 10 and '0' not in str(y)
                     and tuple(sorted(int(c) for c in str(y) if c != '1')) in bset)
        print(f"node {s}: all 7-smooth y with digit product {s}: {pre}   ({time.time()-t0:.0f}s)", flush=True)
        out[s] = {'preimages': pre, 'Lbound': Lb, 'bags': recs}
        json.dump(out, open(f'proof6_{sys.argv[1]}.json', 'w'))
