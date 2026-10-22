"""Dry-run plan of one bag in prove_cum's cost model: every step (side, level j, G, XL, bit length of the top lattice
constant C, cases) in the order prove_bag would run them.  usage (from work/): python3 bag_plan.py D1_D2_..."""
import sys, math
import prove_cum, cbottom2
from twoadic import tmax

def plan(bag):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    snaps = []; orig = prove_cum.matveev_cum
    def spy(k_, tm_, top, bot):
        XL = orig(k_, tm_, top, bot); snaps.append((dict(top), dict(bot), XL)); return XL
    prove_cum.matveev_cum = spy
    try: L, res = prove_cum.prove_bag(bag, log=lambda x: None, dry=True)
    finally: prove_cum.matveev_cum = orig
    rows = []
    for (t0, b0, XL), (t1, b1, _) in zip(snaps, snaps[1:]):
        G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
        U = int(2.0960 * XL + 3); C = U * 10**G // 80 + 1
        dt = [j for j in t1 if t0.get(j) != t1[j]]; db = [j for j in b1 if b0.get(j) != b1[j]]
        if dt:
            j = dt[0]; n = (prove_cum.ncomp(prove_cum.deff_top(j, t0)) if j else 1); no = max(1, len(prove_cum.ordered(bag, j)))
            rows.append(('top', j, G, XL, C.bit_length(), n * no, no, prove_cum.deff_top(j, t0) if j else []))
        else:
            j = db[0]; heff = cbottom2.heff_of(j, k, b0); no = len(prove_cum.ordered(bag, k - j))
            rows.append(('bottom', j, Gb, XL, None, prove_cum.nbottom(heff) * no, no, heff))
    return L, res, rows

if __name__ == '__main__':
    bag = tuple(int(x) for x in sys.argv[1].split('_'))
    L, res, rows = plan(bag)
    print(f"bag {bag}: k={len(bag)} tmax={tmax(tuple(sorted(bag)))} final L <= {L}, cost {res['stats']['cost']:.3e}")
    for side, j, G, XL, cb, cases, no, extra in rows:
        print(f"  {side:6s} j={j:2d} G={G:3d} XL={XL:.3e} Cbits={cb} cases={cases:.3e} orders={no} ranges/heff={extra}")
