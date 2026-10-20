"""prove_iter2 with the C inner loops (csrc/top, csrc/bottom).  Identical proof chain logic, facts and
checkpoint files (ckpt/, bagres/); only the case loops run in C, and the cost ratio of a bottom solve to a
top case is re-measured for the C code (~4 instead of 25)."""
import sys, json, time
from fractions import Fraction
import prove_iter2
from prove_iter import ordered
from solver import valid_points
import ctop, cbottom

def prove_top(j, bag, G, top, known, stats, skip=0, progress=lambda d: None):
    lat = top.lat
    if j == 0:
        if lat.kills_homogeneous(): return True
        stats['enum'] += 1
        return not valid_points(lat, 0, Fraction(1))
    return ctop.prove_top_c(j, bag, G, top, known, stats, ordered(bag, j), skip, progress)

def prove_bottom(j, bag, G, ctx, known, stats, skip=0, progress=lambda d: None):
    k = len(bag)
    return cbottom.prove_bottom_c(j, bag, G, ctx.tm, ctx.U, ctx.W, known, stats, ordered(bag, k - j), skip, progress)

prove_iter2.prove_top = prove_top
prove_iter2.prove_bottom = prove_bottom
prove_iter2.BOTTOM_COST = 4.0
prove_iter2.SOURCE = 'prove_iter3 run (C inner loops)'
prove_bag = prove_iter2.prove_bag

if __name__ == "__main__":
    from bags import bags_of
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
