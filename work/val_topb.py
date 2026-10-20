"""Validate csrc/topb (base 6, 2-dim lattices, several cofactor targets) against exact Python enumeration:
for a 2-dim lattice, level-2 failure <=> some lattice point within R of the target."""
import math, random, itertools
from fractions import Fraction
import flint
from basesolver import Base, cofactors, Top2
from solver import close_points, _gs
import prove_cum_b as P
random.seed(2)
base = Base(6)
for bag, XL in [((2, 2, 4, 4, 5), 2 * 10**15), ((2,) * 9, 10**15), ((2, 4, 5, 5), 10**12)]:
    Bs, pb = cofactors(bag, base)
    Xc = XL * math.log(6) / math.log(5) + 2
    G = math.ceil((math.log(Xc) + 7 * math.log(10)) / math.log(6)) + 2
    tb = P.TopB(base, XL, G); T = tb.T
    gs2, mu = _gs(tb.B)
    Bi = flint.fmpq_mat(flint.fmpz_mat(tb.B)).inv()
    Binv = [[Fraction(int(Bi[i, j].p), int(Bi[i, j].q)) for j in range(2)] for i in range(2)]
    j = 3; bounds = [40, 45, 50]
    orders = P.ordered(bag, j)[:6]
    Ws = [5 * B for B in Bs]
    surv = set(); cases = 0
    for f in P.run_topb(tb, base, Ws, j, bounds, orders):
        if f[0] == 'SURV': surv.add((int(f[1]),) + tuple(map(int, f[2:2 + j - 1])) + (int(f[-4]), int(f[-2])))
        if f[0] == 'DONE': cases = int(f[1])
    # exact reference on a random sample of (configuration, B) + all C survivors
    def comps(bounds):
        out = []
        def rec(i, s, cur):
            if i == len(bounds): out.append(tuple(cur)); return
            for g in range(1, bounds[i] - s): rec(i + 1, s + g, cur + [g])
        rec(0, 0, []); return out
    allc = [(o,) + c for o in range(len(orders)) for c in comps(bounds)]
    sample = random.sample(allc, 400)
    mism = 0; nclose = 0
    for cfg in sample + [s[:-1] for s in surv]:
        o, gaps = cfg[0], cfg[1:]
        s = sum(gaps); N = 6**s; acc = 0
        for x, g in zip(orders[o], gaps): acc += g; N += 5 * (x - 1) * 6**(s - acc)
        for ti, B in enumerate(Bs):
            flint.ctx.prec = T.prec
            kappa = flint.arb(5 * B).log() - (flint.arb(N) / flint.arb(6)**s).log()
            m = T.rnd(kappa)
            pts = close_points(tb.B, gs2, mu, Binv, [0, -m], T.R2)
            has = bool(pts); nclose += has
            if has != ((o,) + tuple(gaps) + (ti,) in surv): mism += 1
    print(f"bag {bag} G={G} |Bs|={len(Bs)}: C cases {cases}, C survivors {len(surv)}, checked {len(sample) + len(surv)} configs x {len(Bs)} targets, "
          f"with close points {nclose}, mismatches {mism}")
