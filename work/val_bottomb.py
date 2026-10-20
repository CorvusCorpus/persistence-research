"""C bottomb vs Python (solve_c over all cofactors) on identical configuration spaces, bases 3, 4, 6."""
import math, random, time
from basesolver import Base, cofactors, solve_c
import prove_cum_b as P, cbottomb
random.seed(4)
LN = math.log
for b, bag, j, botv in [(6, (2, 2, 4, 5), 1, 30), (6, (2, 4, 4, 5, 5), 2, 34), (6, (2, 2, 2, 2, 2, 2), 2, 36), (6, (3, 3, 5), 0, 31),
                        (4, (2, 3, 3), 1, 40), (3, (2, 2, 2, 2), 1, 60), (6, (2, 4), 2, 33)]:
    base = Base(b); k = len(bag)
    Bs, pb = cofactors(bag, base)
    XL = 10**15 if b == 6 else 10**12
    Xc = XL * LN(b) / LN(base.q) + 2
    G = math.ceil((LN(Xc) + 7 * LN(10)) / LN(b)) + 2
    Tp = max([pb[p] for p in pb] + [0]); Gb = max(G, Tp + 6)
    bot = {i: botv for i in range(j + 1, k + 1)}
    heff = P.heff_of(j, k, bot)
    orders = P.ordered(bag, k - j) if k - j else [()]
    st = {'bottom': 0}
    t0 = time.time(); okc = cbottomb.bottom_test_c(j, bag, base, pb, int(Xc), Gb, bot, st, orders, heff); tc = time.time() - t0
    st2 = {'bottom': 0}
    t0 = time.time(); okp = P.bottom_test_py(j, bag, base, Bs, XL, Gb, bot, st2); tp = time.time() - t0
    print(f"base {b} bag {bag} j={j} Gb={Gb} |Bs|={len(Bs)}: C {okc} ({st['bottom']} configs, {tc:.2f}s)  Python {okp} ({st2['bottom']} (config,B) solves, {tp:.1f}s)")
