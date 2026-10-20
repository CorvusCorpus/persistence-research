import math, time
from prove_iter import matveev_partial, ordered
from fastbottom import BottomContext
import prove_iter2, cbottom
from twoadic import tmax
for bag, j, kn in [((2, 3, 7, 8, 9), 2, 40), ((2, 2, 3, 4, 7, 9), 3, 30), ((3, 6, 7, 9), 3, 52), ((2, 4, 7, 7, 8), 1, 20)]:
    k = len(bag); tm = tmax(bag)
    known = {i: kn for i in range(j + 1, k)}
    XL = matveev_partial(k, tm, {i: 52 for i in range(k) if i != j}); G = math.ceil(2 * math.log10(XL)) + 12
    U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
    ctx = BottomContext(G, tm, U, W)
    orders = ordered(bag, k - j)
    st = {'top': 0, 'bottom': 0, 'enum': 0}
    t0 = time.time(); okp = prove_iter2.prove_bottom(j, bag, G, ctx, known, st); tp = time.time() - t0
    st2 = {'bottom': 0}
    t0 = time.time(); okc = cbottom.prove_bottom_c(j, bag, G, tm, U, W, known, st2, orders); tc = time.time() - t0
    print(f"bag {bag} j={j} G={G}: Python {okp} ({st['bottom']} solves, {tp:.1f}s)  C {okc} ({st2['bottom']} solves, {tc:.2f}s)  speedup {tp/tc:.1f}x")
