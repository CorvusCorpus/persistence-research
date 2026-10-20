"""Batched base-10 tests (topc3/bottomc3) vs unbatched (topc/bottomc2): identical counts and verdicts."""
import math
import prove_cum
from prove_iter import FastTop, matveev_partial
from twoadic import tmax
import ctop2, ctop3, cbottom3, cbottom4
for bag, j, tv, bv in [((2, 3, 7, 8, 9), 2, 52, 52), ((2, 2, 3, 3, 7, 9), 3, 45, 52), ((2,) * 7, 3, 52, 52)]:
    k = len(bag); tm = tmax(bag)
    XL = matveev_partial(k, tm, {i: 52 for i in range(k - 1)}); G = math.ceil(2 * math.log10(XL)) + 12
    top = FastTop(XL, tm, G)
    known = {i: tv - i for i in range(j)}
    orders = prove_cum.ordered(bag, j)
    s1 = {'top': 0, 'enum': 0}; o1 = ctop2.prove_top_c(j, bag, G, top, known, s1, orders, cum=True)
    s2 = {'top': 0, 'enum': 0}; o2 = ctop3.prove_top_c(j, bag, G, top, known, s2, orders)
    U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
    bot = {i: bv for i in range(j + 1, k)}; bot[k] = 2
    ob = prove_cum.ordered(bag, k - j)
    s3 = {'bottom': 0}; o3 = cbottom3.prove_bottom_cum(j, bag, G, tm, U, W, bot, s3, ob)
    s4 = {'bottom': 0}; o4 = cbottom4.prove_bottom_cum(j, bag, G, tm, U, W, bot, s4, ob)
    print(f"{bag} j={j}: top {o1}/{o2} cases {s1['top']}/{s2['top']}; bottom {o3}/{o4} solves {s3['bottom']}/{s4['bottom']}")
