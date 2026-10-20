import math, itertools, random
from prove_iter import FastTop, matveev_partial, ordered
import ctop2
random.seed(5)
def comps(bounds):
    """prefixes g_0..g_{j-1} with s_{i+1} <= bounds[i]-1"""
    out = []
    def rec(i, s, cur):
        if i == len(bounds): out.append(tuple(cur)); return
        for g in range(1, bounds[i] - s):
            rec(i + 1, s + g, cur + [g])
    rec(0, 0, []); return out
k, tm = 9, 14
known = {i: 52 for i in range(8)}
XL = matveev_partial(k, tm, known); G = math.ceil(2 * math.log10(XL)) + 12
top = FastTop(XL, tm, G)
for j, bounds, bag in [(1, [52], (2, 3)), (2, [52, 52], (2, 3, 7)), (3, [40, 52, 45], (2, 3, 7)), (4, [30, 35, 40, 52], (2, 2, 9, 9))]:
    orders = ordered(bag, j)
    Ms = {}
    for f in ctop2.run_top(top, j, bounds, orders, debug=2, cum=True):
        if f[0] == 'M': Ms[(int(f[1]),) + tuple(map(int, f[2:-1]))] = int(f[-1])
    ref = {(o,) + c for o in range(len(orders)) for c in comps(bounds)}
    sample = random.sample(sorted(ref), min(2000, len(ref)))
    bad = sum(1 for key in sample if top.target(*ctop2.config_N(orders[key[0]], key[1:])) != Ms.get(key))
    print(f"j={j} bounds={bounds}: C cases {len(Ms)}, reference {len(ref)}, same set {set(Ms) == ref}, m mismatches {bad}/{len(sample)}")
