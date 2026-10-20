import sys, math, time, itertools, random
from fractions import Fraction
from prove_iter import FastTop, matveev_partial, ordered
from solver import valid_points, close_points, _attach
import ctop2 as ctop
random.seed(7)
def check(k, tm, j, bag, known_val=52, sample_m=3000):
    known = {i: known_val for i in range(k - 1)}
    XL = matveev_partial(k, tm, known); G = math.ceil(2 * math.log10(XL)) + 12
    top = FastTop(XL, tm, G); lat = top.lat; _attach(lat)
    orders = ordered(bag, j)
    ranges = [known_val] * j
    # 1) m values agree on a sample
    Ms = {}
    for f in ctop.run_top(top, j, ranges, orders, debug=2):
        if f[0] == 'M':
            Ms[(int(f[1]),) + tuple(map(int, f[2:-1]))] = int(f[-1])
    keys = random.sample(list(Ms), min(sample_m, len(Ms)))
    bad = 0
    for key in keys:
        N, s = ctop.config_N(orders[key[0]], key[1:])
        if top.target(N, s) != Ms[key]: bad += 1
    # 2) level-1 failure sets agree (full enumeration), level-2 kills have no close lattice points
    F1 = {}
    for f in ctop.run_top(top, j, ranges, orders, debug=1):
        if f[0] == 'F1':
            F1[(int(f[1]),) + tuple(map(int, f[2:-1]))] = int(f[-1])
    py_fail = set()
    for o, digs in enumerate(orders):
        for gaps in itertools.product(*[range(1, r) for r in ranges]):
            N, s = ctop.config_N(digs, gaps)
            if not top.kills(top.target(N, s)): py_fail.add((o,) + gaps)
    l2bad = 0; l2killed = 0; l2surv = 0
    for key, sv in F1.items():
        N, s = ctop.config_N(orders[key[0]], key[1:])
        m = top.target(N, s)
        pts = close_points(lat.B, lat.gs2, lat._mu, lat._Binv, [0] * (lat.dim - 1) + [-m], lat.R2)
        if sv == 0:
            l2killed += 1
            if pts: l2bad += 1
        else:
            l2surv += 1
    print(f"k={k} tm={tm} j={j} G={G} dim={lat.dim}: cases {len(Ms)}, m mismatches {bad}/{len(keys)}, "
          f"L1 fail C {len(F1)} vs Py {len(py_fail)} equal={set(F1) == py_fail}, "
          f"L2 killed {l2killed} (with close points: {l2bad}), L2 survivors {l2surv}", flush=True)
check(9, 14, 2, (2, 3, 4, 6, 7, 8, 9))
check(7, 20, 2, (2, 2, 3, 7, 9))
check(5, 0, 2, (3, 3, 7, 9, 9))
check(9, 30, 3, (2, 3, 7, 8), known_val=30)
