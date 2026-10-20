import random, math, time, itertools
from prove_iter import matveev_partial, ordered
from fastbottom import BottomContext
import cbottom
random.seed(11)
k, tm = 9, 14
known = {i: 52 for i in range(8)}
XL = matveev_partial(k, tm, known); G = math.ceil(2 * math.log10(XL)) + 12
U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
ctx = BottomContext(G, tm, U, W); Mm = ctx.Mmax
# 1) planted solutions (as in fastbottom.__main__) + random realistic queries
qs = []; truth = []
for trial in range(400):
    t = random.randint(0, tm); up = random.randint(2, U); w = random.randint(0, W)
    p = random.randint(0, 3 * G); M = p + G
    if M > Mm: continue
    Fm = pow(2, t, 10**M) * pow(3, up, 10**M) * pow(7, w, 10**M) % 10**M
    c = random.choice([9, 18, 27, 45, 54, 63, 72])
    Flow = (Fm - c * 10**p) % 10**M
    if Flow % 5 == 0 or Flow == 0: continue
    qs.append((Flow, c, p, M)); truth.append(('planted', t))
for trial in range(1500):
    nb = random.randint(1, 4)
    pos = [random.randint(0, 1)]
    for _ in range(nb - 1): pos.append(pos[-1] + random.randint(1, 51))
    cs = [9 * (random.choice([2, 3, 4, 6, 7, 8, 9]) - 1) for _ in range(nb)]
    Flow = -1 + sum(c * 10**p for c, p in zip(cs[1:], pos[:-1])) if nb > 1 else -1
    p = pos[-1] if nb == 1 else pos[-1]
    M = p + G
    if M > Mm: continue
    qs.append((Flow, cs[0], p, M)); truth.append(('random', None))
t0 = time.time(); res_c = cbottom.selftest(G, tm, U, W, Mm, qs); tc = time.time() - t0
t0 = time.time(); res_py = [ctx.solutions(*q) for q in qs]; tp = time.time() - t0
agree = sum((rc is not None) == bool(rp) and (rc is None or rc in [s[0] for s in rp]) for rc, rp in zip(res_c, res_py))
planted_ok = sum(1 for (kind, t), rc in zip(truth, res_c) if kind == 'planted' and rc == t)
nplanted = sum(1 for kind, _ in truth if kind == 'planted')
print(f"{len(qs)} queries: C/Python agree on {agree}; planted recovered by C {planted_ok}/{nplanted}; "
      f"random with solutions: {sum(1 for (kd,_), rc in zip(truth, res_c) if kd=='random' and rc is not None)}; "
      f"time C {tc:.2f}s (incl. setup) vs Python {tp:.2f}s")
