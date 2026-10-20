"""(1) enumeration/verdict agreement with a Python reference; (2) planted: genuine numbers from the exhaustive
list must make the cumulative bottom test fail exactly when they satisfy its hypothesis."""
import math, json, random, itertools
from prove_iter import ordered
from fastbottom import BottomContext
from twoadic import tmax, v2
import cbottom2
random.seed(3)
# (1) reference comparison at proof-like parameters
def ref_configs(j, k, heff):
    b = k - j; out = []
    def rec(r, lo, pos):
        if r < 0: out.append(tuple(pos)); return
        for a in range(lo, heff[r]):
            rec(r - 1, a + 1, [a] + pos)
    # positions for r = b-1 (bottom) .. 0 (top); build from the bottom
    def rec2(r, lo, acc):
        if r < 0: out.append(tuple(acc)); return
        for a in range(lo, heff[r]): rec2(r - 1, a + 1, [a] + acc)
    rec2(b - 1, 0, [])
    return out          # tuples (a_{j+1}, ..., a_k)
XL = 86191814072450236416; G = math.ceil(2 * math.log10(XL)) + 12
U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
for bag, j, bot in [((2, 3, 7), 1, {2: 52, 3: 2}), ((2, 3, 7, 9), 1, {2: 52, 3: 52, 4: 2}), ((2, 4, 7, 8, 9), 2, {3: 60, 4: 52, 5: 2})]:
    k = len(bag); tm = tmax(bag)
    heff = cbottom2.heff_of(j, k, bot)
    orders = ordered(bag, k - j)
    st = {'bottom': 0}
    okc = cbottom2.prove_bottom_cum(j, bag, G, tm, U, W, bot, st, orders)
    ctx = BottomContext(G, tm, U, W); ctx.Mmax = max(ctx.Mmax, heff[0])
    ctx.__init__(G, tm, U, W) if ctx.Mmax == 4 * G + 8 else None
    cfgs = ref_configs(j, k, heff)
    n = 0; okp = True
    for o in orders:
        cs = [9 * (d - 1) for d in o]
        for pos in cfgs:
            Flow = -1 + sum(c * 10**a for c, a in zip(cs[1:], pos[1:]))
            p = pos[0]; M = max(G, p + 1)
            n += 1
            if ctx.solutions(Flow, cs[0], p, M): okp = False
    print(f"bag {bag} j={j} heff={heff}: C {okc} ({st['bottom']} solves)  Python {okp} ({n} solves)")
# (2) planted genuine solutions with a small box
Z = json.load(open('Z.json'))
XLs = 400; Gs = 12; Us = int(2.0960 * XLs + 3); Ws = int(1.1834 * XLs + 2)
tested = fails = expect = 0
for y in Z:
    s = str(y)
    if y < 10 or y % 5 == 0 or '0' in s or len(s) > XLs: continue
    L = len(s)
    pos = [(L - 1 - i, int(ch)) for i, ch in enumerate(s) if ch != '1']
    k = len(pos)
    if k < 2 or k > 5: continue
    bag = tuple(sorted(d for _, d in pos)); tm = tmax(bag)
    a = [L] + [p for p, _ in pos]                    # a[0] = L, a[i] = position of digit i
    for j in range(k - 1, -1, -1):
        hyp = a[j + 1] < Gs + 3 and a[j] >= max(Gs, a[j + 1] + 1)     # genuine y satisfies the test's hypothesis
        if a[j + 1] >= Gs + 3: break
        bot = {i: Gs + 3 for i in range(j + 1, k + 1)}
        orders = [tuple(d for _, d in pos[j:])]
        st = {'bottom': 0}
        ok = cbottom2.prove_bottom_cum(j, bag, Gs, tm, Us, Ws, bot, st, orders)
        tested += 1; expect += hyp
        if hyp and ok: print("MISSED planted", y, j); 
        if not ok: fails += 1
    if tested > 400: break
print(f"planted: {tested} bottom tests on genuine numbers; hypothesis satisfied in {expect}; tests failing {fails} "
      f"(must be >= {expect}: every genuine solution satisfying the hypothesis is detected)")
