"""Planted test for the cumulative bottom test: genuine zero-free 7-smooth y from the exhaustive list whose
digit j sits at height >= M = max(G, a_{j+1}+1) while digits j+1..k sit low must make the test FAIL."""
import json
from twoadic import v2
import cbottom3 as cbottom2
Z = json.load(open('Z.json'))
Gs = 14
tested = detected = 0
for y in Z:
    s = str(y)
    if y % 5 == 0 or '0' in s or len(s) > 60: continue
    L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch) for ch in s if ch != '1']
    k = len(pos); a = [L] + pos
    for j in range(k - 1, -1, -1):
        if k - j > 12: break
        if a[j + 1] < 17 and a[j] >= max(Gs, a[j + 1] + 1) and v2(y) <= Gs - 8:
            bag = tuple(sorted(digs)); tm = min(v2(y) + 3, Gs - 8)   # any upper bound >= the genuine t works here
            XLs = 60; U = int(2.0960 * XLs + 3); W = int(1.1834 * XLs + 2)
            heff = [a[j + 1 + r] + 1 for r in range(k - j)]
            st = {'bottom': 0}
            orders = [tuple(digs[j:])]
            bot = {j + 1 + r: heff[r] for r in range(k - j)}
            ok = cbottom2.prove_bottom_cum(j, bag, Gs, tm, U, W, bot, st, orders)
            tested += 1; detected += (not ok)
            if ok: print("MISSED", y, j, heff)
            break
    if tested >= 80: break
print(f"planted genuine configurations: {detected}/{tested} detected")
