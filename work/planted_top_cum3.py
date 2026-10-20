"""Planted test for the cumulative top test: genuine y (exhaustive list) with digits 1..j inside the top
D_{j-1} < G positions and digit j+1 at depth D_j >= G must make the top test of level j FAIL."""
import json, math
from twoadic import v2
from prove_iter import FastTop
import ctop3 as ctop2
Z = json.load(open('Z.json'))
XLs = 80
G = math.ceil(2 * math.log10(XLs)) + 12
tested = detected = 0
tops = {}
for y in Z:
    s = str(y)
    if y % 5 == 0 or '0' in s or len(s) > XLs: continue
    L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch) for ch in s if ch != '1']
    k = len(pos)
    for j in range(1, k):
        Dprev = L - pos[j - 1]; Dj = L - pos[j]
        if Dprev < G and Dj >= G:
            tm = v2(y) + 2
            if tm not in tops: tops[tm] = FastTop(XLs, tm, G)
            top = tops[tm]
            bounds = {i: (L - pos[i]) + 1 for i in range(j)}          # D_i < bounds[i]: genuine prefix included
            st = {'top': 0, 'enum': 0}
            ok = ctop2.prove_top_c(j, tuple(digs), G, top, bounds, st, [tuple(digs[:j])])
            tested += 1; detected += (not ok)
            if ok: print("MISSED", y, j)
            break
    if tested >= 150: break
print(f"G={G}: planted genuine top configurations: {detected}/{tested} detected")
