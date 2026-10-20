"""Planted test for the base-b cumulative top test (csrc/topb + exact Python re-check): genuine zero-free
y = B q^c in base b with digits 1..j inside the top D_{j-1} < G positions and digit j+1 at depth D_j >= G."""
import math, random
import flint
from basesolver import Base, cofactors
import prove_cum_b as P
random.seed(12)
for b in (6, 4, 3):
    base = Base(b); q = base.q
    XL = 60; Xc = XL * math.log(b) / math.log(q) + 2
    G = math.ceil((math.log(Xc) + 7 * math.log(10)) / math.log(b)) + 2
    tested = detected = 0; seen = set()
    for c in range(0, int(Xc)):
        for e in range(0, 12):
            for p in ([1] + base.bdiv):
                B = p**e if p > 1 else 1
                if (B, c) in seen or B % b == 0: continue
                seen.add((B, c))
                y = B * q**c
                s = flint.fmpz(y).str(b)
                if '0' in s or len(s) > XL or len(s) < G + 1: continue
                L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch, b) for ch in s if ch != '1']
                k = len(pos)
                for j in range(1, k):
                    if L - pos[j - 1] < G and L - pos[j] >= G:
                        bag = tuple(sorted(digs))
                        Bs = [B]                                 # the genuine cofactor
                        top = {i: (L - pos[i]) + 1 for i in range(j)}
                        st = {'top': 0, 'enum': 0}
                        orders_full = P.ordered(bag, j)
                        # restrict to the genuine ordering for speed: monkeypatch ordered
                        old = P.ordered
                        P.ordered = lambda bg, r, d=tuple(digs[:j]): [d]
                        try:
                            ok = P.top_test(j, bag, base, Bs, XL, G, top, st)
                        finally:
                            P.ordered = old
                        tested += 1; detected += (not ok)
                        if ok: print("MISSED", b, B, c, s, j)
                        break
                if tested >= 40: break
            if tested >= 40: break
        if tested >= 40: break
    print(f"base {b} (G={G}): planted genuine top configurations {detected}/{tested} detected")
