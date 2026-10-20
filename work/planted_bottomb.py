"""Planted test for csrc/bottomb: choose (B, c), take the last M base-b digits of y = B q^c (resampled until
they are zero-free) as a bottom configuration; the cumulative bottom test must find a solution (FAIL)."""
import random
from basesolver import Base
import prove_cum_b as P, cbottomb
random.seed(9)
for b, pbv, Gb, Xc in [(6, {2: 8, 3: 0}, 15, 10**9), (6, {2: 0, 3: 6}, 13, 10**8), (4, {2: 5}, 16, 10**8), (3, {3: 0}, 24, 10**9)]:
    base = Base(b); q = base.q
    tested = detected = 0
    while tested < 40:
        # cofactor
        B = 1
        for p in base.bdiv:
            if pbv.get(p, 0) and random.random() < 0.7: B *= p**random.randint(1, pbv[p])
        if B % b == 0: continue
        c = random.randint(0, Xc)
        M = Gb
        y = B * pow(q, c, b**(M + 2)) % b**M
        ds = []; t = y
        for a in range(M):
            ds.append(t % b); t //= b
        if 0 in ds: continue
        pos = [a for a in range(M - 1, -1, -1) if ds[a] != 1]          # top to bottom
        digs = tuple(ds[a] for a in pos)
        if not digs: continue
        k = len(digs); j = 0
        heff = [a + 1 for a in pos]                                    # contains the genuine configuration
        st = {'bottom': 0}
        ok = cbottomb.bottom_test_c(j, tuple(sorted(digs)), base, pbv, Xc, Gb, {}, st, [digs], heff)
        tested += 1; detected += (not ok)
        if ok: print("MISSED", b, B, c, digs, pos)
    print(f"base {b} (Gb={Gb}, Xc={Xc}): planted {detected}/{tested} detected")
