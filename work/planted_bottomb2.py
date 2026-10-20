"""Planted test for csrc/bottomb: for random (B, c) take the last M base-b digits of y = B q^c (resampled until
zero-free) as a bottom configuration; the C solver must report exactly this exponent c as its candidate."""
import random, os
from basesolver import Base
import prove_cum_b as P, cbottomb
random.seed(9)
for b, pbv, Gb, Xc in [(6, {2: 8, 3: 0}, 15, 10**9), (6, {2: 0, 3: 6}, 13, 10**8), (4, {2: 5}, 16, 10**8), (3, {3: 0}, 24, 10**9)]:
    base = Base(b); q = base.q
    tested = found = 0
    while tested < 40:
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
        pos = [a for a in range(M - 1, -1, -1) if ds[a] != 1]
        digs = tuple(ds[a] for a in pos)
        if not digs: continue
        heff = [a + 1 for a in pos]
        out = list(cbottomb.run(base, pbv, Xc, Gb, len(digs), 0, heff, [digs]))
        nb = len(digs)
        hit = False
        for f in out:
            if f[0] == 'FAIL': print("unexpected FAIL", f)
            if f[0] == 'CAND' and int(f[3 + nb + 3]) == c: hit = True
            if f[0] == 'WEAK':
                c0, N = int(f[3 + nb + 3]), int(f[3 + nb + 5])
                if c >= c0 and (c - c0) % N == 0: hit = True
        cands = [f for f in out if f[0] in ('CAND', 'WEAK')]
        tested += 1; found += hit
        if not hit: print("MISSED", b, B, c, digs, pos, [f[-3:] for f in cands][:3])
    print(f"base {b} (Gb={Gb}, Xc={Xc}): planted exponent reported as candidate {found}/{tested}")
