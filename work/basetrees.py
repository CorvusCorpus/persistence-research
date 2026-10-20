"""Nonzero-final-digit trees in base b, from an exhaustive list of zero-free smooth numbers.

In base b every digit product is S_b-smooth (S_b = primes < b).  A node of a nonzero-target tree is
a digit product whose f_b-orbit stays zero-free and ends at a nonzero digit.  We list every
S_b-smooth y < b^D not divisible by b (b | y would end in 0), keep the zero-free ones, and follow
their orbits."""
import sys, json, math, itertools
from collections import defaultdict
import flint

def primes_below(b):
    return [p for p in range(2, b) if all(p % q for q in range(2, p))]

def smooth_not_div(b, D):
    """All S_b-smooth y < b^D with b not dividing y."""
    P = primes_below(b); bound = b**D
    out = []
    def rec(i, cur):
        if i == len(P):
            if cur % b: out.append(cur)
            return
        x = cur
        while x < bound:
            rec(i + 1, x)
            x *= P[i]
    rec(0, 1)
    return out

def digits(y, b):
    return [int(ch, 36) for ch in flint.fmpz(y).str(b)]

def fprod(y, b):
    r = 1
    for d in digits(y, b): r *= d
    return r

def bags_of(s, b):
    """Multisets of digits in 2..b-1 with product s."""
    D = list(range(b - 1, 1, -1)); out = []
    def rec(x, start, cur):
        if x == 1: out.append(tuple(cur)); return
        for i in range(start, len(D)):
            d = D[i]
            if x % d == 0: rec(x // d, i, cur + [d])
    rec(s, 0, [])
    return out

def trees(b, D):
    ys = smooth_not_div(b, D)
    zf = [y for y in ys if '0' not in flint.fmpz(y).str(b)]
    nodes = defaultdict(dict)                      # target -> {node: depth}
    for y in zf:
        x, depth, ok = y, 0, True
        while x >= b:
            x = fprod(x, b); depth += 1
            if x == 0: ok = False; break
            if x >= b and '0' in flint.fmpz(x).str(b): ok = False; break
        if ok and x != 0:
            nodes[x][y] = depth
    return len(ys), len(zf), nodes

if __name__ == "__main__":
    out = {}
    for b, D in ((3, 3000), (4, 3000), (5, 800), (6, 800), (7, 200), (8, 200), (9, 200)):
        n, nz, nodes = trees(b, D)
        info = {}
        line = []
        for d in sorted(nodes):
            ns = nodes[d]
            H = max(ns.values())
            km = {s: max((len(t) for t in bags_of(s, b)), default=0) for s in ns}
            big = sorted((k, s) for s, k in km.items() if k > 9)
            info[d] = {'nodes': sorted(ns), 'height': H, 'kmax': max(km.values()), 'big': big[:8], 'nbig': len(big)}
            line.append(f"δ={d}: {len(ns)} nodes, height {H} (max persistence {H+1}), largest digit set {max(km.values())}, nodes >9 digits: {len(big)}")
        print(f"base {b} (searched below {b}^{D}; {n} smooth numbers, {nz} zero-free):")
        for l in line: print("   ", l)
        out[b] = {str(k): v for k, v in info.items()}
    json.dump(out, open('basetrees.json', 'w'))
