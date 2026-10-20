"""Run after zf7: exact verification of survivors + structure of the zero-free 7-smooth set."""
import sys, math
from collections import Counter
import verify as V

K, fname, D = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
zl, pers = V.final(K, fname)
zs = set(zl)

# chains: zero-free 7-smooth m whose digit product is again zero-free 7-smooth (and >= 10)
edges = [(m, V.P(m)) for m in zl if m > 9 and V.P(m) in zs and V.P(m) > 9]
print(f"edges m -> P(m) staying inside the set (both >= 10): {len(edges)}")

T = {x for x in zl if pers[x] == 10}
print("persistence-10 members:", sorted(T))
pre = [m for m in zl if V.P(m) in T]
print("members whose digit product is a persistence-10 member (would give persistence 12):", pre)

# empirical vs heuristic count by digit length
by_len = Counter(len(str(x)) for x in zl)
L = math.log(10)
cA = 1 / (6 * math.log(2) * math.log(3) * math.log(7)); cB = 1 / (6 * math.log(5) * math.log(3) * math.log(7))
f = lambda m: (m * L) ** 3 * (cA + cB)
for lo, hi in ((10, 20), (21, 40), (41, 60), (61, 100), (101, D)):
    obs = sum(v for k, v in by_len.items() if lo <= k <= hi)
    exp = sum((f(n) - f(n - 1)) * 0.9 ** (n - 2) for n in range(lo, hi + 1))
    print(f"digits {lo:>3}-{hi:<5}: found {obs:>5}   naive-heuristic {exp:9.1f}")
print("longest zero-free 7-smooth numbers:")
for x in zl[-3:]:
    e = {p: 0 for p in (2, 3, 5, 7)}
    y = x
    for p in e:
        while y % p == 0: y //= p; e[p] += 1
    print(f"  {len(str(x))} digits = 2^{e[2]} 3^{e[3]} 5^{e[5]} 7^{e[7]}, persistence {pers[x]}")
