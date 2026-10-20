import random, sys
from fastbottom import BottomContext
from solver import solve_bottom
random.seed(5)
G, tm = 80, 15; U, W = 10**33, 6*10**32
ctx = BottomContext(G, tm, U, W)
agree = n = 0
cases = []
for d in [2,3,4,6,7,8,9]:
    for ak in (0, 1):
        cases.append(([d], ak, []))
        for d2 in [2,3,4,6,7,8,9]:
            for g in (1, 2, 3, 7, 40, 79):
                cases.append(([d, d2], ak, [g]))
for digs, ak, gaps in cases:                       # digs bottom-up
    pos = [ak]
    for g in gaps: pos.append(pos[-1] + g)
    if len(digs) == 1: Flow, c, p = -1, 9*(digs[0]-1), ak
    else: Flow, c, p = -1 + 9*(digs[0]-1)*10**pos[0], 9*(digs[1]-1), pos[1]
    M = p + G
    fast = sorted(ctx.solutions(Flow, c, p, M))
    r = ((10**M - 1)//9 + sum((dd-1)*10**pp for dd, pp in zip(digs, pos))) % 10**M
    slow = sorted((t, u+2, w) for (t, u, w) in solve_bottom(r, M, tm, U-2, W))
    n += 1; agree += (fast == slow)
print(f"fast vs slow agree on {agree}/{n} small-position configurations")
ok = tot = 0
for trial in range(300):
    t = random.randint(0, tm); up = random.randint(2, U); w = random.randint(0, W)
    p = random.choice([0, 1, 2, 3]); M = p + G
    Fm = pow(2, t, 10**M) * pow(3, up, 10**M) * pow(7, w, 10**M) % 10**M
    c = random.choice([9, 18, 27, 45, 54, 63, 72]); Flow = (Fm - c*10**p) % 10**M
    if Flow % 5 == 0: continue
    tot += 1; ok += (t, up, w) in ctx.solutions(Flow, c, p, M)
print(f"planted small-p recovered {ok}/{tot}")
