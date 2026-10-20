"""Replay prove_iter2's deterministic proof chain for a bag (all tests succeed, as they did in the
recorded runs) and report the largest modulus exponent M used by any bottom solve vs Mmax = 4G+8."""
import sys, math, json, glob, os
from prove_iter import matveev_partial, ordered
from twoadic import tmax
LIMIT = 5000; BOTTOM_COST = 25.0
def prod(it):
    r = 1
    for x in it: r *= x
    return r
def replay(bag):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag); known = {}
    worst = (-10**9, None); steps = []
    while True:
        XL = matveev_partial(k, tm, known); G = math.ceil(2 * math.log10(XL)) + 12
        if len(known) == k and sum(v - 1 for v in known.values()) + 1 <= LIMIT:
            return worst, [known[j] for j in range(k)], steps
        best = None
        for j in range(k):
            if j in known and known[j] <= G + 4: continue
            if all(i in known for i in range(j)):
                c = prod(known[i] - 1 for i in range(j)) * max(1, len(ordered(bag, j)))
                if best is None or c < best[0]: best = (c, j, 'top')
            if all(i in known for i in range(j + 1, k)):
                b = k - j
                c = 2 * prod(known[i] - 1 for i in range(j + 1, k)) * len(ordered(bag, b)) * BOTTOM_COST
                if best is None or c < best[0]: best = (c, j, 'bottom')
        _, j, side = best
        if side == 'bottom':
            b = k - j
            if b == 1: Mx = 1 + G
            else: Mx = 1 + sum(known[i] - 1 for i in range(j + 1, k)) + G
            Mmax = 4 * G + 8
            if Mx - Mmax > worst[0]: worst = (Mx - Mmax, (j, G, Mx, Mmax, dict(known)))
        steps.append((j, side, G))
        known[j] = min(known.get(j, 10**18), G)
if __name__ == "__main__":
    bags = {}
    for fn in glob.glob('bagres/*.json'):
        r = json.load(open(fn)); bags[tuple(r['bag'])] = r['bounds']
    for fn in glob.glob('proofJ_*.json'):
        d = json.load(open(fn))
        for s, v in d.items():
            for rec in v['bags']: bags[tuple(sorted(rec[0]))] = rec[2]
    print(len(bags), "proved bags")
    bad = []; mism = 0
    for bag, bounds in sorted(bags.items(), key=lambda x: (len(x[0]), x[0])):
        worst, fb, steps = replay(bag)
        if fb != bounds: mism += 1; print("bound mismatch", bag, fb, bounds)
        if worst[0] > 0: bad.append((bag, worst))
    print("chain replays matching recorded bounds:", len(bags) - mism, "/", len(bags))
    print("bags with some bottom solve at M > Mmax:", len(bad))
    for b, w in bad[:30]: print("  ", b, "excess", w[0], "j,G,M,Mmax =", w[1][:4])
