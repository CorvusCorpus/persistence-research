"""Heuristic: expected number of 7-smooth m with digit product exactly x, calibrated on exact data."""
import json, math
from math import comb, log
from collections import Counter
import verify as V

Z = json.load(open('Z.json')); S = set(Z)
pers = {x: V.persistence(x) for x in Z}
LN10 = log(10)

# exact counts of family members per digit length and last-digit class (A1: 2^a.., a>=1; A0: 3^b7^c; B: 5^a.., a>=1)
cnt = Counter()
for tag, a, b, c, x in V.family(10**146):
    cls = 'B' if tag == 'B' else ('A1' if a >= 1 else 'A0')
    cnt[(cls, len(str(x)))] += 1
l2, l3, l5, l7 = (log(p) for p in (2, 3, 5, 7))
def members(cls, D):                     # number of class members with exactly D digits
    if D <= 145: return cnt[(cls, D)]
    f = {'A1': lambda L: L**3/(6*l2*l3*l7), 'A0': lambda L: L**2/(2*l3*l7), 'B': lambda L: L**3/(6*l5*l3*l7)}[cls]
    return f(D*LN10) - f((D-1)*LN10)
frac_last = {'A1': 4/10, 'A0': 4/10, 'B': 1/10}
def cls_of(e): return 'B' if e == 5 else ('A1' if e % 2 == 0 else 'A0')
def logrho(log_s, D, e):                 # log of chance that a specific integer s (D digits, last digit e) is 7-smooth
    c = cls_of(e); n = members(c, D)
    return log(n) - log_s - log(LN10) - log(frac_last[c]) if n > 0 else -math.inf

DIG = {2:(1,0,0,0),3:(0,1,0,0),4:(2,0,0,0),5:(0,0,1,0),6:(1,1,0,0),7:(0,0,0,1),8:(3,0,0,0),9:(0,2,0,0)}
def fac(x):
    e = []
    for p in (2,3,5,7):
        k = 0
        while x % p == 0: x //= p; k += 1
        e.append(k)
    assert x == 1; return tuple(e)

def expected_preimages(x, Dmax=400, split=None):
    ex = fac(x); K = sum(ex)
    divs = [(i,j,k,l) for i in range(ex[0]+1) for j in range(ex[1]+1) for k in range(ex[2]+1) for l in range(ex[3]+1)]
    # W[k][y] = ordered k-tuples of non-1 digits with exponent vector y
    W = [{(0,0,0,0): 1}]
    for k in range(1, K+1):
        nxt = Counter()
        for y, w in W[-1].items():
            for d, v in DIG.items():
                z = tuple(y[t]+v[t] for t in range(4))
                if all(z[t] <= ex[t] for t in range(4)): nxt[z] += w
        W.append(nxt)
    def M(n, y):                         # strings of length n over 1..9 with exponent vector y
        return sum(comb(n, k) * W[k].get(y, 0) for k in range(min(n, K)+1))
    sub = lambda y, d: tuple(y[t]-DIG.get(d,(0,0,0,0))[t] for t in range(4))
    ok = lambda y: all(v >= 0 for v in y)
    tot, beyond = 0.0, 0.0
    for D in range(2, Dmax+1):
        part = 0.0
        if D == 2:
            for f1 in range(1,10):
                for f2 in range(1,10):
                    if f1*f2 == x: part += math.exp(logrho(log(10*f1+f2), 2, f2))
        else:
            for f1 in range(1,10):
                y1 = sub(ex, f1)
                if not ok(y1): continue
                for f2 in range(1,10):
                    y2 = sub(y1, f2)
                    if not ok(y2): continue
                    for e in range(1,10):
                        y3 = sub(y2, e)
                        if not ok(y3): continue
                        m = M(D-3, y3)
                        if m: part += math.exp(log(m) + logrho(log(10*f1+f2+0.5) + (D-2)*LN10, D, e))
        tot += part
        if split and D > split: beyond += part
        if D > K + 12 and part < 1e-18 * max(tot, 1e-300): break
    return tot, beyond

if __name__ == "__main__":
    obs = Counter(pers[x] for x in Z if x >= 10)
    pred = Counter()
    per_node = {}
    for x in Z:
        if x >= 10 and pers[x] >= 3:
            e, _ = expected_preimages(x)
            per_node[x] = e
            pred[pers[x]+1] += e
    print("level  observed  predicted")
    for p in range(4, 12):
        print(f"{p:>5} {obs.get(p,0):>9} {pred.get(p,0):>10.2f}")
    json.dump({str(k): v for k, v in per_node.items()}, open('pernode.json','w'))
