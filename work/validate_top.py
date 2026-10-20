import json, math
from fractions import Fraction
import flint
from solver import TopLattice
flint.ctx.prec = 2000
Z = json.load(open('Z.json'))
def fac(x):
    e=[]
    for p in (2,3,5,7):
        k=0
        while x%p==0: x//=p; k+=1
        e.append(k)
    return e
checked = 0; worst = 0
for y in Z:
    s = str(y); L = len(s)
    if L < 3 or '5' in s: continue
    t,u,_,w = fac(y)
    pos = [(L-1-i, int(ch)) for i,ch in enumerate(s) if ch != '1']   # (position, digit), descending
    k = len(pos)
    if k == 0 or k > 6: continue
    gaps = [L - pos[0][0]] + [pos[i][0]-pos[i+1][0] for i in range(k-1)]
    assert pos[-1][0] <= 1, ("F1 violated", y)
    for j in range(k):
        z = Fraction(0)
        s_acc = 0
        for i in range(j):
            s_acc += gaps[i]
            z += Fraction(9*(pos[i][1]-1), 10**s_acc)
        lam = (t*flint.arb(2).log() + (u+2)*flint.arb(3).log() + w*flint.arb(7).log()
               - L*flint.arb(10).log() - (flint.arb(1) + flint.arb(z.numerator)/z.denominator).log())
        g = gaps[j]
        assert lam > 0 and lam < 80 * flint.arb(10)**(-g), (y, j, lam)
        checked += 1
print("F1 and F3 hold on", checked, "linear forms from real data")

# radius check: the true exponent vector must lie within R of the target when G = actual gap
bad = 0; n = 0
for y in Z:
    s = str(y); L = len(s)
    if L < 6 or '5' in s: continue
    t,u,_,w = fac(y)
    pos = [(L-1-i, int(ch)) for i,ch in enumerate(s) if ch != '1']
    k = len(pos)
    if k == 0 or k > 6: continue
    gaps = [L - pos[0][0]] + [pos[i][0]-pos[i+1][0] for i in range(k-1)]
    j = max(range(k), key=lambda i: gaps[i]); G = gaps[j]
    if G < 3: continue
    lat = TopLattice(XL=L+5, T=max(t,1), G=G)
    z = Fraction(0); s_acc = 0
    for i in range(j):
        s_acc += gaps[i]; z += Fraction(9*(pos[i][1]-1), 10**s_acc)
    m = lat.target_m(1 + z)
    lt = lat.S_t // max(lat.T,1); lw = lat.S_w // ((lat.S_w // round(lat.S_u/(int(1.1834*lat.XL+2)))) or 1) if False else None
    # rebuild unreduced vector v(x)-y
    U = int(2.0960*lat.XL+3); W = int(1.1834*lat.XL+2); S = U
    lw = max(1, round(S/W)); lt = max(1, S//max(lat.T,1))
    last = t*lat.rC[2] + (u+2)*lat.rC[3] + w*lat.rC[7] - L*lat.rC[10] + m
    vec2 = (lt*t)**2 + (u+2)**2 + (lw*w)**2 + last**2
    n += 1
    if vec2 > lat.R2: bad += 1; print("RADIUS VIOLATION", y)
    if lat.kills(m): bad += 1; print("FALSE KILL", y)
print(f"radius/false-kill check on {n} real solutions: {bad} problems")
