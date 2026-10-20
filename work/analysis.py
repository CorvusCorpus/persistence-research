import math
from math import prod
from functools import lru_cache

def P(x): return prod(int(d) for d in str(x))
def chain(x):
    out=[x]
    while x>9: x=P(x); out.append(x)
    return out

# 1. targets of record holders (OEIS A003001)
rec=[0,10,25,39,77,679,6788,68889,2677889,26888999,3778888999,277777788888899]
print("record holders -> final digit:", [(len(chain(r))-1, chain(r)[-1]) for r in rec])

# 2. lifting lemma check: zero-free residue classes of 2^n mod 10^k
for k in range(1,9):
    per=4*5**(k-1); M=10**k
    x=pow(2,k,M); z=0
    for n in range(k,k+per):
        if '0' not in str(x).zfill(k): z+=1
        x=x*2%M
    print(f"k={k}: zero-free classes {z} of {per}  (>= 4^k = {4**k}: {z>=4**k})")

# 3. factor the two persistence-10 7-smooth numbers, min number of non-1 digits of a preimage
def fac(n):
    e={}
    for p in (2,3,5,7):
        while n%p==0: n//=p; e[p]=e.get(p,0)+1
    assert n==1; return e
digs={2:(1,0,0),3:(0,1,0),4:(2,0,0),6:(1,1,0),7:(0,0,1),8:(3,0,0),9:(0,2,0)}
@lru_cache(None)
def mind(a,b,c):
    if a==b==c==0: return 0
    best=10**9
    for d,(x,y,z) in digs.items():
        if x<=a and y<=b and z<=c: best=min(best,1+mind(a-x,b-y,c-z))
    return best
for T in (4996238671872, 937638166841712):
    e=fac(T); print(T, e, "min non-1 digits in a preimage:", mind(e.get(2,0),e.get(3,0),e.get(7,0)), "persistence", len(chain(T))-1)

# 4. heuristic: expected count of zero-free members of the two families with more than D digits
l=math.log
cA=1/(6*l(2)*l(3)*l(7)); cB=1/(6*l(5)*l(3)*l(7))
def dens(n):  # approx number of family members with exactly n digits
    L=math.log(10); f=lambda m:(m*L)**3*(cA+cB)
    return f(n)-f(n-1)
def expect(D,upto=None):
    s=0.0;n=D+1
    while True:
        t=dens(n)*0.9**(n-2)
        s+=t; n+=1
        if t<1e-300*max(s,1e-300) or n>D+20000: break
    return s
for D in (20,40,60,100,145,200,300):
    print(f"expected zero-free family members with > {D} digits: {expect(D):.3g}")
for D in (5000,30000):
    # log10 of leading term sum, computed in logs to avoid underflow
    lg=math.log10(dens(D+1))+(D-1)*math.log10(0.9)-math.log10(1-0.9)
    print(f"expected zero-free family members with > {D} digits: ~10^{lg:.0f}")
