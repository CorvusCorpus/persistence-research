import json, math, verify as V
from tree import expected_preimages
Z = json.load(open('Z.json'))
def chain(x):
    out=[x]
    while x>9: x=V.P(x); out.append(x)
    return out
# lemma check on data: for every n in Z with persistence p>=6
bad = 0
for n in Z:
    c = chain(n); p = len(c)-1
    if p < 6: continue
    tgt = c[-1]
    assert tgt % 2 == 0, (n, tgt)
    last_nonzero = p-2 if tgt == 0 else p-1          # indices whose element must be 0/5-free and even-digit
    for i in range(0, last_nonzero+0):
        s = str(c[i])
        if '0' in s or '5' in s or not set(s) & set('2468'): bad += 1
print("lemma violations in data:", bad)
T1, T2 = 4996238671872, 937638166841712
for T in (T1, T2):
    print(T, "P(T) =", V.P(T), "T in E:", all(ch not in str(y) for y in (T, V.P(T)) for ch in '05'))
    for cut in (15, 20, 30, 50, 100, 5000):
        tot, beyond = expected_preimages(T, Dmax=6000, split=cut)
        print(f"  expected 7-smooth preimages: total {tot:.4f}, with more than {cut} digits {beyond:.3g}")
# heuristic tail for E beyond 10^5000 (digits avoid 0,5 in m and in P(m))
L=math.log(10); cA=1/(6*math.log(2)*math.log(3)*math.log(7))
f=lambda m:(m*L)**3*cA
mu=math.log10(math.factorial(9)/5)/8
D=5000
lg=math.log10(f(D+1)-f(D))+(D-1)*math.log10(0.8)+mu*D*math.log10(0.8)
print(f"E-elements expected beyond 10^5000: ~10^{lg:.0f}")
