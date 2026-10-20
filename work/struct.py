import json, verify as V
def _chain(x):
    out=[x]
    while x>9: x=V.P(x); out.append(x)
    return out
V.chain=_chain
zf,_ = V.small(144)            # exact; the 10^5000 run added nothing above 10^145
Z = sorted(zf | set(range(1,10)))
json.dump(Z, open('Z.json','w'))
S = set(Z); P = V.P; pers = {x: V.persistence(x) for x in Z}
digits = lambda x: set(str(x))
# Z' : even, no digit 0 or 5
Zp = [x for x in Z if x % 2 == 0 and '5' not in str(x)]
# E : m in Z' (m>=10) whose digit product is again in Z' and >= 10
E = [m for m in Zp if m >= 10 and P(m) >= 10 and P(m) in set(Zp)]
print("Z:", len(Z), "largest", len(str(Z[-1])), "digits")
print("Z' (even, no 0/5):", len(Zp), "largest", Zp[-1], len(str(Zp[-1])), "digits")
print("E (m and P(m) both in Z'):", len(E), "largest", E[-1], len(str(E[-1])), "digits")
# sanity: every element of Z with persistence >= 3 that is not in a chain ending in an odd digit
for p in range(3, 11):
    xs = [x for x in Z if pers[x] == p]
    tgt = {}
    for x in xs:
        c = V.chain(x) if hasattr(V,'chain') else None
    print(p, len(xs))
