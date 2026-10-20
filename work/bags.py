import json
DIG = {2:(1,0,0),3:(0,1,0),4:(2,0,0),6:(1,1,0),7:(0,0,1),8:(3,0,0),9:(0,2,0)}
ORDER = [9,8,7,6,4,3,2]
def fac(x):
    e = []
    for p in (2,3,5,7):
        k = 0
        while x % p == 0: x //= p; k += 1
        e.append(k)
    assert x == 1
    return e
def bags_of(s):
    a,b,c,d = fac(s)
    if c: return []          # node divisible by 5: no 5-free preimage (handled separately)
    out = []
    def rec(e, start, cur):
        if e == (0,0,0): out.append(tuple(cur)); return
        for i in range(start, len(ORDER)):
            dg = ORDER[i]; v = DIG[dg]
            if all(e[j] >= v[j] for j in range(3)):
                rec(tuple(e[j]-v[j] for j in range(3)), i, cur+[dg])
    rec((a,b,d), 0, [])
    return out
if __name__ == "__main__":
    nodes = json.load(open('even_nodes.json'))
    for dlt in ('2','4','6','8'):
        rows = []
        for s in nodes[dlt]:
            if s < 10 and str(s) != dlt and s != int(dlt): pass
            bs = bags_of(s)
            kmax = max((len(b) for b in bs), default=0)
            rows.append((kmax, s, len(bs)))
        rows.sort()
        small = [r for r in rows if r[0] <= 5]; big = [r for r in rows if r[0] > 5]
        print(f"target {dlt}: {len(rows)} nodes; max bag size <=5: {len(small)}; >5: {len(big)}")
        print("   big nodes (kmax, node, #bags):", big)
