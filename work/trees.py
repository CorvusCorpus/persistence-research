import json, verify as V
from collections import defaultdict
Z = json.load(open('Z.json'))            # all zero-free 7-smooth numbers < 10^5000 (plus 1..9)
P = V.P
def target(x):
    while x > 9: x = P(x)
    return x
def fac(x):
    e = []
    for p in (2,3,5,7):
        k = 0
        while x % p == 0: x //= p; k += 1
        e.append(k)
    return tuple(e)
nodes = defaultdict(list)
for y in Z:
    if y >= 1:
        nodes[target(y)].append(y)
for d in (2,4,6,8):
    ns = sorted(nodes[d])
    # parent relation y -> P(y)
    print(f"target {d}: {len(ns)} nodes; largest {ns[-1]} = 2^{fac(ns[-1])[0]} 3^{fac(ns[-1])[1]} 5^{fac(ns[-1])[2]} 7^{fac(ns[-1])[3]}")
    if d == 4: print("   nodes:", ns)
json.dump({d: sorted(nodes[d]) for d in (2,4,6,8)}, open('even_nodes.json','w'))
