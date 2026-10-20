import sys, json
import flint
from basesolver import Base, prove_bag
from basetrees import trees, bags_of, fprod, smooth_not_div
b, D, KMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
base = Base(b)
n, nz, nodes = trees(b, D)
zf_list = [y for y in smooth_not_div(b, D) if '0' not in flint.fmpz(y).str(b)]
res = {}
for d in sorted(nodes):
    for s in sorted(nodes[d], key=lambda x: nodes[d][x]):
        bs = bags_of(s, b)
        km = max((len(t) for t in bs), default=0)
        if km > KMAX:
            res[s] = {'target': d, 'depth': nodes[d][s], 'open': True, 'kmax': km}
            print(f"base {b}: node {s} (target {d}) skipped: digit sets up to {km}", flush=True)
            continue
        Lb = max(prove_bag(t, base, D - 1, log=lambda x: None) for t in bs)
        pre = sorted(y for y in zf_list if y >= b and fprod(y, b) == s)
        res[s] = {'target': d, 'depth': nodes[d][s], 'open': False, 'preimages': pre, 'Lbound': Lb}
        print(f"base {b}: node {s} (target {d}): all smooth y with digit product {s}: {pre} (L <= {Lb})", flush=True)
json.dump({str(k): v for k, v in res.items()}, open(f'base{b}_partial.json', 'w'))
for d in sorted({v['target'] for v in res.values()}):
    ns = {s: v for s, v in res.items() if v['target'] == d}
    openn = [s for s, v in ns.items() if v['open']]
    if not openn:
        H = max(v['depth'] for v in ns.values())
        print(f"THEOREM base {b}, final digit {d}: tree complete ({len(ns)} nodes), max persistence {H+1}")
    else:
        print(f"base {b}, final digit {d}: open nodes {openn}")
