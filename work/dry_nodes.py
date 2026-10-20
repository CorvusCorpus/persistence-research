import json, math, sys, time
import prove_cum
from bags import bags_of
nodes = json.load(open('even_nodes.json'))
rows = []
for d, ns in nodes.items():
    for s, depth in (ns.items() if isinstance(ns, dict) else [(x, None) for x in ns]):
        s = int(s)
        bs = bags_of(s)
        if not bs: continue
        km = max(len(b) for b in bs)
        if km < 10 or km > 16: continue
        t0 = time.time(); tot = 0
        for b in bs:
            L, res = prove_cum.prove_bag(b, log=lambda x: None, dry=True)
            tot += res['stats']['cost']
        rows.append((tot, s, d, depth, km, len(bs)))
        print(f"node {s} (target {d}, depth {depth}): kmax {km}, {len(bs)} bags, est. cost {tot:.3g}  [{time.time()-t0:.0f}s]", flush=True)
rows.sort()
print("cheapest first:", [(s, km, f"{c:.2g}") for c, s, d, dp, km, nb in rows])
