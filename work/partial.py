import sys, json, time
from prove_iter2 import prove_bag
from bags import bags_of
s = int(sys.argv[1]); kmax = int(sys.argv[2])
bs = [b for b in bags_of(s) if len(b) <= kmax]
print(f"node {s}: proving {len(bs)} of {len(bags_of(s))} digit sets (those with <= {kmax} non-1 digits)", flush=True)
recs = []
for b in sorted(bs, key=len):
    Lmax, bounds = prove_bag(b, log=lambda x: print(x, flush=True))
    recs.append([list(b), Lmax, bounds])
    json.dump({'node': s, 'kmax': kmax, 'proved_bags': recs}, open(f'partial_{s}.json', 'w'))
print(f"node {s}: all {len(bs)} digit sets with <= {kmax} non-1 digits have no solution beyond the exhaustive list", flush=True)
