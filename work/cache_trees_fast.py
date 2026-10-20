"""Cache nonzero-target trees and the zero-free smooth list for a base; the zero-free test first rejects
numbers with a 0 among their lowest 24 base-b digits (cheap), then checks the full expansion."""
import sys, json, time
import flint
from collections import defaultdict
from basetrees import smooth_not_div, fprod
b, D = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
ys = smooth_not_div(b, D)
zf = []
for y in ys:
    x = y; ok = True
    for n in range(24):
        if x == 0: break
        if x % b == 0: ok = False; break
        x //= b
    if ok and '0' not in flint.fmpz(y).str(b): zf.append(y)
nodes = defaultdict(dict)
for y in zf:
    x, depth, ok = y, 0, True
    while x >= b:
        x = fprod(x, b); depth += 1
        if x == 0: ok = False; break
        if x >= b and '0' in flint.fmpz(x).str(b): ok = False; break
    if ok and x != 0: nodes[x][y] = depth
json.dump({'b': b, 'D': D, 'n_smooth': len(ys), 'n_zerofree': len(zf),
           'nodes': {str(d): {str(s): dp for s, dp in v.items()} for d, v in nodes.items()},
           'zerofree': [str(y) for y in zf]}, open(f'trees_b{b}_D{D}.json', 'w'))
print(f"base {b}: {len(ys)} smooth, {len(zf)} zero-free, {sum(len(v) for v in nodes.values())} nodes; {time.time()-t0:.0f}s")
