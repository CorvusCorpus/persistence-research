"""Cache nonzero-target trees and the zero-free smooth list for a base (one-time cost)."""
import sys, json, time
import flint
from basetrees import trees, smooth_not_div
b, D = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
n, nz, nodes = trees(b, D)
zf = [y for y in smooth_not_div(b, D) if '0' not in flint.fmpz(y).str(b)]
json.dump({'b': b, 'D': D, 'n_smooth': n, 'n_zerofree': nz, 'nodes': {str(d): {str(s): dp for s, dp in v.items()} for d, v in nodes.items()},
           'zerofree': [str(y) for y in zf]}, open(f'trees_b{b}_D{D}.json', 'w'))
print(f"base {b}: {n} smooth, {nz} zero-free, {sum(len(v) for v in nodes.values())} nodes; {time.time()-t0:.0f}s")
