import math, time, subprocess, os, tempfile
from basesolver import Base, cofactors
import prove_cum_b as P
base = Base(6)
def job(tb, Ws, j, bounds, orders, debug):
    d = tb.data
    lines = [f"prec {tb.prec}"] + [f"{kk} {d[kk]}" for kk in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
    lines += [f"base {base.b}", f"ntg {len(Ws)}"] + [f"W{i} {w}" for i, w in enumerate(Ws)]
    lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, bounds)), f"nord {len(orders)}"]
    lines += [" ".join(str((base.b - 1) * (x - 1)) for x in o) for o in orders] + ["skip 0", "maxsurv 0", f"debug {debug}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    return path
for bag, XL, j, bounds in [((2, 2, 4, 4, 5), 2 * 10**15, 3, [40, 45, 50]), ((2,) * 9, 10**15, 4, [30, 32, 34, 36]), ((2, 4, 5, 5), 10**12, 3, [27, 28, 29]), ((2,)*12, 5 * 10**42, 3, [65, 66, 67])]:
    Bs, pb = cofactors(bag, base)
    Xc = XL * math.log(6) / math.log(5) + 2
    G = math.ceil((math.log(Xc) + 7 * math.log(10)) / math.log(6)) + 2
    tb = P.TopB(base, XL, G)
    orders = P.ordered(bag, j)
    Ws = [5 * B for B in Bs]
    path = job(tb, Ws, j, bounds, orders, 0); path2 = job(tb, Ws, j, bounds, orders, 2)
    t0 = time.time(); a = subprocess.run(['csrc/topb', path], capture_output=True, text=True).stdout; ta = time.time() - t0
    t0 = time.time(); bb = subprocess.run(['csrc/topb_fast', path], capture_output=True, text=True).stdout; tbb = time.time() - t0
    r2 = subprocess.run(['csrc/topb_fast', path2], capture_output=True, text=True)
    sa = sorted(l for l in a.splitlines() if l.startswith('SURV')); sb = sorted(l for l in bb.splitlines() if l.startswith('SURV'))
    # compare survivor (config, target) sets; m may differ by the allowed +-1 at exact ties only
    ka = {tuple(l.split()[:-1]) for l in sa}; kb = {tuple(l.split()[:-1]) for l in sb}
    print(f"bag {bag} j={j} G={G}: old {len(sa)} survivors {ta:.2f}s, fast {len(sb)} survivors {tbb:.2f}s, same set {ka == kb}, "
          f"debug cross-check rc={r2.returncode} {r2.stderr.strip()[:60]}")
    os.unlink(path); os.unlink(path2)
