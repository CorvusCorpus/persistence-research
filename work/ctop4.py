"""Base-10 cumulative top test on the GPU (csrc/topgpu): one run per step, all orderings from `skip`, every first gap.
Statement and survivor handling identical to ctop3.prove_top_c (topgpu prints the same SURV runs and counters as topc3:
val_topgpu.py). Resumable at ordering boundaries: progress reports (ordering, 1), which ctop3 also accepts."""
import os, subprocess, tempfile
from fractions import Fraction
from ctop2 import lattice_job, config_N, candidates, could_be_solution
from solver import valid_points

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', os.environ.get('TOPGPU_BIN', 'topgpu'))
CPU = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'topc3')

def run_top(top, j, ranges, orders, skip, g0lo=1):
    d = lattice_job(top)
    prec = max(top.prec, d['C'].bit_length() + 64)
    lines = [f"prec {prec}"] + [f"{k} {d[k]}" for k in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
    lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, ranges)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (x - 1)) for x in o) for o in orders]
    lines += [f"skip {skip}", "maxsurv 10000000", "debug 0", f"g0lo {g0lo}", f"g0hi {ranges[0] - 1}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1, env=dict(os.environ, TP_CPU=CPU))
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("topgpu failed")
    finally:
        os.unlink(path)

def prove_top_c(j, bag, G, top, known, stats, orders, skip=0, progress=lambda d: None):
    lat = top.lat
    ranges = [known[i] for i in range(j)]
    cache = {}
    ok = True
    o0, g0 = (skip, 1) if isinstance(skip, int) else skip    # an (ordering, first gap) checkpoint: finish that ordering from
    runs = []                                                # its first gap g0 alone, then the rest, as ctop3 would resume
    if g0 > 1: runs.append(([orders[o0]], 0, g0, o0)); o0 += 1
    if o0 < len(orders): runs.append((orders, o0, 1, 0))
    for sel, sk, glo, base in runs:
        progress((base + sk if base == 0 else base, glo))
        if not _run(j, bag, top, ranges, sel, sk, glo, base, stats, cache, lat, progress): ok = False
    return ok

def _run(j, bag, top, ranges, sel, sk, glo, base, stats, cache, lat, progress):
    ok = True; last = 0; orders = sel
    for f in run_top(top, j, ranges, sel, sk, glo):
        if f[0] == 'SURV':                                   # as ctop3.prove_top_c
            o = int(f[1]); pre = list(map(int, f[2:2 + j - 1])); lo, hi, mC = int(f[-3]), int(f[-2]), int(f[-1])
            stats['enum'] += 1
            for gl in (lo, hi):
                N, s = config_N(orders[o], pre + [gl])
                if top.target(N, s) != mC: raise RuntimeError("C/Python target mismatch")
            if mC not in cache: cache[mC] = candidates(lat, mC)
            if cache[mC]:
                for gl in range(lo, hi + 1):
                    N, s = config_N(orders[o], pre + [gl])
                    for pt in valid_points(lat, top.target(N, s), Fraction(N, 10**s)):
                        if could_be_solution(bag, pt): ok = False
                        else: stats['resolved'] = stats.get('resolved', 0) + 1
        elif f[0] == 'PROG':                                 # cases counted per finished ordering, as checkpoints record
            cur = int(f[2]); stats['top'] += cur - last; last = cur
            progress((base + int(f[1]), 1))                  # (base: index of a single-ordering run's ordering)
        elif f[0] == 'DONE':
            stats['top'] += int(f[1]) - last; last = int(f[1])
    return ok
