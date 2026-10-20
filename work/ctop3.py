"""Base-10 cumulative top test with fine-grained resumability (csrc/topc3): one C run per (ordering, first gap).
Statement and survivor handling identical to ctop2.prove_top_c(cum=True)."""
import os, subprocess, tempfile
from fractions import Fraction
import ctop2
from ctop2 import lattice_job, config_N, candidates, could_be_solution
from solver import valid_points

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'topc3')

def run_top(top, j, ranges, orders, g0lo, g0hi):
    d = lattice_job(top)
    prec = max(top.prec, d['C'].bit_length() + 64)
    lines = [f"prec {prec}"] + [f"{k} {d[k]}" for k in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
    lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, ranges)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (x - 1)) for x in o) for o in orders]
    lines += ["skip 0", "maxsurv 10000000", "debug 0", f"g0lo {g0lo}", f"g0hi {g0hi}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("C topc3 failed")
    finally:
        os.unlink(path)

def count_g0(ranges, v):
    """number of prefixes (g_0 = v, g_1, ..., g_{j-1}) with g_0+...+g_i <= ranges[i]-1"""
    if v > ranges[0] - 1: return 0
    cnt = {v: 1}
    for bd in ranges[1:]:
        nc = {}
        for s, c in cnt.items():
            for s2 in range(s + 1, bd): nc[s2] = nc.get(s2, 0) + c
        cnt = nc
    return sum(cnt.values())

def groups(counts_lo, vmax, count, target=200000):
    """consecutive value groups [a, b] with about `target` cases each"""
    out = []; a = counts_lo; acc = 0
    for v in range(counts_lo, vmax + 1):
        acc += count(v)
        if acc >= target:
            out.append((a, v)); a = v + 1; acc = 0
    if a <= vmax: out.append((a, vmax))
    return out

def prove_top_c(j, bag, G, top, known, stats, orders, skip=0, progress=lambda d: None):
    lat = top.lat
    ranges = [known[i] for i in range(j)]
    cache = {}
    ok = True
    o0, i0 = (skip, 1) if isinstance(skip, int) else skip
    for oo in range(o0, len(orders)):
        for g0, g1 in groups(i0 if oo == o0 else 1, ranges[0] - 1, lambda v: count_g0(ranges, v)):
            progress((oo, g0))
            for f in run_top(top, j, ranges, [orders[oo]], g0, g1):
                if f[0] == 'SURV':
                    pre = list(map(int, f[2:2 + j - 1])); lo, hi, mC = int(f[-3]), int(f[-2]), int(f[-1])
                    stats['enum'] += 1
                    for gl in (lo, hi):
                        N, s = config_N(orders[oo], pre + [gl])
                        if top.target(N, s) != mC: raise RuntimeError("C/Python target mismatch")
                    if mC not in cache: cache[mC] = candidates(lat, mC)
                    if cache[mC]:
                        for gl in range(lo, hi + 1):
                            N, s = config_N(orders[oo], pre + [gl])
                            for pt in valid_points(lat, top.target(N, s), Fraction(N, 10**s)):
                                if could_be_solution(bag, pt): ok = False
                                else: stats['resolved'] = stats.get('resolved', 0) + 1
                elif f[0] == 'DONE':
                    stats['top'] += int(f[1])
    return ok
