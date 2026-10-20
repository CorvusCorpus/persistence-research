"""Cumulative bottom test for base 10 with candidate resolution (csrc/bottomc2).  Every exponent vector
(t, u', w) in the box compatible with a configuration is reported; it is discarded when the corresponding
y = 2^t 3^(u'-2) 7^w cannot be a solution: one of its digits at heights M..M+K-1 is 0 or not among the
digits left after removing the configuration from the bag (with multiplicity); short y are checked exactly."""
import os, subprocess, tempfile
from collections import Counter
import cbottom2
from cbottom import consts
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'bottomc2')

def run(G, tm, U, W, j, k, heff, orders, skip=0):
    Mmax = max(4 * G + 8, heff[0], G)
    c = consts(Mmax)
    lines = ["cum", f"G {G}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
    lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (d - 1)) for d in o) for o in orders] + [f"skip {skip}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        rc = p.wait()
        if rc: raise RuntimeError(f"C bottom (cum, all candidates) failed rc={rc}")
    finally:
        os.unlink(path)

def candidate_possible(bag, cfg, M, t, up, w, K=40):
    import math
    if up < 2: return True
    if (t * math.log10(2) + (up - 2) * math.log10(3) + w * math.log10(7)) < M + K + 12:
        y = 2**t * 3**(up - 2) * 7**w
        ds = str(y)
        return '0' not in ds and sorted(int(ch) for ch in ds if ch != '1') == sorted(bag)
    mod = 10**(M + K)
    yy = pow(2, t, mod) * pow(3, up - 2, mod) * pow(7, w, mod) % mod // 10**M
    rem = Counter(bag) - Counter(cfg)
    for _ in range(K):
        d = yy % 10; yy //= 10
        if d == 0: return False
        if d != 1:
            rem[d] -= 1
            if rem[d] < 0: return False
    return True

def prove_bottom_cum(j, bag, G, tm, U, W, bot, stats, orders, skip=0, progress=lambda d: None):
    k = len(bag)
    heff = cbottom2.heff_of(j, k, bot)
    ok = None; possible = False
    for f in run(G, tm, U, W, j, k, heff, orders, skip):
        if f[0] == 'PROG': progress(int(f[1]))
        elif f[0] == 'CAND':
            o = int(f[1]); i = f.index('p')
            M = int(f[i + 3]); t = int(f[i + 5]); up = int(f[i + 7]); w = int(f[i + 9])
            if candidate_possible(bag, orders[o], M, t, up, w):
                possible = True; stats['bottom_fail'] = f
            else:
                stats['resolved_bottom'] = stats.get('resolved_bottom', 0) + 1
        elif f[0] == 'FAIL': ok = False; stats['bottom_fail'] = f
        elif f[0] == 'DONE': ok = True; stats['bottom'] += int(f[1])
    if ok is None: raise RuntimeError("C bottom gave no verdict")
    return ok and not possible
