"""Base-10 cumulative bottom test with candidate resolution and fine-grained resumability (csrc/bottomc3):
one C run per (ordering, height of the lowest digit).  Statement identical to cbottom3.prove_bottom_cum."""
import os, subprocess, tempfile
import cbottom2
from cbottom import consts
from cbottom3 import candidate_possible

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'bottomc3')

def run(G, tm, U, W, j, k, heff, orders, aklo, akhi):
    Mmax = max(4 * G + 8, heff[0], G)
    c = consts(Mmax)
    lines = ["cum", f"G {G}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
    lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (d - 1)) for d in o) for o in orders] + ["skip 0", f"aklo {aklo}", f"akhi {akhi}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("C bottomc3 failed")
    finally:
        os.unlink(path)

def count_ak(heff, v):
    """number of height sequences a_{j+1} > ... > a_k = v with a < heff[r]"""
    b = len(heff)
    if v >= heff[b - 1]: return 0
    cnt = {v: 1}
    for r in range(b - 2, -1, -1):
        nc = {}
        for a0, c in cnt.items():
            for a in range(a0 + 1, heff[r]): nc[a] = nc.get(a, 0) + c
        cnt = nc
    return sum(cnt.values())

def prove_bottom_cum(j, bag, G, tm, U, W, bot, stats, orders, skip=0, progress=lambda d: None):
    k = len(bag)
    heff = cbottom2.heff_of(j, k, bot)
    nb = k - j
    o0, i0 = (skip, 0) if isinstance(skip, int) else skip
    from ctop3 import groups
    for oo in range(o0, len(orders)):
        for ak, ak1 in groups(i0 if oo == o0 else 0, heff[nb - 1] - 1, lambda v: count_ak(heff, v), 50000):
            progress((oo, ak))
            ok = None
            for f in run(G, tm, U, W, j, k, heff, [orders[oo]], ak, ak1):
                if f[0] == 'CAND':
                    i = f.index('p')
                    M = int(f[i + 3]); t = int(f[i + 5]); up = int(f[i + 7]); w = int(f[i + 9])
                    if candidate_possible(bag, orders[oo], M, t, up, w):
                        stats['bottom_fail'] = f; return False
                    stats['resolved_bottom'] = stats.get('resolved_bottom', 0) + 1
                elif f[0] == 'FAIL': stats['bottom_fail'] = f; return False
                elif f[0] == 'DONE': ok = True; stats['bottom'] += int(f[1])
            if ok is None: raise RuntimeError("C bottom gave no verdict")
    return True
