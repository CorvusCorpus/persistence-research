"""Base-10 cumulative bottom test on the GPU (csrc/bottomgpu): one run per step, all orderings from `skip`, every height.
Statement and candidate resolution identical to cbottom4.prove_bottom_cum (the GPU binary makes the same decisions as
csrc/bottomc3 and prints the same CAND lines; NOTEBOOK and val_bottomgpu.py). Resumable at ordering boundaries: progress
reports (ordering, 0), which cbottom4 also accepts. Jobs the kernel does not cover are run by bottomc3 inside bottomgpu."""
import os, subprocess, tempfile
import cbottom2
from cbottom import consts
from cbottom3 import candidate_possible

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', os.environ.get('BOTTOMGPU_BIN', 'bottomgpu'))
CPU = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'bottomc3')

def run(G, tm, U, W, j, k, heff, orders, skip):
    Mmax = max(4 * G + 8, heff[0], G)
    c = consts(Mmax)
    lines = ["cum", f"G {G}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
    lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (d - 1)) for d in o) for o in orders] + [f"skip {skip}", "aklo 0", f"akhi {heff[-1] - 1}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1, env=dict(os.environ, BG_CPU=CPU))
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("bottomgpu failed")
    finally:
        os.unlink(path)

def prove_bottom_cum(j, bag, G, tm, U, W, bot, stats, orders, skip=0, progress=lambda d: None):
    k = len(bag)
    heff = cbottom2.heff_of(j, k, bot)
    o0 = skip if isinstance(skip, int) else skip[0]          # a (ordering, height) checkpoint restarts its ordering
    progress((o0, 0))
    done = None; last = 0                                    # solves are added per finished ordering (PROG), as checkpoints
    for f in run(G, tm, U, W, j, k, heff, orders, o0):       # record (ordering, 0) with the stats up to that ordering
        if f[0] == 'CAND':
            o = int(f[1]); i = f.index('p')
            M = int(f[i + 3]); t = int(f[i + 5]); up = int(f[i + 7]); w = int(f[i + 9])
            if candidate_possible(bag, orders[o], M, t, up, w):
                stats['bottom_fail'] = f; return False
            stats['resolved_bottom'] = stats.get('resolved_bottom', 0) + 1
        elif f[0] == 'PROG':
            cur = int(f[2]); stats['bottom'] += cur - last; last = cur
            progress((int(f[1]), 0))
        elif f[0] == 'DONE': done = int(f[1])
    if done is None: raise RuntimeError("bottomgpu gave no verdict")
    stats['bottom'] += done - last
    return True
