"""Driver for the cumulative-height mode of the C bottom solver (csrc/bottomc)."""
import subprocess, os, tempfile
from cbottom import consts

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'csrc', 'bottomc')

def heff_of(j, k, bot):
    """effective bounds a_i < heff[i] for i = j+1..k from the proved height bounds bot[i] (a_i strictly decreasing)."""
    out = []
    for i in range(j + 1, k + 1):
        out.append(min(bot[ip] - (i - ip) for ip in range(j + 1, i + 1)))
    return out                                    # out[r] <-> digit index j+1+r

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
        if rc: raise RuntimeError(f"C bottom (cum) failed rc={rc}")
    finally:
        os.unlink(path)

def prove_bottom_cum(j, bag, G, tm, U, W, bot, stats, orders, skip=0, progress=lambda d: None):
    k = len(bag)
    heff = heff_of(j, k, bot)
    ok = None
    for f in run(G, tm, U, W, j, k, heff, orders, skip):
        if f[0] == 'PROG': progress(int(f[1]))
        elif f[0] == 'FAIL': ok = False; stats['bottom_fail'] = f
        elif f[0] == 'DONE': ok = True; stats['bottom'] += int(f[1])
    if ok is None: raise RuntimeError("C bottom gave no verdict")
    return ok
