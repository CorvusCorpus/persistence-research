"""Python driver for the C bottom-side solver (csrc/bottom)."""
import subprocess, os, tempfile
from solver import plog

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'csrc', 'bottom')
_consts = {}

def consts(Mm):
    if Mm not in _consts:
        _consts[Mm] = dict(L5_3=plog(pow(3, 4, 5**(Mm + 1)), 5, Mm), L5_7=plog(pow(7, 4, 5**(Mm + 1)), 5, Mm),
                           L5_16=plog(16, 5, Mm), L2_3=plog(-3 % 2**(Mm + 2), 2, Mm), L2_7=plog(-7 % 2**(Mm + 2), 2, Mm))
    return _consts[Mm]

def header(mode, G, tm, U, W, Mmax):
    c = consts(Mmax)
    return [mode, f"G {G}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]

def _run(lines):
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        rc = p.wait()
        if rc: raise RuntimeError(f"C bottom failed rc={rc}")
    finally:
        os.unlink(path)

def selftest(G, tm, U, W, Mmax, queries):
    lines = header("selftest", G, tm, U, W, Mmax) + [f"n {len(queries)}"] + [f"{F} {c} {p} {M}" for F, c, p, M in queries]
    return [None if f[1] == '-' else int(f[1]) for f in _run(lines) if f[0] == 'R']

def mneed(j, k, known, G):
    b = k - j
    return 1 + G if b == 1 else 1 + sum(known[i] - 1 for i in range(j + 1, k)) + G

def prove_bottom_c(j, bag, G, tm, U, W, known, stats, orders, skip=0, progress=lambda d: None):
    k = len(bag)
    Mmax = max(4 * G + 8, mneed(j, k, known, G))
    lines = header("proof", G, tm, U, W, Mmax) + [f"k {k}", f"j {j}",
             "known " + " ".join(str(known.get(i, -1)) for i in range(k)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (d - 1)) for d in o) for o in orders]
    lines += [f"skip {skip}"]
    ok = None
    for f in _run(lines):
        if f[0] == 'PROG': progress(int(f[1]))
        elif f[0] == 'FAIL': ok = False; stats['bottom_fail'] = f
        elif f[0] == 'DONE': ok = True; stats['bottom'] += int(f[1])
    if ok is None: raise RuntimeError("C bottom gave no verdict")
    return ok
