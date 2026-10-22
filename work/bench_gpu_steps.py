"""Throughput of the GPU solvers on whole steps of a bag's bootstrapping plan (all orderings, every first gap / height).
usage (from work/): python3 bench_gpu_steps.py BAG [NSTEPS=6]   (the NSTEPS costliest steps; GPU binaries only, plus the
CPU binaries on one slice of each step for the per-case rate)"""
import sys, os, math, time, subprocess, tempfile
import prove_cum, cbottom2, prof_steps
from prove_iter import FastTop
from twoadic import tmax
from cbottom import consts
import val_topgpu

W = os.path.dirname(os.path.abspath(__file__))

def timed(exe, text, env=None):
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, text.encode()); os.close(fd)
    try:
        t0 = time.time(); r = subprocess.run([exe, path], capture_output=True, text=True, env=env); dt = time.time() - t0
        if r.returncode: raise RuntimeError(r.stderr[-300:])
        return dt, r.stdout.splitlines()[-1]
    finally:
        os.unlink(path)

def main(bag, nsteps):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    st = prof_steps.plan_states(bag)
    for i in sorted(sorted(range(len(st)), key=lambda i: -st[i][5])[:nsteps]):
        side, j, XL, top, bot, cost = st[i]
        G = math.ceil(2 * math.log10(XL)) + 12
        if side == 'top':
            if j == 0: continue
            T = FastTop(XL, tm, G); ranges = prove_cum.deff_top(j, top); orders = prove_cum.ordered(bag, j)
            text, d, prec = val_topgpu.job_text(T, j, ranges, orders, 1, ranges[0] - 1)
            dt, done = timed(os.path.join(W, 'csrc', 'topgpu'), text, env=dict(os.environ, TP_CPU=os.path.join(W, 'csrc', 'topc3')))
            n = int(done.split()[1])
            print(f"step {i:2d} top    j={j} G={G}: {n:.3e} cases in {dt:.1f} s = {n / dt:.3e} cases/s  [{done}]", flush=True)
        else:
            Gb = max(G, tm + 8); U = int(2.0960 * XL + 3); Wb = int(1.1834 * XL + 2)
            heff = cbottom2.heff_of(j, k, bot); orders = prove_cum.ordered(bag, k - j); Mmax = max(4 * Gb + 8, heff[0], Gb)
            c = consts(Mmax)
            lines = ["cum", f"G {Gb}", f"tm {tm}", f"U {U}", f"W {Wb}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
            lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(orders)}"]
            lines += [" ".join(str(9 * (dd - 1)) for dd in o) for o in orders] + ["skip 0", "aklo 0", f"akhi {heff[-1] - 1}"]
            dt, done = timed(os.path.join(W, 'csrc', 'bottomgpu'), "\n".join(lines) + "\n", env=dict(os.environ, BG_CPU=os.path.join(W, 'csrc', 'bottomc3')))
            n = int(done.split()[1])
            print(f"step {i:2d} bottom j={j} G={Gb}: {n:.3e} solves in {dt:.1f} s = {n / dt:.3e} solves/s  [{done}]", flush=True)

if __name__ == '__main__':
    main(tuple(int(x) for x in sys.argv[1].split('_')), int(sys.argv[2]) if len(sys.argv) > 2 else 6)
