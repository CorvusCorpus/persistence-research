"""Throughput of the GPU step processes when N run at once, on a cost-weighted sample of one node's own steps.

prove_cum_gpu.py runs one topgpu/bottomgpu process per step in each of its workers, so the pool's speed is the aggregate rate
of N such processes sharing the GPU. This script draws steps of the node's dry-run plans (bag_plan/prof_steps) with
probability proportional to their cost, cuts each into a slice of about SLICE units (a subset of the orderings, every height),
writes the job files exactly as ctop4.run_top / cbottom5.run do (top jobs with the "fp" section), and runs them N at a time for
a fixed window. The steady rate counts each job's units in proportion to its overlap with [60 s, window]. The units are
prove_cum's: one per top case, BOTTOM_COST = 4 per bottom solve. Nothing is written to ckptC/, bagresC/ or any proof file.
usage (from work/): python3 pool_bench.py NODE WINDOW_S N [N ...]   env SLICE (units per job, default 1.5e9), NJOBS (300),
SEED (1), BENCHDIR (/root/scratch/poolbench)"""
import sys, os, math, time, random, subprocess, json, threading
import prove_cum, cbottom2, prof_steps, val_topgpu
from prove_iter import FastTop
from twoadic import tmax
from cbottom import consts
from bags import bags_of

W = os.path.dirname(os.path.abspath(__file__))
TOP, BOT = os.path.join(W, 'csrc', 'topgpu'), os.path.join(W, 'csrc', 'bottomgpu')
ENV = dict(os.environ, TP_CPU=os.path.join(W, 'csrc', 'topc3'), BG_CPU=os.path.join(W, 'csrc', 'bottomc3'))

def make_jobs(node, slice_units, njobs, seed, outdir):
    steps = []
    for b in bags_of(node):
        bag = tuple(sorted(b))
        for side, j, XL, top, bot, cost in prof_steps.plan_states(bag):
            units = cost * (prove_cum.BOTTOM_COST if side == 'bottom' else 1)
            if side == 'top' and j == 0: continue
            steps.append((units, side, j, XL, top, bot, bag))
    total = sum(s[0] for s in steps)
    rnd = random.Random(seed); pick = rnd.choices(steps, weights=[s[0] for s in steps], k=njobs)
    pick.sort(key=lambda s: (s[3], tmax(s[6])))          # one lattice per (XL, tm) in turn; the run order is shuffled below
    os.makedirs(outdir, exist_ok=True); cache = {}; jobs = []
    for i, (units, side, j, XL, top, bot, bag) in enumerate(pick):
        k = len(bag); tm = tmax(bag); G = math.ceil(2 * math.log10(XL)) + 12
        if side == 'top':
            orders = prove_cum.ordered(bag, j)
            per = units / len(orders); no = max(1, min(len(orders), round(slice_units / per))); o0 = rnd.randrange(len(orders) - no + 1)
            key = (XL, tm, G)
            if key not in cache: cache.clear(); cache[key] = FastTop(XL, tm, G)
            T = cache[key]; ranges = prove_cum.deff_top(j, top)
            text, d, prec = val_topgpu.job_text(T, j, ranges, orders[o0:o0 + no], 1, ranges[0] - 1)
            text += val_topgpu.fp_section(T); exe = TOP
        else:
            Gb = max(G, tm + 8); U = int(2.0960 * XL + 3); Wb = int(1.1834 * XL + 2)
            heff = cbottom2.heff_of(j, k, bot); orders = prove_cum.ordered(bag, k - j)
            per = units / len(orders); no = max(1, min(len(orders), round(slice_units / per))); o0 = rnd.randrange(len(orders) - no + 1)
            Mmax = max(4 * Gb + 8, heff[0], Gb); c = consts(Mmax)
            lines = ["cum", f"G {Gb}", f"tm {tm}", f"U {U}", f"W {Wb}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
            lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {no}"]
            lines += [" ".join(str(9 * (dd - 1)) for dd in o) for o in orders[o0:o0 + no]] + ["skip 0", "aklo 0", f"akhi {heff[-1] - 1}"]
            text = "\n".join(lines) + "\n"; exe = BOT; G = Gb
        path = os.path.join(outdir, f"job{i:04d}_{side}_G{G}.job"); open(path, 'w').write(text)
        jobs.append({'path': path, 'exe': exe, 'side': side, 'G': G, 'est': per * no})
    rnd.shuffle(jobs)
    json.dump({'node': node, 'total_units': total, 'jobs': jobs}, open(os.path.join(outdir, 'jobs.json'), 'w'))
    return total, jobs

def run_window(jobs, n, window):
    lock = threading.Lock(); nxt = [0]; recs = []; t0 = time.time()
    def worker():
        while True:
            with lock:
                if nxt[0] >= len(jobs) or time.time() - t0 > window: return
                jb = jobs[nxt[0]]; nxt[0] += 1
            s = time.time() - t0
            r = subprocess.run([jb['exe'], jb['path']], capture_output=True, text=True, env=ENV)
            e = time.time() - t0
            done = [l.split() for l in r.stdout.splitlines() if l.startswith('DONE')]
            if r.returncode or not done: raise RuntimeError(f"{jb['path']}: rc {r.returncode} {r.stderr[-200:]}")
            u = int(done[-1][1]) * (prove_cum.BOTTOM_COST if jb['side'] == 'bottom' else 1)
            with lock: recs.append((s, e, u, jb['side'], jb['G']))
    ths = [threading.Thread(target=worker) for _ in range(n)]
    for t in ths: t.start()
    for t in ths: t.join()
    a, b = 60.0, float(window)
    units = sum(u * max(0.0, min(e, b) - max(s, a)) / max(e - s, 1e-9) for s, e, u, _, _ in recs)
    return units / (b - a), recs

if __name__ == '__main__':
    node, window = int(sys.argv[1]), float(sys.argv[2]); ns = [int(x) for x in sys.argv[3:]]
    out = os.environ.get('BENCHDIR', '/root/scratch/poolbench')
    total, jobs = make_jobs(node, float(os.environ.get('SLICE', '1.5e9')), int(os.environ.get('NJOBS', '300')), int(os.environ.get('SEED', '1')), out)
    print(f"node {node}: {total:.3e} units in its plans; {len(jobs)} sampled slices, mean {sum(j['est'] for j in jobs) / len(jobs):.2e} units,"
          f" {sum(j['side'] == 'top' for j in jobs)} top / {sum(j['side'] == 'bottom' for j in jobs)} bottom", flush=True)
    for n in ns:
        rate, recs = run_window(jobs, n, window)
        busy = sum(e - s for s, e, *_ in recs)
        print(f"N = {n:2d}: {len(recs)} slices in the window, steady rate {rate:.3e} units/s -> whole node {total / rate / 3600:.1f} h"
              f" of shared GPU time (mean slice {busy / max(1, len(recs)):.1f} s)", flush=True)
