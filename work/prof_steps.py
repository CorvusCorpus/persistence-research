"""Time the C case loops (csrc/topc3, csrc/bottomc3) on slices of the steps that prove_bag would run for one bag.

For each selected step of the dry-run plan (bag_plan.py) one slice is run exactly as ctop3/cbottom4 would run it: one
ordering and one group of first gaps (top, about 200,000 cases) or lowest heights (bottom, about 50,000 solves).
Prints the wall time per case / per solve and an estimate of the whole step.
usage (from work/): python3 prof_steps.py BAG [NSTEPS=6] [--keep DIR]   (the NSTEPS most costly steps of the plan)
With --keep DIR the job files are written to DIR (for perf), named step<i>_<side>.job."""
import sys, os, time, math
import prove_cum, cbottom2, ctop3, cbottom4
from prove_iter import FastTop
from twoadic import tmax

def plan_states(bag):
    snaps = []; orig = prove_cum.matveev_cum
    def spy(k_, tm_, top, bot):
        XL = orig(k_, tm_, top, bot); snaps.append((dict(top), dict(bot), XL)); return XL
    prove_cum.matveev_cum = spy
    try: prove_cum.prove_bag(bag, log=lambda x: None, dry=True)
    finally: prove_cum.matveev_cum = orig
    k = len(bag); out = []
    for (t0, b0, XL), (t1, b1, _) in zip(snaps, snaps[1:]):
        dt = [j for j in t1 if t0.get(j) != t1[j]]; db = [j for j in b1 if b0.get(j) != b1[j]]
        if dt:
            j = dt[0]; cost = (prove_cum.ncomp(prove_cum.deff_top(j, t0)) if j else 1) * max(1, len(prove_cum.ordered(bag, j)))
            out.append(('top', j, XL, t0, b0, cost))
        else:
            j = db[0]; cost = prove_cum.nbottom(cbottom2.heff_of(j, k, b0)) * len(prove_cum.ordered(bag, k - j))
            out.append(('bottom', j, XL, t0, b0, cost))
    return out

def slice_top(bag, j, XL, top, keep=None, name=''):
    tm = tmax(bag); G = math.ceil(2 * math.log10(XL)) + 12
    T = FastTop(XL, tm, G)
    ranges = prove_cum.deff_top(j, top); orders = prove_cum.ordered(bag, j); oo = len(orders) // 2
    grp = ctop3.groups(1, ranges[0] - 1, lambda v: ctop3.count_g0(ranges, v))
    g0, g1 = grp[len(grp) // 2]
    if keep:   # same job text as ctop3.run_top
        d = ctop3.lattice_job(T); prec = max(T.prec, d['C'].bit_length() + 64)
        lines = [f"prec {prec}"] + [f"{kk} {d[kk]}" for kk in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
        lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, ranges)), "nord 1", " ".join(str(9 * (x - 1)) for x in orders[oo])]
        lines += ["skip 0", "maxsurv 10000000", "debug 0", f"g0lo {g0}", f"g0hi {g1}"]
        open(os.path.join(keep, name + '.job'), 'w').write("\n".join(lines) + "\n")
    t0 = time.time(); done = None; nsurv = 0
    for f in ctop3.run_top(T, j, ranges, [orders[oo]], g0, g1):
        if f[0] == 'DONE': done = list(map(int, f[1:]))
        elif f[0] == 'SURV': nsurv += 1
    dt = time.time() - t0
    return dict(G=G, Cbits=T.C.bit_length(), prec=max(T.prec, T.C.bit_length() + 64), cases=done[0], f1=done[1], f2=done[2],
                surv=nsurv, wall=dt, per=dt / max(1, done[0]), slice=f"order {oo}/{len(orders)} g0 [{g0},{g1}]")

def slice_bottom(bag, j, XL, bot, keep=None, name=''):
    k = len(bag); tm = tmax(bag); G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
    U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
    heff = cbottom2.heff_of(j, k, bot); orders = prove_cum.ordered(bag, k - j); oo = len(orders) // 2; nb = k - j
    grp = ctop3.groups(0, heff[nb - 1] - 1, lambda v: cbottom4.count_ak(heff, v), 50000)
    a0, a1 = grp[len(grp) // 2]
    if keep:   # same job text as cbottom4.run
        from cbottom import consts
        Mmax = max(4 * Gb + 8, heff[0], Gb); c = consts(Mmax)
        lines = ["cum", f"G {Gb}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{kk} {c[kk]}" for kk in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]
        lines += [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), "nord 1", " ".join(str(9 * (d - 1)) for d in orders[oo]), "skip 0", f"aklo {a0}", f"akhi {a1}"]
        open(os.path.join(keep, name + '.job'), 'w').write("\n".join(lines) + "\n")
    t0 = time.time(); done = None; ncand = 0
    for f in cbottom4.run(Gb, tm, U, W, j, k, heff, [orders[oo]], a0, a1):
        if f[0] == 'DONE': done = list(map(int, f[1:]))
        elif f[0] == 'CAND': ncand += 1
    dt = time.time() - t0
    return dict(G=Gb, M5bits=round(Gb * math.log2(5)), solves=done[0], direct5=done[1], direct2=done[2], cand=ncand, wall=dt,
                per=dt / max(1, done[0]), slice=f"order {oo}/{len(orders)} ak [{a0},{a1}]")

if __name__ == '__main__':
    args = sys.argv[1:]; keep = None
    if '--keep' in args: i = args.index('--keep'); keep = args[i + 1]; del args[i:i + 2]; os.makedirs(keep, exist_ok=True)
    bag = tuple(sorted(int(x) for x in args[0].split('_'))); nsteps = int(args[1]) if len(args) > 1 else 6
    st = plan_states(bag)
    idx = sorted(range(len(st)), key=lambda i: -st[i][5])[:nsteps]
    tot_est = {'top': 0.0, 'bottom': 0.0}
    for i in sorted(idx):
        side, j, XL, top, bot, cost = st[i]
        if side == 'top' and j == 0: continue
        if side == 'top':
            r = slice_top(bag, j, XL, top, keep, f"step{i}_top")
            est = cost * r['per']
            print(f"step {i:2d} top    j={j} G={r['G']} Cbits={r['Cbits']} prec={r['prec']}: {r['cases']} cases in {r['wall']:.2f}s = "
                  f"{1e6 * r['per']:.2f} us/case (f1 {r['f1']}, f2 {r['f2']}, runs {r['surv']}; {r['slice']}); step {cost:.3e} cases -> {est:.0f} s", flush=True)
        else:
            r = slice_bottom(bag, j, XL, bot, keep, f"step{i}_bottom")
            est = cost * r['per']
            print(f"step {i:2d} bottom j={j} G={r['G']} (5^G: {r['M5bits']} bits): {r['solves']} solves in {r['wall']:.2f}s = "
                  f"{1e6 * r['per']:.2f} us/solve (direct5 {r['direct5']}, direct2 {r['direct2']}, cand {r['cand']}; {r['slice']}); "
                  f"step {cost:.3e} solves -> {est:.0f} s", flush=True)
        tot_est[side] += est
    print(f"estimated wall time of these steps on one core: top {tot_est['top']:.0f} s, bottom {tot_est['bottom']:.0f} s")
