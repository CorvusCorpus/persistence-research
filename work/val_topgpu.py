"""Validation of the GPU top filter (csrc/topgpu, or topgpu_emu: the same kernel compiled for the CPU).
  K  accuracy: in dump mode (TP_DUMP=1) the kernel returns psi1 = frac(-Theta lambda), psi2 = frac(-Theta2 lambda) for every
     case of a job (the dumped cases must be exactly the job's cases); a sample, plus a sample of the cases whose first digit
     lies below height 5 (small-height table), is recomputed with arb at prec + 128 bits and every error must be <= 4 units of
     2^-64 (the bound the margins assume, docs/GPU_EVEN_TARGETS.md section 3). Jobs: top steps of bootstrapping plans.
  C  whole top steps (slices of the large ones: the middle and the first group of first gaps) of the plans of several bags:
     output byte-identical to csrc/topc3 (SURV runs, PROG, DONE cases f1 f2 nsurv).
  P  planted genuine numbers (work/Z.json, as planted_top_cum3.py): the top step of level j over the configuration space
     that contains a genuine number's first j digits must report it (a SURV run covering it), byte-identical to topc3.
usage (from work/): python3 val_topgpu.py EXE [K|C|P ...] [--sample N]"""
import sys, os, math, random, subprocess, tempfile, json
from fractions import Fraction
import flint

def run(exe, text, env=None):
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, text.encode()); os.close(fd)
    try:
        r = subprocess.run([exe, path], capture_output=True, text=True, env=env)
        if r.returncode: raise RuntimeError(f"{exe} failed ({r.returncode}): {r.stderr[-400:]}")
        return r.stdout
    finally:
        os.unlink(path)

def job_text(T, j, ranges, orders, g0lo, g0hi):
    """the job file ctop3.run_top writes"""
    import ctop2
    d = ctop2.lattice_job(T); prec = max(T.prec, d['C'].bit_length() + 64)
    lines = [f"prec {prec}"] + [f"{k} {d[k]}" for k in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd')]
    lines += [f"j {j}", "cum 1", "ranges " + " ".join(map(str, ranges)), f"nord {len(orders)}"]
    lines += [" ".join(str(9 * (x - 1)) for x in o) for o in orders]
    lines += ["skip 0", "maxsurv 10000000", "debug 0", f"g0lo {g0lo}", f"g0hi {g0hi}"]
    return "\n".join(lines) + "\n", d, prec

def plan_jobs(bag, maxcases=400_000, groups=('mid', 'first')):
    """the top steps of the bag's plan: whole if at most maxcases cases, else slices (the middle ordering, and the middle and/or
    the first group of first gaps: the first group holds the small heights, which take the kernel's small-height table)"""
    import prove_cum, ctop3
    from prove_iter import FastTop
    from twoadic import tmax
    import prof_steps
    out = []; tm = tmax(bag); cache = {}
    for side, j, XL, top, bot, cost in prof_steps.plan_states(bag):
        if side != 'top' or j == 0: continue
        G = math.ceil(2 * math.log10(XL)) + 12
        if (XL, G) not in cache: cache.clear(); cache[(XL, G)] = FastTop(XL, tm, G)
        T = cache[(XL, G)]
        ranges = prove_cum.deff_top(j, top); orders = prove_cum.ordered(bag, j)
        if cost <= maxcases: sl = [(orders, 1, ranges[0] - 1, 'whole')]
        else:
            grp = ctop3.groups(1, ranges[0] - 1, lambda v: ctop3.count_g0(ranges, v))
            sel = [orders[len(orders) // 2]]
            sl = [(sel, *grp[len(grp) // 2 if g == 'mid' else 0], f'slice-{g}') for g in groups]
        for sel, g0, g1, what in sl:
            text, d, prec = job_text(T, j, ranges, sel, g0, g1)
            out.append((f"{'_'.join(map(str, bag))} j={j} G={G} {what}", text, d, prec, sel, T))
    return out

def part_k(exe, sample, seed=1):
    rnd = random.Random(seed); ok = True
    for bag in [(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 3, 4, 6, 7, 8, 9)]:
        for name, text, d, prec, sel, _T in plan_jobs(bag, maxcases=60_000):
            outp = run(exe, text, env=dict(os.environ, TP_DUMP='1'))
            rows = [l.split() for l in outp.splitlines() if l.startswith('PSI')]
            if not rows: continue
            # the dumped cases must be exactly the job's cases (a wrapped field once made the kernel see other cases)
            import ctop4
            j_ = int(text.split("\nj ")[1].split()[0]); rng_ = list(map(int, text.split("\nranges ")[1].split("\n")[0].split()))
            g0lo_ = int(text.split("\ng0lo ")[1].split()[0]); g0hi_ = int(text.split("\ng0hi ")[1].split()[0])
            import top_filter_ref
            kv = {'j': j_, 'ranges': rng_, 'orders': sel, 'skip': 0, 'g0lo': g0lo_, 'g0hi': g0hi_}
            want = sorted((o, tuple(g)) for o, g in top_filter_ref.cases(kv))   # every case (small heights: the table)
            got = sorted((int(f[1]), tuple(map(int, f[2:2 + j_]))) for f in rows)
            if got != want:
                ok = False; print(f"K {name}: dumped cases differ from the job's cases ({len(got)} vs {len(want)}) -> WRONG CASES", flush=True)
                continue
            flint.ctx.prec = prec + 128
            C, qa, qb, a2, D = d['C'], d['qa'], d['qb'], d['a2'], d['D']
            Th = flint.arb(C) * qa / qb; Th2 = flint.arb(C) * a2 / D
            j = int(text.split("\nj ")[1].split()[0]); worst = 0; wpath = {}
            small = [f for f in rows if int(f[2]) < 5]   # first digit below height 5: lambda from the small-height table
            chk = rnd.sample(rows, min(sample, len(rows))) + rnd.sample(small, min(sample // 4, len(small)))
            bypath = {}                                  # psi1's path (0 table, 1 direct, 2 recurrence start, 3 recurrence step)
            for f in rows: bypath.setdefault(int(f[4 + j]) if len(f) > 4 + j else -1, []).append(f)
            for pth, fs in sorted(bypath.items()): chk += rnd.sample(fs, min(sample // 4, len(fs)))
            for f in chk:
                o = int(f[1]); gaps = list(map(int, f[2:2 + j])); p1, p2 = int(f[2 + j]), int(f[3 + j])
                digs = sel[o]; s = sum(gaps); N = 10**s; acc = 0
                for x, g in zip(digs, gaps): acc += g; N += 9 * (x - 1) * 10**(s - acc)
                lam = flint.arb(N).log() - s * flint.arb(10).log()
                pth = int(f[4 + j]) if len(f) > 4 + j else -1
                for which, (th, pg) in enumerate(((Th, p1), (Th2, p2))):
                    x = -th * lam; fr = x - x.mid().floor().unique_fmpz()
                    t = int((fr * 2**64).mid().floor().unique_fmpz()) % 2**64
                    e = min((pg - t) % 2**64, (t - pg) % 2**64)
                    worst = max(worst, e)
                    if which == 0: wpath[pth] = max(wpath.get(pth, 0), e)
            good = worst <= 4; ok &= good
            paths = " ".join(f"{p}:{len(bypath[p])}/{wpath.get(p, '-')}" for p in sorted(bypath))
            print(f"K {name}: {len(rows)} cases dumped ({len(small)} with a digit below height 5), {len(chk)} checked, worst error {worst} units of 2^-64 -> {'OK' if good else 'TOO LARGE'}  [psi1 path:cases/worst {paths}]", flush=True)
    return ok

def part_c(exe, cpu):
    ok = True; n = 0
    for bag in [(2, 2, 3), (2, 3, 7, 7), (2, 2, 2, 3, 3, 7), (2, 3, 4, 6, 7, 8, 9), (2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3)]:
        for name, text, d, prec, sel, _T in plan_jobs(bag):
            a = run(cpu, text); b = run(exe, text, env=dict(os.environ, TP_CPU=cpu)); n += 1
            same = a == b; ok &= same
            print(f"C {name}: {a.splitlines()[-1]} -> {'IDENTICAL' if same else 'DIFFERENT'}", flush=True)
    print(f"C: {n} jobs")
    return ok

def part_p(exe, cpu, nmax=150, fp=False):
    """planted_top_cum3.py's configurations: genuine y with digits 1..j inside the top D_{j-1} < G positions and digit j+1 at
    depth D_j >= G; the top job of level j (the genuine ordering, bounds D_i < D_i(y) + 1) must keep it. fp: with the fp
    section (the genuine run must still be printed; the output is checked with fp_compare instead of byte identity)"""
    from twoadic import v2
    from prove_iter import FastTop
    Z = json.load(open('Z.json'))
    XLs = 80; G = math.ceil(2 * math.log10(XLs)) + 12
    tested = detected = 0; ok = True; tops = {}
    for y in Z:
        s = str(y)
        if y % 5 == 0 or '0' in s or '5' in s or len(s) > XLs: continue
        L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch) for ch in s if ch != '1']
        k = len(pos)
        for j in range(1, k):
            Dprev = L - pos[j - 1]; Dj = L - pos[j]
            if Dprev < G and Dj >= G:
                tm = v2(y) + 2
                if tm not in tops: tops[tm] = FastTop(XLs, tm, G)
                T = tops[tm]
                bounds = [(L - pos[i]) + 1 for i in range(j)]
                text, d, prec = job_text(T, j, bounds, [tuple(digs[:j])], 1, bounds[0] - 1)
                if fp: text += fp_section(T)
                a = run(cpu, text); b = run(exe, text, env=dict(os.environ, TP_CPU=cpu))
                gaps = [L - pos[0]] + [pos[i - 1] - pos[i] for i in range(1, j)]
                hit = False
                for l in b.splitlines():
                    if not l.startswith('SURV'): continue
                    f = l.split(); pre = list(map(int, f[2:2 + j - 1])); lo, hi = int(f[1 + j]), int(f[2 + j])
                    if pre == gaps[:j - 1] and lo <= gaps[j - 1] <= hi: hit = True
                same = fp_compare(a, b, T.lat)[0] if fp else a == b
                ok &= same and hit; tested += 1; detected += hit
                if not (same and hit): print(f"P y={y} j={j}: {'IDENTICAL' if same else 'DIFFERENT'}, genuine {'kept' if hit else 'MISSED'}")
                break
        if tested >= nmax: break
    print(f"P{' (fp)' if fp else ''} G={G}: {detected}/{tested} planted genuine configurations kept, outputs {'consistent with' if fp else 'identical to'} topc3: {ok}")
    return ok

def fp_section(T):
    """the job section with which topgpu resolves survivor runs whose exact lattice enumeration is empty (as ctop4.run_top)"""
    lat = T.lat; R2 = Fraction(lat.R2)
    return "\n".join([f"fp {lat.dim}"] + [" ".join(map(str, r)) for r in lat.B] + [f"R2 {R2.numerator} {R2.denominator}"]) + "\n"

def py_nonempty(lat, m):
    """Python's exact enumeration (solver.close_points, the base of ctop2.candidates) around (0, .., 0, -m): non-empty?"""
    from solver import close_points, _attach
    _attach(lat)
    return bool(close_points(lat.B, lat.gs2, lat._mu, lat._Binv, [0] * (lat.dim - 1) + [-m], lat.R2))

def fp_compare(a, b, lat):
    """a: topc3's output, b: topgpu's on the same job with the fp section. DONE and the first PROG fields identical; the printed
    SURV runs are exactly topc3's runs with a non-empty Python enumeration, in order; the resolved counts (fifth PROG field)
    are, per ordering, the runs with an empty one. Returns (ok, runs, resolved)."""
    la, lb = a.splitlines(), b.splitlines()
    pa = [l.split() for l in la if l.startswith('PROG')]; pb = [l.split() for l in lb if l.startswith('PROG')]
    if pb and all(len(p) == 5 for p in pb):                  # handed to topc3 (exec): no fp mode, the output must be topc3's
        return a == b, sum(l.startswith('SURV') for l in la), 0
    ok = [l for l in la if l.startswith('DONE')] == [l for l in lb if l.startswith('DONE')]
    ok &= [p[:5] for p in pa] == [p[:5] for p in pb]
    runs = [l for l in la if l.startswith('SURV')]; keep = []; res = {}
    for l in runs:
        f = l.split()
        if py_nonempty(lat, int(f[-1])): keep.append(l)
        else: res[int(f[1])] = res.get(int(f[1]), 0) + 1
    ok &= keep == [l for l in lb if l.startswith('SURV')]
    ok &= all(len(p) > 5 and int(p[5]) == res.get(int(p[1]) - 1, 0) for p in pb)
    return ok, len(runs), len(runs) - len(keep)

def part_f(exe, cpu, nsyn=200, seed=7):
    """fp mode (survivor runs resolved by topgpu's exact lattice enumeration):
    F1 the C enumeration against Python's close_points on m values: random ones, and ones planted next to genuine lattice
       vectors (t lt, u', w lw, t rC2 + u' rC3 + w rC7 - L rC10) of the box, at distances inside, on and just beyond R;
    F2 whole steps and slices (as part C): fp_compare against topc3 on every survivor run;
    F3 planted genuine configurations (part P) with the fp section: the genuine run must still be printed."""
    from prove_iter import FastTop
    from solver import _attach
    rnd = random.Random(seed); ok = True; tot = agree = nonempty = 0
    lats = []
    for bag in [(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (3, 3, 7, 7, 9)]:          # dim 4, and dim 3 (no even digit: t = 0)
        for name, text, d, prec, sel, T in plan_jobs(bag, maxcases=0, groups=('mid',)):
            if all(T.lat is not x[1].lat for x in lats): lats.append((name, T, text))
    for name, T, text in lats:
        lat = T.lat; _attach(lat); n = lat.dim; R2 = Fraction(lat.R2); ms = []
        for _ in range(nsyn // 4): ms.append(rnd.randrange(-3 * T.C, 3 * T.C))
        while len(ms) < nsyn:
            t = rnd.randrange(0, lat.T + 1) if n == 4 else 0
            up, w, L = rnd.randrange(0, lat._U + 1), rnd.randrange(0, lat._W + 1), rnd.randrange(1, lat.XL + 1)
            last = t * lat.rC[2] + up * lat.rC[3] + w * lat.rC[7] - L * lat.rC[10]
            rem = R2 - (t * lat._lt) ** 2 - up ** 2 - (w * lat._lw) ** 2
            if rem < 0: continue
            s = math.isqrt(rem.numerator // rem.denominator)
            for dl in (rnd.randrange(-s, s + 1), s, -s, s + 1, -s - 1): ms.append(-last + dl)
        job = text + fp_section(T) + f"fptest {len(ms)}\n" + "\n".join(map(str, ms)) + "\n"
        out = run(exe, job, env=dict(os.environ, TP_CPU=cpu))
        got = {int(l.split()[1]): int(l.split()[2]) for l in out.splitlines() if l.startswith('FP')}
        bad = 0
        for i, m in enumerate(ms):
            want = py_nonempty(lat, m); tot += 1; nonempty += want
            if got.get(i) != int(want): bad += 1
            else: agree += 1
        ok &= bad == 0
        print(f"F1 {name} dim={n}: {len(ms)} values, C and Python agree on {len(ms) - bad}"
              f" (non-empty {sum(py_nonempty(lat, m) for m in ms)})", flush=True)
    print(f"F1: {agree}/{tot} agree ({nonempty} non-empty)")
    nj = nr = nres = 0
    for bag in [(2, 2, 3), (2, 3, 7, 7), (2, 3, 4, 6, 7, 8, 9), (2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3)]:
        for name, text, d, prec, sel, T in plan_jobs(bag):
            job = text + fp_section(T)
            a = run(cpu, job); b = run(exe, job, env=dict(os.environ, TP_CPU=cpu))
            good, runs, res = fp_compare(a, b, T.lat); ok &= good; nj += 1; nr += runs; nres += res
            print(f"F2 {name}: {runs} survivor runs, {res} resolved by topgpu -> {'CONSISTENT' if good else 'INCONSISTENT'}", flush=True)
    print(f"F2: {nj} jobs, {nr} survivor runs, {nres} resolved")
    ok &= part_p(exe, cpu, fp=True)
    return ok

if __name__ == '__main__':
    args = sys.argv[1:]; exe = args[0]; parts = [a for a in args[1:] if a in ('K', 'C', 'P', 'F')] or ['K', 'C', 'P', 'F']
    sample = int(args[args.index('--sample') + 1]) if '--sample' in args else 400
    cpu = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'topc3')
    ok = True
    if 'K' in parts: ok &= part_k(exe, sample)
    if 'C' in parts: ok &= part_c(exe, cpu)
    if 'P' in parts: ok &= part_p(exe, cpu)
    if 'F' in parts: ok &= part_f(exe, cpu)
    print("TOPGPU VALIDATION PASSED" if ok else "TOPGPU VALIDATION FAILED")
    sys.exit(0 if ok else 1)
