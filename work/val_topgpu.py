"""Validation of the GPU top filter (csrc/topgpu, or topgpu_emu: the same kernel compiled for the CPU).
  K  accuracy: in dump mode (TP_DUMP=1) the kernel returns psi1 = frac(-Theta lambda), psi2 = frac(-Theta2 lambda) for every
     case of a job; a sample is recomputed with arb at prec + 128 bits and every error must be <= 4 units of 2^-64 (the bound
     the margins assume, docs/GPU_EVEN_TARGETS.md section 3). Jobs: top steps of bootstrapping plans, G = 51 .. 383.
  C  whole top steps (slices of the large ones) of the plans of several bags: output byte-identical to csrc/topc3 (SURV
     runs, PROG, DONE cases f1 f2 nsurv).
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

def plan_jobs(bag, maxcases=400_000):
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
        if cost <= maxcases: sel, g0, g1, what = orders, 1, ranges[0] - 1, 'whole'
        else:
            grp = ctop3.groups(1, ranges[0] - 1, lambda v: ctop3.count_g0(ranges, v))
            sel = [orders[len(orders) // 2]]; (g0, g1) = grp[len(grp) // 2]; what = 'slice'
        text, d, prec = job_text(T, j, ranges, sel, g0, g1)
        out.append((f"{'_'.join(map(str, bag))} j={j} G={G} {what}", text, d, prec, sel))
    return out

def part_k(exe, sample, seed=1):
    rnd = random.Random(seed); ok = True
    for bag in [(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 3, 4, 6, 7, 8, 9)]:
        for name, text, d, prec, sel in plan_jobs(bag, maxcases=60_000):
            outp = run(exe, text, env=dict(os.environ, TP_DUMP='1'))
            rows = [l.split() for l in outp.splitlines() if l.startswith('PSI')]
            if not rows: continue
            # the dumped cases must be exactly the job's cases (a wrapped field once made the kernel see other cases)
            import ctop4
            j_ = int(text.split("\nj ")[1].split()[0]); rng_ = list(map(int, text.split("\nranges ")[1].split("\n")[0].split()))
            g0lo_ = int(text.split("\ng0lo ")[1].split()[0]); g0hi_ = int(text.split("\ng0hi ")[1].split()[0])
            import top_filter_ref
            kv = {'j': j_, 'ranges': rng_, 'orders': sel, 'skip': 0, 'g0lo': g0lo_, 'g0hi': g0hi_}
            def kernel_case(g):                          # the kernel hands heights below Smin = 5 to the host (no dump)
                if j_ >= 3 and g[0] + g[1] < 5: return False
                return sum(g) >= 5
            want = sorted((o, tuple(g)) for o, g in top_filter_ref.cases(kv) if kernel_case(g))
            got = sorted((int(f[1]), tuple(map(int, f[2:2 + j_]))) for f in rows)
            if got != want:
                ok = False; print(f"K {name}: dumped cases differ from the job's cases ({len(got)} vs {len(want)}) -> WRONG CASES", flush=True)
                continue
            flint.ctx.prec = prec + 128
            C, qa, qb, a2, D = d['C'], d['qa'], d['qb'], d['a2'], d['D']
            Th = flint.arb(C) * qa / qb; Th2 = flint.arb(C) * a2 / D
            j = int(text.split("\nj ")[1].split()[0]); worst = 0
            for f in rnd.sample(rows, min(sample, len(rows))):
                o = int(f[1]); gaps = list(map(int, f[2:2 + j])); p1, p2 = int(f[2 + j]), int(f[3 + j])
                digs = sel[o]; s = sum(gaps); N = 10**s; acc = 0
                for x, g in zip(digs, gaps): acc += g; N += 9 * (x - 1) * 10**(s - acc)
                lam = flint.arb(N).log() - s * flint.arb(10).log()
                for th, pg in ((Th, p1), (Th2, p2)):
                    x = -th * lam; fr = x - x.mid().floor().unique_fmpz()
                    t = int((fr * 2**64).mid().floor().unique_fmpz()) % 2**64
                    e = min((pg - t) % 2**64, (t - pg) % 2**64)
                    worst = max(worst, e)
            good = worst <= 4; ok &= good
            print(f"K {name}: {len(rows)} cases dumped, {min(sample, len(rows))} checked, worst error {worst} units of 2^-64 -> {'OK' if good else 'TOO LARGE'}", flush=True)
    return ok

def part_c(exe, cpu):
    ok = True; n = 0
    for bag in [(2, 2, 3), (2, 3, 7, 7), (2, 2, 2, 3, 3, 7), (2, 3, 4, 6, 7, 8, 9), (2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3)]:
        for name, text, d, prec, sel in plan_jobs(bag):
            a = run(cpu, text); b = run(exe, text, env=dict(os.environ, TP_CPU=cpu)); n += 1
            same = a == b; ok &= same
            print(f"C {name}: {a.splitlines()[-1]} -> {'IDENTICAL' if same else 'DIFFERENT'}", flush=True)
    print(f"C: {n} jobs")
    return ok

def part_p(exe, cpu, nmax=150):
    """planted_top_cum3.py's configurations: genuine y with digits 1..j inside the top D_{j-1} < G positions and digit j+1 at
    depth D_j >= G; the top job of level j (the genuine ordering, bounds D_i < D_i(y) + 1) must keep it"""
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
                a = run(cpu, text); b = run(exe, text, env=dict(os.environ, TP_CPU=cpu))
                gaps = [L - pos[0]] + [pos[i - 1] - pos[i] for i in range(1, j)]
                hit = False
                for l in b.splitlines():
                    if not l.startswith('SURV'): continue
                    f = l.split(); pre = list(map(int, f[2:2 + j - 1])); lo, hi = int(f[1 + j]), int(f[2 + j])
                    if pre == gaps[:j - 1] and lo <= gaps[j - 1] <= hi: hit = True
                same = a == b; ok &= same and hit; tested += 1; detected += hit
                if not (same and hit): print(f"P y={y} j={j}: {'IDENTICAL' if same else 'DIFFERENT'}, genuine {'kept' if hit else 'MISSED'}")
                break
        if tested >= nmax: break
    print(f"P G={G}: {detected}/{tested} planted genuine configurations kept, outputs identical to topc3: {ok}")
    return ok

if __name__ == '__main__':
    args = sys.argv[1:]; exe = args[0]; parts = [a for a in args[1:] if a in ('K', 'C', 'P')] or ['K', 'C', 'P']
    sample = int(args[args.index('--sample') + 1]) if '--sample' in args else 400
    cpu = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'topc3')
    ok = True
    if 'K' in parts: ok &= part_k(exe, sample)
    if 'C' in parts: ok &= part_c(exe, cpu)
    if 'P' in parts: ok &= part_p(exe, cpu)
    print("TOPGPU VALIDATION PASSED" if ok else "TOPGPU VALIDATION FAILED")
    sys.exit(0 if ok else 1)
