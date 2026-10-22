"""Validation of the GPU bottom solver (csrc/bottomgpu, or bottomgpu_emu: the same kernel compiled for the CPU).
  K  kernel tests: queries (Flow, c, p) with arbitrary Flow at production scale (G = 14 .. 360, boxes U, W ~ 2 XL with
     G = ceil(2 log10 XL) + 12). Planted: F = 2^t 3^u' 7^w mod 10^G with (t, u', w) in the box, split as c 10^p + Flow
     (all p, including p = 0 and small p, which take the direct 5-adic and 2-adic paths); random: arbitrary Flow.
     For every query the kernel's candidates must equal bottom_ref.py's exact list, and every planted point must be found.
  C  configuration spaces: whole steps of the bootstrapping plan of small bags (and slices of large ones), and the steps with
     G > 200 of the largest bag of node 221184 (size classes up to G = 360), byte-identical output (CAND lines, PROG, DONE
     counters) to csrc/bottomc3.
  P  genuine numbers: jobs whose configuration space contains the digits of a genuine zero-free 7-smooth number (work/Z.json);
     byte-identical to bottomc3, and the genuine solution must be among the CAND lines.
usage (from work/): python3 val_bottomgpu.py EXE [K|C|P ...] [--n NQ] [--seed S]"""
import sys, os, math, random, subprocess, tempfile, json
import bottom_ref
from cbottom import consts

DIGC = [9, 18, 27, 45, 54, 63, 72]

def header(mode, G, tm, U, W, Mmax):
    c = consts(Mmax)
    return [mode, f"G {G}", f"tm {tm}", f"U {U}", f"W {W}", f"Mmax {Mmax}"] + [f"{k} {c[k]}" for k in ('L5_3', 'L5_7', 'L5_16', 'L2_3', 'L2_7')]

def run(exe, lines, env=None):
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        r = subprocess.run([exe, path], capture_output=True, text=True, env=env)
        if r.returncode: raise RuntimeError(f"{exe} failed ({r.returncode}): {r.stderr[-500:]}")
        return r.stdout
    finally:
        os.unlink(path)

def part_k(exe, nq, seed):
    rnd = random.Random(seed); ok = True
    for G in (14, 25, 40, 51, 60, 87, 90, 123, 125, 160, 165, 197, 200, 234, 240, 271, 280, 308, 320, 346, 360):
        XL = 10 ** ((G - 12) // 2) + rnd.randrange(10 ** max(1, (G - 12) // 2 - 1))
        U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2); tm = min(38, G - 8); Mmax = max(4 * G + 8, G)
        cs = consts(Mmax)
        ctx = bottom_ref.Ctx(G, tm, U, W, Mmax, (cs['L5_3'], cs['L5_7'], cs['L5_16'], cs['L2_3'], cs['L2_7']))
        qs, planted = [], {}
        while len(qs) < nq:
            kind = rnd.random()
            c = rnd.choice(DIGC)
            p = rnd.choice([0, 1, 2, 3]) if rnd.random() < 0.3 else rnd.randrange(0, G)
            if kind < 0.7:                                    # planted
                t = rnd.randrange(0, min(tm, 12) + 1); up = rnd.randrange(2, U + 1); w = rnd.randrange(0, W + 1)
                F = pow(2, t, 10**G) * pow(3, up, 10**G) * pow(7, w, 10**G) % 10**G
                Flow = F - c * 10**p
                if Flow % 5 == 0: continue                    # (p = 0 only) never in real jobs: Flow = c - 1 mod 5, c != 1 mod 5
                try:
                    ctx.mt(G, t)                              # skip t where bottom3.c itself would stop (lattice too fine)
                except AssertionError:
                    continue
                planted[len(qs)] = (t, up, w)
            else:                                             # arbitrary Flow (never divisible by 5, as in real jobs)
                Flow = rnd.randrange(-10**G, 10**G)
                if Flow % 5 == 0: continue
            qs.append((Flow, c, p))
        out = run(exe, header("ktest", G, tm, U, W, Mmax) + [f"n {len(qs)}"] + [f"{a} {b} {c}" for a, b, c in qs])
        got = {}; deferred = 0
        for line in out.splitlines():
            f = line.split()
            if f[0] == 'Q': got.setdefault(int(f[1]), []).append((int(f[3]), int(f[5]), int(f[7])))
            elif f[0] == 'H': deferred += 1
            elif f[0] == 'CAND':                              # bottom3.c's decision for a deferred query (rep_o = -1 - index)
                i = -1 - int(f[1]); k_ = f.index('t'); got.setdefault(i, []).append((int(f[k_ + 1]), int(f[k_ + 3]), int(f[k_ + 5])))
        bad = 0; found = 0
        for i, (Flow, c, p) in enumerate(qs):
            lw = None
            try:
                lw = ctx.low(Flow)
            except (ValueError, ZeroDivisionError):
                lw = None                                     # Flow divisible by 5: no solve can pass
            exp = []
            if lw is not None:
                try:
                    res = ctx.beta(lw, c, p, G)
                    if res is not None: exp = ctx.candidates(*res)
                except AssertionError:
                    exp = None                                # bottom3.c would stop here (lattice too fine)
            if exp is None: continue
            g = sorted(got.get(i, []))
            if g != sorted(exp): bad += 1; ok = False; print(f"  MISMATCH G={G} q{i} Flow={Flow} c={c} p={p}: kernel {g} ref {sorted(exp)}")
            if i in planted:
                if planted[i] in g: found += 1
                else: bad += 1; ok = False; print(f"  PLANTED NOT FOUND G={G} q{i} {planted[i]}")
        print(f"K G={G:3d} XL=1e{math.log10(XL):.0f}: {len(qs)} queries, planted {found}/{len(planted)} found, mismatches {bad}, deferred {deferred}  [{out.splitlines()[-1]}]", flush=True)
    return ok

def plan_jobs(bag, maxsolves=3_000_000):
    """bottom steps of the dry-run plan of `bag`, as job texts (whole steps when small, else one ordering / ak group)"""
    import prove_cum, cbottom2, ctop3, cbottom4
    from twoadic import tmax
    import prof_steps
    out = []
    for side, j, XL, top, bot, cost in prof_steps.plan_states(bag):
        if side != 'bottom': continue
        k = len(bag); tm = tmax(bag); G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
        U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
        heff = cbottom2.heff_of(j, k, bot); orders = prove_cum.ordered(bag, k - j); nb = k - j
        Mmax = max(4 * Gb + 8, heff[0], Gb)
        if cost <= maxsolves: sel, a0, a1 = orders, 0, heff[nb - 1] - 1
        else:
            grp = ctop3.groups(0, heff[nb - 1] - 1, lambda v: cbottom4.count_ak(heff, v), 50000)
            sel = [orders[len(orders) // 2]]; a0, a1 = grp[len(grp) // 2]
        lines = header("cum", Gb, tm, U, W, Mmax) + [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(sel)}"]
        lines += [" ".join(str(9 * (d - 1)) for d in o) for o in sel] + ["skip 0", f"aklo {a0}", f"akhi {a1}"]
        out.append((f"{'_'.join(map(str, bag))} j={j} G={Gb} {'whole' if cost <= maxsolves else 'slice'}", lines))
    return out

def part_c(exe, cpu):
    ok = True; n = 0
    bags = [((2, 2, 3), 0), ((2, 3, 7, 7), 0), ((2, 2, 2, 3, 3, 7), 0), ((2, 3, 4, 6, 7, 8, 9), 0), ((2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), 0),
            ((2,) * 13 + (3, 3, 3), 200)]                 # the largest bag of node 221184: only its steps with G > 200
    for bag, gmin in bags:
        for name, lines in plan_jobs(bag):
            G = int(lines[1].split()[1]); h0 = int([l for l in lines if l.startswith('heff')][0].split()[1])
            if G <= gmin: continue
            a = run(cpu, lines); b = run(exe, lines); n += 1
            same = a == b
            ok &= same
            path = 'kernel' if (h0 <= G and G <= 360) else 'bottomc3 (exec)'
            print(f"C {name} [{path}]: {a.splitlines()[-1]} -> {'IDENTICAL' if same else 'DIFFERENT'}", flush=True)
            if not same:
                print("   cpu:", a[-300:].replace("\n", " | ")); print("   gpu:", b[-300:].replace("\n", " | "))
    print(f"C: {n} jobs")
    return ok

def part_p(exe, cpu, nmax=60):
    """the lowest b digits of a genuine y below G and its next non-1 digit at height >= G: then 9y = F mod 10^G for the
    configuration of those b digits, so (t, u', w) of 9y is a solution in any box containing it. Zero-free numbers have
    short runs of 1s, so G is taken per number in [max(a_(j+1) + 1, tm + 8), a_j] (small G; the kernel is the same)."""
    Z = json.load(open('Z.json')); ok = True; tested = 0
    from twoadic import v2
    for y in Z:
        s = str(y)
        if y % 5 == 0 or '0' in s or '5' in s or len(s) < 4: continue      # bags never hold a 5 (Lemma 0)
        L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; digs = [int(ch) for ch in s if ch != '1']
        k = len(pos); a = [L] + pos
        t = v2(9 * y); up = 0; z = 9 * y // 2**t
        while z % 3 == 0: z //= 3; up += 1
        w = 0
        while z % 7 == 0: z //= 7; w += 1
        tm = t                                            # any bound >= t; the smallest allows the smallest G (>= tm + 8)
        for j in range(k - 1, -1, -1):
            b = k - j
            if b > 12: break
            Glo, Ghi = max(a[j + 1] + 1, tm + 8), min(a[j], 60)
            G = next((g for g in range(Ghi, Glo - 1, -1) if math.comb(g, b) * g <= 6_000_000), None)
            if G is None: continue
            XL = max(up, w) + 5; U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
            heff = [G - r for r in range(b)]
            Mmax = max(4 * G + 8, G)
            lines = header("cum", G, tm, U, W, Mmax) + [f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), "nord 1",
                                                       " ".join(str(9 * (d - 1)) for d in digs[j:]), "skip 0", "aklo 0", f"akhi {heff[b - 1] - 1}"]
            try:
                A = run(cpu, lines)
            except RuntimeError:
                continue                                  # bottomc3 stops (e.g. lattice too fine for this tiny G): not a test case
            B = run(exe, lines)
            genuine = f" t {t} up {up} w {w}"
            hit = any(l.startswith('CAND') and l.endswith(genuine) for l in B.splitlines())
            same = A == B; ok &= same and hit; tested += 1
            nc = sum(l.startswith('CAND') for l in B.splitlines())
            print(f"P y={y} j={j} b={b} G={G}: {'IDENTICAL' if same else 'DIFFERENT'}, {nc} CAND, genuine {'found' if hit else 'MISSING'}  [{B.splitlines()[-1]}]", flush=True)
            break
        if tested >= nmax: break
    print(f"P: {tested} genuine configurations")
    return ok

if __name__ == '__main__':
    args = sys.argv[1:]; exe = args[0]; parts = [a for a in args[1:] if a in ('K', 'C', 'P')] or ['K', 'C', 'P']
    nq = int(args[args.index('--n') + 1]) if '--n' in args else 300
    seed = int(args[args.index('--seed') + 1]) if '--seed' in args else 1
    cpu = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'csrc', 'bottomc3')
    ok = True
    if 'K' in parts: ok &= part_k(exe, nq, seed)
    if 'C' in parts: ok &= part_c(exe, cpu)
    if 'P' in parts: ok &= part_p(exe, cpu)
    print("BOTTOMGPU VALIDATION PASSED" if ok else "BOTTOMGPU VALIDATION FAILED")
    sys.exit(0 if ok else 1)
