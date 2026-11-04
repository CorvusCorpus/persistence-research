"""Fast MITM top step: job export for csrc/mitm_bottom (bottom-window table) and csrc/mitm_top (top probe; or csrc/mitmgpu,
env MITM_TOP), and the exact decision of the candidate pairs with mitm_ref (docs/MITM_STEPS.md section 2).
usage (from work/): python3 mitm_fast.py test             (C bottom table against mitm_ref.bottom_entries)
                    python3 mitm_fast.py test_top         (mitm_top's pairs against brute force over all pairs)
                    python3 mitm_fast.py test_gpu EXE     (csrc/mitmgpu or mitmgpu_emu against mitm_top and mpmath; MG_QUICK=1:
                                                           the six small lattices only)
                    python3 mitm_fast.py test_exact EXE   (the probe's exact decisions, csrc/mitm_exact.h, against mitm_ref)
env MITM_EXACT=py: decide the candidates here with mitm_ref.exact_pairs instead of in the probe (default 'c')"""
import os, sys, subprocess, tempfile
_here = sys.path.pop(0)           # work/struct.py would shadow the standard struct that numpy needs (as in prove_cum_gpu.py)
import struct
import numpy as np
sys.path.insert(0, _here)
import mpmath as mp
import prove_cum, mitm_ref

W = os.path.dirname(os.path.abspath(__file__))
BOTTOM = os.path.join(W, 'csrc', 'mitm_bottom')
TOP = os.environ.get('MITM_TOP', os.path.join(W, 'csrc', 'mitm_top'))     # or csrc/mitmgpu (same job, same output)
VALS = (2, 3, 4, 6, 7, 8, 9)
ENTRY = np.dtype([('key', '<u8'), ('q', '<u8', 3), ('hcode', '<u8'), ('dcode', '<u8'), ('t', '<u4'), ('b', '<u4'), ('sig', '<u8')])

def frac128(x):
    """floor(frac(x) 2^128) for an mpmath real x"""
    f = x - mp.floor(x)
    return int(mp.floor(f * mp.mpf(2) ** 128)) % (1 << 128)

def bottom_constants(lat):
    Mi = lat.Minv
    def row(vec): return [sum(vec[r] * Mi[r, i] for r in range(3)) for i in range(3)]
    tau0 = [mp.mpf(lat.U + 2) / (2 * lat.U), mp.mpf(1) / 2, mp.mpf(1) / 2]
    du = [mp.mpf(-1) / lat.U, 0, -lat.l3 / lat.eps]
    dw = [0, mp.mpf(-1) / lat.W, -lat.l7 / lat.eps]
    dt = [0, 0, -lat.l2 / lat.eps]
    neg = lambda v: [-x for x in v]
    return {'c0': [frac128(x) for x in row(neg(tau0))], 'Au': [frac128(x) for x in row(neg(du))],
            'Aw': [frac128(x) for x in row(neg(dw))], 'At': [frac128(x) for x in row(neg(dt))]}

def write_bottom_job(lat, bag, bot, bmax, path):
    c = bottom_constants(lat)
    cnt = [bag.count(v) for v in VALS]
    lines = [f"Gp {lat.Gp}", f"tm {lat.tm}", f"k {len(bag)}", f"bmax {bmax}", "cnt " + " ".join(map(str, cnt)),
             f"phi {lat.phi}", f"e7 {lat.e7}", f"b1 {lat.b1[0]} {lat.b1[1]}", f"b2 {lat.b2[0]} {lat.b2[1]}",
             "ncell " + " ".join(map(str, lat.ncell))]
    lines += [f"bot {i} {v}" for i, v in sorted(bot.items())]
    lines += [f"{key} " + " ".join(map(str, c[key])) for key in ('c0', 'Au', 'Aw', 'At')]
    open(path, 'w').write("\n".join(lines) + "\n")

def run_bottom(lat, bag, bot, bmax, workdir):
    job = os.path.join(workdir, 'bottom.job'); out = os.path.join(workdir, 'bottom.bin')
    write_bottom_job(lat, bag, bot, bmax, job)
    r = subprocess.run([BOTTOM, job, out], capture_output=True, text=True)
    if r.returncode: raise RuntimeError(f"mitm_bottom failed: {r.stderr[-300:]}")
    with open(out, 'rb') as fh:
        n = struct.unpack('<Q', fh.read(8))[0]
    return np.memmap(out, dtype=ENTRY, mode='r', offset=8, shape=(n,)) if n else np.zeros(0, dtype=ENTRY)

def decode(e):
    """mitm_bottom's window codes: hcode the heights as a bit mask (from the top down), dcode 3-bit VALS indices in that order"""
    b = int(e['b']); m = int(e['hcode'])
    hs = [h for h in range(63, -1, -1) if m >> h & 1]
    if len(hs) != b: raise ValueError("window code")
    ds = tuple(VALS[(int(e['dcode']) >> (3 * r)) & 7] for r in range(b))
    return hs, ds, int(e['t'])

def write_top_job(lat, bag, j, ranges, table, path, exact=False, cands=True):
    """exact: the probe decides the candidates itself (csrc/mitm_exact.h) and prints PAIR lines; cands: CAND lines too"""
    th = [lat.Minv[2, i] / lat.eps for i in range(3)]
    mag = max(int(mp.log10(abs(t))) + 1 if t else 1 for t in th)
    prec = int(mag * 3.33) + 200
    hw = [int(mp.ceil((h + lat.margin) * mp.mpf(2) ** 64)) + 16 for h in lat.hw]
    lines = [f"prec {prec}", f"j {j}", "cnt " + " ".join(str(bag.count(v)) for v in VALS), "ranges " + " ".join(map(str, ranges))]
    lines += [f"theta{i} " + mp.nstr(th[i], mag + 45, strip_zeros=False) for i in range(3)]
    lines += ["hw " + " ".join(map(str, hw)), "ncell " + " ".join(map(str, lat.ncell)), f"table {table}"]
    if exact:
        lines += ["exact 1", f"xcands {int(cands)}", f"xGs {lat.Gs}", f"xGp {lat.Gp}", f"xdps {lat.dps}", f"xU {lat.U}", f"xW {lat.W}",
                  f"xb1 {lat.b1[0]} {lat.b1[1]}", f"xb2 {lat.b2[0]} {lat.b2[1]}", f"xe7 {lat.e7}", f"xphi {lat.phi}",
                  "xTr " + " ".join(str(lat.Tr[i][k]) for i in range(3) for k in range(3))]
    open(path, 'w').write("\n".join(lines) + "\n")

def run_top(lat, bag, j, ranges, workdir, threads=None, exe=None, exact=False, cands=True):
    """(candidates [(o, depths, idx)], exact pairs [(o, depths, idx, (t, u, w, L))] when exact, DONE fields)"""
    job = os.path.join(workdir, 'top.job')
    write_top_job(lat, bag, j, ranges, os.path.join(workdir, 'bottom.bin'), job, exact, cands)
    r = subprocess.run([exe or TOP, job] + ([str(threads)] if threads else []), capture_output=True, text=True)
    if r.returncode: raise RuntimeError(f"{exe or TOP} failed: {r.stderr[-300:]}")
    cl = []; pl = []; done = None
    for line in r.stdout.splitlines():
        f = line.split()
        if f[0] == 'CAND': cl.append((int(f[1]), [int(x) for x in f[2:2 + j]], int(f[2 + j])))
        elif f[0] == 'PAIR': pl.append((int(f[1]), [int(x) for x in f[2:2 + j]], int(f[2 + j]), tuple(int(x) for x in f[3 + j:7 + j])))
        elif f[0] == 'DONE': done = list(map(int, f[1:]))
    if exact and (done is None or len(done) < 4 or done[3] != len(pl)): raise RuntimeError("exact decision: PAIR count mismatch")
    return (cl, pl, done) if exact else (cl, done)

def entry_dict(lat, e):
    hs, ds, t = decode(e)
    F = -1 + sum(9 * (d - 1) * 10 ** h for d, h in zip(ds, hs))
    dl = mitm_ref.dlog3(F * pow(2, -t, lat.M5) % lat.M5, lat.Gp)
    return {'heights': hs, 'digits': ds, 'F': F, 't': t, 'vB': lat.v_rep(dl)}

EXACT = os.environ.get('MITM_EXACT', 'c')       # 'c': the probe decides the candidates (csrc/mitm_exact.h); 'py': mitm_ref here

def fast_pairs(lat, bag, j, ranges, bot, workdir, threads=None):
    """exactly solvable pairs ((ordering, depths), (heights, digits, t)) -> lattice points; and the counters"""
    tab = run_bottom(lat, bag, bot, len(bag) - j, workdir)
    orders = prove_cum.ordered(bag, j); out = {}
    if EXACT == 'c':
        _, pl, done = run_top(lat, bag, j, ranges, workdir, threads, exact=True, cands=False)
        for o, ds, idx, pt in pl:
            hs, dg, t = decode(tab[idx])
            out.setdefault(((orders[o], tuple(ds)), (tuple(hs), dg, t)), []).append(pt)
    else:
        cands, done = run_top(lat, bag, j, ranges, workdir, threads)
        for o, ds, idx in cands:
            e = entry_dict(lat, tab[idx])
            pts = mitm_ref.exact_pairs(lat, mitm_ref.lam_T(orders[o], ds), e)
            if pts: out[((orders[o], tuple(ds)), (tuple(e['heights']), e['digits'], e['t']))] = pts
    return out, {'entries': len(tab), 'cases': done[0], 'probes': done[1], 'candidates': done[2]}

def test_exact(exe):
    """the probe's own decisions (exact 1, PAIR lines) against mitm_ref.exact_pairs on every one of its candidates"""
    six = (2, 2, 3, 3, 7, 7); wide = (2,) * 12 + (3,); ok = True    # wide: windows of up to 11 digits
    cases = [(six, 3000, 12, 6, 4, 2, [6, 7]), (six, 3000, 12, 9, 7, 2, [6, 7]), (six, 3000, 12, 8, 6, 3, [7, 8, 9]),
             (six, 3000, 12, 10, 1, 2, [9, 10]), (six, 10 ** 9, 20, 16, 9, 2, [8, 9]), (six, 10 ** 12, 30, 22, 12, 3, [20, 21, 22]),
             (six, 10 ** 60, 55, 118, 12, 2, [117, 118]),
             (six, 10 ** 90, 55, 185, 2, 2, [184, 185]),               # log10 XL > Gp + 70 (mitm_ref's precision fix)
             (wide, 10 ** 12, 40, 18, 11, 1, [18])]
    for bag, XL, tm, Gs, Gp, j, ranges in cases:
        k = len(bag)
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        with tempfile.TemporaryDirectory() as d:
            tab = run_bottom(lat, bag, {k: 2}, k - j, d)
            cl, pl, done = run_top(lat, bag, j, ranges, d, exe=exe, exact=True, cands=True)
            orders = prove_cum.ordered(bag, j); ref = set()
            for o, ds, idx in cl:
                for pt in mitm_ref.exact_pairs(lat, mitm_ref.lam_T(orders[o], ds), entry_dict(lat, tab[idx])):
                    ref.add((o, tuple(ds), idx, pt))
        got = {(o, tuple(ds), idx, pt) for o, ds, idx, pt in pl}
        good = got == ref and len(cl) == done[2]
        ok &= good
        print(f"XL 10^{len(str(XL)) - 1} Gs {Gs} Gp {Gp} j {j}: {len(cl)} candidates; solutions: probe {len(got)}, mitm_ref {len(ref)} "
              f"-> {'IDENTICAL' if good else 'DIFFERENT'}", flush=True)
    print("MITM_EXACT TEST", "PASSED" if ok else "FAILED")
    return ok

def test_top():
    """the fast probe finds exactly the reference's exactly-solvable pairs"""
    from collections import Counter
    bag = (2, 2, 3, 3, 7, 7); k = len(bag); ok = True
    for XL, tm, Gs, Gp, j, ranges in ((3000, 12, 9, 7, 2, [6, 7]), (3000, 12, 8, 6, 3, [7, 8, 9]), (10 ** 9, 20, 16, 9, 2, [8, 9]),
                                      (3000, 12, 10, 1, 2, [9, 10])):
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        with tempfile.TemporaryDirectory() as d:
            fast, st = fast_pairs(lat, bag, j, ranges, {k: 2}, d)
        tops = list(mitm_ref.top_configs(bag, j, ranges)); ents = mitm_ref.bottom_entries(lat, bag, {k: 2}, k - j)
        full = Counter(bag); ref = {}
        for o, ds in tops:
            lam = mitm_ref.lam_T(o, ds)
            for _, e in ents:
                if not Counter(o) + Counter(e['digits']) <= full: continue
                pts = mitm_ref.exact_pairs(lat, lam, e)
                if pts: ref[((o, tuple(ds)), (tuple(e['heights']), tuple(e['digits']), e['t']))] = pts
        same = set(fast) == set(ref)
        ok &= same
        print(f"top Gs {Gs} Gp {Gp} j {j}: {st}; exactly solvable pairs: fast {len(fast)}, all pairs (brute force) {len(ref)} -> "
              f"{'IDENTICAL' if same else 'DIFFERENT'}", flush=True)
    print("MITM_TOP TEST", "PASSED" if ok else "FAILED")
    return ok

def test_gpu(exe):
    """csrc/mitmgpu (or mitmgpu_emu) against csrc/mitm_top: the same cases, a superset of the candidates (extras only within
    the extra margin of 64 units); and its dumped P_i (MG_DUMP=1) against mpmath at 60 digits beyond |Theta_i|, every case"""
    import random
    bag = (2, 2, 3, 3, 7, 7); k = len(bag); ok = True
    cases = [(3000, 12, 9, 7, 2, [6, 7]), (3000, 12, 8, 6, 3, [7, 8, 9]), (10 ** 9, 20, 16, 9, 2, [8, 9]),
             (10 ** 9, 20, 13, 9, 3, [12, 13, 14]), (10 ** 12, 30, 22, 12, 3, [20, 21, 22]), (10 ** 20, 40, 34, 14, 2, [33, 34]),
             (10 ** 60, 55, 118, 12, 2, [117, 118]), (10 ** 100, 55, 212, 10, 2, [211, 212]), (10 ** 130, 55, 272, 8, 1, [272]),
             (10 ** 180, 55, 372, 8, 1, [372])]                     # the last: beyond the size classes (mitm_top fallback)
    if os.environ.get('MG_QUICK'): cases = cases[:6]
    for XL, tm, Gs, Gp, j, ranges in cases:
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        with tempfile.TemporaryDirectory() as d:
            run_bottom(lat, bag, {k: 2}, k - j, d)
            cpu, dc = run_top(lat, bag, j, ranges, d, exe=os.path.join(W, 'csrc', 'mitm_top'))
            gpu, dg = run_top(lat, bag, j, ranges, d, exe=exe)
            env = dict(os.environ, MG_DUMP='1')
            r = subprocess.run([exe, os.path.join(d, 'top.job')], capture_output=True, text=True, env=env)
        sc = {(o, tuple(ds), i) for o, ds, i in cpu}; sg = {(o, tuple(ds), i) for o, ds, i in gpu}
        same = dc[0] == dg[0] and sc <= sg
        # dump: every case's P_i against mpmath
        th = [lat.Minv[2, i] / lat.eps for i in range(3)]
        orders = prove_cum.ordered(bag, j); worst = 0; npt = 0
        lines = [l.split() for l in r.stdout.splitlines() if l.startswith('PT ')]
        rnd = random.Random(1); sample = lines if len(lines) <= 3000 else rnd.sample(lines, 3000)
        with mp.workdps(int(max(mp.log10(abs(t)) for t in th)) + 60):
            for f in sample:
                o = int(f[1]); ds = [int(x) for x in f[2:2 + j]]; P = [int(x) for x in f[2 + j:5 + j]]
                lam = mp.log1p(sum(9 * (dd - 1) * mp.mpf(10) ** (-D) for dd, D in zip(orders[o], ds)))
                for i in range(3):
                    x = th[i] * lam; want = int(mp.floor((x - mp.floor(x)) * mp.mpf(2) ** 64))
                    dd = (P[i] - want) % (1 << 64); worst = max(worst, min(dd, (1 << 64) - dd))
                npt += 1
        fallback = not r.stdout.startswith('MGDUMP')    # the job went to mitm_top (no dump)
        good = same and r.returncode == 0 and (fallback or (len(lines) == dg[0] and worst <= 4))
        ok &= good
        print(f"Gs {Gs} Gp {Gp} j {j}: cases cpu {dc[0]} gpu {dg[0]}; candidates cpu {len(sc)} gpu {len(sg)} "
              f"(extra {len(sg - sc)}, missing {len(sc - sg)}); " + ("CPU fallback (beyond the size classes)" if fallback else
              f"dump {len(lines)} cases, {npt} checked, worst |dP| {worst} units of 2^-64") + f" -> {'OK' if good else 'FAIL'}", flush=True)
    print("MITM_GPU TEST", "PASSED" if ok else "FAILED")
    return ok

def test():
    ok = True
    six = (2, 2, 3, 3, 7, 7); wide = (2,) * 12 + (3,)                # wide: windows of up to 11 digits (codes, b > 10)
    for bag, bmax, XL, tm, Gs, Gp in ((six, 4, 3000, 12, 9, 1), (six, 4, 3000, 12, 9, 7), (six, 4, 10 ** 12, 30, 20, 12),
                                      (six, 4, 10 ** 20, 40, 34, 18), (wide, 12, 10 ** 12, 40, 22, 11)):
        k = len(bag)
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        bot = {k: 2}
        ref = mitm_ref.bottom_entries(lat, bag, bot, bmax)
        with tempfile.TemporaryDirectory() as d:
            tab = run_bottom(lat, bag, bot, bmax, d)
        refd = {}
        for key, e in ref:
            refd[(tuple(e['heights']), tuple(e['digits']), e['t'])] = (key, [int(mp.floor((q - mp.floor(q)) * mp.mpf(2) ** 64)) for q in e['Q']])
        got = {}
        for e in tab:
            hs, ds, t = decode(e)
            got[(tuple(hs), ds, t)] = (int(e['key']), [int(x) for x in e['q']])
        same_set = set(refd) == set(got)
        worst = 0; keydiff = 0
        for kk in refd.keys() & got.keys():
            for a, b in zip(refd[kk][1], got[kk][1]):
                dd = (a - b) % (1 << 64); worst = max(worst, min(dd, (1 << 64) - dd))
            rk = refd[kk][0]; nc = lat.ncell; ck = (rk[0] * nc[1] + rk[1]) * nc[2] + rk[2]
            keydiff += ck != got[kk][0]
        sorted_ok = bool(np.all(np.diff(tab['key'].astype(np.uint64)) >= 0)) if len(tab) > 1 else True
        good = same_set and worst <= 4 and sorted_ok
        ok &= good
        print(f"Gs {Gs} Gp {Gp}: reference {len(refd)} entries, C {len(got)}; same set {same_set}; worst |dQ| {worst} units of 2^-64; "
              f"cell keys differing {keydiff}; sorted {sorted_ok} -> {'OK' if good else 'FAIL'}", flush=True)
    # large windows (the reference would be slow): every discrete log of the table checked against Pohlig-Hellman in C
    for XL, tm, Gs, Gp in ((10 ** 20, 40, 34, 14), (10 ** 30, 55, 50, 20), (10 ** 40, 55, 70, 27)):
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, 'bottom.job'); write_bottom_job(lat, six, {6: 2}, 4, job)
            r = subprocess.run([BOTTOM, job, os.path.join(d, 'bottom.bin')], capture_output=True, text=True,
                               env=dict(os.environ, MB_CHECK='1'))
        good = r.returncode == 0 and r.stdout.startswith('DONE')
        ok &= good
        print(f"Gp {Gp}: MB_CHECK=1 {r.stdout.strip() or r.stderr.strip()[-200:]} -> {'OK' if good else 'FAIL'}", flush=True)
    print("MITM_BOTTOM TEST", "PASSED" if ok else "FAILED")
    return ok

if __name__ == '__main__':
    if sys.argv[1:] == ['test']: sys.exit(0 if test() else 1)
    if sys.argv[1:] == ['test_top']: sys.exit(0 if test_top() else 1)
    if sys.argv[1:2] == ['test_gpu']: sys.exit(0 if test_gpu(os.path.abspath(sys.argv[2])) else 1)
    if sys.argv[1:2] == ['test_exact']: sys.exit(0 if test_exact(os.path.abspath(sys.argv[2])) else 1)
