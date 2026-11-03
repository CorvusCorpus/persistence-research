"""Fast MITM top step: job export for csrc/mitm_bottom (bottom-window table) and csrc/mitm_top (top probe), and the exact
decision of the candidate pairs with mitm_ref (docs/MITM_STEPS.md section 2).
usage (from work/): python3 mitm_fast.py test   (C bottom table against mitm_ref.bottom_entries)"""
import os, sys, subprocess, tempfile
_here = sys.path.pop(0)           # work/struct.py would shadow the standard struct that numpy needs (as in prove_cum_gpu.py)
import struct
import numpy as np
sys.path.insert(0, _here)
import mpmath as mp
import prove_cum, mitm_ref

W = os.path.dirname(os.path.abspath(__file__))
BOTTOM = os.path.join(W, 'csrc', 'mitm_bottom')
TOP = os.path.join(W, 'csrc', 'mitm_top')
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
        tab = np.frombuffer(fh.read(), dtype=ENTRY, count=n)
    return tab

def decode(e):
    b = int(e['b']); hs = [(int(e['hcode']) >> (6 * r)) & 63 for r in range(b)]; ds = tuple((int(e['dcode']) >> (4 * r)) & 15 for r in range(b))
    return hs, ds, int(e['t'])

def write_top_job(lat, bag, j, ranges, table, path):
    th = [lat.Minv[2, i] / lat.eps for i in range(3)]
    mag = max(int(mp.log10(abs(t))) + 1 if t else 1 for t in th)
    prec = int(mag * 3.33) + 200
    hw = [int(mp.ceil((h + lat.margin) * mp.mpf(2) ** 64)) + 16 for h in lat.hw]
    lines = [f"prec {prec}", f"j {j}", "cnt " + " ".join(str(bag.count(v)) for v in VALS), "ranges " + " ".join(map(str, ranges))]
    lines += [f"theta{i} " + mp.nstr(th[i], mag + 45, strip_zeros=False) for i in range(3)]
    lines += ["hw " + " ".join(map(str, hw)), "ncell " + " ".join(map(str, lat.ncell)), f"table {table}"]
    open(path, 'w').write("\n".join(lines) + "\n")

def run_top(lat, bag, j, ranges, workdir, threads=None):
    job = os.path.join(workdir, 'top.job')
    write_top_job(lat, bag, j, ranges, os.path.join(workdir, 'bottom.bin'), job)
    r = subprocess.run([TOP, job] + ([str(threads)] if threads else []), capture_output=True, text=True)
    if r.returncode: raise RuntimeError(f"mitm_top failed: {r.stderr[-300:]}")
    cands = []; done = None
    for line in r.stdout.splitlines():
        f = line.split()
        if f[0] == 'CAND': cands.append((int(f[1]), [int(x) for x in f[2:2 + j]], int(f[2 + j])))
        elif f[0] == 'DONE': done = list(map(int, f[1:]))
    return cands, done

def entry_dict(lat, e):
    hs, ds, t = decode(e)
    F = -1 + sum(9 * (d - 1) * 10 ** h for d, h in zip(ds, hs))
    dl = mitm_ref.dlog3(F * pow(2, -t, lat.M5) % lat.M5, lat.Gp)
    return {'heights': hs, 'digits': ds, 'F': F, 't': t, 'vB': lat.v_rep(dl)}

def fast_pairs(lat, bag, j, ranges, bot, workdir, threads=None):
    """exactly solvable pairs ((ordering, depths), (heights, digits, t)) -> lattice points; and the counters"""
    tab = run_bottom(lat, bag, bot, len(bag) - j, workdir)
    cands, done = run_top(lat, bag, j, ranges, workdir, threads)
    orders = prove_cum.ordered(bag, j); out = {}
    for o, ds, idx in cands:
        e = entry_dict(lat, tab[idx])
        pts = mitm_ref.exact_pairs(lat, mitm_ref.lam_T(orders[o], ds), e)
        if pts: out[((orders[o], tuple(ds)), (tuple(e['heights']), e['digits'], e['t']))] = pts
    return out, {'entries': len(tab), 'cases': done[0], 'probes': done[1], 'candidates': done[2]}

def test_top():
    """the fast probe finds exactly the reference's exactly-solvable pairs"""
    from collections import Counter
    bag = (2, 2, 3, 3, 7, 7); k = len(bag); ok = True
    for XL, tm, Gs, Gp, j, ranges in ((3000, 12, 9, 7, 2, [6, 7]), (3000, 12, 8, 6, 3, [7, 8, 9]), (10 ** 9, 20, 16, 9, 2, [8, 9])):
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

def test():
    bag = (2, 2, 3, 3, 7, 7); k = len(bag)
    ok = True
    for XL, tm, Gs, Gp in ((3000, 12, 9, 7), (10 ** 12, 30, 20, 12), (10 ** 20, 40, 34, 18)):
        lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
        bot = {k: 2}
        ref = mitm_ref.bottom_entries(lat, bag, bot, k - 2)
        with tempfile.TemporaryDirectory() as d:
            tab = run_bottom(lat, bag, bot, k - 2, d)
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
    print("MITM_BOTTOM TEST", "PASSED" if ok else "FAILED")
    return ok

if __name__ == '__main__':
    if sys.argv[1:] == ['test']: sys.exit(0 if test() else 1)
    if sys.argv[1:] == ['test_top']: sys.exit(0 if test_top() else 1)
