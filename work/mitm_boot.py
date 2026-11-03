"""Bootstrapping with MITM top steps (reference implementation, mitm_ref.py), for validation on small bags.

prove_cum's loop with the two-chain bound, standard steps on the CPU path (ctop3/topc3, cbottom4/bottomc3), and a MITM top
step for top levels j >= 1 whenever its threshold Gs is below the standard G by more than 4. A MITM fact D_j < Gs holds for
solutions with L >= S_j = max(top[j-1], Gs) + Gp; with S = max S_j the proof bounds L by max(XL, S - 1), where XL is the
two-chain bound computed from all facts (valid for L >= S).
usage (from work/): PC_BOOTSTRAP=twochain python3 mitm_boot.py BAG ... [--gp N]   (BAG as 2_3_7_7)
Checks per bag: the bound <= 5000, and every genuine solution of the bag (Z.json) satisfies every fact that applies to it."""
import sys, os, math, json, time, tempfile
from fractions import Fraction
import prove_cum, ctop2, ctop3, cbottom2, cbottom4
from prove_iter import FastTop
from solver import valid_points
from twoadic import tmax
import mitm_ref

def mitm_top_step(bag, j, top, bot, XL, tm, Gs, Gp, stats):
    k = len(bag)
    lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
    ranges = prove_cum.deff_top(j, top)
    tops = list(mitm_ref.top_configs(bag, j, ranges))
    ents = mitm_ref.bottom_entries(lat, bag, bot, k - j)
    found, probes = mitm_ref.mitm_join(lat, tops, ents, bag)
    S = max(top.get(j - 1, 0), Gs) + Gp
    ok = True; npts = 0
    for ti, e, lam in found:
        for pt in mitm_ref.exact_pairs(lat, lam, e):
            npts += 1
            if pt[3] < S: continue
            if ctop2.could_be_solution(bag, pt): ok = False
    stats['mitm_top'] = stats.get('mitm_top', 0) + len(tops); stats['mitm_bot'] = stats.get('mitm_bot', 0) + len(ents)
    stats['mitm_pairs'] = stats.get('mitm_pairs', 0) + len(found); stats['mitm_points'] = stats.get('mitm_points', 0) + npts
    return ok, S, lat

def prove_bag_mitm(bag, Gp=8, log=print):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    top, bot = {}, {k: 2}; S = 0
    stats = {'top': 0, 'bottom': 0, 'enum': 0, 'proofs': 0, 'cost': 0, 'mitm_steps': 0}
    while True:
        XL = prove_cum.matveev_cum(k, tm, top, bot)
        G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
        if max(XL, S - 1) <= prove_cum.LIMIT:
            return max(XL, S - 1), top, bot, S, stats
        best = None
        for j in range(k):
            if all(i in top for i in range(j)):
                Gs = max(4, math.ceil(2 * math.log10(XL) + 6 - 0.699 * Gp)) if j else G
                use = j >= 1 and Gs < G - 4
                Gj = Gs if use else G
                if not (j in top and top[j] <= Gj + 4):
                    c = (prove_cum.ncomp(prove_cum.deff_top(j, top)) if j else 1) * max(1, len(prove_cum.ordered(bag, j)))
                    if best is None or c < best[0]: best = (c, j, 'mtop' if use else 'top', Gj)
            if all(i in bot for i in range(j + 1, k + 1)) and not (j in bot and bot[j] <= Gb + 4):
                c = prove_cum.nbottom(cbottom2.heff_of(j, k, bot)) * len(prove_cum.ordered(bag, k - j)) * prove_cum.BOTTOM_COST
                if best is None or c < best[0]: best = (c, j, 'bottom', Gb)
        cost, j, side, Gj = best
        stats['cost'] += cost; stats['proofs'] += 1
        if side == 'mtop':
            ok, Sj, lat = mitm_top_step(bag, j, top, bot, XL, tm, Gj, Gp, stats)
            stats['mitm_steps'] += 1; S = max(S, Sj)
            log(f"    mitm top j={j} XL 10^{math.log10(XL):.1f} G {G} -> Gs {Gj} (Gp {Gp}, det {float(lat.det3):.3g}, cells {lat.ncell}) "
                f"{'ok' if ok else 'FAILED'}")
        elif side == 'top':
            topl = FastTop(XL, tm, G)
            if j == 0:
                ok = topl.lat.kills_homogeneous() or not [pt for pt in valid_points(topl.lat, 0, Fraction(1)) if ctop2.could_be_solution(bag, pt)]
            else:
                D = prove_cum.deff_top(j, top)
                ok = ctop3.prove_top_c(j, bag, G, topl, {i: D[i] for i in range(j)}, stats, prove_cum.ordered(bag, j))
        else:
            U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
            ok = cbottom4.prove_bottom_cum(j, bag, Gj, tm, U, W, bot, stats, prove_cum.ordered(bag, k - j))
        if not ok: raise RuntimeError(f"bag {bag}: level {j} ({side}) not provable at G={Gj}")
        if side in ('top', 'mtop'): top[j] = min(top.get(j, 10 ** 18), Gj)
        else: bot[j] = min(bot.get(j, 10 ** 18), Gj)

def check_genuine(bag, top, bot, S, Z):
    """every genuine solution of the bag satisfies every fact that applies to it (MITM facts only for L >= S)"""
    viol = 0; n = 0
    for y in Z:
        s = str(y)
        if y < 10 or '0' in s or tuple(sorted(int(c) for c in s if c != '1')) != bag: continue
        n += 1; L = len(s); pos = [L - 1 - i for i, ch in enumerate(s) if ch != '1']; a = [L] + pos
        for j, T in top.items():
            if not L - pos[j] < T and L >= S: viol += 1
        for i, H in bot.items():
            if not a[i] < H: viol += 1
    return n, viol

if __name__ == '__main__':
    args = sys.argv[1:]; Gp = 8
    if '--gp' in args: i = args.index('--gp'); Gp = int(args[i + 1]); del args[i:i + 2]
    tmp = tempfile.mkdtemp(); prove_cum.CK = os.path.join(tmp, 'c'); prove_cum.RES = os.path.join(tmp, 'r')
    os.makedirs(prove_cum.CK); os.makedirs(prove_cum.RES)
    Z = json.load(open('Z.json')); allok = True
    for arg in args:
        bag = tuple(sorted(int(x) for x in arg.split('_')))
        t0 = time.time()
        Lb, top, bot, S, st = prove_bag_mitm(bag, Gp)
        n, viol = check_genuine(bag, top, bot, S, Z)
        ok = Lb <= 5000 and viol == 0
        allok &= ok
        print(f"{arg}: L <= {Lb} (S {S}); top {dict(sorted(top.items()))} bot {dict(sorted(bot.items()))}; {st['mitm_steps']} MITM steps "
              f"({st.get('mitm_top', 0)} top configs, {st.get('mitm_bot', 0)} bottom entries, {st.get('mitm_pairs', 0)} pairs probed, "
              f"{st.get('mitm_points', 0)} lattice points); genuine solutions {n}, fact violations {viol}; {time.time() - t0:.0f}s "
              f"-> {'OK' if ok else 'FAIL'}", flush=True)
    print("MITM BOOTSTRAP CHECK", "PASSED" if allok else "FAILED")
    sys.exit(0 if allok else 1)
