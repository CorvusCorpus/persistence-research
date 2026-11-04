"""prove_cum's bootstrapping with MITM top steps (docs/MITM_STEPS.md section 2), selected by PC_MITM=1 (two chains required).

Rounds as in plan_mitm.plan: each round re-proves one side's levels up to its frontier at one threshold, with standard steps
(top at G, bottom at max(G, tm + 8)) or MITM top steps at Gs < G; the round with the smallest estimated cost per digit of
log10 XL gained is taken. A MITM top step at level j with window Gp (the smallest that keeps the expected false pairs within
the step's cost, with tables of at most MITM_TABLE_CAP entries) proves D_j < Gs for every solution with
L >= S_j = max(top[j-1], Gs) + Gp; the bag's bound is max(XL, S - 1), S = max S_j over all steps so far, with XL the
two-chain bound of all facts: a solution with L >= S satisfies every fact, so L <= XL. Every step's box uses max(XL, S - 1).
Level 0 takes the homogeneous test at G (in a MITM round with Gs >= G - 8, at Gs first).
Records (bagresC) carry 'bootstrap': 'twochain', 'mitm': 1 and 'lmin': S; checkpoints carry S.
usage: as prove_cum.prove_bag (prove_cum dispatches here when PC_MITM=1)."""
import os, math, json, time, tempfile
from fractions import Fraction
import prove_cum, cbottom2, ctop2
from prove_iter import FastTop
from solver import valid_points
from twoadic import tmax
import plan_mitm as PL
import mitm_ref, mitm_fast

CAP = float(os.environ.get('MITM_TABLE_CAP', '2.5e7'))    # 64-byte entries; mitm_bottom's peak is about 2.5x the table
PL.CAP = CAP; PL.GPMAX = 27                               # mitm_bottom's 64-bit arithmetic needs Gp <= 27
# Steps with tables above MITM_BIG entries run one at a time across processes (an flock on MITM_LOCK): several workers
# building 2e7-entry tables at once would exhaust the host's and the GPU's memory.
BIG = float(os.environ.get('MITM_BIG', '2e6'))
LOCK = os.environ.get('MITM_LOCK', os.path.join(tempfile.gettempdir(), 'mitm_big_table.lock'))

class big_step:
    def __init__(self, nb): self.on = nb > BIG
    def __enter__(self):
        if self.on:
            import fcntl
            self.fh = open(LOCK, 'a'); fcntl.flock(self.fh, fcntl.LOCK_EX)
    def __exit__(self, *a):
        if self.on:
            import fcntl
            fcntl.flock(self.fh, fcntl.LOCK_UN); self.fh.close()

def choose_round(bag, XL, top, bot, tm):
    k = len(bag)
    G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
    opts = []
    for side in ('top', 'bot'):
        grid = [(G if side == 'top' else Gb, 'std')]
        if side == 'top': grid += [(g, 'mitm') for g in range(max(4, G - 60), G)]
        for Gs, kind in grid:
            r = (PL.round_top if side == 'top' else PL.round_bot)(bag, XL, top, bot, Gs, kind)
            if r is None: continue
            cost, facts, lv = r
            nt, nbt = (facts, bot) if side == 'top' else (top, facts)
            XLn = prove_cum.twochain_cum(k, tm, nt, nbt)
            gain = math.log10(XL) - math.log10(XLn)
            if gain <= 0.05: continue
            opts.append((cost / gain, cost, side, kind, Gs))
    return min(opts, key=lambda o: o[0])

def window_for(bag, XL, top, bot, j, Gs):
    k = len(bag)
    nt = prove_cum.ncomp(prove_cum.deff_top(j, top)) * PL.nord(bag, j)
    for Gp in range(1, PL.GPMAX + 1):
        nb = PL.NB_window(bag, bot, Gp, k - j)
        if nb > CAP: break
        if Gs + PL.L5 * Gp >= PL.need(XL, nt, PL.BC * nb): return Gp, nb
    return None, 0

def mitm_step(bag, j, top, bot, XL, tm, Gs, Gp, nb, stats):
    lat = mitm_ref.MitmLattice(XL, tm, Gs, Gp)
    ranges = prove_cum.deff_top(j, top)
    with big_step(nb), tempfile.TemporaryDirectory(prefix="mitm_") as workdir:   # the step's table and jobs
        pairs, st = mitm_fast.fast_pairs(lat, bag, j, ranges, bot, workdir)
    S = max(top.get(j - 1, 0), Gs) + Gp
    ok = True
    for key, pts in pairs.items():
        for pt in pts:
            if pt[3] >= S and ctop2.could_be_solution(bag, pt): ok = False
    stats['mitm_top'] = stats.get('mitm_top', 0) + st['cases']; stats['mitm_bot'] = stats.get('mitm_bot', 0) + st['entries']
    stats['mitm_cand'] = stats.get('mitm_cand', 0) + st['candidates']; stats['mitm_pairs'] = stats.get('mitm_pairs', 0) + len(pairs)
    return ok, S

def prove_bag(bag, log=print, limit=prove_cum.LIMIT, dry=False):
    assert prove_cum.BOOTSTRAP == 'twochain', "PC_MITM needs PC_BOOTSTRAP=twochain"
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    key = '_'.join(map(str, bag))
    rpath = os.path.join(prove_cum.RES, key + '.json'); cpath = os.path.join(prove_cum.CK, key + '.json')
    if not dry and os.path.exists(rpath):
        r = json.load(open(rpath)); log(f"  bag {bag}: proved earlier (cumulative) -> L <= {r['Lmax']}"); return r['Lmax'], r
    top, bot, S = {}, {k: 2}, 0
    stats = {'top': 0, 'bottom': 0, 'enum': 0, 'proofs': 0, 'cost': 0, 'mitm_steps': 0}
    if not dry and os.path.exists(cpath):
        c = json.load(open(cpath))
        top = {int(a): v for a, v in c['top'].items()}; bot = {int(a): v for a, v in c['bot'].items()}
        stats = c['stats']; S = c.get('lmin', 0)
        log(f"  bag {bag}: resuming top={top} bot={bot} lmin={S}")
    t0 = time.time()
    def save():
        if not dry: prove_cum._save(cpath, {'top': top, 'bot': bot, 'stats': stats, 'partial': None, 'lmin': S})
    while True:
        XL = max(prove_cum.matveev_cum(k, tm, top, bot), S - 1)
        if max(XL, S - 1) <= limit:
            res = {'bag': list(bag), 'Lmax': max(XL, S - 1), 'top': top, 'bot': bot, 'stats': stats,
                   'bootstrap': 'twochain', 'mitm': 1, 'lmin': S}
            if prove_cum.GB_FLOOR: res['gb_floor'] = 1
            if not dry: prove_cum._save(rpath, res)
            log(f"  bag {bag}: k={k} tmax={tm} top {dict(sorted(top.items()))} bot {dict(sorted(bot.items()))} lmin {S} -> every solution "
                f"has L <= {res['Lmax']} (top cases {stats['top']}, bottom solves {stats['bottom']}, MITM steps {stats['mitm_steps']}, "
                f"{time.time()-t0:.0f}s)")
            return res['Lmax'], res
        _, cost, side, kind, Gs = choose_round(bag, XL, top, bot, tm)
        G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
        before = (dict(top), dict(bot))
        if side == 'top':
            J = 0
            while J < k and J in top: J += 1
            if J >= k: J = k - 1
            for j in range(0, J + 1):
                if j in top and top[j] <= Gs + 4: continue
                XL = max(prove_cum.matveev_cum(k, tm, top, bot), S - 1)
                nt = (prove_cum.ncomp(prove_cum.deff_top(j, top)) if j else 1) * max(1, len(prove_cum.ordered(bag, j)))
                stats['cost'] += nt; stats['proofs'] += 1
                ok = True
                Gp, nb = window_for(bag, XL, top, bot, j, Gs) if (j and kind == 'mitm') else (None, 0)
                if j == 0:                                   # homogeneous test at G; at Gs first if Gs >= G - 8 (below that the
                    gs0 = [Gs, G] if kind == 'mitm' and Gs >= G - 8 else [G]   # box holds too many lattice points to list)
                    Gj = gs0[0]
                    if not dry:
                        for g in gs0:
                            topl = FastTop(XL, tm, g); Gj = g
                            ok = topl.lat.kills_homogeneous() or not [pt for pt in valid_points(topl.lat, 0, Fraction(1)) if ctop2.could_be_solution(bag, pt)]
                            if ok: break
                elif Gp is not None:                         # MITM top step
                    Gj = Gs; stats['cost'] += PL.BC * nb
                    if not dry: ok, Sj = mitm_step(bag, j, top, bot, XL, tm, Gs, Gp, nb, stats)
                    else: Sj = max(top.get(j - 1, 0), Gs) + Gp
                    S = max(S, Sj); stats['mitm_steps'] += 1
                    log(f"    bag {bag}: MITM top j={j} XL 10^{math.log10(XL):.1f} G {G} -> {Gs} (window {Gp}, table {nb:.3g}, "
                        f"top configurations {nt:.3g}) {'ok' if ok else 'FAILED'}")
                else:                                        # standard top step at G
                    Gj = G
                    if not dry:
                        topl = FastTop(XL, tm, G); D = prove_cum.deff_top(j, top)
                        ok = prove_cum.ctop3.prove_top_c(j, bag, G, topl, {i: D[i] for i in range(j)}, stats, prove_cum.ordered(bag, j))
                if not ok: raise RuntimeError(f"bag {bag}: top level {j} not provable at G={Gj}")
                top[j] = min(top.get(j, 10 ** 18), Gj); save()
        else:
            I = k
            while I > 0 and I in bot: I -= 1
            if I <= 0: I = 1
            for i in range(k, I - 1, -1):
                if i in bot and bot[i] <= Gb + 4: continue
                XL = max(prove_cum.matveev_cum(k, tm, top, bot), S - 1)
                G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
                g = prove_cum.gb_step(i, k, bot, Gb)
                nb = prove_cum.nbottom(cbottom2.heff_of(i, k, bot)) * len(prove_cum.ordered(bag, k - i))
                stats['cost'] += nb * prove_cum.BOTTOM_COST; stats['proofs'] += 1
                if not dry:
                    U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
                    ok = prove_cum.cbottom4.prove_bottom_cum(i, bag, g, tm, U, W, bot, stats, prove_cum.ordered(bag, k - i))
                    if not ok: raise RuntimeError(f"bag {bag}: bottom level {i} not provable at G={g}")
                bot[i] = min(bot.get(i, 10 ** 18), g); save()
        if (top, bot) == before: raise RuntimeError(f"bag {bag}: round {side} {kind} at {Gs} changed no fact")
