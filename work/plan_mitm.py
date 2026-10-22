"""Cost model of the meet-in-the-middle steps and the deep-digit split (docs/MITM_STEPS.md); dry runs only.

Facts and bound are prove_cum's with the two-chain bootstrap (prove_cum.twochain_cum). A plan is a sequence of rounds; each
round re-proves one side's levels up to its frontier at one threshold, with one kind of step:
  standard  top levels at G = ceil(2 log10 XL) + 12 (bottom at max(G, tm + 8)), cost = configurations (top 1, bottom 4 units);
  MITM      threshold Gs < G; each level pairs its configurations with a table of windows of the other side (bottom: lowest
            digits at heights < Gp; top: top digits at depths < Gt), 5-adic congruence only.  False pairs
            N_self * N_table * 108 * 10^(2 log10 XL - precision), each CFP units, kept at most equal to the step's cost.
            Tables are charged once per round (4 units per bottom-window entry, 1 per top-window entry);
  deep      threshold G; the digits beyond the linear boundary (top: depth >= (G+5)/2; bottom: height >= G/2) are matched in
            two halves per configuration of the others (deep_cost).
Each round takes the side and kind with the smallest cost per digit of log10 XL gained.
usage (from work/): python3 plan_mitm.py BAG ...   (BAG as 2_2_3_7)   env MODES (default 0,1,2,3 = standard, MITM, deep,
MITM or deep per round), CFP (1000), SHOW=1 (print the rounds)"""
import sys, os, math
import prove_cum, cbottom2
from twoadic import tmax

INF = float('inf')
CFP = float(os.environ.get('CFP', '1000'))
BC = prove_cum.BOTTOM_COST
L5 = math.log10(5)

def nord(bag, r):
    return max(1, len(prove_cum.ordered(bag, r))) if r else 1

def ranges_window(top, a, Gw):
    """depth bounds D_i < r_i (i < a) of the top a digits inside the window [1, Gw)"""
    out = []
    for i in range(a):
        r = Gw - (a - 1 - i)
        for ip in range(i, a):
            if ip in top: r = min(r, top[ip] - (ip - i))
        out.append(r)
    return out

def NT_window(bag, top, Gw, amax):
    """top windows: the top a digits (a <= amax) at depths < Gw, digit a+1 at depth >= Gw"""
    k = len(bag); tot = 0
    for a in range(0, amax + 1):
        if a < k and top.get(a, INF) <= Gw: continue          # D_a < top[a] <= Gw: digit a+1 lies inside the window
        if a == 0: tot += 1; continue
        rg = ranges_window(top, a, Gw)
        if min(rg) <= 0: continue
        tot += prove_cum.ncomp(rg) * nord(bag, a)
    return tot

def heff_window(bot, k, b, Gw):
    """height bounds a_i < h of the lowest b digits i = k-b+1..k inside the window [0, Gw)"""
    out = []
    for i in range(k - b + 1, k + 1):
        h = Gw - (i - (k - b + 1))
        for ip in range(k - b + 1, i + 1):
            if ip in bot: h = min(h, bot[ip] - (i - ip))
        out.append(h)
    return out

def NB_window(bag, bot, Gw, bmax):
    """bottom windows: the lowest b digits (b <= bmax) at heights < Gw, digit k-b at height >= Gw"""
    k = len(bag); tot = 0
    for b in range(0, bmax + 1):
        if b < k and bot.get(k - b, INF) <= Gw: continue      # digit k-b lies inside the window
        if b == 0: tot += 1; continue
        he = heff_window(bot, k, b, Gw)
        if min(he) <= 0: continue
        tot += prove_cum.nbottom(he) * nord(bag, b)
    return tot

def count_inc(lo, his):
    """strictly increasing x_0 < x_1 < ... with lo <= x_i < his[i]"""
    cnt = {lo - 1: 1}
    for hi in his:
        nc = {}
        for last, c in cnt.items():
            for x in range(last + 1, hi): nc[x] = nc.get(x, 0) + c
        cnt = nc
        if not cnt: return 0
    return sum(cnt.values())

def deep_cost(bounds, lo0, h, width, nvals):
    """positions x_0 < x_1 < ... (x_i < bounds[i]) split into the a below h, enumerated, and the b at >= h, whose
    contribution is linear given the first a and is matched in two halves (n1 + n2 instead of n1 * n2), after a setup of
    width * nvals vectors per configuration of the first a"""
    tot = 0
    for a in range(len(bounds) + 1):
        nA = count_inc(lo0, [min(h, r) for r in bounds[:a]])
        if nA == 0: continue
        rest = bounds[a:]; b = len(rest)
        if b == 0: tot += nA; continue
        if b == 1: tot += nA * count_inc(h, rest); continue
        b1 = (b + 1) // 2
        tot += nA * (width * nvals + count_inc(h, rest[:b1]) + count_inc(h, rest[b1:]))
    return tot

def need(XL, ns, npart):
    """digits of precision (self + partner) that keep the false pairs' cost at most the step's own cost"""
    ns, npart = max(1, ns), max(1, npart)
    return 2 * math.log10(XL) + math.log10(108 * CFP) + math.log10(ns * npart / (ns + npart))

def round_top(bag, XL, top, bot, Gs, kind):
    k = len(bag); top = dict(top); cost = 0.0; tables = {}; lv = []
    J = 0
    while J < k and J in top: J += 1                           # frontier: first level without a fact
    if J >= k: J = k - 1
    for j in range(0, J + 1):
        if j in top and top[j] <= Gs + 4: continue
        nt = (prove_cum.ncomp(prove_cum.deff_top(j, top)) if j else 1) * nord(bag, j)
        if kind == 'deep' and j >= 2:
            h = math.ceil((Gs + 5) / 2)
            nt = min(nt, deep_cost(prove_cum.deff_top(j, top), 1, h, Gs - h, len(set(bag))) * nord(bag, j))
        if j and kind == 'mitm':
            got = None
            for Gp in range(1, 80):
                nb = NB_window(bag, bot, Gp, k - j)
                if Gs + L5 * Gp >= need(XL, nt, BC * nb): got = (Gp, nb); break
                if nb > 1e14: break
            if got is None: return None
            Gp, nb = got
            if tables.get(Gp, 0) < nb: cost += BC * (nb - tables.get(Gp, 0)); tables[Gp] = nb
        cost += nt; top[j] = min(top.get(j, INF), Gs); lv.append(j)
    return cost, top, lv

def round_bot(bag, XL, top, bot, Gs, kind):
    k = len(bag); bot = dict(bot); cost = 0.0; tables = {}; lv = []
    I = k
    while I > 0 and I in bot: I -= 1                           # frontier: highest level without a fact
    if I <= 0: I = 1
    for i in range(k, I - 1, -1):
        if i in bot and bot[i] <= Gs + 4: continue
        nb = prove_cum.nbottom(cbottom2.heff_of(i, k, bot)) * nord(bag, k - i) if i < k else 1
        if kind == 'deep' and k - i >= 2:
            h = math.ceil(Gs / 2)
            nb = min(nb, deep_cost(list(reversed(cbottom2.heff_of(i, k, bot))), 0, h, Gs - h, len(set(bag))) * nord(bag, k - i))
        if kind == 'mitm':
            got = None
            for Gt in range(1, 80):
                ntw = NT_window(bag, top, Gt, i)
                if L5 * Gs + Gt >= need(XL, BC * nb, ntw): got = (Gt, ntw); break
                if ntw > 1e14: break
            if got is None: return None
            Gt, ntw = got
            if tables.get(Gt, 0) < ntw: cost += ntw - tables.get(Gt, 0); tables[Gt] = ntw
        cost += BC * nb; bot[i] = min(bot.get(i, INF), Gs); lv.append(i)
    return cost, bot, lv

def plan(bag, mode, verbose=False):
    """mode 0 standard, 1 standard or MITM, 2 standard or deep, 3 any of the three (chosen per round)"""
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag)
    top, bot = {}, {k: 2}; total = 0.0; rounds = 0
    while True:
        XL = prove_cum.twochain_cum(k, tm, top, bot)
        if XL <= prove_cum.LIMIT: break
        G = math.ceil(2 * math.log10(XL)) + 12; Gb = max(G, tm + 8)
        opts = []
        for side in ('top', 'bot'):
            grid = [(G if side == 'top' else Gb, 'std')]
            if mode in (1, 3): grid += [(g, 'mitm') for g in range(max(4, G - 60), G)]
            if mode in (2, 3): grid += [(G if side == 'top' else Gb, 'deep')]
            for Gs, kind in grid:
                r = (round_top if side == 'top' else round_bot)(bag, XL, top, bot, Gs, kind)
                if r is None: continue
                cost, facts, lv = r
                nt, nbt = (facts, bot) if side == 'top' else (top, facts)
                XLn = prove_cum.twochain_cum(k, tm, nt, nbt)
                gain = math.log10(XL) - math.log10(XLn)
                if gain <= 0.05: continue
                opts.append((cost / gain, cost, side, kind, Gs, nt, nbt, lv, XLn))
        if not opts: raise RuntimeError(f"no round possible at XL 10^{math.log10(XL):.1f}")
        _, cost, side, kind, Gs, top, bot, lv, XLn = min(opts, key=lambda o: o[0])
        total += cost; rounds += 1
        if verbose:
            print(f"    round {rounds:2d} {side} {kind:4s} G {G} -> {Gs}, levels {lv[0]}..{lv[-1]}: cost {cost:.2e}, "
                  f"log10 XL {math.log10(XL):.1f} -> {math.log10(XLn):.1f}", flush=True)
    return total, rounds

if __name__ == '__main__':
    names = ['standard', 'MITM', 'deep', 'MITM or deep']
    for arg in sys.argv[1:]:
        bag = tuple(sorted(int(x) for x in arg.split('_'))); out = []
        for m in [int(x) for x in os.environ.get('MODES', '0,1,2,3').split(',')]:
            tot, nr = plan(bag, m, verbose=bool(os.environ.get('SHOW')))
            out.append(f"{names[m]} {tot:.3e} ({nr} rounds)")
        print(f"{arg}: " + ";  ".join(out), flush=True)
