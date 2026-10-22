"""The GPU step drivers resume like the CPU ones: ctop4 vs ctop3 and cbottom5 vs cbottom4 started from the same
artificial checkpoints (ordering, first gap / lowest height) on steps of a bag's plan; verdicts and counters must agree.
usage (from work/): python3 val_resume_gpu.py [BAG]"""
import sys, math
import prove_cum, cbottom2, prof_steps, ctop3, ctop4, cbottom4, cbottom5
from prove_iter import FastTop
from twoadic import tmax

def main(bag):
    bag = tuple(sorted(bag)); k = len(bag); tm = tmax(bag); ok = True; n = 0
    for side, j, XL, top, bot, cost in prof_steps.plan_states(bag):
        G = math.ceil(2 * math.log10(XL)) + 12
        if side == 'top' and j >= 1:
            T = FastTop(XL, tm, G); orders = prove_cum.ordered(bag, j); known = dict(enumerate(prove_cum.deff_top(j, top)))
            rng0 = known[0]
            for skip in [0, (len(orders) // 2, 1), (len(orders) // 2, max(2, rng0 // 3)), (len(orders) - 1, max(2, rng0 // 2))]:
                if not isinstance(skip, int) and skip[1] >= rng0: continue
                sa = {'top': 0, 'enum': 0}; sb = {'top': 0, 'enum': 0}; pa = []; pb = []
                ra = ctop3.prove_top_c(j, bag, G, T, known, sa, orders, skip, lambda d: pa.append(d))
                rb = ctop4.prove_top_c(j, bag, G, T, known, sb, orders, skip, lambda d: pb.append(d))
                same = ra == rb and sa == sb; ok &= same; n += 1
                print(f"top    j={j} G={G} skip={skip}: ctop3 {ra} {sa} | ctop4 {rb} {sb} -> {'SAME' if same else 'DIFFERENT'}", flush=True)
        elif side == 'bottom':
            Gb = max(G, tm + 8); U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
            orders = prove_cum.ordered(bag, k - j); heff = cbottom2.heff_of(j, k, bot)
            for skip in [0, (len(orders) // 2, 0), (len(orders) // 2, 1), (len(orders) - 1, 1)]:
                if not isinstance(skip, int) and skip[1] >= heff[-1]: continue
                sa = {'bottom': 0}; sb = {'bottom': 0}
                ra = cbottom4.prove_bottom_cum(j, bag, Gb, tm, U, W, bot, sa, orders, skip)
                rb = cbottom5.prove_bottom_cum(j, bag, Gb, tm, U, W, bot, sb, orders, skip)
                same = ra == rb and sa == sb; ok &= same; n += 1
                print(f"bottom j={j} G={Gb} skip={skip}: cbottom4 {ra} {sa} | cbottom5 {rb} {sb} -> {'SAME' if same else 'DIFFERENT'}", flush=True)
    print(f"{n} resumed steps: {'ALL SAME' if ok else 'DIFFERENCES'}")
    return ok

if __name__ == '__main__':
    sys.exit(0 if main(tuple(int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '2_2_2_3_3_7').split('_'))) else 1)
