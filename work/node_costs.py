"""Estimated cost of even-target nodes in prove_cum's own cost model (dry runs), split into top and bottom.

For every bag of bags.bags_of(node), prove_cum.prove_bag(bag, dry=True) replays the bootstrapping loop (each test assumed
to succeed at its threshold) and adds up the cost of the chosen steps: top cases, and bottom solves times BOTTOM_COST.
The split by side is recovered from the facts before and after each step (every step lowers exactly one fact), by
recomputing that step's cost with prove_cum's formula; the two parts must add up to res['stats']['cost'].
Bags already proved in bagresC/ are reported separately (a dry run ignores checkpoints).
usage (from work/): python3 node_costs.py [--json OUT] [NODE ...]   (default: the open nodes of docs/RESULTS.md 2.3)"""
import sys, os, json, time
import prove_cum, cbottom2
from bags import bags_of
from twoadic import tmax

OPEN = [93312, 314928, 1229312, 11239424, 1161216, 24192, 41472, 122472, 221184, 131712, 33718272]

def bag_cost(b):
    snaps = []
    orig = prove_cum.matveev_cum
    def spy(k, tm, top, bot):                          # called once per step, with the facts before the step
        snaps.append((dict(top), dict(bot)))
        return orig(k, tm, top, bot)
    prove_cum.matveev_cum = spy
    try:
        L, res = prove_cum.prove_bag(b, log=lambda x: None, dry=True)
    finally:
        prove_cum.matveev_cum = orig
    bag = tuple(sorted(b)); k = len(bag)
    split = {'top': 0.0, 'bottom': 0.0}; steps = []
    for (t0, b0), (t1, b1) in zip(snaps, snaps[1:]):
        dt = [j for j in t1 if t0.get(j) != t1[j]]; db = [j for j in b1 if b0.get(j) != b1[j]]
        assert len(dt) + len(db) == 1, (t0, t1, b0, b1)
        if dt:
            j = dt[0]
            c = (prove_cum.ncomp(prove_cum.deff_top(j, t0)) if j else 1) * max(1, len(prove_cum.ordered(bag, j)))
            split['top'] += c; steps.append(('top', j, c))
        else:
            j = db[0]
            c = prove_cum.nbottom(cbottom2.heff_of(j, k, b0)) * len(prove_cum.ordered(bag, k - j)) * prove_cum.BOTTOM_COST
            split['bottom'] += c; steps.append(('bottom', j, c))
    tot = res['stats']['cost']
    assert abs(split['top'] + split['bottom'] - tot) <= 1e-9 * max(1.0, tot), (split, tot)
    return L, tot, split, steps

def main(argv):
    out = None
    if argv[:1] == ['--json']: out = argv[1]; argv = argv[2:]
    nodes = [int(x) for x in argv] or OPEN
    allrows = {}
    for s in nodes:
        t0 = time.time()
        bs = bags_of(s)
        row = {'node': s, 'bags': len(bs), 'kmax': max(len(b) for b in bs), 'total': 0.0, 'top': 0.0, 'bottom': 0.0,
               'done_bags': 0, 'done_cost': 0.0, 'no_estimate': [], 'per_bag': []}
        for b in bs:
            key = '_'.join(map(str, sorted(b)))
            try:                                       # prove_cum.matveev_cum asserts when the a-priori bound leaves its float range
                prove_cum.matveev_cum(len(b), tmax(tuple(sorted(b))), {}, {len(b): 2})
            except AssertionError:
                row['no_estimate'].append(key); continue
            L, tot, split, steps = bag_cost(b)
            done = os.path.exists(os.path.join(prove_cum.RES, key + '.json'))
            row['total'] += tot; row['top'] += split['top']; row['bottom'] += split['bottom']
            if done: row['done_bags'] += 1; row['done_cost'] += tot
            row['per_bag'].append({'bag': key, 'k': len(b), 'cost': tot, 'top': split['top'], 'bottom': split['bottom'],
                                   'done': done, 'steps': len(steps)})
        row['per_bag'].sort(key=lambda r: -r['cost'])
        allrows[s] = row
        big = ', '.join(f"{r['bag']} {r['cost']:.2e}" for r in row['per_bag'][:3])
        ne = row['no_estimate']
        kne = max((len(x.split('_')) for x in ne), default=0)
        print(f"node {s}: {row['bags']} bags, kmax {row['kmax']}, cost {row['total']:.3e} (top {row['top']:.3e}, bottom "
              f"{row['bottom']:.3e}); proved bags {row['done_bags']} ({row['done_cost']:.2e}); largest: {big}"
              + (f"; NO ESTIMATE for {len(ne)} bags (k up to {kne}: Matveev bound overflows prove_cum)" if ne else "")
              + f"  [{time.time() - t0:.0f}s]", flush=True)
    order = sorted(allrows, key=lambda x: (len(allrows[x]['no_estimate']) > 0, allrows[x]['total']))
    print("cheapest first:", ', '.join(f"{s} {allrows[s]['total']:.2e}" + ('+?' if allrows[s]['no_estimate'] else '') for s in order))
    if out:
        with open(out + '.tmp', 'w') as fh: json.dump(allrows, fh, indent=1)
        os.replace(out + '.tmp', out)

if __name__ == '__main__':
    main(sys.argv[1:])
