"""End-to-end check: with the Python cost model (BOTTOM_COST = 25) the C-backed prover must follow the
same proof chain as the recorded Python runs: same final bounds and the same numbers of top cases and
bottom solves."""
import sys, json, os, time, tempfile, re
import prove_iter3, prove_iter2
tmp = tempfile.mkdtemp()
prove_iter2.BAGRES = os.path.join(tmp, 'r'); prove_iter2.CKPT = os.path.join(tmp, 'c')
os.makedirs(prove_iter2.BAGRES); os.makedirs(prove_iter2.CKPT)
prove_iter2.BOTTOM_COST = 25.0
pat = re.compile(r"^\s+bag \(([\d, ]+)\): k=(\d+) tmax=(\d+) bounds \[([\d, ]+)\] -> every solution has L <= (\d+) \(top (\d+), bottom (\d+), enum (\d+), proofs (\d+), (\d+)s\)")
rec = {}
for fn in ['runJ_A.log', 'runJ_B.log']:
    for line in open(fn):
        m = pat.match(line)
        if m:
            bag = tuple(sorted(int(x) for x in m.group(1).replace(' ', '').split(',') if x))
            rec[bag] = dict(bounds=[int(x) for x in m.group(4).replace(' ', '').split(',')], top=int(m.group(6)),
                            bottom=int(m.group(7)), enum=int(m.group(8)), secs=int(m.group(10)))
for bag in [tuple(map(int, a.split(','))) for a in sys.argv[1:]]:
    out = []
    t0 = time.time()
    Lmax, bounds = prove_iter3.prove_bag(bag, log=out.append)
    dt = time.time() - t0
    st = json.load(open(os.path.join(prove_iter2.BAGRES, '_'.join(map(str, sorted(bag))) + '.json')))['stats']
    r = rec[tuple(sorted(bag))]
    print(f"{bag}: bounds equal {bounds == r['bounds']}; top cases {st['top']} vs {r['top']}; bottom solves {st['bottom']} vs {r['bottom']}; "
          f"survivor runs {st['enum']} vs Python enum {r['enum']}; time {dt:.0f}s vs {r['secs']}s", flush=True)
