"""End-to-end runs of the cumulative prover on bags already proved by the gap chain; also checks that every
genuine solution of the bag (from the exhaustive list) satisfies every proved fact D_j < top[j], a_i < bot[i]."""
import sys, os, json, time, tempfile, re
import prove_cum
tmp = tempfile.mkdtemp(); prove_cum.CK = os.path.join(tmp, 'c'); prove_cum.RES = os.path.join(tmp, 'r')
os.makedirs(prove_cum.CK); os.makedirs(prove_cum.RES); prove_cum.USE_GAP_RESULTS = False
Z = json.load(open('Z.json'))
pat = re.compile(r"^\s+bag \(([\d, ]+)\): k=(\d+).* \((top \d+, bottom \d+), enum \d+, proofs \d+, (\d+)s\)")
rec = {}
for fn in ['runJ_A.log', 'runJ_B.log', 'runL_A.log', 'runL_B.log']:
    if not os.path.exists(fn): continue
    for line in open(fn):
        m = pat.match(line)
        if m: rec[tuple(sorted(int(x) for x in m.group(1).replace(' ', '').split(',') if x))] = (m.group(3), int(m.group(4)), fn)
for arg in sys.argv[1:]:
    bag = tuple(sorted(map(int, arg.split(','))))
    t0 = time.time()
    L, res = prove_cum.prove_bag(bag, log=lambda x: None)
    dt = time.time() - t0
    sols = [y for y in Z if y >= 10 and '0' not in str(y) and tuple(sorted(int(c) for c in str(y) if c != '1')) == bag]
    viol = 0
    for y in sols:
        s = str(y); Ly = len(s); pos = [Ly - 1 - i for i, ch in enumerate(s) if ch != '1']
        a = [Ly] + pos
        for j, T in res['top'].items():
            if not (Ly - pos[int(j)] < T): viol += 1
        for i, H in res['bot'].items():
            if not (a[int(i)] < H): viol += 1
    st = res['stats']
    print(f"{bag}: L <= {L}; {st['proofs']} proofs, top cases {st['top']}, bottom solves {st['bottom']}, {dt:.0f}s "
          f"| gap chain: {rec.get(bag, ('?', '?', ''))[:2]} | genuine solutions {len(sols)}, fact violations {viol}", flush=True)
