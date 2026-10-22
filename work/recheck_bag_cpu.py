"""Independent re-proof of one bag on the CPU path (cbottom4 -> bottomc3, ctop3 -> topc3), in scratch directories, compared
with the bag's record in bagresC/ (proved by prove_cum_gpu.py): facts, length bound and every step counter must agree.
usage (from work/): python3 recheck_bag_cpu.py BAG [SCRATCH]"""
import sys, os, json, time

def main(bag, scratch):
    import prove_cum
    key = '_'.join(map(str, sorted(bag)))
    rec = json.load(open(os.path.join('bagresC', key + '.json')))
    os.makedirs(os.path.join(scratch, 'ckptC'), exist_ok=True); os.makedirs(os.path.join(scratch, 'bagresC'), exist_ok=True)
    prove_cum.CK = os.path.join(scratch, 'ckptC'); prove_cum.RES = os.path.join(scratch, 'bagresC')
    prove_cum.USE_GAP_RESULTS = False
    for d in (prove_cum.CK, prove_cum.RES):
        for f in os.listdir(d): os.unlink(os.path.join(d, f))
    t0 = time.time()
    L, res = prove_cum.prove_bag(tuple(bag), log=print)
    print(f"CPU re-proof: {time.time() - t0:.0f} s")
    bad = []
    for k in ('Lmax',):
        if res[k] != rec[k]: bad.append((k, res[k], rec[k]))
    for k in ('top', 'bot'):
        a = {int(x): v for x, v in res[k].items()}; b = {int(x): v for x, v in rec[k].items()}
        if a != b: bad.append((k, a, b))
    for k in ('top', 'bottom', 'enum', 'proofs', 'cost', 'resolved', 'resolved_bottom'):
        if res['stats'].get(k) != rec['stats'].get(k): bad.append(('stats.' + k, res['stats'].get(k), rec['stats'].get(k)))
    print("record:", {k: rec['stats'].get(k) for k in ('top', 'bottom', 'enum', 'proofs', 'cost')})
    print("CPU:   ", {k: res['stats'].get(k) for k in ('top', 'bottom', 'enum', 'proofs', 'cost')})
    print("IDENTICAL" if not bad else f"DIFFERENT: {bad}")
    return not bad

if __name__ == '__main__':
    sys.exit(0 if main([int(x) for x in sys.argv[1].split('_')], sys.argv[2] if len(sys.argv) > 2 else '/root/scratch/recheck') else 1)
