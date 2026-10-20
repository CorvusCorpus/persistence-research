"""Assemble the proof certificate: for each proved node, every digit multiset with its proved facts
(gap-chain bounds from proofJ_*.json, or cumulative top/bottom facts from proofC_*.json) and the
resulting length bound, plus the complete list of 7-smooth preimages (from the exhaustive list
below 10^5000)."""
import json, glob, re
from bags import bags_of

def load():
    res = {}
    for f in sorted(glob.glob('proofJ_*.json')):
        for k, v in json.load(open(f)).items():
            res[int(k)] = dict(v, method='gap chain (prove_iter2/3)')
    for f in sorted(glob.glob('proofC_*.json')):
        for k, v in json.load(open(f)).items():
            if int(k) not in res: res[int(k)] = dict(v, method='cumulative (prove_cum)')
    return res

if __name__ == "__main__":
    res = load()
    nodes = json.load(open('even_nodes.json'))
    out = {'method': 'gap chain: case split on runs of >=G ones; cumulative: bounds on depths D_j = L - a_(j+1) (top, '
                     'LLL lattice tests) and heights a_i (bottom, exact 10-adic logarithms); Matveev bound; 2-adic bound via '
                     'digit-placement recursion; exhaustive list of zero-free 7-smooth numbers below 10^5000',
           'nodes': {}}
    for d in ('2', '4', '6', '8'):
        for s in nodes[d]:
            if s in res:
                out['nodes'][s] = {'target': int(d), **res[s]}
    json.dump(out, open('certificate.json', 'w'), indent=1)
    proved = {s: v['preimages'] for s, v in out['nodes'].items()}
    print(f"{len(proved)} nodes with complete preimage sets")
    for d in ('2', '4', '6', '8'):
        level = [int(d)]; depth = 0
        while level and all(x in proved for x in level):
            level = sorted(set(y for x in level for y in proved[x])); depth += 1
        print(f"target {d}: levels 0..{depth} exactly determined; blocking nodes: "
              f"{[x for x in level if x not in proved]}")
