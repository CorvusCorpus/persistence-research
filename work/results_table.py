"""Readable summary of all proved nodes (proofJ_*.json gap chain, proofC_*.json cumulative) plus tree-level consequences."""
import json, glob
import verify as V
from bags import bags_of

res = {}
for f in glob.glob('proofJ_*.json') + glob.glob('proofC_*.json'):
    for k, v in json.load(open(f)).items():
        res[int(k)] = v
nodes = json.load(open('even_nodes.json'))
lines = ["# Complete 7-smooth preimage sets for even-target nodes", "",
         "Each row: every integer y = 2^a 3^b 7^c whose digit product equals the node, over ALL sizes.",
         "", "| target | depth | node | max non-1 digits | digit sets | proved length bound | all 7-smooth preimages |",
         "|---|---|---|---|---|---|---|"]
def depth(x):
    k = 0
    while x > 9: x = V.P(x); k += 1
    return k
for d in ('2', '4', '6', '8'):
    for s in sorted(nodes[d], key=lambda x: (depth(x), x)):
        if s not in res: continue
        bs = bags_of(s)
        lines.append(f"| {d} | {depth(s)} | {s} | {max(len(b) for b in bs)} | {len(bs)} | "
                     f"{res[s]['Lbound']} | {', '.join(map(str, res[s]['preimages'])) or 'none'} |")
lines.append("")
proved = {s: v['preimages'] for s, v in res.items()}
lines.append("## Tree levels determined exactly")
for d in ('2', '4', '6', '8'):
    level = [int(d)]; dep = 0
    while level and all(x in proved for x in level):
        level = sorted(set(y for x in level for y in proved[x])); dep += 1
    lines.append(f"- target {d}: every node at depth <= {dep} is known; open nodes at depth {dep}: "
                 f"{[x for x in level if x not in proved]}")
open('results.md', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines[-6:]))
print(len(res), "nodes in table")
