"""Summary of the small-base results (bases 3, 4, 6) from base{b}_cum.json and the cached trees."""
import json, os
import flint
lines = ["# Multiplicative persistence in small bases: complete antecedent trees", "",
         "For a base b and a nonzero final digit d, the tree T_d consists of the digit products whose base-b",
         "digit-product orbit ends at d. Each node row lists EVERY integer y (all sizes) that is a product of",
         "primes below b, has no 0 digit, and has digit product equal to the node; when every node of T_d has",
         "been treated, T_d is complete and every positive integer whose orbit ends at d has persistence at",
         "most height(T_d) + 1.", "",
         "These statements concern nonzero final digits only. As in base 10, numbers whose orbit reaches 0",
         "(where the longest persistences live) are not covered.", ""]
for b, D in ((3, 3000), (4, 3000), (6, 800)):
    fn = f'base{b}_cum.json'
    if not os.path.exists(fn): continue
    res = {int(s): v for s, v in json.load(open(fn)).items()}
    tf = f'trees_b{b}_D{D}.json'
    if os.path.exists(tf):
        nodes = {int(d): {int(s): dp for s, dp in v.items()} for d, v in json.load(open(tf))['nodes'].items()}
    else:
        from basetrees import trees
        nodes = trees(b, D)[2]
    lines += [f"## Base {b}", "", "| final digit | depth | node (base 10) | node (base %d) | largest digit set | proved length bound | all preimages |" % b, "|---|---|---|---|---|---|---|"]
    for d in sorted(nodes):
        for s in sorted(nodes[d], key=lambda x: (nodes[d][x], x)):
            v = res.get(s)
            pre = (', '.join(f"{y} ({flint.fmpz(y).str(b)}_{b})" for y in v['preimages']) or 'none') if v else 'OPEN'
            lines.append(f"| {d} | {nodes[d][s]} | {s} | {flint.fmpz(s).str(b)} | {v['kmax'] if v else '?'} | {v['Lbound'] if v else '-'} | {pre} |")
    lines.append("")
    for d in sorted(nodes):
        ns = nodes[d]
        if all(s in res for s in ns):
            H = max(ns.values())
            lines.append(f"- **Theorem (base {b}, final digit {d}).** The tree has exactly {len(ns)} node(s) and height {H}; "
                         f"every positive integer whose base-{b} digit-product orbit ends at {d} has persistence at most {H + 1}, and this is attained.")
        else:
            lines.append(f"- base {b}, final digit {d}: open nodes {[s for s in ns if s not in res]}")
    lines.append("")
open('smallbases.md', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(l for l in lines if l.startswith('- ')))
