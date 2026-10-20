# Exhaustive check of all small 7-smooth numbers (below 10^XMAX) in the two forms: which lie in G?
import json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
XMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 131
LIM = 10 ** XMAX
def inG(s):
    if '0' in s: return False
    return not ('5' in s and any(ch in s for ch in '2468'))
members = set(); count = 0
# form (i): 2^a 3^b 7^c
x7 = 1
while x7 < LIM:
    x3 = x7
    while x3 < LIM:
        x = x3
        while x < LIM:
            count += 1
            if x >= 10 and inG(str(x)): members.add(x)
            x *= 2
        x3 *= 3
    x7 *= 7
# form (iii): 3^a 5^b 7^c, b >= 1
x7 = 1
while x7 < LIM:
    x5 = x7 * 5
    while x5 < LIM:
        x = x5
        while x < LIM:
            count += 1
            if x >= 10 and inG(str(x)): members.add(x)
            x *= 3
        x5 *= 5
    x7 *= 7
Z = json.load(open(os.path.join(HERE, '..', 'Z.json')))
known = {y for y in Z if y >= 10 and y < LIM and inG(str(y))}
print(f"numbers checked below 10^{XMAX}: {count}; G-members: {len(members)}; largest has {len(str(max(members)))} digits")
print("matches the known list:", members == known, f"(known {len(known)})")
assert members == known
json.dump(sorted(members), open(os.path.join(HERE, 'G_small.json'), 'w'))
