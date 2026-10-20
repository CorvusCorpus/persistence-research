"""bottomc2 (all candidates reported) vs bottomc on proof-scale spaces: identical solve counts and verdicts."""
import math
from prove_iter import ordered
from twoadic import tmax
import cbottom2 as old
import cbottom3 as new
XL = 86191814072450236416; G = math.ceil(2 * math.log10(XL)) + 12
U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
import os
for bag, j, bot in [((2, 3, 7), 1, {2: 52, 3: 2}), ((2, 3, 7, 9), 1, {2: 52, 3: 52, 4: 2}), ((2, 4, 7, 8, 9), 2, {3: 60, 4: 52, 5: 2})]:
    k = len(bag); tm = tmax(bag); orders = ordered(bag, k - j)
    s1 = {'bottom': 0}; s2 = {'bottom': 0}
    o1 = old.prove_bottom_cum(j, bag, G, tm, U, W, bot, s1, orders)
    o2 = new.prove_bottom_cum(j, bag, G, tm, U, W, bot, s2, orders)
    print(bag, j, "bottomc:", o1, s1['bottom'], " bottomc2:", o2, s2['bottom'], s2.get('resolved_bottom', 0))
