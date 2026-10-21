# Parameters for the tree sieve (tsieve.c), prototype. Exact integer arithmetic only.
#
# Idea. For x = 2^a 3^b 7^c with a >= K, the last K digits of x are those of 2^e, e = a + b mu3 + c mu7 mod P_K
# (lemma L1), and exactly 4^K of the P_K = 4*5^(K-1) residues give K digits outside {0,5} (lemma L2). The points of
# the region with a given residue r mod P_k form a coset of the lattice L_k = {(a,b,c): a + b mu3 + c mu7 = 0 mod P_k}.
# tsieve walks the tree of admissible residues down to depth k (4^k leaves, never touching a point), then enumerates
# the points of each leaf's coset inside the region. Within a coset, digits k+1..2k are h0 + 2^k*Y with Y affine in
# the lattice coordinates mod 5^k (lemma L3), so the walk along a lattice row is one addition per point.
#
# usage: tparams.py K OUTFILE     (writes the constants tsieve.c reads)
import sys
from flint import fmpz_mat

def dlog2_mod5(target, m):
    """x mod 4*5^(m-1) with 2^x == target (mod 5^m); 2 is a primitive root mod 5^m."""
    x = next(i for i in range(4) if pow(2, i, 5) == target % 5)
    per = 4
    for i in range(1, m):
        mod = 5 ** (i + 1)
        x = next(x + j * per for j in range(5) if pow(2, x + j * per, mod) == target % mod)
        per *= 5
    assert pow(2, x, 5 ** m) == target % 5 ** m
    return x

K = int(sys.argv[1]); out = sys.argv[2]
assert 2 <= K <= 19                       # the k-digit block h0 + 2^k Y must fit in 64 bits (10^19 < 2^64)
M1, M2 = 5 ** K, 5 ** (2 * K)
P = lambda j: 4 * 5 ** (j - 1)
PK, P2K = P(K), P(2 * K)
mu3, mu7 = dlog2_mod5(3, 2 * K), dlog2_mod5(7, 2 * K)        # level-2K logs; reduced mod P_K they are the level-K logs

# reduced basis of L_K, in the norm that makes the region a/3 + b/2 + c <= N roughly isotropic (weights 2, 3, 6)
Bs = fmpz_mat([[2 * PK, 0, 0], [-2 * (mu3 % PK), 3, 0], [-2 * (mu7 % PK), 0, 6]]).lll()
V = []
for i in range(3):
    r = [int(Bs[i, j]) for j in range(3)]
    assert r[0] % 2 == 0 and r[1] % 3 == 0 and r[2] % 6 == 0
    V.append((r[0] // 2, r[1] // 3, r[2] // 6))
# walk direction V[2]: the shortest basis vector whose step changes J by a unit mod 5 (so Y runs through all residues)
e0 = lambda v: v[0] + v[1] * mu3 + v[2] * mu7
cand = sorted((sum(x * x * wt for x, wt in zip(v, (4, 9, 36))), i) for i, v in enumerate(V) if (e0(v) % P2K) // PK % 5)
assert cand, "no basis vector leaves L_(K+1)"
iw = cand[0][1]; V = [v for i, v in enumerate(V) if i != iw] + [V[iw]]
det = (V[0][0] * (V[1][1] * V[2][2] - V[1][2] * V[2][1]) - V[0][1] * (V[1][0] * V[2][2] - V[1][2] * V[2][0])
       + V[0][2] * (V[1][0] * V[2][1] - V[1][1] * V[2][0]))
assert abs(det) == PK, (det, PK)
e = lambda v: v[0] + v[1] * mu3 + v[2] * mu7
for v in V: assert e(v) % PK == 0
eps = [(e(v) % P2K) // PK for v in V]                        # J(p + v) = J(p) + eps (mod 5^K)
G = pow(2, PK, M2); assert (G - 1) % M1 == 0
w = (G - 1) // M1
weps = [w * x % M1 for x in eps]

with open(out, "w") as f:
    W = lambda name, *vals: f.write(name + " " + " ".join(str(v) for v in vals) + "\n")
    W("K", K); W("M1", M1); W("PK", PK)
    for i, v in enumerate(V): W(f"V{i}", *v)
    W("WEPS", *weps)
    # roots: residues r = 0..3 mod P_1 = 4 (every last digit 2, 4, 8, 6 is admissible); y = 2^(r - 2K) mod 5^(2K)
    for r in range(4):
        y = pow(2, r - 2 * K, M2) if r >= 2 * K else pow(pow(2, 2 * K - r, M2), -1, M2)
        W(f"ROOT{r}", y % M1, y // M1)
    # level j -> j+1 (j = 1..K-1): children r + t P_j, y_t = y * G_j^t mod 5^(2K); digit j+1 test constants
    for j in range(1, K):
        Gj = pow(2, P(j), M2)
        u = ((pow(2, P(j), 5 ** (j + 1)) - 1) // 5 ** j) % 5; assert u != 0
        Cj = pow(2, 2 * K - j - 1, 5 ** (j + 1))
        Gt = [pow(Gj, t, M2) for t in range(5)]
        W(f"LEV{j}", P(j), 5 ** j, Cj, u, *[x for g in Gt for x in (g % M1, g // M1)])
print(f"K={K}: P_K={PK}, basis {V}, eps {eps}")
