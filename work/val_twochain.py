"""Checks of the two-chain bootstrapping bound prove_cum.twochain_cum (PC_BOOTSTRAP=twochain, docs/BOOTSTRAP_TWO_CHAINS.md).

C  constants: YU_C2, YU_K, YU_HMIN, C4, PI, PI3 re-evaluated with mpmath at 60 digits from the formulas of Yu (2013) and
   Matveev (2000); prove_cum's floats must not be below them by more than 1e-12 (relative).
R  regression: with the default bootstrap prove_cum.matveev_cum equals a frozen copy of the original function on random
   fact states (records written so far stay reproducible).
D  dominance: twochain_cum <= matveev_cum on random fact states and on every state of the dry-run plans of BAGS.
G  genuine numbers: for every zero-free 7-smooth y in Z.json without digit 5 and with k >= 1 non-1 digits, exactly:
   F_i is a 5-adic unit and 5^(a_i) | 9y - F_i (a_0 = L); ord_5(9y - F_i) is below the Yu bound with the exact height of
   2^t/F_i and B = max(u', w, 3), and that height is below the chain's estimate from a_(i+1); 0 < Lam_j < 80*10^(-D_j);
   N_j < 9*10^(D_(j-1)); D_j is below the Matveev bound (4 logarithms) with the exact heights; a_k <= 1; and twochain_cum
   with every partial set of exact facts (top levels < J, bottom levels >= I) is >= L (k <= 20; tm = the number's own t).
usage (from work/): python3 val_twochain.py [BAG ...]   (BAG as 2_2_3_7; default: three bags of proved nodes)"""
import sys, os, json, math, random
import mpmath as mp
import prove_cum

mp.mp.dps = 60
LN = math.log
bad = 0
def fail(msg):
    global bad
    bad += 1
    if bad <= 20: print("FAIL", msg)

# ---- C: constants
n, p = 3, 5
c2 = mp.mpf(557) * (p - 1) / (p - 2)
C2 = c2 * mp.mpf(7) ** n * mp.mpf(n + 1) ** (n + 1) / (mp.factorial(n - 1) * 2 * mp.log(p) ** 3) \
    * max(mp.mpf(p), (mp.e / n) ** n * mp.log(p) ** (n + 1)) * max(mp.log(mp.e ** 4 * (n + 1)), 1, mp.log(p))
YK = (n + 1) * C2 * mp.log(3) * mp.log(7)
HMIN = max(n / mp.log(p), 1) / 25
C4 = mp.mpf('1.4') * mp.mpf(30) ** 7 * mp.mpf(4) ** mp.mpf('4.5')
PI = mp.log(2) * mp.log(3) * mp.log(7) * mp.log(10)
PI3 = mp.log(3) * mp.log(7) * mp.log(10)
for name, flt, exact in (('YU_C2', prove_cum.YU_C2, C2), ('YU_K', prove_cum.YU_K, YK), ('C4', prove_cum.C4, C4),
                         ('PI', prove_cum.PI, PI), ('PI3', prove_cum.PI3, PI3)):
    rel = (mp.mpf(flt) - exact) / exact
    print(f"C {name:6s} = {mp.nstr(exact, 15):>22s}  float {flt:.15e}  rel {mp.nstr(rel, 3)}")
    if rel < -1e-12: fail(f"constant {name} below its exact value")
if abs(mp.mpf(prove_cum.YU_HMIN) - HMIN) > 1e-15: fail("YU_HMIN")
print(f"C YU_HMIN = {prove_cum.YU_HMIN:.6f};  C2 = {mp.nstr(C2, 6)} (Theorem 2), per-level factor at log B = 170: "
      f"{mp.nstr(YK * 170 * mp.log(10), 4)}")

# ---- R: regression of the default bound against the original code
def matveev_orig(k, tm, top, bot):
    INF = float('inf'); C5 = prove_cum.C5
    def f(B):
        beta = 1 + LN(B)
        Lb = bot[0] - 1 if 0 in bot else INF
        s = 0.0
        for j in range(k):
            m = (LN(80) + (prove_cum.C4 * beta * prove_cum.PI if j == 0 else C5 * beta * prove_cum.PI * (s * LN(10) + LN(9)))) / LN(10)
            if j in top: m = min(m, top[j] - 1)
            s = m
            if (j + 1) in bot: Lb = min(Lb, s + bot[j + 1] - 1)
        return max(2.0960 * Lb + 2.2, tm + 1), Lb
    B = 1e300
    for _ in range(300):
        nb, _ = f(B)
        if nb >= B * (1 - 1e-9): break
        B = nb
    B = max(B * 1.01, 10.0)
    nb, Lb = f(B)
    assert nb <= B
    return max(10, math.ceil(Lb * 1.01))

rnd = random.Random(1)
def random_state():
    k = rnd.randint(2, 13); tm = rnd.randint(0, 60)
    J = rnd.randint(1, k); I = rnd.randint(0, k)
    top = {j: rnd.randint(20, 600) for j in range(J)}
    bot = {i: rnd.randint(20, 600) for i in range(I, k)}; bot[k] = 2
    return k, tm, top, bot
nR = nD = 0
for _ in range(3000):
    k, tm, top, bot = random_state()
    try: a = matveev_orig(k, tm, top, bot)
    except AssertionError: continue
    prove_cum.BOOTSTRAP = 'matveev'; b = prove_cum.matveev_cum(k, tm, top, bot)
    prove_cum.BOOTSTRAP = 'twochain'; c = prove_cum.matveev_cum(k, tm, top, bot)
    nR += 1; nD += 1
    if a != b: fail(f"regression {k} {tm} {top} {bot}: {a} != {b}")
    if c > b: fail(f"dominance {k} {tm} {top} {bot}: twochain {c} > matveev {b}")
prove_cum.BOOTSTRAP = 'matveev'
print(f"R default bound == original code on {nR} random states;  D twochain <= matveev on them")

bags = [tuple(sorted(int(x) for x in a.split('_'))) for a in sys.argv[1:]] or \
       [(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 7), (2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 7), (3, 3, 3, 3, 4, 4, 4, 6, 8)]
for bag in bags:
    snaps = []; orig = prove_cum.matveev_cum
    def spy(k_, tm_, top, bot):
        snaps.append((k_, tm_, dict(top), dict(bot))); return orig(k_, tm_, top, bot)
    prove_cum.matveev_cum = spy
    try: prove_cum.prove_bag(bag, log=lambda x: None, dry=True)
    finally: prove_cum.matveev_cum = orig
    worse = 0
    for k, tm, top, bot in snaps:
        if prove_cum.twochain_cum(k, tm, top, bot) > matveev_orig(k, tm, top, bot): worse += 1
    if worse: fail(f"dominance on the plan of {bag}: {worse} states")
    print(f"D plan of {bag}: twochain <= matveev on all {len(snaps)} states")

# ---- G: genuine numbers
def v(x, q):
    c = 0
    while x % q == 0: x //= q; c += 1
    return c
def h_rat(a, b):
    g = math.gcd(a, b); return max(math.log(abs(a) // g), math.log(abs(b) // g))
Z = json.load(open('Z.json'))
KMAX = 20                        # bag sizes of the prover
ng = nchk = 0
for y in Z:
    s = str(y)
    if y < 10 or '0' in s or '5' in s or set(s) == {'1'} or len(s) - s.count('1') > KMAX: continue
    t = v(y, 2); r = y >> t; u = v(r, 3); r //= 3 ** u; w = v(r, 7); r //= 7 ** w
    assert r == 1, y
    up = u + 2; L = len(s)
    dig = [(L - 1 - i, int(ch)) for i, ch in enumerate(s) if ch != '1']          # (a_m, d_m), a_1 > ... > a_k
    k = len(dig); a = [L] + [x for x, _ in dig]; c = [None] + [9 * (d - 1) for _, d in dig]
    tm = t; ng += 1               # the bounds increase with tm; the number's own t is the strictest choice (and tmax(bag)
    if a[k] > 1: fail(f"{y}: a_k = {a[k]} > 1")   # of large bags needs gigabytes)
    B = max(up, w, 3); lB = max(math.log(B), math.log(5))
    for i in range(k):
        F = -1 + sum(c[m] * 10 ** a[m] for m in range(i + 1, k + 1))
        if F % 5 == 0: fail(f"{y}: F_{i} divisible by 5")
        X = 9 * y - F
        o5 = v(X, 5)
        if o5 < a[i]: fail(f"{y}: ord_5(9y - F_{i}) = {o5} < a_{i} = {a[i]}")
        hF = h_rat(2 ** t, F)
        if not hF <= max(math.log(80) + a[i + 1] * math.log(10), tm * math.log(2)) + math.log(2):
            fail(f"{y}: height of 2^t/F_{i} above the chain estimate")
        if not o5 < prove_cum.YU_K * max(hF, prove_cum.YU_HMIN) * lB: fail(f"{y}: Yu bound violated at i = {i}")
        nchk += 1
    lg = {q: mp.log(q) for q in (2, 3, 7, 10)}
    for j in range(k):
        Nfull = 10 ** L + sum(c[m] * 10 ** a[m] for m in range(1, j + 1))       # 10^L (1 + z_j)
        D = L - a[j + 1]
        lam = t * lg[2] + up * lg[3] + w * lg[7] - L * lg[10] - (mp.log(Nfull) - L * lg[10])
        if not (0 < lam < 80 * mp.mpf(10) ** (-D)): fail(f"{y}: Lam_{j} outside (0, 80*10^-{D})")
        Bm = max(t, up, w, L, 3); beta = 1 + mp.log(Bm)
        if j == 0: mb = (mp.log(80) + C4 * beta * PI) / mp.log(10)
        else:
            Nj = Nfull // 10 ** a[j]; assert Nj * 10 ** a[j] == Nfull
            if not Nj < 9 * 10 ** (L - a[j]): fail(f"{y}: N_{j} >= 9*10^D_(j-1)")
            A = max(mp.mpf(h_rat(Nj, 2 ** t)), abs(mp.log(Nj) - t * lg[2]), mp.mpf('0.16'))
            mb = (mp.log(80) + C4 * beta * PI3 * A) / mp.log(10)
        if not D < mb: fail(f"{y}: Matveev bound violated at j = {j}")
    for J in range(k + 1):
        for I in range(J, k + 1):
            top = {j: L - a[j + 1] + 1 for j in range(J)}
            bot = {i: a[i] + 1 for i in range(I, k)}; bot[k] = 2
            if prove_cum.twochain_cum(k, tm, top, bot) < L: fail(f"{y}: twochain bound below L with facts J={J} I={I}")
print(f"G {ng} genuine numbers, {nchk} bottom-chain inequalities and every top inequality checked")
print("ALL OK" if not bad else f"{bad} FAILURES")
sys.exit(1 if bad else 0)
