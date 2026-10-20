import random, math, flint
from fractions import Fraction
from solver import TopLattice, matveev_bound, valid_points, solve_bottom
random.seed(7)
for k, tm in ((2, 6), (3, 9), (4, 15), (5, 13)):
    XL, B = matveev_bound(k, tm)
    G = math.ceil(2 * math.log10(XL)) + 12
    lat = TopLattice(XL, tm, G)
    U = int(2.0960 * XL + 3); W = int(1.1834 * XL + 2)
    flint.ctx.prec = 4000
    miss = 0
    for trial in range(40):
        while True:
            t = random.randint(0, tm); L = random.randint(2, XL); up = random.randint(2, min(U, int(2.09*L)))
            w = int(((L*flint.arb(10).log() - up*flint.arb(3).log() - t*flint.arb(2).log()) / flint.arb(7).log()).mid().floor().unique_fmpz())
            if 0 <= w <= W: break
        # choose 1+z = exp(t log2 + u' log3 + w log7 - L log10 - eps) with tiny eps>0, as an exact rational
        eps = flint.arb(random.randint(1, 79)) * flint.arb(10)**(-G)
        val = (t*flint.arb(2).log() + up*flint.arb(3).log() + w*flint.arb(7).log() - L*flint.arb(10).log() - eps).exp()
        # rational approximation of val with denominator 10^(G+60) (perturbation far below eps)
        den = 10**(G + 120)
        num = int((val * den).mid().floor().unique_fmpz())
        opz = Fraction(num, den)
        m = lat.target_m(opz)
        found = valid_points(lat, m, opz)
        if not any(p[:4] == (t, up, w, L) for p in found):
            miss += 1
            print("MISSED planted point", k, (t, up, w, L), found)
        if lat.kills(m): print("QUICK TEST WRONGLY KILLS", k)
    print(f"k={k}: top-side planted points detected {40-miss}/40 (G={G}, XL={XL:.2e})")
    # bottom side
    Umax, Wmax = int(2.0960*XL + 3), int(1.1834*XL + 2)
    bad = 0
    for trial in range(40):
        t = random.randint(0, tm); u = random.randint(0, Umax); w = random.randint(0, Wmax)
        M = G + random.randint(0, G)
        r = pow(2, t, 10**M) * pow(3, u, 10**M) * pow(7, w, 10**M) % 10**M
        if (t, u, w) not in solve_bottom(r, M, tm, Umax, Wmax): bad += 1
    print(f"k={k}: bottom-side planted solutions recovered {40-bad}/40")
