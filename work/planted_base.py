import random, math, flint
from fractions import Fraction
from basesolver import Base, Top2, solve_c, matveev_partial, lam_ok
random.seed(5)
for b in (3, 4, 6):
    base = Base(b); q = base.q
    for k in (1, 3, 5):
        XL, _ = matveev_partial(k, base, 6, {})
        G = math.ceil((math.log(XL * math.log(b) / math.log(q) + 2) + 7 * math.log(10)) / math.log(b)) + 2
        T = Top2(base, XL, G)
        hit = 0; tot = 0
        for _ in range(25):
            B = random.choice([1, 2, 3]) if b == 6 else (random.choice([1, 2]) if b == 4 else 1)
            # pick L, then c so that log((b-1)B q^c / b^L) lies in a sane range; 1+z = that ratio * exp(-eps)
            flint.ctx.prec = T.prec + 200
            L = random.randint(2, XL)
            c = int(((L * flint.arb(b).log() - flint.arb((b - 1) * B).log()) / flint.arb(q).log()).mid().floor().unique_fmpz())
            if not (0 <= c <= T.Xc): continue
            eps = flint.arb(random.randint(1, base.eps0 - 1)) * flint.arb(b) ** (-G) / 2
            val = (flint.arb((b - 1) * B).log() + c * flint.arb(q).log() - L * flint.arb(b).log() - eps).exp()
            den = b ** (G + 80)
            opz = Fraction(int((val * den).mid().floor().unique_fmpz()), den)
            flint.ctx.prec = T.prec
            kappa = flint.arb((b - 1) * B).log() - (flint.arb(opz.numerator) / opz.denominator).log()
            surv = T.survivors(kappa, lambda cc, LL: lam_ok(base, B, cc, LL, opz, G))
            tot += 1; hit += (c, L) in surv
        okb = 0; totb = 0
        for _ in range(25):
            B = random.choice([1, 2, 3]) if b == 6 else (random.choice([1, 2]) if b == 4 else 1)
            c = random.randint(0, T.Xc); M = G + random.randint(0, 40)
            r = B * pow(q, c, b ** M) % b ** M
            totb += 1; okb += c in solve_c(base, B, r, M, T.Xc)
        print(f"base {b} k={k}: X_L={XL:.1e} G={G}: top planted {hit}/{tot}, bottom planted {okb}/{totb}")
