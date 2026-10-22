"""Reference check of the top-filter mathematics (docs/GPU_EVEN_TARGETS.md section 3) on a csrc/topc3 job.

For every case of the job (same enumeration as topc3): m = round(-C log(1+z)) as topc3 computes it, topc3's exact level-1 and
level-2 decisions, and the filter's quantities psi1 = frac(-Theta lambda), psi2 = frac(-Theta2 lambda) with
Theta = C qa/qb, Theta2 = C a2/D, lambda = log(1+z) (arb, high precision). Checks:
  (1) |dist(m theta, Z) - dist(psi1, Z)| <= 0.6 |theta|  (m enters only through a negligible term);
  (2) the filter's kill rule, with margin eps, never kills a case that topc3 keeps (level 1 and level 2);
and reports how many cases the filter passes to the exact path. usage (from work/): python3 top_filter_ref.py JOB [MAXCASES]"""
import sys
from fractions import Fraction
import flint

def read_job(path):
    tok = open(path).read().split(); i = 0; kv = {}
    def nxt():
        nonlocal i; i += 1; return tok[i - 1]
    for name in ('prec', 'C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd', 'j', 'cum'):
        assert nxt() == name; kv[name] = int(nxt())
    assert nxt() == 'ranges'; kv['ranges'] = [int(nxt()) for _ in range(kv['j'])]
    assert nxt() == 'nord'; nord = int(nxt())
    kv['orders'] = [[int(nxt()) for _ in range(kv['j'])] for _ in range(nord)]
    for name in ('skip', 'maxsurv', 'debug', 'g0lo', 'g0hi'):
        assert nxt() == name; kv[name] = int(nxt())
    return kv

def cases(kv):
    """(ordering, gaps) in topc3's order (cum mode)"""
    j, rng = kv['j'], kv['ranges']
    for o in range(kv['skip'], len(kv['orders'])):
        def rec(i, acc, g):
            if i == j - 1:
                ghi = rng[j - 1] - 1 - acc; glo = 1
                if j == 1: glo = max(glo, kv['g0lo']); ghi = min(ghi, kv['g0hi'])
                for gl in range(glo, ghi + 1): yield o, g + [gl]
                return
            lo = max(1, kv['g0lo']) if i == 0 else 1; hi = rng[i] - 1 - acc
            if i == 0: hi = min(hi, kv['g0hi'])
            for gi in range(lo, hi + 1): yield from rec(i + 1, acc + gi, g + [gi])
        yield from rec(0, 0, [])

def main(path, maxc=20000):
    kv = read_job(path)
    flint.ctx.prec = kv['prec'] + 128
    C, qa, qb, D, a2, P, Q = (kv[k] for k in ('C', 'qa', 'qb', 'D', 'a2', 'P', 'Q'))
    theta = Fraction(qa, qb); theta2 = Fraction(a2, D); mu = Fraction(P, Q)
    g2 = Fraction(kv['g2n'], kv['g2d']); g3 = Fraction(kv['g3n'], kv['g3d']); R2 = Fraction(kv['Rn'], kv['Rd'])
    q = R2 / g3; rho = (flint.arb(q.numerator) / q.denominator).sqrt()
    Th = flint.arb(C) * flint.arb(qa) / qb; Th2 = flint.arb(C) * flint.arb(a2) / D
    eta = 0.6 * abs(float(theta))
    eps = 2.0 ** -55                                         # the filter's error budget (GPU target), plus eta
    n = n1 = n2 = pass1 = pass2 = 0; worst = 0.0
    for o, g in cases(kv):
        n += 1
        if n > maxc: break
        cs = kv['orders'][o]; s = sum(g); N = 10**s; acc = 0
        for c, gi in zip(cs, g): acc += gi; N += c * 10**(s - acc)
        lam = flint.arb(N).log() - s * flint.arb(10).log()  # log(1+z)
        v = -flint.arb(C) * lam
        m = int(v.mid().floor().unique_fmpz()); m += 1 if (v - m) > 0.5 else 0
        # exact level 1 (topc3)
        r = (-m * qa) % qb; d = min(r, qb - r); keep1 = not (d * d * kv['K1'] > kv['K2'])
        # filter level 1
        p1 = -Th * lam; f1 = float((p1 - p1.mid().floor().unique_fmpz()).mid())   # in [0, 1)
        dist1 = min(f1, 1 - f1)
        worst = max(worst, abs(float(Fraction(d, qb)) - dist1))
        kill1 = dist1 > float(rho.upper()) + eps + eta
        if kill1 and keep1: raise SystemExit(f"UNSOUND level 1 at case {n}: o={o} g={g}")
        if keep1: n1 += 1
        if not kill1: pass1 += 1
        if not keep1: continue
        # exact level 2 (topc3's loop over x3), and the filter's level 2
        a3 = kv['a3']; base3 = m * a3; fl = (-base3) // D
        R3a = kv['g3n'] * kv['Rd']; R3b = kv['Rn'] * kv['g3d'] * D * D
        L2a = kv['g2n'] * kv['g3d'] * kv['Rd']; L2b = Q * Q * kv['g3n'] * kv['g2d'] * kv['Rd']; L2c = kv['Rn'] * kv['g2d'] * kv['g3d'] * Q * Q * D * D
        Kl = int((flint.arb(kv['Rn'] * kv['g3d']) / (kv['Rd'] * kv['g3n'])).sqrt().floor().unique_fmpz()) + 1
        surv = False
        for dx in range(-Kl, Kl + 2):
            x3 = fl + dx; E3 = x3 * D + base3
            if E3 * E3 * R3a > R3b: continue
            N2 = -m * a2 * Q - P * E3; r2 = N2 % (Q * D); d2 = min(r2, Q * D - r2)
            if d2 * d2 * L2a + E3 * E3 * L2b <= L2c: surv = True; break
        if surv: n2 += 1
        sgn = f1 if f1 < 0.5 else f1 - 1                    # signed distance of psi1 to Z (= x3 + m theta for the nearest x3)
        p2 = -Th2 * lam; f2 = float((p2 - p2.mid().floor().unique_fmpz()).mid())
        y = -f2 - float(mu) * sgn; y -= round(y); e2 = (eps + eta) * (1 + abs(float(mu)))
        lhs = max(0.0, abs(y) - e2) ** 2 * float(g2) + max(0.0, abs(sgn) - eps - eta) ** 2 * float(g3)
        kill2 = lhs > float(R2)
        if kill2 and surv: raise SystemExit(f"UNSOUND level 2 at case {n}: o={o} g={g}")
        if not kill2: pass2 += 1
    n = min(n, maxc)
    print(f"{path}: {n} cases; theta = 2^{float(abs(theta)).hex()} ({abs(float(theta)):.3g}); worst |dist(m theta) - dist(psi1)| = {worst:.3g}")
    print(f"  topc3 keeps {n1} at level 1, {n2} at level 2; the filter passes {pass1} (level 1), {pass2} (level 2) to the exact path; no unsound kill")

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20000)
