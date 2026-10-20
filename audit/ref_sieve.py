#!/usr/bin/env python3
"""Independent reference implementation of the chain sieve's digit tests, for audits.

Plain Python integers, no floating point, nothing imported from work/gsieve. Every decision is computed from
its definition (the decimal digits of x), not from the lemmas the fast program relies on, and can be compared
with the fast program's dump of the same region. It is slow (roughly 10^6 points per second per process) and
meant for spot regions.

Definitions (docs/CHAIN_SIEVE.md). Mode G (chain numbers): an even candidate must avoid the digits 0 and 5;
x = 3^a 7^c (odd, no factor 5) must avoid 0 and must not contain both a 5 and an even digit; x = 3^a 5^b 7^c
(b >= 1) must have only odd digits. Mode E (zero-free, used for validation): avoid 0. Regions: exponent triples
with digit length m(x) <= N (and > N0 for a shell) and x >= 10^130 (smaller x are enumerated exhaustively by
work/gsieve/tiny.py):
    main    x = 2^a 3^b 7^c,  a >= 26,   2a + 3b + 6c <= 6N   (> 6N0)
    direct  x = 2^a 3^b 7^c,  a <  26,   2a + 3b + 6c <= 6N   (> 6N0)
    odd     x = 3^a 5^b 7^c,  b >=  1,   a + 2b + 2c  <= 2N   (> 2N0)
Decisions reproduced (fast-program dump tag in brackets):
    main    S1: last 13 digits allowed;  S2: last 26 digits allowed                          [S2 a b c]
    direct  last 54 digits allowed                                                           [D4 0 a b c]
    odd     G mode, b <= 39: last 54 digits odd;  E mode, every b: last 54 digits nonzero    [D4 1 a b c]
            G mode, b >= 40: last 40 digits odd                                              [T4 1 a b c]
For b >= 40 the fast program does not test points one by one; it visits only the residue class that lemma L4
allows (the "odd thread"). Here the thread is derived independently (`lemmas`), every candidate is checked
digit by digit, and the candidate set is compared with the fast program's.

usage:
    ref_sieve.py lemmas [GCONST_H]
        check every constant of gconst.h from its defining property, and lemmas L1-L4 (L4 by an independent
        derivation of the odd thread)
    ref_sieve.py main|direct|odd|thread G|E N c0 c1 [n0=M] [--exe BIN] [--threads T] [--jobs J] [--full]
        reference decisions on the region (c0 <= c < c1); `thread` is the b >= 40 part of `odd` alone.
        --exe BIN: also run BIN on the same region with `dump` and require equal point counts, survivor sets
        and counters. --full: also decide chain membership (mode G) or zero-freeness (mode E) of every survivor
        from all its digits, and compare with the fast program's MEMBER lines.
"""
import multiprocessing, os, random, re, subprocess, sys
from decimal import Decimal, getcontext

sys.set_int_max_str_digits(0)
T130 = 10 ** 130                         # the fast program handles x >= 10^130 (XMIN = 9 * NL + 4 with NL = 14)

# ------------------------------------------------------------------ digit tests (definitions)
def allowed_even(s, dm):                 # digits of an even candidate
    return dm == 'P' or ('0' not in s and (dm == 'E' or '5' not in s))
def allowed_odd_no5(s, dm):              # digits of 3^b 7^c
    if dm == 'P': return True
    if '0' in s: return False
    return dm == 'E' or not ('5' in s and any(d in s for d in '2468'))
def allowed_odd5(s, dm):                 # digits of 3^a 5^b 7^c, b >= 1
    if dm == 'P': return True
    if '0' in s: return False
    return dm == 'E' or not any(d in s for d in '2468')

def smallest_exp(base, mult, target=T130):
    """smallest e >= 0 with mult * base^e >= target (exact integers)"""
    e, v = 0, mult
    while v < target: v *= base; e += 1
    return e

# ------------------------------------------------------------------ the odd thread, derived independently
def odd_thread(K=40):
    """The unique K-digit residue r with 5^K | r and only odd digits, built digit by digit: if r works for k
    digits, the lifts r + j 10^k divisible by 5^(k+1) are j = j0 and j0 + 5, one odd and one even, so exactly
    one survives. Returns r and (s, t) with r = s * 5^t (mod 2^K), s = +-1, t mod 2^(K-2)."""
    r = 5                                                          # k = 1: the only odd digit divisible by 5
    for k in range(1, K):
        lifts = [j for j in range(10) if (r + j * 10 ** k) % 5 ** (k + 1) == 0 and j % 2 == 1]
        assert len(lifts) == 1, (k, lifts)
        r += lifts[0] * 10 ** k
    s = 1 if r % 4 == 1 else -1
    return r, s, dlog5(s * r, K)

def dlog5(u, K):
    """t mod 2^(K-2) with 5^t = u (mod 2^K), for u = 1 (mod 4); 5 generates that subgroup, of order 2^(K-2)"""
    u %= 2 ** K; assert u % 4 == 1
    t = 0
    for i in range(K - 2):                                         # fix bit i using the modulus 2^(i+3)
        if pow(5, t, 2 ** (i + 3)) != u % 2 ** (i + 3): t += 2 ** i
    assert pow(5, t, 2 ** K) == u
    return t

THREAD_R, THREAD_S, THREAD_T = odd_thread()
F3, F7 = dlog5(-3, 40), dlog5(-7, 40)                              # 3 = -5^F3, 7 = -5^F7 (mod 2^40)

# ------------------------------------------------------------------ reference decisions per c
def main_c(N, N0, c, dm, full):
    pts = s1 = 0; S2 = set(); members = set()
    if 6 * N - 6 * c - 52 < 0: return pts, s1, S2, members, 0
    m26 = 10 ** 26; m7 = 7 ** c; mult = m7
    for b in range((6 * N - 6 * c - 52) // 3 + 1):
        ahi = (6 * N - 3 * b - 6 * c) // 2
        alo = max(26, smallest_exp(2, mult) if mult < T130 else 0)
        if N0 >= 0 and 6 * N0 - 3 * b - 6 * c >= 0: alo = max(alo, (6 * N0 - 3 * b - 6 * c) // 2 + 1)
        if alo <= ahi:
            x = pow(2, alo, m26) * pow(3, b, m26) * pow(7, c, m26) % m26
            for a in range(alo, ahi + 1):
                pts += 1
                s = '%026d' % x                                    # x >= 10^130: these are real digits
                if allowed_even(s[13:], dm):
                    s1 += 1
                    if allowed_even(s, dm):
                        S2.add((a, b, c))
                        if full and allowed_even(str(2 ** a * 3 ** b * m7), dm): members.add((0, a, b, c))
                x = 2 * x % m26
        if mult < T130: mult *= 3
    return pts, s1, S2, members, 0

def direct_c(N, N0, c, dm, full):
    pts = 0; D = set(); members = set(); m54 = 10 ** 54; m7 = 7 ** c
    for a in range(26):
        if 6 * N - 2 * a - 6 * c < 0: continue
        bmax = (6 * N - 2 * a - 6 * c) // 3
        blo = smallest_exp(3, 2 ** a * m7)
        if N0 >= 0 and 6 * N0 - 2 * a - 6 * c >= 0: blo = max(blo, (6 * N0 - 2 * a - 6 * c) // 3 + 1)
        test = allowed_odd_no5 if a == 0 else allowed_even
        x = pow(2, a, m54) * pow(3, blo, m54) * pow(7, c, m54) % m54
        for b in range(blo, bmax + 1):
            pts += 1
            if test('%054d' % x, dm):
                D.add((0, a, b, c))
                if full and test(str(2 ** a * 3 ** b * m7), dm): members.add((0, a, b, c))
            x = 3 * x % m54
    return pts, 0, D, members, 0

def odd_c(N, N0, c, dm, full, part):
    pts = s1 = 0; D = set(); members = set(); m54 = 10 ** 54; m7 = 7 ** c
    if part in ('all', 'direct'):                                  # b <= 39 (G) or every b (E, P), point by point
        for b in range(1, (39 if dm == 'G' else N) + 1):
            amax = 2 * N - 2 * b - 2 * c
            if amax < 0: break
            m57 = 5 ** b * m7
            alo = smallest_exp(3, m57)
            if N0 >= 0 and 2 * N0 - 2 * b - 2 * c >= 0: alo = max(alo, 2 * N0 - 2 * b - 2 * c + 1)
            x = pow(3, alo, m54) * pow(5, b, m54) * pow(7, c, m54) % m54
            for a in range(alo, amax + 1):
                pts += 1
                if allowed_odd5('%054d' % x, dm):
                    D.add((1, a, b, c))
                    if full and allowed_odd5(str(3 ** a * m57), dm): members.add((1, a, b, c))
                x = 3 * x % m54
    if part in ('all', 'thread') and dm == 'G':                   # b >= 40: only the odd thread can survive
        for a in range(0, 2 * N - 2 * c - 80 + 1):                 # some b >= 40 fits: a + 80 + 2c <= 2N
            if (a + c) % 2 == 0: continue                          # x = (-1)^(a+c) 5^(b + a F3 + c F7) (mod 2^40)
            s1 += 1
            b0 = (THREAD_T - a * F3 - c * F7) % 2 ** 38            # the only b mod 2^38 (b <= N < 2^38)
            if b0 < 40 or a + 2 * b0 + 2 * c > 2 * N: continue
            if N0 >= 0 and a + 2 * b0 + 2 * c <= 2 * N0: continue
            if 3 ** a * 5 ** b0 * m7 < T130: continue
            r = 3 ** a * pow(5, b0, 10 ** 40) * m7 % 10 ** 40
            assert r == THREAD_R and allowed_odd5('%040d' % r, 'G'), (a, b0, c)     # direct check of the candidate
            D.add((1, a, b0, c))
            if full and allowed_odd5(str(3 ** a * 5 ** b0 * m7), dm): members.add((1, a, b0, c))
    return pts, s1, D, members, 0

def work(args):
    mode, N, N0, c, dm, full = args
    if mode == 'main': return main_c(N, N0, c, dm, full)
    if mode == 'direct': return direct_c(N, N0, c, dm, full)
    return odd_c(N, N0, c, dm, full, 'thread' if mode == 'thread' else 'all')

# ------------------------------------------------------------------ comparison with the fast program
TAG = {'main': ('S2',), 'direct': ('D4',), 'odd': ('D4', 'T4'), 'thread': ('T4',)}
def run_fast(exe, mode, dm, N, N0, c0, c1, threads):
    cmd = [exe, 'odd' if mode == 'thread' else mode, dm, str(N), str(c0), str(c1), 'dump', f'threads={threads}']
    if N0 >= 0: cmd.append(f'n0={N0}')
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    sets = {t: set() for t in ('S2', 'D4', 'T4')}; members = set(); done = None
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'S2': sets['S2'].add(tuple(map(int, p[1:4])))
        elif p[0] in ('D4', 'T4'): sets[p[0]].add(tuple(map(int, p[1:5])))
        elif p[0] == 'MEMBER':
            d = dict(kv.split('=') for kv in p[1:])
            members.add(tuple(int(d[k]) for k in ('form', 'a', 'b', 'c')))
        elif p[0] == 'DONE': done = dict(kv.split('=', 1) for kv in p[1:])
    return sets, members, done

def region_cmd(argv):
    mode, dm, N, c0, c1 = argv[0], argv[1], int(argv[2]), int(argv[3]), int(argv[4])
    N0, exe, threads, jobs, full = -1, None, 1, 1, False
    rest = argv[5:]
    while rest:
        a = rest.pop(0)
        if a.startswith('n0='): N0 = int(a[3:])
        elif a == '--exe': exe = rest.pop(0)
        elif a == '--threads': threads = int(rest.pop(0))
        elif a == '--jobs': jobs = int(rest.pop(0))
        elif a == '--full': full = True
        else: sys.exit(f'unknown argument {a}')
    tasks = [(mode, N, N0, c, dm, full) for c in range(c0, c1)]
    if jobs > 1:
        with multiprocessing.Pool(jobs) as pool: res = pool.map(work, tasks, chunksize=1)
    else: res = list(map(work, tasks))
    pts = sum(r[0] for r in res); s1 = sum(r[1] for r in res)
    surv = set().union(*(r[2] for r in res)) if res else set()
    mem = set().union(*(r[3] for r in res)) if res else set()
    print(f"reference {mode} {dm} N={N} N0={N0} c=[{c0},{c1}): points={pts} s1={s1} survivors={len(surv)}"
          + (f" members={len(mem)} {sorted(mem)}" if full else ""), flush=True)
    if not exe: return 0
    sets, fmem, done = run_fast(exe, mode, dm, N, N0, c0, c1, threads)
    fast = set().union(*(sets[t] for t in TAG[mode]))
    ok = True
    def check(name, ref, fst):
        nonlocal ok
        good = ref == fst; ok &= good
        if isinstance(ref, set):
            detail = f"(reference {len(ref)}, fast {len(fst)}" + ("" if good else
                     f"; only reference {sorted(ref - fst)[:5]}, only fast {sorted(fst - ref)[:5]}") + ")"
        else: detail = f"(reference {ref}, fast {fst})"
        print(f"  {name}: {'OK' if good else 'MISMATCH'} {detail}", flush=True)
    check('survivor set', surv, fast)
    if mode != 'thread': check('points', pts, int(done['points']))      # odd: the b <= 39 points (G mode)
    if mode == 'main':
        check('s1 = stage-1 survivors', s1, int(done['s1']))
        check('s2 = stage-2 survivors', len(surv), int(done['s2']))
    if mode in ('odd', 'thread') and dm == 'G': check('s1 = thread steps', s1, int(done['s1']))
    if mode in ('direct', 'odd'): check('s4 = survivors passed to GMP', len(surv), int(done['s4']))
    if full: check('members', mem, fmem)
    print(f"  fast program: {exe} variant={done.get('variant')} threads={done.get('threads')}")
    print("REFERENCE AGREES" if ok else "REFERENCE DISAGREES", flush=True)
    return 0 if ok else 1

# ------------------------------------------------------------------ constants and lemmas
def lemmas_cmd(argv):
    path = argv[0] if argv else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'work', 'gsieve', 'gconst.h')
    C = {m.group(1): int(m.group(2)) for m in re.finditer(r'#define (\w+) (\d+)ULL', open(path).read())}
    ok = True
    def check(name, good):
        nonlocal ok
        ok &= bool(good)
        print(f"  {name}: {'OK' if good else 'FAIL'}", flush=True)
    rnd = random.Random(20261020)
    print(f"constants of {path}, each from its defining property:")
    for K in (12, 13):
        p = f'K{K}_'; P1, Q, P2, Q2 = 4 * 5 ** (K - 1), 5 ** K, 4 * 5 ** (2 * K - 1), 5 ** (2 * K)
        check(f'K={K}: P1, Q, P2, Q2, 10^K, 2^K', (C[p + 'P1'], C[p + 'Q'], C[p + 'P2'], C[p + 'Q2'], C[p + 'TEN'], C[p + 'TWO'])
              == (P1, Q, P2, Q2, 10 ** K, 2 ** K))
        check(f'K={K}: 2^MU3 = 3 and 2^MU7 = 7 (mod 5^{2 * K}), 0 <= MU < P2',
              pow(2, C[p + 'MU3'], Q2) == 3 and pow(2, C[p + 'MU7'], Q2) == 7 and C[p + 'MU3'] < P2 and C[p + 'MU7'] < P2)
        g = pow(2, P1, Q2)
        check(f'K={K}: W = (2^P1 mod 5^{2 * K} - 1) / 5^{K} mod 5^{K}', (g - 1) % Q == 0 and C[p + 'W'] == (g - 1) // Q % Q)
        check(f'K={K}: Y0 * 2^{2 * K} = 1 (mod 5^{2 * K})', C[p + 'Y0'] * 2 ** (2 * K) % Q2 == 1)
        check(f'K={K}: X1 = 2^P1 mod 10^{K}', C[p + 'X1'] == pow(2, P1, 10 ** K))
    getcontext().prec = 80
    for pr in (2, 3, 5, 7):
        v = int((Decimal(pr).ln() / Decimal(10).ln() * Decimal(2) ** 64).to_integral_value())
        check(f'L{pr} = round(2^64 log10 {pr})', C[f'L{pr}'] == v)
    check('OF3: 5^OF3 = -3 (mod 2^40), OF3 < 2^38', pow(5, C['OF3'], 2 ** 40) == (-3) % 2 ** 40 and C['OF3'] < 2 ** 38)
    check('OF7: 5^OF7 = -7 (mod 2^40), OF7 < 2^38', pow(5, C['OF7'], 2 ** 40) == (-7) % 2 ** 40 and C['OF7'] < 2 ** 38)
    print("lemma L4 (odd thread), derived independently by digit lifting:")
    check(f'unique 40-digit all-odd residue divisible by 5^40: ...{THREAD_R}', len(str(THREAD_R)) == 40
          and THREAD_R % 5 ** 40 == 0 and allowed_odd5(str(THREAD_R), 'G')
          and str(THREAD_R).endswith('9313975971735999337397515773773193359375'))
    check(f'r = s 5^t (mod 2^40) with s = {THREAD_S}, t = {THREAD_T}; gconst.h OTHREAD = {C["OTHREAD"]}',
          THREAD_S == -1 and THREAD_T == C['OTHREAD'] and pow(5, THREAD_T, 2 ** 40) == (-THREAD_R) % 2 ** 40)
    check('F3, F7 equal OF3, OF7', (F3, F7) == (C['OF3'], C['OF7']))
    hit = miss = 0
    for _ in range(300):                                           # planted thread points, and nearby non-thread points
        a = rnd.randrange(0, 3000); c = rnd.randrange(0, 3000); c += (a + c + 1) % 2
        b = (THREAD_T - a * F3 - c * F7) % 2 ** 38 + 2 ** 38 * rnd.randrange(0, 3)
        b = b if b >= 40 else b + 2 ** 38
        hit += allowed_odd5('%040d' % (3 ** a * pow(5, b, 10 ** 40) * 7 ** c % 10 ** 40), 'G')
        miss += allowed_odd5('%040d' % (3 ** a * pow(5, b + 2 ** rnd.randrange(0, 38), 10 ** 40) * 7 ** c % 10 ** 40), 'G')
    check(f'planted thread points have 40 odd digits ({hit}/300); shifted exponents never do ({miss}/300)', hit == 300 and miss == 0)
    print("lemma L1 (suffix collapse) and L3 (affine lift), random tests against direct computation:")
    for K in (12, 13):
        p = f'K{K}_'; P1, Q, Q2 = 4 * 5 ** (K - 1), 5 ** K, 5 ** (2 * K)
        bad = 0
        for _ in range(1000):                                      # L1 for K and 2K digits
            a = rnd.randrange(2 * K, 2000); b = rnd.randrange(0, 2000); c = rnd.randrange(0, 2000)
            for k, mu3, mu7 in ((K, C[p + 'MU3'] % P1, C[p + 'MU7'] % P1), (2 * K, C[p + 'MU3'], C[p + 'MU7'])):
                per = 4 * 5 ** (k - 1); e = k + (a + b * mu3 + c * mu7 - k) % per
                bad += (2 ** a * 3 ** b * 7 ** c) % 10 ** k != pow(2, e, 10 ** k)
        check(f'K={K}: L1, 2^a 3^b 7^c = 2^e (mod 10^K and 10^2K), 1000 random (a, b, c)', bad == 0)
        bad = 0
        for _ in range(1000):                                      # L3: digits K+1..2K of 2^(e + P1 j), e >= 2K
            e = rnd.randrange(2 * K, P1 + 2 * K); j = rnd.randrange(0, Q)
            y = pow(2, e - 2 * K, Q2); alo, ahi = y % Q, y // Q
            hi = (2 ** K * alo // Q) + 2 ** K * ((ahi + alo * C[p + 'W'] * j) % Q)
            bad += pow(2, e + P1 * j, 10 ** (2 * K)) // 10 ** K != hi
        check(f'K={K}: L3, digits {K + 1}-{2 * K} = h0 + 2^K ((A_hi + A_lo W j) mod 5^K), 1000 random (e, j)', bad == 0)
    print("lemma L2 (exact branching), exhaustive for k <= 9:")
    for k in range(1, 10):
        per = 4 * 5 ** (k - 1); x = pow(2, k * per, 10 ** k); n = 0     # 2^E for E = k*per, k*per + 1, ...: E >= k
        for _ in range(per):
            n += allowed_even('%0*d' % (k, x), 'G'); x = 2 * x % 10 ** k
        check(f'k={k}: {n} admissible endings of 2^E (E >= k) among {per} residues; 4^k = {4 ** k}', n == 4 ** k)
    print("ALL CONSTANTS AND LEMMAS OK" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1

if __name__ == '__main__':
    if len(sys.argv) < 2: sys.exit(__doc__)
    if sys.argv[1] == 'lemmas': sys.exit(lemmas_cmd(sys.argv[2:]))
    if sys.argv[1] in ('main', 'direct', 'odd', 'thread') and len(sys.argv) >= 7: sys.exit(region_cmd(sys.argv[1:]))
    sys.exit(__doc__)
