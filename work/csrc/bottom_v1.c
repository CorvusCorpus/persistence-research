/* Bottom-side case loop of the base-10 prover: port of fastbottom.BottomContext.solutions and
   prove_iter2.prove_bottom.  Decides, for every bottom configuration with a run of >= G ones above it,
   whether 2^t 3^u' 7^w = F (mod 10^M) has a solution with t <= tm, 2 <= u' <= U, 0 <= w <= W.
   Exact integer arithmetic throughout (GMP).  Mathematics as in fastbottom.py:
     2-adic: t = v2(F); sigma = sign class of F/2^t mod 4; u' L2_3 + w L2_7 = log(sigma F/2^t) mod 2^n
     5-adic: 3u'+w = dlog_2(F 2^-t mod 5) mod 4;  u' L5_3 + w L5_7 = log(F^4) - t log 16 mod 5^M
   (L2_3 = log(-3), L2_7 = log(-7) 2-adically; L5_x = log(x^4) 5-adically), giving u' = alpha w + beta
   mod Q = 2^(n-2) 5^(M-1) and parity/mod-4 conditions; each class w = w0 mod 4 is a coset of a 2-dim
   lattice whose shortest vector exceeds the box diameter, so the box holds at most one point; with a
   Gauss-reduced basis that point has coordinates within 0.79 / 0.58 of those of the box centre.
   Precision: all p-adic data are carried modulo 5^Mmax, 2^(Mmax+2) with Mmax >= every M used (checked). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

static char buf[1 << 20];
static void rdz(FILE *f, const char *name, mpz_t x) {
    char key[64];
    if (fscanf(f, "%63s %s", key, buf) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    if (mpz_set_str(x, buf, 10)) { fprintf(stderr, "bad int %s\n", name); exit(2); }
}
static long rdl(FILE *f, const char *name) {
    char key[64]; long v;
    if (fscanf(f, "%63s %ld", key, &v) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s (got %s)\n", name, key); exit(2); }
    return v;
}
static void die(const char *m) { fprintf(stderr, "%s\n", m); exit(3); }

static long Mmax, tm;
static mpz_t U, W, diam2, L5_3, L5_7, L5_16, L2_3, L2_7;
static mpz_t *pow5, *pow2, *pow10;          /* pow5[i] = 5^i (i <= Mmax+12), pow2[i] = 2^i (i <= Mmax+44), pow10[i] */
static mpz_t T1, T2, T3, T4, T5, T6;
#define NINV 2000
static mpz_t inv5n[NINV], inv2n[NINV];   /* inverse of (n with p-part removed) mod 5^(Mmax+12) / 2^(Mmax+44) */

static long v2z(const mpz_t x) { if (mpz_sgn(x) == 0) return 1000000000L; return (long)mpz_scan1(x, 0); }
static long vpl(long n, long p) { long e = 0; while (n % p == 0) { n /= p; e++; } return e; }
static long floorlog(long n, long p) { long lg = 0, q = p; while (q <= n) { q *= p; lg++; } return lg; }
static mpz_ptr P(unsigned long p, long i) {
    if (p == 5) { if (i > Mmax + 12) die("pow5 range"); return pow5[i]; }
    if (i > Mmax + 44) die("pow2 range");
    return pow2[i];
}

/* log(1+x) mod p^prec, x integer with v_p(x) >= vx >= 1 (>= 2 if p == 2); exact port of _log1p_series */
static void log1p_series(mpz_t out, const mpz_t x, unsigned long p, long prec, long vx) {
    if (vx < (p > 2 ? 1 : 2)) die("log1p_series: vx too small");
    if (!mpz_divisible_p(x, P(p, vx))) die("log1p_series: valuation");
    mpz_t xp, term, nninv; mpz_inits(xp, term, nninv, NULL);
    mpz_set_ui(out, 0);
    /* xp = x^n mod p^(prec + emax) maintained incrementally; emax = floorlog(nmax) <= 12 */
    long ext = 12;
    mpz_ptr modx = P(p, prec + ext);
    mpz_set_ui(xp, 1);
    for (long n = 1;; n++) {
        long lg = floorlog(n, p);
        if (n * vx - lg >= prec + 2 && n > 1) break;
        long e = vpl(n, p); long nn = n; for (long i = 0; i < e; i++) nn /= p;
        if (e > ext) die("log1p_series: e too large");
        mpz_mul(xp, xp, x); mpz_mod(xp, xp, modx);                 /* x^n mod p^(prec+ext) */
        mpz_divexact(term, xp, P(p, e));                           /* exact: v(x^n) >= n vx >= e */
        if (n >= NINV) die("series too long");
        mpz_mul(term, term, p == 5 ? inv5n[n] : inv2n[n]); mpz_mod(term, term, P(p, prec));
        if (n & 1) mpz_add(out, out, term); else mpz_sub(out, out, term);
        mpz_mod(out, out, P(p, prec));
    }
    mpz_clears(xp, term, nninv, NULL);
}

/* plog: log_p(x) mod p^M for x = 1 mod p (p odd) or x = 1 mod 4 (p = 2); exact port of solver.plog */
static void plog(mpz_t out, const mpz_t x, unsigned long p, long M) {
    long s = 0; while ((s + 1) * (s + 1) <= M) s++;            /* isqrt(M) */
    s += 2;
    long N = M + s;
    mpz_t modN, X, z, total, term, zn, nninv; mpz_inits(modN, X, z, total, term, zn, nninv, NULL);
    mpz_ui_pow_ui(modN, p, N);
    mpz_t ps; mpz_init(ps); mpz_ui_pow_ui(ps, p, s);
    mpz_mod(X, x, modN); mpz_powm(X, X, ps, modN);
    mpz_sub_ui(z, X, 1); mpz_mod(z, z, modN);
    long vz = s + (p > 2 ? 1 : 2);
    mpz_t pv; mpz_init(pv); mpz_ui_pow_ui(pv, p, vz);
    if (!mpz_divisible_p(z, pv)) die("plog: z valuation");
    mpz_set_ui(total, 0);
    mpz_t modext; mpz_init(modext); mpz_ui_pow_ui(modext, p, N + 12);
    mpz_set_ui(zn, 1);
    for (long n = 1;; n++) {
        long lg = floorlog(n, p);
        if (n * vz - lg >= N + 5 && n > 2) break;
        long e = vpl(n, p); long nn = n; for (long i = 0; i < e; i++) nn /= p;
        if (e > 12) die("plog: e too large");
        mpz_mul(zn, zn, z); mpz_mod(zn, zn, modext);                 /* z^n mod p^(N+12) */
        mpz_t pe; mpz_init(pe); mpz_ui_pow_ui(pe, p, e);
        mpz_divexact(term, zn, pe);
        mpz_set_ui(nninv, nn); if (!mpz_invert(nninv, nninv, modN)) die("plog invert");
        mpz_mul(term, term, nninv); mpz_mod(term, term, modN);
        if (n & 1) mpz_add(total, total, term); else mpz_sub(total, total, term);
        mpz_mod(total, total, modN);
        mpz_clear(pe);
    }
    if (!mpz_divisible_p(total, ps)) die("plog: total valuation");
    mpz_divexact(total, total, ps);
    mpz_t modM; mpz_init(modM); mpz_ui_pow_ui(modM, p, M);
    mpz_mod(out, total, modM);
    mpz_clears(modN, X, z, total, term, zn, nninv, ps, pv, modM, modext, NULL);
}

/* ---------------- mt(M,t) cache ---------------- */
typedef struct { int ok; mpz_t q2, q5, Q, Ainv, A5inv, alpha2, alpha, e2, e5, b1x, b1y, b2x, b2y, det; } MT;
static MT *mtc; static long mtM, mtT;
static MT *mt(long M, long t) {
    if (M > Mmax || t > tm || M < 1) die("mt range");
    MT *r = &mtc[M * (mtT + 1) + t];
    if (r->ok) return r;
    mpz_inits(r->q2, r->q5, r->Q, r->Ainv, r->A5inv, r->alpha2, r->alpha, r->e2, r->e5, r->b1x, r->b1y, r->b2x, r->b2y, r->det, NULL);
    long n = M - t;
    if (n < 6) die("mt: n too small");
    mpz_set(r->q2, P(2, n - 2)); mpz_set(r->q5, P(5, M - 1)); mpz_mul(r->Q, r->q2, r->q5);
    mpz_t A, Bc, A5, B5, t1; mpz_inits(A, Bc, A5, B5, t1, NULL);
    mpz_fdiv_r_2exp(A, L2_3, n); mpz_fdiv_q_2exp(A, A, 2);
    mpz_fdiv_r_2exp(Bc, L2_7, n); mpz_fdiv_q_2exp(Bc, Bc, 2);
    mpz_mod(A5, L5_3, P(5, M)); mpz_divexact_ui(A5, A5, 5);      /* L5_3 = 0 mod 5 */
    mpz_mod(B5, L5_7, P(5, M)); mpz_divexact_ui(B5, B5, 5);
    if (!mpz_invert(r->Ainv, A, r->q2)) die("Ainv");
    if (!mpz_invert(r->A5inv, A5, r->q5)) die("A5inv");
    mpz_t alpha5; mpz_init(alpha5);
    mpz_mul(r->alpha2, Bc, r->Ainv); mpz_neg(r->alpha2, r->alpha2); mpz_mod(r->alpha2, r->alpha2, r->q2);
    mpz_mul(alpha5, B5, r->A5inv); mpz_neg(alpha5, alpha5); mpz_mod(alpha5, alpha5, r->q5);
    if (!mpz_invert(t1, r->q5, r->q2)) die("crt1");
    mpz_mul(r->e2, r->q5, t1);
    if (!mpz_invert(t1, r->q2, r->q5)) die("crt2");
    mpz_mul(r->e5, r->q2, t1);
    mpz_mul(r->alpha, r->alpha2, r->e2); mpz_addmul(r->alpha, alpha5, r->e5); mpz_mod(r->alpha, r->alpha, r->Q);
    /* Gauss reduction of b1 = (4 alpha mod Q, 4), b2 = (Q, 0) */
    mpz_mul_ui(r->b1x, r->alpha, 4); mpz_mod(r->b1x, r->b1x, r->Q); mpz_set_ui(r->b1y, 4);
    mpz_set(r->b2x, r->Q); mpz_set_ui(r->b2y, 0);
    mpz_t n1, n2, num, den, k, tmp; mpz_inits(n1, n2, num, den, k, tmp, NULL);
    for (;;) {
        mpz_mul(n1, r->b1x, r->b1x); mpz_addmul(n1, r->b1y, r->b1y);
        mpz_mul(n2, r->b2x, r->b2x); mpz_addmul(n2, r->b2y, r->b2y);
        if (mpz_cmp(n1, n2) > 0) { mpz_swap(r->b1x, r->b2x); mpz_swap(r->b1y, r->b2y); mpz_swap(n1, n2); }
        mpz_mul(num, r->b1x, r->b2x); mpz_addmul(num, r->b1y, r->b2y);
        mpz_set(den, n1);
        /* k = floor((2 num + den) / (2 den)) */
        mpz_mul_ui(tmp, num, 2); mpz_add(tmp, tmp, den); mpz_mul_ui(den, den, 2); mpz_fdiv_q(k, tmp, den);
        if (mpz_sgn(k) == 0) break;
        mpz_submul(r->b2x, k, r->b1x); mpz_submul(r->b2y, k, r->b1y);
    }
    if (mpz_cmp(n1, diam2) <= 0) die("lattice too fine");
    mpz_mul(r->det, r->b1x, r->b2y); mpz_submul(r->det, r->b1y, r->b2x);
    mpz_clears(A, Bc, A5, B5, t1, alpha5, n1, n2, num, den, k, tmp, NULL);
    r->ok = 1;
    return r;
}

/* ---------------- low(Flow) ---------------- */
static long tl_; static int sig_; static mpz_t l2low, l5low, inv5, inv2, oddlow;
static void low(const mpz_t Flow) {
    tl_ = v2z(Flow);
    mpz_tdiv_q_2exp(oddlow, Flow, tl_);           /* exact (Flow divisible by 2^tl), keeps sign */
    mpz_fdiv_r_2exp(T1, oddlow, 2);              /* odd mod 4 in {1,3} */
    sig_ = mpz_cmp_ui(T1, 1) == 0 ? 1 : -1;
    mpz_mul_si(T1, oddlow, sig_); mpz_fdiv_r_2exp(T1, T1, Mmax + 2);
    plog(l2low, T1, 2, Mmax);
    mpz_mod(T1, Flow, P(5, Mmax + 1)); mpz_powm_ui(T1, T1, 4, P(5, Mmax + 1));
    plog(l5low, T1, 5, Mmax);
    mpz_mod(T1, Flow, P(5, Mmax)); if (!mpz_invert(inv5, T1, P(5, Mmax))) die("inv5");
    mpz_fdiv_r_2exp(T1, oddlow, Mmax + 2); if (!mpz_invert(inv2, T1, P(2, Mmax + 2))) die("inv2");
}

static const int DL5[5] = {-1, 0, 1, 3, 2};          /* dlog base 2 mod 5: 1->0, 2->1, 4->2, 3->3 */
static long nsolve = 0, ndirect5 = 0, ndirect2 = 0;
static mpz_t F, x5, l5F, x2, l2F, D5, D2, beta5, beta2, beta, off0, cx2, cy2, ci, cj, up, ww, tq;

/* returns 1 and fills (t, up, ww) if a solution exists */
static int solutions(const mpz_t Flow, long c, long p, long M, long *tout) {
    nsolve++;
    if (M > Mmax) die("M > Mmax");
    mpz_mul_ui(F, pow10[p], (unsigned long)c); mpz_add(F, F, Flow);
    long t = v2z(F);
    if (t >= M || t > tm || mpz_divisible_ui_p(F, 5)) return 0;
    long n = M - t;
    /* 5-adic log of F^4 */
    if (p >= 1) {
        mpz_mul_ui(x5, pow10[p], (unsigned long)c); mpz_mul(x5, x5, inv5); mpz_mod(x5, x5, P(5, M + 10));
        log1p_series(T2, x5, 5, M, p);
        mpz_mul_ui(T2, T2, 4); mpz_add(T2, T2, l5low); mpz_mod(l5F, T2, P(5, M));
    } else {
        ndirect5++;
        mpz_mod(T2, F, P(5, M + 1)); mpz_powm_ui(T2, T2, 4, P(5, M + 1));
        plog(l5F, T2, 5, M);
    }
    /* 2-adic */
    int sigma;
    long vx2 = p + vpl(c, 2) - tl_;
    if (tl_ == t && vx2 >= 2) {
        mpz_mul_ui(T2, pow10[p], (unsigned long)c); mpz_fdiv_q_2exp(T2, T2, tl_);     /* c 10^p >> tl (exact, >= 0) */
        mpz_mul(x2, T2, inv2); mpz_fdiv_r_2exp(x2, x2, n + 40);
        log1p_series(T3, x2, 2, n, vx2);
        mpz_add(T3, T3, l2low); mpz_fdiv_r_2exp(l2F, T3, n);
        sigma = sig_;
    } else {
        ndirect2++;
        mpz_tdiv_q_2exp(T2, F, t); mpz_fdiv_r_2exp(T2, T2, n + 2);                  /* (F >> t) mod 2^(n+2) */
        mpz_fdiv_r_2exp(T3, T2, 2); sigma = mpz_cmp_ui(T3, 1) == 0 ? 1 : -1;
        mpz_mul_si(T2, T2, sigma); mpz_fdiv_r_2exp(T2, T2, n + 2);
        plog(l2F, T2, 2, n);
    }
    MT *r = mt(M, t);
    mpz_mul_ui(T2, L5_16, (unsigned long)t); mpz_sub(T2, l5F, T2); mpz_mod(T2, T2, P(5, M));
    if (!mpz_divisible_ui_p(T2, 5)) die("D5 not divisible");
    mpz_divexact_ui(D5, T2, 5);
    mpz_fdiv_r_2exp(T2, l2F, n); mpz_fdiv_q_2exp(D2, T2, 2);
    mpz_mul(beta5, D5, r->A5inv); mpz_mod(beta5, beta5, r->q5);
    mpz_mul(beta2, D2, r->Ainv); mpz_mod(beta2, beta2, r->q2);
    mpz_mul(beta, beta2, r->e2); mpz_addmul(beta, beta5, r->e5); mpz_mod(beta, beta, r->Q);
    /* rho5 = F * 2^-t mod 5 */
    long Fm5 = (long)mpz_fdiv_ui(F, 5);
    long inv2t = 1; for (long i = 0; i < t; i++) inv2t = (inv2t * 3) % 5;          /* 2^-1 = 3 mod 5 */
    long rho5 = (Fm5 * inv2t) % 5;
    long a2m4 = (long)mpz_fdiv_ui(r->alpha2, 4), b2m4 = (long)mpz_fdiv_ui(beta2, 4);
    for (long w0 = 0; w0 < 4; w0++) {
        long u0 = (a2m4 * w0 + b2m4) % 4;
        if ((u0 + w0) % 2 != (sigma == 1 ? 0 : 1)) continue;
        if ((3 * u0 + w0) % 4 != DL5[rho5]) continue;
        mpz_mul_ui(off0, r->alpha, (unsigned long)w0); mpz_add(off0, off0, beta); mpz_mod(off0, off0, r->Q);
        /* 2 cx = U - 2 off0 ; 2 cy = W - 2 w0 ; ci = (cx b2y - cy b2x)/det ; cj = (b1x cy - b1y cx)/det */
        mpz_mul_ui(cx2, off0, 2); mpz_sub(cx2, U, cx2);
        mpz_set_si(cy2, -2 * w0); mpz_add(cy2, cy2, W);
        mpz_mul_ui(tq, r->det, 2);
        mpz_mul(T4, cx2, r->b2y); mpz_submul(T4, cy2, r->b2x); mpz_fdiv_q(ci, T4, tq);
        mpz_mul(T4, r->b1x, cy2); mpz_submul(T4, r->b1y, cx2); mpz_fdiv_q(cj, T4, tq);
        for (long di = 0; di <= 1; di++) for (long dj = 0; dj <= 1; dj++) {
            mpz_add_ui(T4, ci, di); mpz_add_ui(T5, cj, dj);
            mpz_set(up, off0); mpz_addmul(up, T4, r->b1x); mpz_addmul(up, T5, r->b2x);
            mpz_set_ui(ww, w0); mpz_addmul(ww, T4, r->b1y); mpz_addmul(ww, T5, r->b2y);
            if (mpz_cmp_ui(up, 2) >= 0 && mpz_cmp(up, U) <= 0 && mpz_sgn(ww) >= 0 && mpz_cmp(ww, W) <= 0) {
                *tout = t; return 1;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    FILE *f = argc > 1 ? fopen(argv[1], "r") : stdin;
    if (!f) { perror("open"); return 2; }
    mpz_inits(U, W, diam2, L5_3, L5_7, L5_16, L2_3, L2_7, T1, T2, T3, T4, T5, T6, NULL);
    mpz_inits(l2low, l5low, inv5, inv2, oddlow, F, x5, l5F, x2, l2F, D5, D2, beta5, beta2, beta, off0, cx2, cy2, ci, cj, up, ww, tq, NULL);
    char mode[32]; if (fscanf(f, "%31s", mode) != 1) return 2;
    long G = rdl(f, "G"); tm = rdl(f, "tm"); rdz(f, "U", U); rdz(f, "W", W); Mmax = rdl(f, "Mmax");
    rdz(f, "L5_3", L5_3); rdz(f, "L5_7", L5_7); rdz(f, "L5_16", L5_16); rdz(f, "L2_3", L2_3); rdz(f, "L2_7", L2_7);
    mpz_mul(diam2, U, U); mpz_addmul(diam2, W, W);
    long np = Mmax + 50;
    pow5 = malloc(sizeof(mpz_t) * (np)); pow2 = malloc(sizeof(mpz_t) * (np)); pow10 = malloc(sizeof(mpz_t) * (np));
    for (long i = 0; i < np; i++) { mpz_init(pow5[i]); mpz_ui_pow_ui(pow5[i], 5, i); mpz_init(pow2[i]); mpz_ui_pow_ui(pow2[i], 2, i); mpz_init(pow10[i]); mpz_ui_pow_ui(pow10[i], 10, i); }
    for (long n = 1; n < NINV; n++) {
        long nn = n; while (nn % 5 == 0) nn /= 5;
        mpz_init_set_ui(inv5n[n], nn); if (!mpz_invert(inv5n[n], inv5n[n], pow5[Mmax + 12])) die("inv5n");
        nn = n; while (nn % 2 == 0) nn /= 2;
        mpz_init_set_ui(inv2n[n], nn); if (!mpz_invert(inv2n[n], inv2n[n], pow2[Mmax + 44])) die("inv2n");
    }
    mtT = tm; mtM = Mmax;
    mtc = calloc((Mmax + 1) * (tm + 1), sizeof(MT));
    /* self-check of the passed logarithm constants against the C plog */
    {
        mpz_t chk; mpz_init(chk);
        mpz_set_ui(T1, 81); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_3)) die("L5_3 mismatch");
        mpz_set_ui(T1, 2401); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_7)) die("L5_7 mismatch");
        mpz_set_ui(T1, 16); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_16)) die("L5_16 mismatch");
        mpz_set_si(T1, -3); mpz_fdiv_r_2exp(T1, T1, Mmax + 2); plog(chk, T1, 2, Mmax); if (mpz_cmp(chk, L2_3)) die("L2_3 mismatch");
        mpz_set_si(T1, -7); mpz_fdiv_r_2exp(T1, T1, Mmax + 2); plog(chk, T1, 2, Mmax); if (mpz_cmp(chk, L2_7)) die("L2_7 mismatch");
        mpz_clear(chk);
    }
    mpz_t Flow; mpz_init(Flow);
    if (!strcmp(mode, "selftest")) {
        /* lines: Flow c p M  -> prints "R t" or "R -" */
        long nq = rdl(f, "n");
        for (long q = 0; q < nq; q++) {
            long c, p, M;
            if (fscanf(f, "%s %ld %ld %ld", buf, &c, &p, &M) != 4) return 2;
            mpz_set_str(Flow, buf, 10);
            low(Flow);
            long t;
            if (solutions(Flow, c, p, M, &t)) printf("R %ld\n", t); else printf("R -\n");
        }
        return 0;
    }
    /* proof mode */
    long k = rdl(f, "k"), j = rdl(f, "j");
    long b = k - j;
    long *known = malloc(sizeof(long) * (k + 1));
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, "known")) die("known");
    for (long i = 0; i < k; i++) if (fscanf(f, "%ld", &known[i]) != 1) die("known vals");     /* -1 if unknown */
    long nord = rdl(f, "nord");
    long *cs = malloc(sizeof(long) * (nord * b + 1));
    for (long o = 0; o < nord; o++) for (long i = 0; i < b; i++) if (fscanf(f, "%ld", &cs[o * b + i]) != 1) die("ord");
    long skip = rdl(f, "skip");
    for (long i = j + 1; i < k; i++) if (known[i] < 2) die("needed gap bound missing");
    /* largest M used */
    long Mneed = (b == 1) ? 1 + G : 1 + G;
    if (b >= 2) { long sgap = 0; for (long i = j + 1; i < k; i++) sgap += known[i] - 1; Mneed = 1 + sgap + G; }
    if (Mneed > Mmax) die("Mmax too small for this configuration space");
    long nl = b >= 2 ? b - 2 : 0;                     /* lower gaps g_{k-1} .. g_{j+2} (bottom-up) */
    long *lr = malloc(sizeof(long) * (nl + 1)), *lg = malloc(sizeof(long) * (nl + 1)), *pos = malloc(sizeof(long) * (b + 1));
    for (long i = 0; i < nl; i++) lr[i] = known[k - 1 - i];
    for (long o = skip; o < nord; o++) {
        long *c = cs + o * b;                          /* c[0] = top digit d_{j+1} ... c[b-1] = d_k */
        for (long ak = 0; ak <= 1; ak++) {
            for (long i = 0; i < nl; i++) lg[i] = 1;
            int done = 0;
            while (!done) {
                if (b == 1) {
                    mpz_set_si(Flow, -1);
                    low(Flow);
                    long t;
                    if (solutions(Flow, c[0], ak, ak + G, &t)) {
                        printf("FAIL %ld %ld %ld\n", o, ak, t); fflush(stdout); return 0;
                    }
                } else {
                    pos[0] = ak;
                    for (long i = 0; i < nl; i++) pos[i + 1] = pos[i] + lg[i];
                    mpz_set_si(Flow, -1);
                    for (long i = 0; i <= nl; i++) mpz_addmul_ui(Flow, pow10[pos[i]], (unsigned long)c[b - 1 - i]);
                    low(Flow);
                    for (long g = 1; g < known[j + 1]; g++) {
                        long p = pos[nl] + g, t;
                        if (solutions(Flow, c[0], p, p + G, &t)) {
                            printf("FAIL %ld %ld", o, ak);
                            for (long i = 0; i < nl; i++) printf(" %ld", lg[i]);
                            printf(" g %ld t %ld\n", g, t); fflush(stdout); return 0;
                        }
                    }
                }
                long i = nl - 1;
                while (i >= 0) { lg[i]++; if (lg[i] < lr[i]) break; lg[i] = 1; i--; }
                if (i < 0) done = 1;
            }
        }
        printf("PROG %ld %ld\n", o + 1, nsolve); fflush(stdout);
    }
    printf("DONE %ld %ld %ld\n", nsolve, ndirect5, ndirect2);
    return 0;
}
