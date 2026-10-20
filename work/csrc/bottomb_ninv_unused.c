/* Cumulative bottom test for bases with one free prime q (port of basesolver.solve_c + the cumulative
   height enumeration of prove_cum_b.bottom_test).
   For each configuration (digits j+1..k at heights a_{j+1} > ... > a_k >= 0, a_i < heff[i]; none if j = k):
     M = max(Gb, a_{j+1}+1)  (Gb if no digits),   r = (b^M - 1)/(b - 1) + sum (d_i - 1) b^{a_i}   (0 < r < b^M)
   decide whether B q^c = r (mod b^M) has a solution with c in [0, Xc] and B an allowed cofactor
   (B = prod_{p|b} p^{beta_p}, beta_p <= pb_p).  Since every pb_p < M v_p(b), beta_p = v_p(r) is forced.
   Then q^c = rho_p (mod p^{n_p}), n_p = M v_p(b) - beta_p, solved with p-adic logarithms exactly as in
   solve_c (residue mod p by brute force for odd p, sign mod 4 for p = 2), and combined by generalized CRT.
   Any doubt (CRT modulus <= Xc, too little 2-adic precision) is reported as a failure. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

static void die(const char *m) { fprintf(stderr, "%s\n", m); exit(3); }
static long rdl(FILE *f, const char *name) {
    char key[64]; long v;
    if (fscanf(f, "%63s %ld", key, &v) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s (got %s)\n", name, key); exit(2); }
    return v;
}
static char buf[1 << 16];
static void rdz(FILE *f, const char *name, mpz_t x) {
    char key[64];
    if (fscanf(f, "%63s %s", key, buf) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    if (mpz_set_str(x, buf, 10)) die("bad int");
}
static long floorlog(long n, long p) { long lg = 0, q = p; while (q <= n) { q *= p; lg++; } return lg; }
static long vpl(long n, long p) { long e = 0; while (n % p == 0) { n /= p; e++; } return e; }
static long vpz(const mpz_t x, unsigned long p) {
    if (mpz_sgn(x) == 0) return 1000000000L;
    if (p == 2) return (long)mpz_scan1(x, 0);
    mpz_t t; mpz_init_set(t, x); long e = 0;
    while (mpz_divisible_ui_p(t, p)) { mpz_divexact_ui(t, t, p); e++; }
    mpz_clear(t); return e;
}
#define NINV 600
static long BIGE = 0;                       /* inverses valid modulo p^BIGE */
static mpz_t ninv2[NINV], ninv3[NINV], ninv5[NINV], bigmod2, bigmod3, bigmod5;
static int ninv_ready = 0;
static void init_ninv(long E) {
    BIGE = E;
    mpz_inits(bigmod2, bigmod3, bigmod5, NULL);
    mpz_ui_pow_ui(bigmod2, 2, E); mpz_ui_pow_ui(bigmod3, 3, E); mpz_ui_pow_ui(bigmod5, 5, E);
    for (long n = 1; n < NINV; n++) {
        long nn;
        nn = n; while (nn % 2 == 0) nn /= 2; mpz_init_set_ui(ninv2[n], nn); if (!mpz_invert(ninv2[n], ninv2[n], bigmod2)) exit(5);
        nn = n; while (nn % 3 == 0) nn /= 3; mpz_init_set_ui(ninv3[n], nn); if (!mpz_invert(ninv3[n], ninv3[n], bigmod3)) exit(5);
        nn = n; while (nn % 5 == 0) nn /= 5; mpz_init_set_ui(ninv5[n], nn); if (!mpz_invert(ninv5[n], ninv5[n], bigmod5)) exit(5);
    }
    ninv_ready = 1;
}
/* exact port of solver.plog: log_p(x) mod p^M for x = 1 mod p (odd p) / x = 1 mod 4 (p = 2) */
static void plog(mpz_t out, const mpz_t x, unsigned long p, long M) {
    long s = 0; while ((s + 1) * (s + 1) <= M) s++;
    s += 2;
    long N = M + s;
    mpz_t modN, X, z, total, term, zn, nninv, ps, pv, modext, pe, modM;
    mpz_inits(modN, X, z, total, term, zn, nninv, ps, pv, modext, pe, modM, NULL);
    mpz_ui_pow_ui(modN, p, N); mpz_ui_pow_ui(ps, p, s);
    mpz_mod(X, x, modN); mpz_powm(X, X, ps, modN);
    mpz_sub_ui(z, X, 1); mpz_mod(z, z, modN);
    long vz = s + (p > 2 ? 1 : 2);
    mpz_ui_pow_ui(pv, p, vz);
    if (!mpz_divisible_p(z, pv)) die("plog: z valuation");
    mpz_set_ui(total, 0); mpz_ui_pow_ui(modext, p, N + 12); mpz_set_ui(zn, 1);
    for (long n = 1;; n++) {
        long lg = floorlog(n, p);
        if (n * vz - lg >= N + 5 && n > 2) break;
        long e = vpl(n, p); long nn = n; for (long i = 0; i < e; i++) nn /= p;
        if (e > 12) die("plog: e");
        mpz_mul(zn, zn, z); mpz_mod(zn, zn, modext);
        mpz_ui_pow_ui(pe, p, e); mpz_divexact(term, zn, pe);
        if (ninv_ready && n < NINV && N <= BIGE && (p == 2 || p == 3 || p == 5)) {
            mpz_mul(term, term, p == 2 ? ninv2[n] : (p == 3 ? ninv3[n] : ninv5[n]));
        } else {
            mpz_set_ui(nninv, nn); if (!mpz_invert(nninv, nninv, modN)) die("plog invert");
            mpz_mul(term, term, nninv);
        }
        mpz_mod(term, term, modN);
        if (n & 1) mpz_add(total, total, term); else mpz_sub(total, total, term);
        mpz_mod(total, total, modN);
    }
    if (!mpz_divisible_p(total, ps)) die("plog: total valuation");
    mpz_divexact(total, total, ps);
    mpz_ui_pow_ui(modM, p, M); mpz_mod(out, total, modM);
    mpz_clears(modN, X, z, total, term, zn, nninv, ps, pv, modext, pe, modM, NULL);
}

static long b, q, nd, P[4], E[4], PB[4], Gb;
static mpz_t Xc;
/* cache of log_p(q) data per (p index, n) */
static int lqok[4][400]; static mpz_t lqv[4][400]; static long lqval[4][400];
static int sq_;           /* sign class of q mod 4 */
static long ordp[4]; static long qpow_res[4][64];   /* q^c mod p for c < ord */

/* returns 1 if some admissible (B, c) exists (or the decision is unsafe), 0 if none */
static long nsolve = 0;
static int solve(const mpz_t r, long M) {
    nsolve++;
    mpz_t rho, Bp, t, mod, lr, c0, Nn, rr, g, l, kk, tmp, pw;
    mpz_inits(rho, Bp, t, mod, lr, c0, Nn, rr, g, l, kk, tmp, pw, NULL);
    int res = 0;
    long beta[4];
    for (long i = 0; i < nd; i++) {
        long e = M * E[i];
        long v = vpz(r, P[i]);
        if (v >= e) { goto none; }                 /* would need beta >= e > pb: impossible */
        if (v > PB[i]) goto none;
        beta[i] = v;
    }
    /* B = prod p^beta; for each p: rho_p = (r / p^beta) * (B / p^beta)^-1 mod p^n */
    mpz_set_ui(c0, 0); mpz_set_ui(Nn, 1);
    for (long i = 0; i < nd; i++) {
        long p = P[i], n = M * E[i] - beta[i];
        if (n <= 0) continue;
        mpz_ui_pow_ui(mod, p, n);
        mpz_ui_pow_ui(pw, p, beta[i]); mpz_divexact(rho, r, pw);
        mpz_set_ui(Bp, 1);
        for (long i2 = 0; i2 < nd; i2++) if (i2 != i) { mpz_ui_pow_ui(t, P[i2], beta[i2]); mpz_mul(Bp, Bp, t); }
        if (!mpz_invert(t, Bp, mod)) die("B part not invertible");
        mpz_mul(rho, rho, t); mpz_mod(rho, rho, mod);
        /* congruences for c */
        long cr[2]; mpz_t cm[2]; mpz_init(cm[0]); mpz_init(cm[1]); int ncg = 0; mpz_t crz[2]; mpz_init(crz[0]); mpz_init(crz[1]);
        if (n >= 400) die("n too large");
        if (p == 2) {
            if (n < 4) { res = 1; goto done; }       /* too little information: unsafe -> fail */
            long sr = mpz_fdiv_ui(rho, 4) == 1 ? 1 : -1;
            if (sq_ == 1 && sr == -1) { mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]); goto none; }
            if (sq_ == -1) { mpz_set_ui(crz[ncg], sr == 1 ? 0 : 1); mpz_set_ui(cm[ncg], 2); ncg++; }
            if (!lqok[i][n]) {
                mpz_init(lqv[i][n]); mpz_set_si(t, sq_ * q); mpz_fdiv_r_2exp(t, t, n + 2); plog(lqv[i][n], t, 2, n);
                lqval[i][n] = vpz(lqv[i][n], 2); lqok[i][n] = 1;
            }
            mpz_mul_si(t, rho, sr); mpz_fdiv_r_2exp(t, t, n + 2); plog(lr, t, 2, n);
            long v = lqval[i][n];
            mpz_fdiv_r_2exp(tmp, lr, n);
            if (vpz(lr, 2) < v && mpz_sgn(tmp) != 0) { mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]); goto none; }
            mpz_ui_pow_ui(cm[ncg], 2, n - v);
            mpz_fdiv_q_2exp(t, lqv[i][n], v); if (!mpz_invert(t, t, cm[ncg])) die("lq inv");
            mpz_fdiv_q_2exp(tmp, lr, v); mpz_mul(tmp, tmp, t); mpz_mod(crz[ncg], tmp, cm[ncg]); ncg++;
        } else {
            long rm = (long)mpz_fdiv_ui(rho, p), found = -1, cnt = 0;
            for (long c = 0; c < ordp[i]; c++) if (qpow_res[i][c] == rm) { found = c; cnt++; }
            if (found < 0) { mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]); goto none; }
            if (cnt != 1) die("residue not unique");
            mpz_set_ui(crz[ncg], found); mpz_set_ui(cm[ncg], ordp[i]); ncg++;
            if (!lqok[i][n]) {
                mpz_init(lqv[i][n]); mpz_ui_pow_ui(tmp, p, n + 1); mpz_set_ui(t, q); mpz_powm_ui(t, t, p - 1, tmp); plog(lqv[i][n], t, p, n);
                lqval[i][n] = vpz(lqv[i][n], p); lqok[i][n] = 1;
            }
            mpz_ui_pow_ui(tmp, p, n + 1); mpz_powm_ui(t, rho, p - 1, tmp); plog(lr, t, p, n);
            long v = lqval[i][n];
            mpz_mod(tmp, lr, mod);
            if (mpz_sgn(tmp) != 0 && vpz(lr, p) < v) { mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]); goto none; }
            mpz_ui_pow_ui(cm[ncg], p, n - v);
            mpz_ui_pow_ui(pw, p, v); mpz_divexact(t, lqv[i][n], pw); if (!mpz_invert(t, t, cm[ncg])) die("lq inv odd");
            mpz_divexact(tmp, lr, pw); mpz_mul(tmp, tmp, t); mpz_mod(crz[ncg], tmp, cm[ncg]); ncg++;
        }
        /* generalized CRT: c = c0 mod Nn combined with crz[i] mod cm[i] */
        for (int ci = 0; ci < ncg; ci++) {
            mpz_gcd(g, Nn, cm[ci]);
            mpz_sub(tmp, crz[ci], c0); 
            if (!mpz_divisible_p(tmp, g)) { mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]); goto none; }
            mpz_divexact(l, Nn, g); mpz_mul(l, l, cm[ci]);                 /* lcm */
            mpz_divexact(kk, tmp, g);                                       /* (r - c0)/g */
            mpz_divexact(t, Nn, g); mpz_divexact(rr, cm[ci], g);
            if (mpz_cmp_ui(rr, 1) == 0) mpz_set_ui(kk, 0);
            else { if (!mpz_invert(t, t, rr)) die("crt inv"); mpz_mul(kk, kk, t); mpz_mod(kk, kk, rr); }
            mpz_mul(tmp, Nn, kk); mpz_add(c0, c0, tmp); mpz_mod(c0, c0, l); mpz_set(Nn, l);
        }
        mpz_clear(cm[0]); mpz_clear(cm[1]); mpz_clear(crz[0]); mpz_clear(crz[1]);
    }
    if (mpz_cmp(Nn, Xc) <= 0) { res = 1; goto done; }                  /* congruence too weak: unsafe -> fail */
    res = mpz_cmp(c0, Xc) <= 0;
    goto done;
none:
    res = 0;
done:
    mpz_clears(rho, Bp, t, mod, lr, c0, Nn, rr, g, l, kk, tmp, pw, NULL);
    return res;
}

static long kk_, jj_, nb_, *heff, *pos, *dm1;
static mpz_t *powb;
static long cur_o;
static int rec(long r, long lo) {
    if (r < 0) {
        long M = nb_ ? (pos[0] + 1 > Gb ? pos[0] + 1 : Gb) : Gb;
        mpz_t y; mpz_init(y);
        mpz_sub_ui(y, powb[M], 1); mpz_divexact_ui(y, y, b - 1);
        for (long i = 0; i < nb_; i++) mpz_addmul_ui(y, powb[pos[i]], (unsigned long)dm1[cur_o * nb_ + i]);
        int s = solve(y, M);
        if (s) {
            printf("FAIL %ld pos", cur_o); for (long i = 0; i < nb_; i++) printf(" %ld", pos[i]); printf(" M %ld\n", M); fflush(stdout);
        }
        mpz_clear(y);
        return s;
    }
    for (long a = lo; a < heff[r]; a++) { pos[r] = a; if (rec(r - 1, a + 1)) return 1; }
    return 0;
}

int main(int argc, char **argv) {
    FILE *f = argc > 1 ? fopen(argv[1], "r") : stdin;
    if (!f) { perror("open"); return 2; }
    mpz_init(Xc);
    b = rdl(f, "b"); q = rdl(f, "q"); nd = rdl(f, "nd");
    for (long i = 0; i < nd; i++) { P[i] = rdl(f, "p"); E[i] = rdl(f, "e"); PB[i] = rdl(f, "pb"); }
    rdz(f, "Xc", Xc); Gb = rdl(f, "Gb");
    kk_ = rdl(f, "k"); jj_ = rdl(f, "j"); nb_ = kk_ - jj_;
    heff = malloc(sizeof(long) * (nb_ + 1)); pos = malloc(sizeof(long) * (nb_ + 1));
    char key[64]; if (fscanf(f, "%63s", key) != 1 || strcmp(key, "heff")) die("heff");
    for (long i = 0; i < nb_; i++) if (fscanf(f, "%ld", &heff[i]) != 1) die("heff v");
    long nord = rdl(f, "nord");
    dm1 = malloc(sizeof(long) * (nord * nb_ + 1));
    for (long o = 0; o < nord; o++) for (long i = 0; i < nb_; i++) if (fscanf(f, "%ld", &dm1[o * nb_ + i]) != 1) die("ord");
    long skip = rdl(f, "skip");
    for (long i = 0; i < nd; i++) if (PB[i] >= Gb * E[i]) die("Gb too small for cofactor bounds");
    long Mmax = Gb; if (nb_ && heff[0] > Mmax) Mmax = heff[0];
    powb = malloc(sizeof(mpz_t) * (Mmax + 2));
    for (long i = 0; i <= Mmax + 1; i++) { mpz_init(powb[i]); mpz_ui_pow_ui(powb[i], b, i); }
    sq_ = (q % 4 == 1) ? 1 : -1;
    init_ninv(Mmax + 40);
    for (long i = 0; i < nd; i++) if (P[i] != 2) {
        long o = 1, x = q % P[i]; while (x != 1) { x = x * q % P[i]; o++; } ordp[i] = o;
        long y = 1; for (long c = 0; c < o; c++) { qpow_res[i][c] = y; y = y * q % P[i]; }
    }
    for (long o = skip; o < nord; o++) {
        cur_o = o;
        if (nb_ == 0) { if (rec(-1, 0)) return 0; }
        else if (rec(nb_ - 1, 0)) return 0;
        printf("PROG %ld %ld\n", o + 1, nsolve); fflush(stdout);
    }
    printf("DONE %ld\n", nsolve);
    return 0;
}
