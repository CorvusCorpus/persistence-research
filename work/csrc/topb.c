/* Top-side case loop of the base-10 prover (port of prove_iter2.prove_top + FastTop).
   For every top configuration (ordered digits d_1..d_j, gaps g_0..g_{j-1} in proved ranges):
     1+z = 1 + sum_i c_i 10^{-(g_0+..+g_i)},  c_i = 9(d_i - 1)
     m   = nearest integer to v = -C log(1+z)        (|v - m| <= 0.5 + 2^-50, rigorous: see below)
   Level-1 kill (as in Python):  z_last = -m qa/qb;  d = dist(z_last, Z);  kill if d^2 gs2[n-1] > R2.
   Level-2 kill (new, exact rational): for every integer x_{n-1} with (x-z_last)^2 gs2[n-1] <= R2,
     c = z_{n-2} - mu (x - z_last);  kill if dist(c,Z)^2 gs2[n-2] + (x-z_last)^2 gs2[n-1] > R2 for all x.
   This is the first two levels of the Fincke-Pohst enumeration used by Python's valid_points, so a
   configuration killed here has no lattice point within R of the target.  Survivors are printed and
   re-checked exactly in Python (valid_points).
   Error bound for v: every MPFR operation is correctly rounded (relative error <= 2^-prec); z is a sum of
   positive terms, log1p preserves relative error for positive arguments, so the computed log(1+z) has
   relative error <= (j+12) 2^-prec and |v_computed - v| <= C*2.2*(j+12)*2^-prec; with
   prec >= bitlen(C) + 64 and j <= 40 this is < 2^-50.  We check prec >= bitlen(C)+64 at start. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

static char buf[1 << 20];
static void rdz(FILE *f, const char *name, mpz_t x) {
    char key[64];
    if (fscanf(f, "%63s %s", key, buf) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    if (mpz_set_str(x, buf, 10)) { fprintf(stderr, "bad int %s\n", name); exit(2); }
}
static long rdl(FILE *f, const char *name) {
    char key[64]; long v;
    if (fscanf(f, "%63s %ld", key, &v) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    return v;
}

static void flush_run(long o, long *g, long j, long lo, long hi, mpz_t mm) {
    printf("SURV %ld", o);
    for (long i = 0; i < j - 1; i++) printf(" %ld", g[i]);
    printf(" %ld %ld ", lo, hi); mpz_out_str(stdout, 10, mm); printf("\n");
}
static mpz_t run_m;

int main(int argc, char **argv) {
    mpz_init(run_m);
    FILE *f = argc > 1 ? fopen(argv[1], "r") : stdin;
    if (!f) { perror("open"); return 2; }
    mpz_t C, qa, qb, K1, K2, D, a2, a3, P, Q, g2n, g2d, g3n, g3d, Rn, Rd;
    mpz_inits(C, qa, qb, K1, K2, D, a2, a3, P, Q, g2n, g2d, g3n, g3d, Rn, Rd, NULL);
    long prec = rdl(f, "prec");
    rdz(f, "C", C);
    rdz(f, "qa", qa); rdz(f, "qb", qb); rdz(f, "K1", K1); rdz(f, "K2", K2);
    rdz(f, "D", D); rdz(f, "a2", a2); rdz(f, "a3", a3); rdz(f, "P", P); rdz(f, "Q", Q);
    rdz(f, "g2n", g2n); rdz(f, "g2d", g2d); rdz(f, "g3n", g3n); rdz(f, "g3d", g3d);
    rdz(f, "Rn", Rn); rdz(f, "Rd", Rd);
    long base = rdl(f, "base");
    long ntg = rdl(f, "ntg");
    mpz_t *Wt = malloc(sizeof(mpz_t) * (ntg + 1));
    for (long i = 0; i < ntg; i++) { mpz_init(Wt[i]); char nm[16]; snprintf(nm, sizeof nm, "W%ld", i); rdz(f, nm, Wt[i]); if (mpz_sgn(Wt[i]) <= 0) { fprintf(stderr, "bad W\n"); return 2; } }
    long j = rdl(f, "j");
    long cum = rdl(f, "cum");
    long *rng = malloc(sizeof(long) * (j + 1));
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, "ranges")) { fprintf(stderr, "bad ranges\n"); return 2; }
    for (long i = 0; i < j; i++) if (fscanf(f, "%ld", &rng[i]) != 1) return 2;
    long nord = rdl(f, "nord");
    long *cs = malloc(sizeof(long) * (nord * j + 1));
    for (long o = 0; o < nord; o++) for (long i = 0; i < j; i++) if (fscanf(f, "%ld", &cs[o * j + i]) != 1) return 2;
    long skip = rdl(f, "skip");
    long maxsurv = rdl(f, "maxsurv");
    long debug = rdl(f, "debug");
    long g0lo = rdl(f, "g0lo"), g0hi = rdl(f, "g0hi");      /* restrict the first gap g_0 to [g0lo, g0hi] */
    if (prec < (long)mpz_sizeinbase(C, 2) + 64) { fprintf(stderr, "prec too small\n"); return 2; }
    if (j < 1 || j > 40) { fprintf(stderr, "bad j\n"); return 2; }
    if (mpz_sgn(Q) <= 0 || mpz_sgn(D) <= 0 || mpz_sgn(qb) <= 0) { fprintf(stderr, "bad denominators\n"); return 2; }

    mpfr_set_default_prec(prec);
    long emax = 0; for (long i = 0; i < j; i++) emax += rng[i];   /* covers both modes (cum: acc+gl <= rng[j-1]-1) */
    mpfr_t *p10 = malloc(sizeof(mpfr_t) * (emax + 2));      /* p10[e] = 10^-e (rounded) */
    for (long e = 0; e <= emax + 1; e++) { mpfr_init2(p10[e], prec); mpfr_ui_pow_ui(p10[e], (unsigned long)base, e, MPFR_RNDN); mpfr_ui_div(p10[e], 1, p10[e], MPFR_RNDN); }
    /* offsets C*log(W_i) (correctly rounded log of an exact integer, then one rounded product), split as
       Ib_i + fb_i (integer part, fractional part as double); RIb_i = (-Ib_i qa) mod qb */
    mpfr_t *offs = malloc(sizeof(mpfr_t) * (ntg + 1));
    mpz_t *Ib = malloc(sizeof(mpz_t) * (ntg + 1)), *RIb = malloc(sizeof(mpz_t) * (ntg + 1));
    double *fb = malloc(sizeof(double) * (ntg + 1));
    mpfr_t fr; mpfr_init2(fr, prec);
    for (long i = 0; i < ntg; i++) {
        mpfr_init2(offs[i], prec); mpfr_set_z(offs[i], Wt[i], MPFR_RNDN); mpfr_log(offs[i], offs[i], MPFR_RNDN); mpfr_mul_z(offs[i], offs[i], C, MPFR_RNDN);
        mpz_init(Ib[i]); mpfr_get_z(Ib[i], offs[i], MPFR_RNDD);
        mpfr_sub_z(fr, offs[i], Ib[i], MPFR_RNDN); fb[i] = mpfr_get_d(fr, MPFR_RNDN);
        mpz_init(RIb[i]); mpz_mul(RIb[i], Ib[i], qa); mpz_neg(RIb[i], RIb[i]); mpz_fdiv_r(RIb[i], RIb[i], qb);
    }
    /* level-1 kill  <=>  d^2 K1 > K2  <=>  d >= Th,  Th = ceil(sqrt(floor(K2/K1) + 1)) */
    mpz_t Th, Q0, qam; mpz_inits(Th, Q0, qam, NULL);
    mpz_fdiv_q(Q0, K2, K1); mpz_add_ui(Q0, Q0, 1);
    mpz_sqrt(Th, Q0); { mpz_t sq; mpz_init(sq); mpz_mul(sq, Th, Th); if (mpz_cmp(sq, Q0) < 0) mpz_add_ui(Th, Th, 1); mpz_clear(sq); }
    mpz_fdiv_r(qam, qa, qb);
    mpz_t Iv, Rv; mpz_inits(Iv, Rv, NULL); double fv;
    mpfr_t zp, A, LA, u, LU, v, tmp;
    mpfr_inits2(prec, zp, A, LA, u, LU, v, tmp, (mpfr_ptr)0);
    mpz_t m, t, r, d, lhs, rhs, x3, E3, N2, QD, r2, d2, base3, KK, fl;
    mpz_inits(m, t, r, d, lhs, rhs, x3, E3, N2, QD, r2, d2, base3, KK, fl, NULL);
    /* level-2 constants */
    mpz_t L2a, L2b, L2c, R3a, R3b;          /* survive iff d2^2*L2a + E3^2*L2b <= L2c ; level-3 pre-test E3^2*R3a <= R3b */
    mpz_inits(L2a, L2b, L2c, R3a, R3b, NULL);
    mpz_mul(L2a, g2n, g3d); mpz_mul(L2a, L2a, Rd);                       /* d2^2 g2n g3d Rd */
    mpz_mul(L2b, Q, Q); mpz_mul(L2b, L2b, g3n); mpz_mul(L2b, L2b, g2d); mpz_mul(L2b, L2b, Rd);  /* E3^2 Q^2 g3n g2d Rd */
    mpz_mul(L2c, Rn, g2d); mpz_mul(L2c, L2c, g3d); mpz_mul(L2c, L2c, Q); mpz_mul(L2c, L2c, Q); mpz_mul(L2c, L2c, D); mpz_mul(L2c, L2c, D);
    mpz_mul(R3a, g3n, Rd); mpz_mul(R3b, Rn, g3d); mpz_mul(R3b, R3b, D); mpz_mul(R3b, R3b, D);  /* E3^2 g3n Rd <= Rn g3d D^2 */
    mpz_mul(QD, Q, D);
    /* KK >= sqrt(R2/gs2[n-1]) : KK = isqrt(floor(Rn*g3d/(Rd*g3n))) + 1 */
    mpz_mul(t, Rn, g3d); mpz_mul(r, Rd, g3n); mpz_fdiv_q(t, t, r); mpz_sqrt(KK, t); mpz_add_ui(KK, KK, 1);
    long Kl = mpz_get_si(KK);
    if (!mpz_fits_slong_p(KK) || Kl > 1000) { fprintf(stderr, "level-2 radius too large\n"); return 2; }

    long *g = malloc(sizeof(long) * (j + 1));
    long *ss = malloc(sizeof(long) * (j + 1));
    if (!cum) { for (long i = 0; i < j; i++) if (rng[i] < 2) { printf("DONE 0 0 0 0\n"); return 0; } }
    unsigned long long cases = 0, f1 = 0, f2 = 0, nsurv = 0;
    for (long o = skip; o < nord; o++) {
        long *c = cs + o * j;
        /* enumerate prefixes g_0..g_{j-2}; old mode: g_i in [1, rng_i); cum mode: s_{i+1} = s_i + g_i <= rng_i - 1 */
        long i = 0; ss[0] = 0;
        if (j >= 2) g[0] = 0;
        for (;;) {
            if (j >= 2) {
                /* advance the prefix odometer to the next complete prefix (depth j-1) */
                int found = 0;
                while (i >= 0) {
                    g[i]++;
                    if (i == 0 && g[0] < g0lo) g[0] = g0lo;
                    long hi = cum ? rng[i] - 1 - ss[i] : rng[i] - 1;
                    if (i == 0 && hi > g0hi) hi = g0hi;
                    if (g[i] > hi) { i--; continue; }
                    ss[i + 1] = ss[i] + g[i];
                    if (i == j - 2) { found = 1; break; }
                    i++; g[i] = 0;
                }
                if (!found) break;
            }
            /* prefix complete: z' = sum_{i<j-1} c_i 10^-acc_i */
            long acc = 0;
            mpfr_set_ui(zp, 0, MPFR_RNDN);
            for (long q = 0; q < j - 1; q++) { acc += g[q]; mpfr_mul_ui(tmp, p10[acc], c[q], MPFR_RNDN); mpfr_add(zp, zp, tmp, MPFR_RNDN); }
            mpfr_add_ui(A, zp, 1, MPFR_RNDN);
            mpfr_log1p(LA, zp, MPFR_RNDN);
            long cl = c[j - 1];
            long ghi = cum ? rng[j - 1] - 1 - acc : rng[j - 1] - 1;
            long glo = 1;
            if (j == 1) { if (glo < g0lo) glo = g0lo; if (ghi > g0hi) ghi = g0hi; }
            for (long gl = glo; gl <= ghi; gl++) {
                long e = acc + gl;
                mpfr_mul_ui(u, p10[e], cl, MPFR_RNDN);
                mpfr_div(u, u, A, MPFR_RNDN);
                mpfr_log1p(LU, u, MPFR_RNDN);
                mpfr_add(v, LA, LU, MPFR_RNDN);
                mpfr_mul_z(v, v, C, MPFR_RNDN);
                cases++;
                mpfr_get_z(Iv, v, MPFR_RNDD);
                mpfr_sub_z(fr, v, Iv, MPFR_RNDN); fv = mpfr_get_d(fr, MPFR_RNDN);
                mpz_mul(Rv, Iv, qa); mpz_fdiv_r(Rv, Rv, qb);
                for (long ti = 0; ti < ntg; ti++) {
                /* m = round(C log W - v) = Ib - Iv + round(fb - fv) (either rounding at a near-tie keeps |m - x| < 0.51) */
                double fx = fb[ti] - fv;
                long dl = fx >= 0.5 ? 1 : (fx < -0.5 ? -1 : 0);
                /* level 1: r = (-m qa) mod qb = RIb + Iv qa - dl qa (mod qb) */
                mpz_add(r, RIb[ti], Rv);
                if (dl == 1) mpz_sub(r, r, qam); else if (dl == -1) mpz_add(r, r, qam);
                while (mpz_cmp(r, qb) >= 0) mpz_sub(r, r, qb);
                while (mpz_sgn(r) < 0) mpz_add(r, r, qb);
                mpz_sub(d, qb, r); if (mpz_cmp(r, d) < 0) mpz_set(d, r);
                if (debug == 2) {                 /* cross-check against direct rounding and the multiplication form */
                    mpfr_sub(tmp, offs[ti], v, MPFR_RNDN); mpz_t mm, rr2; mpz_inits(mm, rr2, NULL); mpfr_get_z(mm, tmp, MPFR_RNDN);
                    mpz_sub(m, Ib[ti], Iv); if (dl == 1) mpz_add_ui(m, m, 1); else if (dl == -1) mpz_sub_ui(m, m, 1);
                    mpz_sub(rr2, mm, m);
                    if (mpz_sgn(rr2) != 0 && !(fabs(fx - 0.5) < 1e-6 || fabs(fx + 0.5) < 1e-6)) { fprintf(stderr, "target mismatch away from a tie\n"); exit(4); }
                    mpz_mul(rr2, m, qa); mpz_neg(rr2, rr2); mpz_fdiv_r(rr2, rr2, qb);
                    if (mpz_cmp(rr2, r) != 0) { fprintf(stderr, "residue mismatch\n"); exit(4); }
                    mpz_mul(lhs, d, d); mpz_mul(lhs, lhs, K1);
                    if ((mpz_cmp(lhs, K2) > 0) != (mpz_cmp(d, Th) >= 0)) { fprintf(stderr, "threshold mismatch\n"); exit(4); }
                    mpz_clears(mm, rr2, NULL);
                }
                if (mpz_cmp(d, Th) >= 0) continue;
                mpz_sub(m, Ib[ti], Iv); if (dl == 1) mpz_add_ui(m, m, 1); else if (dl == -1) mpz_sub_ui(m, m, 1);
                f1++;
                /* level 2 */
                int surv = 0;
                mpz_mul(base3, m, a3);
                mpz_neg(t, base3); mpz_fdiv_q(fl, t, D);
                for (long dx = -Kl; dx <= Kl + 1 && !surv; dx++) {
                    if (dx >= 0) mpz_add_ui(x3, fl, (unsigned long)dx); else mpz_sub_ui(x3, fl, (unsigned long)(-dx));
                    mpz_mul(E3, x3, D); mpz_add(E3, E3, base3);
                    mpz_mul(lhs, E3, E3); mpz_mul(rhs, lhs, R3a);
                    if (mpz_cmp(rhs, R3b) > 0) continue;
                    mpz_mul(N2, m, a2); mpz_mul(N2, N2, Q); mpz_neg(N2, N2);
                    mpz_submul(N2, P, E3);
                    mpz_fdiv_r(r2, N2, QD); mpz_sub(d2, QD, r2); if (mpz_cmp(r2, d2) < 0) mpz_set(d2, r2);
                    mpz_mul(rhs, d2, d2); mpz_mul(rhs, rhs, L2a);
                    mpz_addmul(rhs, lhs, L2b);
                    if (mpz_cmp(rhs, L2c) <= 0) surv = 1;
                }
                if (debug == 1) {
                    printf("F1 %ld", o);
                    for (long q = 0; q < j - 1; q++) printf(" %ld", g[q]);
                    printf(" %ld %d\n", gl, surv);
                }
                if (!surv) continue;
                f2++;
                /* survivors are reported individually with their target index (no run compression) */
                printf("SURV %ld", o);
                for (long q = 0; q < j - 1; q++) printf(" %ld", g[q]);
                printf(" %ld %ld %ld ", gl, gl, ti); mpz_out_str(stdout, 10, m); printf("\n");
                nsurv++;
                }
            }
            if (j < 2) break;
        }
        printf("PROG %ld %llu %llu %llu\n", o + 1, cases, f1, f2);
        fflush(stdout);
    }
    printf("DONE %llu %llu %llu %llu\n", cases, f1, f2, nsurv);
    (void)maxsurv;
    return 0;
}
