/* topgpu: drop-in replacement for csrc/topc3 (cum mode, debug 0), with the GPU filter of topgpu.cu in front of top3.c's exact
 * per-case code. Same job file, same output: SURV runs, PROG per ordering, DONE cases f1 f2 nsurv.
 * Every case the kernel does not certainly kill is decided here by top3.c's own code, copied below unchanged in substance
 * (same MPFR precision, the same operations in the same order, the same exact level-1 and level-2 tests and run merging); a
 * killed case is one top3.c kills (level 1, or level 2 after a level-1 pass, which the kernel counts in f1). Jobs outside the
 * filter (cum 0, debug != 0, |theta| >= 2^-60, more than 40 digits, heights >= 600) are handed to topc3 (exec).
 * usage: topgpu JOBFILE      env TP_BATCH (prefixes per launch, default 2^18), TP_CPU (path of topc3)
 * build: see work/Makefile (topgpu: CUDA; topgpu_emu: the same kernel compiled for the CPU) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "topgpu.h"

static char buf[1 << 20];
static void rdz(FILE *f, const char *name, mpz_t x) {                 /* as top3.c */
    char key[64];
    if (fscanf(f, "%63s %s", key, buf) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    if (mpz_set_str(x, buf, 10)) { fprintf(stderr, "bad int %s\n", name); exit(2); }
}
static long rdl(FILE *f, const char *name) {
    char key[64]; long v;
    if (fscanf(f, "%63s %ld", key, &v) != 2 || strcmp(key, name)) { fprintf(stderr, "bad input at %s\n", name); exit(2); }
    return v;
}
static void die(const char *m) { fprintf(stderr, "%s\n", m); exit(3); }

/* size classes: must match TP_CLASSES in topgpu.cu (L limbs; FL = 32 (L - 1)); larger Theta: topc3 */
#define NCLS 8
static const int CLS[NCLS] = {6, 8, 10, 12, 14, 16, 20, 22};
static const int CVAL[TP_NC] = {9, 18, 27, 45, 54, 63, 72};

static unsigned int *TB; static long long TW, TCAP;
static long long talloc(long long words) { long long o = TW; TW += words; if (TW > TCAP) { TCAP = 2 * TW + 4096; TB = realloc(TB, TCAP * 4); } memset(TB + o, 0, words * 4); return o; }
static void put_mpz(long long off, const mpz_t x, int limbs) {        /* 0 <= x < 2^(32 limbs) */
    size_t cnt = 0; unsigned int tmp[512];
    if (mpz_sgn(x) < 0 || mpz_sizeinbase(x, 2) > (size_t)32 * limbs) die("put_mpz: value out of range");
    memset(tmp, 0, sizeof tmp); mpz_export(tmp, &cnt, -1, 4, 0, 0, x);
    for (int i = 0; i < limbs; i++) TB[off + i] = tmp[i];
}

/* ---------------- top3.c's exact per-case code (state shared with main) ---------------- */
static long j_, prec_; static mpz_t C, qa, qb, K1, K2, D, a2, a3, P_, Q_, run_m;
static mpz_t m, t, r, d, lhs, rhs, x3, E3, N2, QD, r2, d2, base3, fl, L2a, L2b, L2c, R3a, R3b;
static long Kl; static mpfr_t *p10; static mpfr_t zp, A, LA, u, LU, v, tmp;
static unsigned long long cases, f1, f2, nsurv;
static void flush_run(long o, long *g, long j, long lo, long hi, mpz_t mm) {   /* as top3.c */
    printf("SURV %ld", o);
    for (long i = 0; i < j - 1; i++) printf(" %ld", g[i]);
    printf(" %ld %ld ", lo, hi); mpz_out_str(stdout, 10, mm); printf("\n");
}
/* top3.c's prefix computation (zp, A, LA) for gaps g[0..j-2] of ordering coefficients c; returns acc */
static long prefix_exact(const long *c, const long *g) {
    long acc = 0;
    mpfr_set_ui(zp, 0, MPFR_RNDN);
    for (long q = 0; q < j_ - 1; q++) { acc += g[q]; mpfr_mul_ui(tmp, p10[acc], c[q], MPFR_RNDN); mpfr_add(zp, zp, tmp, MPFR_RNDN); }
    mpfr_add_ui(A, zp, 1, MPFR_RNDN);
    mpfr_log1p(LA, zp, MPFR_RNDN);
    return acc;
}
/* top3.c's per-case body for last gap gl (after the prefix): level 1, level 2, run merging; counts f1, f2 (not cases) */
static long run_lo, run_hi;
static void case_exact(long o, const long *g, long acc, long cl, long gl) {
    long e = acc + gl;
    mpfr_mul_ui(u, p10[e], cl, MPFR_RNDN);
    mpfr_div(u, u, A, MPFR_RNDN);
    mpfr_log1p(LU, u, MPFR_RNDN);
    mpfr_add(v, LA, LU, MPFR_RNDN);
    mpfr_mul_z(v, v, C, MPFR_RNDN);
    mpfr_neg(v, v, MPFR_RNDN);
    mpfr_get_z(m, v, MPFR_RNDN);
    mpz_mul(t, m, qa); mpz_neg(t, t); mpz_fdiv_r(r, t, qb);
    mpz_sub(d, qb, r); if (mpz_cmp(r, d) < 0) mpz_set(d, r);
    mpz_mul(lhs, d, d); mpz_mul(lhs, lhs, K1);
    if (mpz_cmp(lhs, K2) > 0) return;
    f1++;
    int surv = 0;
    mpz_mul(base3, m, a3);
    mpz_neg(t, base3); mpz_fdiv_q(fl, t, D);
    for (long dx = -Kl; dx <= Kl + 1 && !surv; dx++) {
        if (dx >= 0) mpz_add_ui(x3, fl, (unsigned long)dx); else mpz_sub_ui(x3, fl, (unsigned long)(-dx));
        mpz_mul(E3, x3, D); mpz_add(E3, E3, base3);
        mpz_mul(lhs, E3, E3); mpz_mul(rhs, lhs, R3a);
        if (mpz_cmp(rhs, R3b) > 0) continue;
        mpz_mul(N2, m, a2); mpz_mul(N2, N2, Q_); mpz_neg(N2, N2);
        mpz_submul(N2, P_, E3);
        mpz_fdiv_r(r2, N2, QD); mpz_sub(d2, QD, r2); if (mpz_cmp(r2, d2) < 0) mpz_set(d2, r2);
        mpz_mul(rhs, d2, d2); mpz_mul(rhs, rhs, L2a);
        mpz_addmul(rhs, lhs, L2b);
        if (mpz_cmp(rhs, L2c) <= 0) surv = 1;
    }
    if (!surv) return;
    f2++;
    if (run_lo > 0 && mpz_cmp(m, run_m) == 0 && run_hi == gl - 1) { run_hi = gl; return; }
    if (run_lo > 0) { flush_run(o, (long *)g, j_, run_lo, run_hi, run_m); nsurv++; }
    run_lo = gl; run_hi = gl; mpz_set(run_m, m);
}

/* ---------------- batches ---------------- */
static int cmp_out(const void *A_, const void *B_) {
    const tp_out *a_ = (const tp_out *)A_, *b_ = (const tp_out *)B_;
    if (a_->prefix != b_->prefix) return a_->prefix < b_->prefix ? -1 : 1;
    if (a_->kind != b_->kind) return a_->kind > b_->kind ? -1 : 1;  /* a whole-prefix record first */
    return a_->gl < b_->gl ? -1 : (a_->gl > b_->gl);
}
static tp_out *OUT; static long long OUTCAP = 1 << 16; static long *CS; static int DUMP;
static void run_batch(const tp_prefix *pf, long n, long o) {
    unsigned long long cn[2] = {0, 0}; long long no;
    if (!OUT) OUT = malloc(sizeof(tp_out) * OUTCAP);
    while ((no = tp_run(pf, (int)n, OUT, OUTCAP, cn)) < 0) { OUTCAP *= 4; OUT = realloc(OUT, sizeof(tp_out) * OUTCAP); }
    cases += cn[0]; f1 += cn[1];
    qsort(OUT, (size_t)no, sizeof(tp_out), cmp_out);
    const long *c = CS + o * j_;
    long g[TP_MAXJ + 1];
    if (DUMP) {                                          /* validation: the kernel's psi1, psi2 for every case */
        for (long long i = 0; i < no; i++) {
            if (OUT[i].kind != 3) continue;
            const tp_prefix *q = &pf[OUT[i].prefix];
            printf("PSI %ld", o); for (long s = 0; s < j_ - 1; s++) printf(" %d", q->g[s]);
            printf(" %u %llu %llu\n", OUT[i].gl, OUT[i].psi1, OUT[i].psi2);
        }
        return;
    }
    for (long long i = 0; i < no;) {
        const unsigned p = OUT[i].prefix; const tp_prefix *q = &pf[p];
        for (long s = 0; s < j_ - 1; s++) g[s] = q->g[s];
        long acc = prefix_exact(c, g);
        run_lo = 0; run_hi = 0;
        if (OUT[i].kind == 2) {                          /* the whole prefix */
            for (long gl = q->glo; gl <= q->ghi; gl++) case_exact(o, g, acc, c[j_ - 1], gl);
            while (i < no && OUT[i].prefix == p) i++;
        } else {
            for (; i < no && OUT[i].prefix == p; i++) case_exact(o, g, acc, c[j_ - 1], (long)OUT[i].gl);
        }
        if (run_lo > 0) { flush_run(o, g, j_, run_lo, run_hi, run_m); nsurv++; run_lo = 0; }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: topgpu JOBFILE\n"); return 2; }
    const char *cpu = getenv("TP_CPU") ? getenv("TP_CPU") : "csrc/topc3";
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 2; }
    mpz_inits(C, qa, qb, K1, K2, D, a2, a3, P_, Q_, run_m, NULL);
    mpz_t g2n, g2d, g3n, g3d, Rn, Rd; mpz_inits(g2n, g2d, g3n, g3d, Rn, Rd, NULL);
    prec_ = rdl(f, "prec");
    rdz(f, "C", C);
    rdz(f, "qa", qa); rdz(f, "qb", qb); rdz(f, "K1", K1); rdz(f, "K2", K2);
    rdz(f, "D", D); rdz(f, "a2", a2); rdz(f, "a3", a3); rdz(f, "P", P_); rdz(f, "Q", Q_);
    rdz(f, "g2n", g2n); rdz(f, "g2d", g2d); rdz(f, "g3n", g3n); rdz(f, "g3d", g3d);
    rdz(f, "Rn", Rn); rdz(f, "Rd", Rd);
    long j = rdl(f, "j"); j_ = j;
    long cum = rdl(f, "cum");
    long *rng = malloc(sizeof(long) * (j + 1));
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, "ranges")) { fprintf(stderr, "bad ranges\n"); return 2; }
    for (long i = 0; i < j; i++) if (fscanf(f, "%ld", &rng[i]) != 1) return 2;
    long nord = rdl(f, "nord");
    long *cs = malloc(sizeof(long) * (nord * j + 1)); CS = cs;
    for (long o = 0; o < nord; o++) for (long i = 0; i < j; i++) if (fscanf(f, "%ld", &cs[o * j + i]) != 1) return 2;
    long skip = rdl(f, "skip");
    long maxsurv = rdl(f, "maxsurv"); (void)maxsurv;
    long debug = rdl(f, "debug");
    long g0lo = rdl(f, "g0lo"), g0hi = rdl(f, "g0hi");
    fclose(f);
    if (prec_ < (long)mpz_sizeinbase(C, 2) + 64) { fprintf(stderr, "prec too small\n"); return 2; }
    if (j < 1 || j > 40) { fprintf(stderr, "bad j\n"); return 2; }
    if (mpz_sgn(Q_) <= 0 || mpz_sgn(D) <= 0 || mpz_sgn(qb) <= 0) { fprintf(stderr, "bad denominators\n"); return 2; }
    long smax = 0; for (long i = 0; i < j; i++) if (rng[i] > smax) smax = rng[i];
    /* filter coverage; anything else is topc3's (exec) */
    int ok = cum == 1 && debug == 0 && j <= TP_MAXJ && smax < 590;
    int ci_[TP_NC * 64]; (void)ci_;
    for (long o = 0; o < nord && ok; o++) for (long i = 0; i < j; i++) { int found = 0; for (int c = 0; c < TP_NC; c++) found |= CVAL[c] == cs[o * j + i]; ok &= found; }
    /* |theta|, |theta2| < 2^-60 (theta = qa/qb, theta2 = a2/D) */
    if (ok) { mpz_t x; mpz_init(x); mpz_abs(x, qa); mpz_mul_2exp(x, x, 60); ok &= mpz_cmp(x, qb) < 0;
              mpz_abs(x, a2); mpz_mul_2exp(x, x, 60); ok &= mpz_cmp(x, D) < 0; mpz_clear(x); }
    if (!ok) { execv(cpu, argv); perror("exec topc3"); return 2; }

    /* top3.c's setup */
    mpfr_set_default_prec(prec_);
    long emax = 0; for (long i = 0; i < j; i++) emax += rng[i];
    p10 = malloc(sizeof(mpfr_t) * (emax + 2));
    for (long e = 0; e <= emax + 1; e++) { mpfr_init2(p10[e], prec_); mpfr_ui_pow_ui(p10[e], 10, e, MPFR_RNDN); mpfr_ui_div(p10[e], 1, p10[e], MPFR_RNDN); }
    mpfr_inits2(prec_, zp, A, LA, u, LU, v, tmp, (mpfr_ptr)0);
    mpz_inits(m, t, r, d, lhs, rhs, x3, E3, N2, QD, r2, d2, base3, fl, L2a, L2b, L2c, R3a, R3b, NULL);
    mpz_mul(L2a, g2n, g3d); mpz_mul(L2a, L2a, Rd);
    mpz_mul(L2b, Q_, Q_); mpz_mul(L2b, L2b, g3n); mpz_mul(L2b, L2b, g2d); mpz_mul(L2b, L2b, Rd);
    mpz_mul(L2c, Rn, g2d); mpz_mul(L2c, L2c, g3d); mpz_mul(L2c, L2c, Q_); mpz_mul(L2c, L2c, Q_); mpz_mul(L2c, L2c, D); mpz_mul(L2c, L2c, D);
    mpz_mul(R3a, g3n, Rd); mpz_mul(R3b, Rn, g3d); mpz_mul(R3b, R3b, D); mpz_mul(R3b, R3b, D);
    mpz_mul(QD, Q_, D);
    { mpz_t KK; mpz_init(KK); mpz_mul(t, Rn, g3d); mpz_mul(r, Rd, g3n); mpz_fdiv_q(t, t, r); mpz_sqrt(KK, t); mpz_add_ui(KK, KK, 1);
      Kl = mpz_get_si(KK); if (!mpz_fits_slong_p(KK) || Kl > 1000) { fprintf(stderr, "level-2 radius too large\n"); return 2; } mpz_clear(KK); }

    /* ---------------- filter tables ---------------- */
    tp_params Pp; memset(&Pp, 0, sizeof Pp);
    mpz_t Th, Th2, x, y; mpz_inits(Th, Th2, x, y, NULL);
    mpz_mul(Th, C, qa); mpz_tdiv_q(x, Th, qb); long tb = (long)mpz_sizeinbase(x, 2);          /* bits of |Theta| (integer part) */
    mpz_mul(Th2, C, a2); mpz_tdiv_q(y, Th2, D); long tb2 = (long)mpz_sizeinbase(y, 2); if (tb2 > tb) tb = tb2;
    int cls = -1; for (int c = 0; c < NCLS; c++) if (32 * (CLS[c] - 1) >= tb + 80) { cls = c; break; }
    if (cls < 0) { execv(cpu, argv); perror("exec topc3"); return 2; }
    const int L = CLS[cls], FL = 32 * (L - 1);
    Pp.j = (int)j; Pp.cls = cls; Pp.LW = L; Pp.LT = L + 1; Pp.FL = FL; Pp.FT = 64; Pp.Smin = 5; Pp.smax = (int)smax;
    for (long S = 0; S < 600; S++) {                     /* |w| <= 72 10^-S: (n + 1) (S log2 10 - log2 72) >= FL + 1 */
        double bits = S * 3.321928094887362 - 6.17;
        Pp.nterm[S] = bits <= 1 ? 0 : (int)((FL + 1) / bits) + 1;
    }
    int nmax = 1; for (long S = Pp.Smin; S <= smax; S++) if (Pp.nterm[S] > nmax) nmax = Pp.nterm[S];
    const long SS = smax + 1;
    Pp.o_cs = talloc(nord * ((j + 3) / 4) + 1);
    for (long o = 0; o < nord; o++) for (long i = 0; i < j; i++) {
        int ci = 0; for (int c = 0; c < TP_NC; c++) if (CVAL[c] == cs[o * j + i]) ci = c;
        TB[Pp.o_cs + o * ((j + 3) / 4) + (i >> 2)] |= (unsigned)ci << (8 * (i & 3));
    }
    Pp.o_lam1 = talloc((long long)TP_NC * SS * L); Pp.o_ia1 = talloc((long long)TP_NC * SS * L); Pp.o_c10 = talloc((long long)TP_NC * SS * L);
    Pp.o_inv = talloc((long long)(nmax + 1) * L);
    {
        mpfr_t a_, b_, s2; mpfr_inits2(FL + 128, a_, b_, s2, (mpfr_ptr)0); mpz_t z; mpz_init(z);
        mpfr_set_ui_2exp(s2, 1, FL, MPFR_RNDN);
        for (int c = 0; c < TP_NC; c++) for (long S = 0; S <= smax; S++) {
            const long long e_ = ((long long)c * SS + S) * L;
            mpfr_ui_pow_ui(a_, 10, S, MPFR_RNDN); mpfr_ui_div(a_, CVAL[c], a_, MPFR_RNDN);            /* c 10^-S */
            mpfr_mul(b_, a_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDD); put_mpz(Pp.o_c10 + e_, z, L);   /* floor */
            mpfr_log1p(b_, a_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_lam1 + e_, z, L);
            mpfr_add_ui(b_, a_, 1, MPFR_RNDN); mpfr_ui_div(b_, 1, b_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_ia1 + e_, z, L);
        }
        for (long n = 1; n <= nmax; n++) { mpz_set_ui(z, 0); mpz_setbit(z, FL); mpz_fdiv_q_ui(z, z, n); put_mpz(Pp.o_inv + n * L, z, L); }
        mpfr_clears(a_, b_, s2, (mpfr_ptr)0); mpz_clear(z);
    }
    /* Theta = C qa / qb, Theta2 = C a2 / D: word 0 = 1 if negative, then the low L + 1 limbs of round(|Theta| 2^64) */
    Pp.o_th = talloc(L + 2); Pp.o_th2 = talloc(L + 2);
    {
        mpz_t n_, dd, q_; mpz_inits(n_, dd, q_, NULL);
        for (int s = 0; s < 2; s++) {
            mpz_set(n_, s ? Th2 : Th); mpz_set(dd, s ? D : qb);   /* numerator C qa or C a2 */
            TB[s ? Pp.o_th2 : Pp.o_th] = mpz_sgn(n_) < 0 ? 1u : 0u;
            mpz_abs(n_, n_); mpz_mul_2exp(n_, n_, 65); mpz_add(n_, n_, dd); mpz_fdiv_q(q_, n_, dd); mpz_fdiv_q_2exp(q_, q_, 1);   /* round */
            mpz_fdiv_r_2exp(q_, q_, 32 * (L + 1));
            put_mpz((s ? Pp.o_th2 : Pp.o_th) + 1, q_, L + 1);
        }
        mpz_clears(n_, dd, q_, NULL);
    }
    /* margins (units 2^-64): eps + eta < 2^-59.5 (topgpu.cu header); take 64 units = 2^-58 */
    {
        const unsigned long long EM = 64;
        mpz_t a_, b_, s_; mpz_inits(a_, b_, s_, NULL);
        mpz_mul(a_, Rn, g3d); mpz_mul_2exp(a_, a_, 128); mpz_mul(b_, Rd, g3n);                /* rho^2 2^128 */
        mpz_fdiv_q(s_, a_, b_); mpz_sqrt(s_, s_);                                            /* floor(rho 2^64) */
        /* the kernel's level 2 takes the one x3 nearest -m theta, which needs rho < 1/4; constants must fit 64 bits */
        int fits = mpz_sizeinbase(s_, 2) <= 61;
        unsigned long long rl = fits ? mpz_get_ui(s_) : 0;
        Pp.rho_hi = rl + 1 + EM; Pp.rho_lo = rl > EM ? rl - EM : 0;
        Pp.ds = EM; Pp.dy = 4 * EM;                                                          /* (2 + |mu|) EM, |mu| <= 1/2, plus slack */
        mpz_mul_2exp(a_, P_, 63); mpz_fdiv_q(s_, a_, Q_);                                    /* mu 2^63 (floor) */
        fits &= mpz_fits_slong_p(s_);
        Pp.mu63 = fits ? mpz_get_si(s_) : 0;
        fits &= Pp.mu63 <= (1LL << 62) + 2 && Pp.mu63 >= -(1LL << 62) - 2;                  /* |mu| <= 1/2 (size-reduced basis) */
        mpz_mul(a_, g2n, Rd); mpz_mul_2exp(a_, a_, 32); mpz_mul(b_, g2d, Rn); mpz_fdiv_q(s_, a_, b_);   /* floor(g2/R2 2^32) */
        fits &= mpz_sizeinbase(s_, 2) <= 60; Pp.g2l = fits ? mpz_get_ui(s_) : 0;
        mpz_mul(a_, g3n, Rd); mpz_mul_2exp(a_, a_, 32); mpz_mul(b_, g3d, Rn); mpz_fdiv_q(s_, a_, b_);
        fits &= mpz_sizeinbase(s_, 2) <= 60; Pp.g3l = fits ? mpz_get_ui(s_) : 0;
        mpz_clears(a_, b_, s_, NULL);
        if (!fits) { execv(cpu, argv); perror("exec topc3"); return 2; }
    }
    Pp.words = TW;
    DUMP = Pp.dump = getenv("TP_DUMP") && atoi(getenv("TP_DUMP")) ? 1 : 0;
    if (DUMP) printf("TPDUMP FL %d L %d theta_bits %ld\n", FL, L, tb);
    if (tp_init(&Pp, TB)) die("tp_init");

    /* ---------------- enumeration (top3.c's order), batches within one ordering ---------------- */
    long batch = getenv("TP_BATCH") ? atol(getenv("TP_BATCH")) : (1L << 18);
    tp_prefix *pf = malloc(sizeof(tp_prefix) * batch); long np = 0;
    long *g = malloc(sizeof(long) * (j + 1)), *ss = malloc(sizeof(long) * (j + 1));
    for (long o = skip; o < nord; o++) {
        long i = 0; ss[0] = 0;
        if (j >= 2) g[0] = 0;
        for (;;) {
            if (j >= 2) {                                /* top3.c's prefix odometer (cum mode) */
                int found = 0;
                while (i >= 0) {
                    g[i]++;
                    if (i == 0 && g[0] < g0lo) g[0] = g0lo;
                    long hi = rng[i] - 1 - ss[i];
                    if (i == 0 && hi > g0hi) hi = g0hi;
                    if (g[i] > hi) { i--; continue; }
                    ss[i + 1] = ss[i] + g[i];
                    if (i == j - 2) { found = 1; break; }
                    i++; g[i] = 0;
                }
                if (!found) break;
            }
            long acc = 0; for (long q = 0; q < j - 1; q++) acc += g[q];
            long ghi = rng[j - 1] - 1 - acc, glo = 1;
            if (j == 1) { if (glo < g0lo) glo = g0lo; if (ghi > g0hi) ghi = g0hi; }
            if (glo <= ghi) {
                tp_prefix *q = &pf[np++]; memset(q, 0, sizeof *q);
                q->o = (unsigned)o; q->glo = (unsigned short)glo; q->ghi = (unsigned short)ghi;
                for (long s = 0; s < j - 1; s++) q->g[s] = (unsigned char)g[s];
                if (np == batch) { run_batch(pf, np, o); np = 0; }
            }
            if (j < 2) break;
        }
        if (np) { run_batch(pf, np, o); np = 0; }
        printf("PROG %ld %llu %llu %llu\n", o + 1, cases, f1, f2);
        fflush(stdout);
    }
    printf("DONE %llu %llu %llu %llu\n", cases, f1, f2, nsurv);
    return 0;
}
