/* bottomgpu: drop-in replacement for csrc/bottomc3 in 'cum' mode, with the per-solve loop on the GPU (bottomgpu.cu).
 * Same job file, same output (CAND lines in bottomc3's order and format, PROG per ordering, DONE nsolve ndirect5 ndirect2).
 * The host part includes bottom3.c unchanged: its parsing helpers, its mt() for the lattice data, and its low()/solutions()
 * for the rare solves the kernel hands back (fraction of a lattice estimate within 2^-60 of an integer, or v2(Flow) too
 * large). Jobs the kernel does not cover (another mode, heff[0] > G so that M varies, G above the largest size class) are
 * handed to bottomc3 itself (exec).
 * usage: bottomgpu JOBFILE      env BG_BATCH (flows per launch, default 2^20), BG_CPU (path of bottomc3)
 * build: see work/Makefile (bottomgpu: CUDA; bottomgpu_emu: the same kernel compiled for the CPU) */
#define main bottom3_main
#include "bottom3.c"
#undef main
#include <unistd.h>
#include "bottomgpu.h"

/* size classes: must match BG_CLASSES in bottomgpu.cu (Gmax, L5, L2, LQ, LH). L5: 5^(G + s) < 2^(32 L5) with >= 8 bits to
 * spare; L2: 32 L2 >= G + s + 2 and tlmax = 32 L2 - G - 8 >= 40; LQ: Q = 2^(G-2) 5^(G-1) < 2^(32 LQ - 20); LH: the coordinates
 * (checked per job by widths_ok); 2 LH <= BG_OUTW */
#define NBCLS 9
static const int CLS[NBCLS][5] = {{60, 6, 5, 10, 7}, {90, 8, 6, 13, 9}, {125, 11, 7, 17, 11}, {165, 13, 9, 21, 13}, {200, 16, 10, 24, 14},
                                  {240, 19, 10, 26, 15}, {280, 22, 11, 30, 16}, {320, 25, 12, 34, 18}, {360, 28, 13, 38, 20}};
static const int CVAL[BG_NC] = {9, 18, 27, 45, 54, 63, 72};

static unsigned int *TB; static long long TW, TCAP;
static long long talloc(long long words) { long long o = TW; TW += words; if (TW > TCAP) { TCAP = 2 * TW + 4096; TB = realloc(TB, TCAP * 4); } memset(TB + o, 0, words * 4); return o; }
static void put_mpz(long long off, const mpz_t x, int limbs) {   /* nonnegative x < 2^(32 limbs), little-endian u32 limbs */
    size_t cnt = 0; unsigned int tmp[256];
    if (mpz_sgn(x) < 0 || mpz_sizeinbase(x, 2) > (size_t)32 * limbs) die("put_mpz: value out of range");
    memset(tmp, 0, sizeof tmp); mpz_export(tmp, &cnt, -1, 4, 0, 0, x);
    for (int i = 0; i < limbs; i++) TB[off + i] = tmp[i];
}
static void put_signed(long long off, const mpz_t x, int limbs) {   /* two's complement in `limbs` limbs */
    mpz_t y; mpz_init(y);
    if (mpz_sgn(x) < 0) { mpz_ui_pow_ui(y, 2, 32 * limbs); mpz_add(y, y, x); if (mpz_sizeinbase(x, 2) >= (size_t)32 * limbs - 1) die("put_signed: range"); }
    else { mpz_set(y, x); if (mpz_sizeinbase(x, 2) >= (size_t)32 * limbs - 1) die("put_signed: range"); }
    put_mpz(off, y, limbs); mpz_clear(y);
}
static void round_div(mpz_t q, const mpz_t a, const mpz_t b) {   /* nearest integer to a / b (b != 0): floor((2a + b) / 2b) */
    mpz_t n, d; mpz_inits(n, d, NULL);
    mpz_set(n, a); mpz_set(d, b);
    if (mpz_sgn(d) < 0) { mpz_neg(n, n); mpz_neg(d, d); }
    mpz_mul_2exp(n, n, 1); mpz_add(n, n, d); mpz_mul_2exp(d, d, 1);
    mpz_fdiv_q(q, n, d); mpz_clears(n, d, NULL);
}

/* ---------------- output of one batch of flows (all of ordering o) in bottom3.c's order ---------------- */
static long B_b, B_G, *B_cs; static bg_out *B_out; static long long B_outcap = 1 << 14; static unsigned long long B_cnt[3];
static int cmp_out(const void *A, const void *B) {   /* flow, p; within a solve host records first, then seq */
    const bg_out *u = (const bg_out *)A, *v = (const bg_out *)B;
    if (u->flow != v->flow) return u->flow < v->flow ? -1 : 1;
    if (u->p != v->p) return u->p < v->p ? -1 : 1;
    if ((u->kind != 1) != (v->kind != 1)) return u->kind != 1 ? -1 : 1;
    return u->seq < v->seq ? -1 : (u->seq > v->seq);
}
static void run_flows(const bg_flow *fls, long n, long o) {
    long long no;
    if (!B_out) B_out = malloc(sizeof(bg_out) * B_outcap);
    while ((no = bg_run(fls, (int)n, B_out, B_outcap, B_cnt)) < 0) { B_outcap *= 4; B_out = realloc(B_out, sizeof(bg_out) * B_outcap); }
    qsort(B_out, (size_t)no, sizeof(bg_out), cmp_out);
    mpz_t vu, vw, Fl; mpz_inits(vu, vw, Fl, NULL);
    const long b = B_b;
    for (long long i = 0; i < no; i++) {
        const bg_out *e = &B_out[i];
        const bg_flow *q = &fls[e->flow];
        long pp[BG_MAXB + 1]; pp[0] = 0; for (long s = 1; s < b; s++) pp[s] = q->pos[s];
        if (e->kind != 1) {                              /* bottom3.c redoes the whole solve and prints its CAND lines */
            mpz_set_si(Fl, -1);
            for (long s = 1; s < b; s++) mpz_addmul_ui(Fl, pow10[pp[s]], (unsigned long)B_cs[o * b + s]);
            low(Fl);
            long t_, ns0 = nsolve, n50 = ndirect5, n20 = ndirect2;
            report_all = 1; rep_o = o; rep_pos = pp; rep_nb = b;
            solutions(Fl, B_cs[o * b], (long)e->p, B_G, &t_);
            fflush(stdout);
            if (e->kind == 3) { B_cnt[0] += (unsigned long long)(nsolve - ns0); B_cnt[1] += (unsigned long long)(ndirect5 - n50); B_cnt[2] += (unsigned long long)(ndirect2 - n20); }
            while (i + 1 < no && B_out[i + 1].flow == e->flow && B_out[i + 1].p == e->p) i++;   /* the solve's other records */
            continue;
        }
        mpz_import(vu, BG_OUTW, -1, 4, 0, 0, e->up); mpz_import(vw, BG_OUTW, -1, 4, 0, 0, e->w);
        printf("CAND %ld pos", o); for (long s = 1; s < b; s++) printf(" %ld", pp[s]);
        printf(" p %u M %ld t %u up ", e->p, B_G, e->t); mpz_out_str(stdout, 10, vu); printf(" w "); mpz_out_str(stdout, 10, vw); printf("\n");
    }
    mpz_clears(vu, vw, Fl, NULL);
}

/* bottom3.c's initialisation (its main), so that low(), solutions() and mt() work */
static void b3_init(void) {
    mpz_mul(diam2, U, U); mpz_addmul(diam2, W, W);
    long np = Mmax + 50;
    pow5 = malloc(sizeof(mpz_t) * np); pow2 = malloc(sizeof(mpz_t) * np); pow10 = malloc(sizeof(mpz_t) * np);
    for (long i = 0; i < np; i++) { mpz_init(pow5[i]); mpz_ui_pow_ui(pow5[i], 5, i); mpz_init(pow2[i]); mpz_ui_pow_ui(pow2[i], 2, i); mpz_init(pow10[i]); mpz_ui_pow_ui(pow10[i], 10, i); }
    for (long n = 1; n < NINV; n++) {
        long nn = n; while (nn % 5 == 0) nn /= 5;
        mpz_init_set_ui(inv5n[n], nn); if (!mpz_invert(inv5n[n], inv5n[n], pow5[Mmax + 12])) die("inv5n");
        nn = n; while (nn % 2 == 0) nn /= 2;
        mpz_init_set_ui(inv2n[n], nn); if (!mpz_invert(inv2n[n], inv2n[n], pow2[Mmax + 44])) die("inv2n");
    }
    mtT = tm; mtM = Mmax; mtc = calloc((Mmax + 1) * (tm + 1), sizeof(MT));
    mpz_t chk; mpz_init(chk);
    mpz_set_ui(T1, 81); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_3)) die("L5_3 mismatch");
    mpz_set_ui(T1, 2401); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_7)) die("L5_7 mismatch");
    mpz_set_ui(T1, 16); plog(chk, T1, 5, Mmax); if (mpz_cmp(chk, L5_16)) die("L5_16 mismatch");
    mpz_set_si(T1, -3); mpz_fdiv_r_2exp(T1, T1, Mmax + 2); plog(chk, T1, 2, Mmax); if (mpz_cmp(chk, L2_3)) die("L2_3 mismatch");
    mpz_set_si(T1, -7); mpz_fdiv_r_2exp(T1, T1, Mmax + 2); plog(chk, T1, 2, Mmax); if (mpz_cmp(chk, L2_7)) die("L2_7 mismatch");
    mpz_clear(chk);
}

static long isqrt_l(long x) { long s = 0; while ((s + 1) * (s + 1) <= x) s++; return s; }

/* would bottom3.c's mt(M, t) succeed? (the same reduction; "lattice too fine" is reported instead of stopping). A wrong
 * answer cannot change a decision: "no" sends the t to bottom3.c's own solutions(), "yes" calls mt() itself. */
static int mt_ok(long M, long t) {
    long n = M - t; if (n < 6) return 0;
    mpz_t q2, q5, Q, A, Bc, A5, B5, Ai, A5i, al2, al5, e2, e5, al, t1, b1x, b1y, b2x, b2y, n1, n2, num, den, kk, tmp;
    mpz_inits(q2, q5, Q, A, Bc, A5, B5, Ai, A5i, al2, al5, e2, e5, al, t1, b1x, b1y, b2x, b2y, n1, n2, num, den, kk, tmp, NULL);
    mpz_ui_pow_ui(q2, 2, n - 2); mpz_ui_pow_ui(q5, 5, M - 1); mpz_mul(Q, q2, q5);
    mpz_fdiv_r_2exp(A, L2_3, n); mpz_fdiv_q_2exp(A, A, 2); mpz_fdiv_r_2exp(Bc, L2_7, n); mpz_fdiv_q_2exp(Bc, Bc, 2);
    mpz_ui_pow_ui(tmp, 5, M); mpz_mod(A5, L5_3, tmp); mpz_divexact_ui(A5, A5, 5); mpz_mod(B5, L5_7, tmp); mpz_divexact_ui(B5, B5, 5);
    int ok = mpz_invert(Ai, A, q2) && mpz_invert(A5i, A5, q5);
    if (ok) {
        mpz_mul(al2, Bc, Ai); mpz_neg(al2, al2); mpz_mod(al2, al2, q2);
        mpz_mul(al5, B5, A5i); mpz_neg(al5, al5); mpz_mod(al5, al5, q5);
        mpz_invert(t1, q5, q2); mpz_mul(e2, q5, t1); mpz_invert(t1, q2, q5); mpz_mul(e5, q2, t1);
        mpz_mul(al, al2, e2); mpz_addmul(al, al5, e5); mpz_mod(al, al, Q);
        mpz_mul_ui(b1x, al, 4); mpz_mod(b1x, b1x, Q); mpz_set_ui(b1y, 4); mpz_set(b2x, Q); mpz_set_ui(b2y, 0);
        for (;;) {
            mpz_mul(n1, b1x, b1x); mpz_addmul(n1, b1y, b1y); mpz_mul(n2, b2x, b2x); mpz_addmul(n2, b2y, b2y);
            if (mpz_cmp(n1, n2) > 0) { mpz_swap(b1x, b2x); mpz_swap(b1y, b2y); mpz_swap(n1, n2); }
            mpz_mul(num, b1x, b2x); mpz_addmul(num, b1y, b2y); mpz_set(den, n1);
            mpz_mul_ui(tmp, num, 2); mpz_add(tmp, tmp, den); mpz_mul_ui(den, den, 2); mpz_fdiv_q(kk, tmp, den);
            if (mpz_sgn(kk) == 0) break;
            mpz_submul(b2x, kk, b1x); mpz_submul(b2y, kk, b1y);
        }
        ok = mpz_cmp(n1, diam2) > 0;
    }
    mpz_clears(q2, q5, Q, A, Bc, A5, B5, Ai, A5i, al2, al5, e2, e5, al, t1, b1x, b1y, b2x, b2y, n1, n2, num, den, kk, tmp, NULL);
    return ok;
}
static long ivp(long n, long p) { long e = 0; while (n % p == 0) { n /= p; e++; } return e; }

/* the kernel's fixed widths cannot overflow for the lattice data of one t: the estimates K - off0 R (2 LQ limbs, off0 < Q),
 * the coordinates ci, cj (LH limbs, signed) and up, ww (2 LH limbs, signed). With |K| the largest over w0:
 * |ci| <= (|K_i| + Q |R_i|) / 2^F + 1 (cj alike), plus 1 for the neighbours; |up| <= Q + (|ci| + 1) |b1x| + (|cj| + 1) |b2x|,
 * |ww| <= 4 + (|ci| + 1) |b1y| + (|cj| + 1) |b2y|. A failure sends the job to bottomc3. */
static int widths_ok(const MT *r, const mpz_t Ri, const mpz_t Rj, const mpz_t Ki, const mpz_t Kj, int F, int LQ, int LH) {
    mpz_t e, ci_, cj_, s, lim; mpz_inits(e, ci_, cj_, s, lim, NULL);
    int ok = 1;
    mpz_ui_pow_ui(lim, 2, 64 * LQ - 2);
    mpz_abs(e, Ri); mpz_mul(e, e, r->Q); mpz_abs(s, Ki); mpz_add(e, e, s); ok &= mpz_cmp(e, lim) < 0;   /* |K_i| + Q |R_i| */
    mpz_fdiv_q_2exp(ci_, e, F); mpz_add_ui(ci_, ci_, 2);
    mpz_abs(e, Rj); mpz_mul(e, e, r->Q); mpz_abs(s, Kj); mpz_add(e, e, s); ok &= mpz_cmp(e, lim) < 0;
    mpz_fdiv_q_2exp(cj_, e, F); mpz_add_ui(cj_, cj_, 2);
    mpz_ui_pow_ui(lim, 2, 32 * LH - 2); ok &= mpz_cmp(ci_, lim) < 0 && mpz_cmp(cj_, lim) < 0;
    mpz_ui_pow_ui(lim, 2, 64 * LH - 2);
    mpz_abs(e, r->b1x); mpz_mul(s, ci_, e); mpz_abs(e, r->b2x); mpz_addmul(s, cj_, e); mpz_add(s, s, r->Q); ok &= mpz_cmp(s, lim) < 0;
    mpz_abs(e, r->b1y); mpz_mul(s, ci_, e); mpz_abs(e, r->b2y); mpz_addmul(s, cj_, e); mpz_add_ui(s, s, 4); ok &= mpz_cmp(s, lim) < 0;
    mpz_clears(e, ci_, cj_, s, lim, NULL);
    return ok;
}

/* coefficient (-1)^(n+1) p^k / n' mod p^Mmod (k >= 0), as a residue; zero if k >= Mmod */
static void series_coef(mpz_t out, long n, unsigned long p, long k, long Mmod) {
    mpz_t m, t; mpz_inits(m, t, NULL); mpz_ui_pow_ui(m, p, Mmod);
    if (k >= Mmod) { mpz_set_ui(out, 0); mpz_clears(m, t, NULL); return; }
    long nn = n; while (nn % (long)p == 0) nn /= p;
    mpz_set_ui(t, nn); if (!mpz_invert(t, t, m)) die("series_coef");
    mpz_ui_pow_ui(out, p, k); mpz_mul(out, out, t); mpz_mod(out, out, m);
    if (!(n & 1)) { mpz_sub(out, m, out); mpz_mod(out, out, m); }
    mpz_clears(m, t, NULL);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: bottomgpu JOBFILE\n"); return 2; }
    const char *cpu = getenv("BG_CPU") ? getenv("BG_CPU") : "csrc/bottomc3";
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 2; }
    mpz_inits(U, W, diam2, L5_3, L5_7, L5_16, L2_3, L2_7, T1, T2, T3, T4, T5, T6, NULL);
    mpz_inits(l2low, l5low, inv5, inv2, oddlow, F, x5, l5F, x2, l2F, D5, D2, beta5, beta2, beta, off0, cx2, cy2, ci, cj, up, ww, tq, NULL);
    char mode[32]; if (fscanf(f, "%31s", mode) != 1) return 2;
    const int ktest = !strcmp(mode, "ktest");          /* kernel test: header, then "n NQ" and NQ lines "Flow c p" */
    if (strcmp(mode, "cum") && !ktest) { fclose(f); execv(cpu, argv); perror("exec bottomc3"); return 2; }
    long G = rdl(f, "G"); tm = rdl(f, "tm"); rdz(f, "U", U); rdz(f, "W", W); Mmax = rdl(f, "Mmax");
    rdz(f, "L5_3", L5_3); rdz(f, "L5_7", L5_7); rdz(f, "L5_16", L5_16); rdz(f, "L2_3", L2_3); rdz(f, "L2_7", L2_7);
    long k = 1, j = 0, b = 1, heff[BG_MAXB + 1], nord = 0, *cs = NULL, skip = 0, nq = 0, *qc = NULL, *qp = NULL;
    mpz_t *qflow = NULL;
    if (!ktest) {
        k = rdl(f, "k"); j = rdl(f, "j"); b = k - j;
        char key2[64]; if (fscanf(f, "%63s", key2) != 1 || strcmp(key2, "heff")) die("heff");
        if (b < 1 || b > BG_MAXB) { fclose(f); execv(cpu, argv); perror("exec bottomc3"); return 2; }
        for (long r = 0; r < b; r++) if (fscanf(f, "%ld", &heff[r]) != 1) die("heff vals");
        nord = rdl(f, "nord");
        cs = malloc(sizeof(long) * (nord * b + 1));
        for (long o = 0; o < nord; o++) for (long i = 0; i < b; i++) if (fscanf(f, "%ld", &cs[o * b + i]) != 1) die("ord");
        skip = rdl(f, "skip"); aklo = rdl(f, "aklo"); akhi = rdl(f, "akhi");
    } else {
        heff[0] = G; nq = rdl(f, "n");
        qflow = malloc(sizeof(mpz_t) * (nq + 1)); qc = malloc(sizeof(long) * (nq + 1)); qp = malloc(sizeof(long) * (nq + 1));
        for (long i = 0; i < nq; i++) {
            if (fscanf(f, "%s %ld %ld", buf, &qc[i], &qp[i]) != 3) die("ktest query");
            mpz_init_set_str(qflow[i], buf, 10);
            if (qp[i] < 0 || qp[i] >= G) die("ktest: need 0 <= p < G");
        }
    }
    fclose(f);
    long Mneed = heff[0] > G ? heff[0] : G;
    if (Mneed > Mmax) die("Mmax too small for this configuration space");
    int cls = -1; for (int c = 0; c < NBCLS; c++) if (G <= CLS[c][0]) { cls = c; break; }
    /* kernel coverage: M = G for every solve (heff[0] <= G), a size class (G <= 360: positions fit 16 bits), n = G - t >= 8 */
    if (heff[0] > G || cls < 0 || tm + 8 > G) {
        if (ktest) die("ktest: G outside the kernel's size classes");
        execv(cpu, argv); perror("exec bottomc3"); return 2;
    }
    b3_init();
    cG = G; cb = b;

    /* ---------------- tables ---------------- */
    bg_params P; memset(&P, 0, sizeof P);
    const int L5 = CLS[cls][1], L2 = CLS[cls][2], LQ = CLS[cls][3], LH = CLS[cls][4];
    P.G = (int)G; P.tm = (int)tm; P.b = (int)b; P.cls = cls; P.L5 = L5; P.L2 = L2; P.LQ = LQ; P.LH = LH;
    P.s5 = (int)isqrt_l(G) + 2; P.vz5 = P.s5 + 1; P.s2 = (int)isqrt_l(G) + 2; P.vz2 = P.s2 + 2;
    P.tlmax = 32 * L2 - (int)G - 8;
    for (int c = 0; c < BG_NC; c++) { P.cval[c] = CVAL[c]; P.v2c[c] = (int)ivp(CVAL[c], 2); }
    mpz_t mN, mG, mQ, R, Rq, x, y, z, w2, m2; mpz_inits(mN, mG, mQ, R, Rq, x, y, z, w2, m2, NULL);
    const long N = G + P.s5;
    mpz_ui_pow_ui(mN, 5, N); mpz_ui_pow_ui(mG, 5, G); mpz_ui_pow_ui(mQ, 5, G - 1);
    mpz_ui_pow_ui(R, 2, 32 * L5); mpz_ui_pow_ui(m2, 2, 32 * L2);
    if (mpz_cmp(mN, R) >= 0) die("size class: 5^N >= R");
    mpz_ui_pow_ui(x, 2, 32);
    mpz_t *mods[3] = {&mN, &mG, &mQ}; unsigned int *mps[3] = {&P.mpN, &P.mpG, &P.mpQ};
    long long *om[3] = {&P.o_mN, &P.o_mG, &P.o_mQ}, *oR2[3] = {&P.o_R2N, &P.o_R2G, &P.o_R2Q};
    for (int i = 0; i < 3; i++) {
        mpz_invert(y, *mods[i], x); mpz_sub(y, x, y); *mps[i] = (unsigned int)mpz_get_ui(y);      /* -m^-1 mod 2^32 */
        *om[i] = talloc(L5); put_mpz(*om[i], *mods[i], L5);
        mpz_mul(y, R, R); mpz_mod(y, y, *mods[i]); *oR2[i] = talloc(L5); put_mpz(*oR2[i], y, L5);
    }
    mpz_mul_ui(y, R, 2); mpz_mod(y, y, mG); P.o_twoRG = talloc(L5); put_mpz(P.o_twoRG, y, L5);
    mpz_set_ui(y, 5); mpz_invert(y, y, R); P.o_inv5 = talloc(L5); put_mpz(P.o_inv5, y, L5);
    mpz_ui_pow_ui(y, 5, P.vz5); mpz_invert(y, y, R); P.o_invz5 = talloc(L5); put_mpz(P.o_invz5, y, L5);
    P.o_T10N = talloc((long long)BG_NC * G * L5); P.o_T10G = talloc((long long)BG_NC * G * L5); P.o_T10_2 = talloc((long long)BG_NC * G * L2);
    P.o_Tc5 = talloc((long long)BG_NC * G * L5); P.o_Tc2 = talloc((long long)BG_NC * G * L2);
    for (int c = 0; c < BG_NC; c++) for (long pos = 0; pos < G; pos++) {
        long long e = (long long)c * G + pos;
        mpz_ui_pow_ui(y, 10, pos); mpz_mul_ui(y, y, CVAL[c]);
        mpz_mod(z, y, mN); put_mpz(P.o_T10N + e * L5, z, L5);
        mpz_mod(z, y, mG); put_mpz(P.o_T10G + e * L5, z, L5);
        mpz_mod(z, y, m2); put_mpz(P.o_T10_2 + e * L2, z, L2);
        mpz_ui_pow_ui(y, 2, pos); mpz_mul_ui(y, y, CVAL[c]); mpz_mul(y, y, R); mpz_mod(z, y, mG); put_mpz(P.o_Tc5 + e * L5, z, L5);
        long codd = CVAL[c]; while (!(codd & 1)) codd >>= 1;
        mpz_ui_pow_ui(y, 5, pos); mpz_mul_ui(y, y, codd); mpz_mod(z, y, m2); put_mpz(P.o_Tc2 + e * L2, z, L2);
    }
    /* series coefficients: log1p in 5^G (v = p = 1 .. G-1) and in 2^G (v = 2 .. G + 4); plog in 5^G and 2^G */
    P.o_n5 = talloc(G + 1); P.o_off5 = talloc(G + 1);
    long tot5 = 0; for (long v = 1; v < G; v++) { long n = 1; while (n * v - floorlog(n, 5) < G) n++; TB[P.o_n5 + v] = (unsigned)(n - 1); TB[P.o_off5 + v] = (unsigned)tot5; tot5 += n - 1; }
    P.o_coef5 = talloc(tot5 * L5);
    for (long v = 1; v < G; v++) for (long n = 1; n <= (long)TB[P.o_n5 + v]; n++) {
        series_coef(z, n, 5, n * v - ivp(n, 5), G); put_mpz(P.o_coef5 + (TB[P.o_off5 + v] + n - 1) * (long long)L5, z, L5);
    }
    const long V2MAX = G + 4;
    P.o_n2 = talloc(V2MAX + 1); P.o_off2 = talloc(V2MAX + 1);
    long tot2 = 0; for (long v = 2; v <= V2MAX; v++) { long n = 1; while (n * v - floorlog(n, 2) < G) n++; TB[P.o_n2 + v] = (unsigned)(n - 1); TB[P.o_off2 + v] = (unsigned)tot2; tot2 += n - 1; }
    P.o_coef2 = talloc(tot2 * L2 + L2);
    for (long v = 2; v <= V2MAX; v++) for (long n = 1; n <= (long)TB[P.o_n2 + v]; n++) {
        series_coef(z, n, 2, n * v - ivp(n, 2), G); put_mpz(P.o_coef2 + (TB[P.o_off2 + v] + n - 1) * (long long)L2, z, L2);
    }
    { long n = 1; while (n * P.vz5 - floorlog(n, 5) - P.s5 < G) n++; P.np5 = (int)(n - 1); }
    { long n = 1; while (n * P.vz2 - floorlog(n, 2) - P.s2 < G) n++; P.np2 = (int)(n - 1); }
    P.o_plog5 = talloc((long long)P.np5 * L5); P.o_plog2 = talloc((long long)P.np2 * L2);
    for (long n = 1; n <= P.np5; n++) { series_coef(z, n, 5, n * P.vz5 - ivp(n, 5) - P.s5, G); put_mpz(P.o_plog5 + (n - 1) * L5, z, L5); }
    for (long n = 1; n <= P.np2; n++) { series_coef(z, n, 2, n * P.vz2 - ivp(n, 2) - P.s2, G); put_mpz(P.o_plog2 + (n - 1) * L2, z, L2); }
    P.o_L16 = talloc((tm + 1) * (long long)L5);
    for (long t = 0; t <= tm; t++) { mpz_mul_ui(y, L5_16, t); mpz_mod(y, y, mG); put_mpz(P.o_L16 + t * L5, y, L5); }
    /* lattice data per t (bottom3.c's mt(G, t)); fixed point with F = bits(Q at t = 0) + 64 fraction bits */
    P.o_Ainv = talloc((tm + 1) * (long long)L2); P.o_A5inv = talloc((tm + 1) * (long long)L5); P.o_q2inv5 = talloc((tm + 1) * (long long)L5);
    P.o_a2m4 = talloc(tm + 1); P.o_alphaW = talloc((tm + 1) * 4LL * LQ); P.o_Q = talloc((tm + 1) * (long long)LQ);
    P.o_b = talloc((tm + 1) * 4LL * LH); P.o_Ri = talloc((tm + 1) * (long long)LQ); P.o_Rj = talloc((tm + 1) * (long long)LQ);
    P.o_Ki = talloc((tm + 1) * 4LL * 2 * LQ); P.o_Kj = talloc((tm + 1) * 4LL * 2 * LQ);
    P.o_U = talloc(LQ); P.o_W = talloc(LQ); put_mpz(P.o_U, U, LQ); put_mpz(P.o_W, W, LQ);
    P.o_tok = talloc(tm + 1);
    if (!mt_ok(G, 0)) { execv(cpu, argv); perror("exec bottomc3"); return 2; }
    { MT *r0 = mt(G, 0); P.F = (int)mpz_sizeinbase(r0->Q, 2) + 64; }
    if (P.F + 64 > 64 * LQ) die("size class: F too large");
    mpz_t sc, num, den; mpz_inits(sc, num, den, NULL); mpz_ui_pow_ui(sc, 2, P.F);
    mpz_t vRi, vRj, vK[8], Kim, Kjm; mpz_inits(vRi, vRj, Kim, Kjm, NULL); for (int i = 0; i < 8; i++) mpz_init(vK[i]);
    for (long t = 0; t <= tm; t++) {
        if (!mt_ok(G, t)) continue;                      /* o_tok[t] = 0: such solves go to bottom3.c */
        TB[P.o_tok + t] = 1;
        MT *r = mt(G, t);
        /* R_i = b2y / det, R_j = -b1y / det; K_i = (U b2y - (W - 2 w0) b2x) / (2 det), K_j = (b1x (W - 2 w0) - b1y U) / (2 det),
         * in fixed point with F fraction bits; first checked against the kernel's widths (else the job goes to bottomc3) */
        mpz_mul(num, r->b2y, sc); round_div(vRi, num, r->det);
        mpz_mul(num, r->b1y, sc); mpz_neg(num, num); round_div(vRj, num, r->det);
        mpz_mul_2exp(den, r->det, 1); mpz_set_ui(Kim, 0); mpz_set_ui(Kjm, 0);
        for (int w0 = 0; w0 < 4; w0++) {
            mpz_sub_ui(w2, W, 2 * w0);
            mpz_mul(num, U, r->b2y); mpz_submul(num, w2, r->b2x); mpz_mul(num, num, sc); round_div(vK[2 * w0], num, den);
            mpz_mul(num, r->b1x, w2); mpz_submul(num, r->b1y, U); mpz_mul(num, num, sc); round_div(vK[2 * w0 + 1], num, den);
            mpz_abs(y, vK[2 * w0]); if (mpz_cmp(y, Kim) > 0) mpz_set(Kim, y);
            mpz_abs(y, vK[2 * w0 + 1]); if (mpz_cmp(y, Kjm) > 0) mpz_set(Kjm, y);
        }
        int wok = widths_ok(r, vRi, vRj, Kim, Kjm, P.F, LQ, LH);
        const mpz_ptr bv[4] = {r->b1x, r->b1y, r->b2x, r->b2y};
        for (int i = 0; i < 4; i++) wok &= mpz_sizeinbase(bv[i], 2) < (size_t)32 * LH - 1;         /* put_signed's ranges */
        wok &= mpz_sizeinbase(vRi, 2) < (size_t)32 * LQ - 1 && mpz_sizeinbase(vRj, 2) < (size_t)32 * LQ - 1;
        if (!wok) { execv(cpu, argv); perror("exec bottomc3"); return 2; }
        put_mpz(P.o_Ainv + t * L2, r->Ainv, L2);
        mpz_mul(y, r->A5inv, R); mpz_mod(y, y, r->q5); put_mpz(P.o_A5inv + t * L5, y, L5);
        mpz_invert(y, r->q2, r->q5); mpz_mul(y, y, R); mpz_mod(y, y, r->q5); put_mpz(P.o_q2inv5 + t * L5, y, L5);
        if (mpz_cmp(r->q5, mQ)) die("q5 != 5^(G-1)");
        TB[P.o_a2m4 + t] = (unsigned)mpz_fdiv_ui(r->alpha2, 4);
        put_mpz(P.o_Q + t * LQ, r->Q, LQ);
        for (int w0 = 0; w0 < 4; w0++) { mpz_mul_ui(y, r->alpha, w0); mpz_mod(y, y, r->Q); put_mpz(P.o_alphaW + (t * 4 + w0) * LQ, y, LQ); }
        put_signed(P.o_b + (t * 4 + 0) * LH, r->b1x, LH); put_signed(P.o_b + (t * 4 + 1) * LH, r->b1y, LH);
        put_signed(P.o_b + (t * 4 + 2) * LH, r->b2x, LH); put_signed(P.o_b + (t * 4 + 3) * LH, r->b2y, LH);
        put_signed(P.o_Ri + t * LQ, vRi, LQ); put_signed(P.o_Rj + t * LQ, vRj, LQ);
        for (int w0 = 0; w0 < 4; w0++) {
            put_signed(P.o_Ki + (t * 4 + w0) * 2LL * LQ, vK[2 * w0], 2 * LQ);
            put_signed(P.o_Kj + (t * 4 + w0) * 2LL * LQ, vK[2 * w0 + 1], 2 * LQ);
        }
    }
    P.words = TW;
    /* orderings as coefficient indices */
    unsigned char *oci = malloc((size_t)nord * b + 1);
    for (long o = 0; o < nord; o++) for (long i = 0; i < b; i++) {
        int ci_ = -1; for (int c = 0; c < BG_NC; c++) if (CVAL[c] == cs[o * b + i]) ci_ = c;
        if (ci_ < 0) die("digit coefficient not in 9 18 27 45 54 63 72");
        oci[o * b + i] = (unsigned char)ci_;
    }
    if (bg_init(&P, TB, oci, (int)nord)) die("bg_init");
    if (ktest) {                                         /* records: "Q i t up w" (kernel), host decisions as bottom3.c's CAND lines */
        bg_query *qs = calloc(nq + 1, sizeof(bg_query)); mpz_t r_; mpz_init(r_);
        for (long i = 0; i < nq; i++) {
            int ci_ = -1; for (int c = 0; c < BG_NC; c++) if (CVAL[c] == qc[i]) ci_ = c;
            if (ci_ < 0) die("ktest: c must be a digit coefficient");
            qs[i].ci = ci_; qs[i].p = (int)qp[i]; qs[i].f5 = (int)mpz_fdiv_ui(qflow[i], 5);
            if (mpz_sgn(qflow[i]) == 0) die("ktest: Flow = 0");
            unsigned int tmp[256]; size_t cnt_; mpz_t *ms[3] = {&mN, &mG, &m2}; unsigned int *dst[3] = {qs[i].f5n, qs[i].f5g, qs[i].f2};
            const int lim[3] = {L5, L5, L2};
            for (int s = 0; s < 3; s++) {
                mpz_mod(r_, qflow[i], *ms[s]); memset(tmp, 0, sizeof tmp); mpz_export(tmp, &cnt_, -1, 4, 0, 0, r_);
                for (int w = 0; w < lim[s]; w++) dst[s][w] = tmp[w];
            }
        }
        long long cap = 1 << 16; bg_out *o2 = malloc(sizeof(bg_out) * cap); unsigned long long cn[3] = {0, 0, 0}; long long no;
        while ((no = bg_runq(qs, (int)nq, o2, cap, cn)) < 0) { cap *= 4; o2 = realloc(o2, sizeof(bg_out) * cap); }
        qsort(o2, (size_t)no, sizeof(bg_out), cmp_out);
        mpz_t vu, vw; mpz_inits(vu, vw, NULL); long pp0[2] = {0, 0};
        for (long long i = 0; i < no; i++) {
            const bg_out *e = &o2[i];
            if (e->kind != 1) {                          /* deferred: bottom3.c decides (rep_o = -1 - query index) */
                low(qflow[e->flow]); long t_; report_all = 1; rep_o = -1 - (long)e->flow; rep_pos = pp0; rep_nb = 1;
                printf("H %u\n", e->flow); solutions(qflow[e->flow], qc[e->flow], (long)e->p, G, &t_);
                while (i + 1 < no && o2[i + 1].flow == e->flow && o2[i + 1].p == e->p) i++;
                continue;
            }
            mpz_import(vu, BG_OUTW, -1, 4, 0, 0, e->up); mpz_import(vw, BG_OUTW, -1, 4, 0, 0, e->w);
            printf("Q %u t %u up ", e->flow, e->t); mpz_out_str(stdout, 10, vu); printf(" w "); mpz_out_str(stdout, 10, vw); printf("\n");
        }
        printf("KDONE %llu %llu %llu device %s\n", cn[0], cn[1], cn[2], bg_device());
        return 0;
    }

    /* ---------------- enumeration (bottom3.c's cum_rec order), batches within one ordering, output ---------------- */
    long batch = getenv("BG_BATCH") ? atol(getenv("BG_BATCH")) : (1L << 20);
    bg_flow *fl = malloc(sizeof(bg_flow) * batch); long nf = 0;
    long pos[BG_MAXB + 1];
    B_b = b; B_G = G; B_cs = cs;
    for (long o = skip; o < nord; o++) {
        /* iterative cum_rec: heights a[r] for r = b-1 (lowest digit) .. 1, each flow with its range of top heights p */
        if (b == 1) {
            long lo = 0, phi = heff[0]; if (lo < aklo) lo = aklo; if (phi > akhi + 1) phi = akhi + 1;
            if (lo < phi) { bg_flow *q = &fl[nf++]; memset(q, 0, sizeof *q); q->o = (unsigned)o; q->lo = (unsigned short)lo; q->phi = (unsigned short)phi; }
        } else {
            long r = b - 1, a[BG_MAXB + 1];
            long lo0 = 0, hi0 = heff[b - 1]; if (lo0 < aklo) lo0 = aklo; if (hi0 > akhi + 1) hi0 = akhi + 1;
            a[b - 1] = lo0 - 1;
            while (r <= b - 1) {
                const long hi = r == b - 1 ? hi0 : heff[r];
                if (++a[r] >= hi) { r++; continue; }
                pos[r] = a[r];
                if (r > 1) { r--; a[r] = pos[r + 1]; continue; }
                const long lo = pos[1] + 1, phi = heff[0];
                if (lo < phi) {
                    bg_flow *q = &fl[nf++]; memset(q, 0, sizeof *q); q->o = (unsigned)o;
                    for (long s2 = 1; s2 < b; s2++) q->pos[s2] = (unsigned short)pos[s2];
                    q->lo = (unsigned short)lo; q->phi = (unsigned short)phi;
                    if (nf == batch) { run_flows(fl, nf, o); nf = 0; }
                }
            }
        }
        if (nf) { run_flows(fl, nf, o); nf = 0; }
        printf("PROG %ld %llu\n", o + 1, B_cnt[0]); fflush(stdout);
    }
    printf("DONE %llu %llu %llu\n", B_cnt[0], B_cnt[1], B_cnt[2]);
    return 0;
}
