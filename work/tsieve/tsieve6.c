/* tsieve6: tree sieve for the main path of the chain sieve, with gsieve's later stages (prototype).
 *
 * Same region, same decisions and same counters as `gsieve_k13 main G N 0 N+1 [n0=N0]`:
 *   points  region points (a >= 26, x >= 10^XMIN as gsieve computes it, 2a + 3b + 6c <= 6N, > 6 N0 if N0 >= 0)
 *   s1      points whose last 13 digits avoid {0, 5}            (gsieve: bitmap)
 *   s2      points whose last 26 digits avoid {0, 5}            (gsieve: affine lift)
 *   s3, s4, gmp, full, members                                  (gsieve's leading-digit table, limb engine, GMP)
 * How the first 26 digits are decided:
 *   1. The residues r mod P_K whose power 2^r has K admissible last digits form a tree with exactly 4 children per
 *      node (lemma L2); its 4^K leaves are the cosets of L_K = {a + b mu3 + c mu7 = 0 mod P_K} that hold the s1 points
 *      (K = 13).
 *   2. In a coset, digits K+1..2K are the block h0 + 2^K Y, Y = Ys + n delta2 (mod 5^K) along a row of the reduced
 *      basis (lemma L3). Its lowest M digits pass iff Y mod 5^M lies in a set T(h0) of exactly 4^M residues, so the
 *      passing points of a row are the set bits of U = delta2^-1 T(h0) in the window [delta2^-1 Ys, + len): one
 *      64-bit read per 64 points. Each hit then has the rest of its block checked.
 *   3. Survivors go through gsieve's stage 3 (TOP table), stage 4 (limb_check) and stage 5 (full_check), included
 *      from gsieve.c, with M = 3^b 7^c mod 10^126 from tables.
 * Exactness: integer arithmetic except the float estimates in mulm1 and floordiv_d, which are corrected exactly.
 *
 * usage: tsieve6 PARAMFILE G|E M N [n0=N0] [threads=T] [dump]   (PARAMFILE from tparams.py 13)
 *   G: chain mode (digits outside {0,5}); E: zero-free mode (validation, positive control). In E mode a node has
 *   5 admissible lifts when their digits are odd and 4 when they are even, so leaf and T counts vary.
 * build: cc -O3 -march=native -DKK=13 -I../gsieve -o tsieve6 tsieve6.c -lgmp -lm -lpthread
 */
#define _GNU_SOURCE                                     /* exp10l in gsieve.c (before any header) */
#include <signal.h>
#define main gsieve_main
#include "gsieve.c"
#undef main

typedef __int128 i128;
static int TK, TM, TSKIP; static u64 TM1, TPK, TMm, TBRM; static i64 TV[3][3]; static u64 TWEPS[3];
static u64 TROOT[4][2]; static u64 TLP[20], TL5[21], TLC[20], TLU[20], TLG[20][5][2], TP2K[20];
static unsigned char TOK4[10000], TOKR[10000]; static int TNFULL, TNREM;
static long double TINVM1;
static u32 *TT; static int TTOK;
static u32 (*T3)[NL], (*T7)[NL];
typedef struct { i128 a[3], b0, br; } ineq;
static ineq RQ3[8], RQ2[32], RQ1[256]; static int n3, n2, n1; static i64 A3[8][3]; static double IC3[8];
#define NDMAX 4096
static __thread u32 *TBUF; static __thread int TN; static __thread u64 *TBMP, *TND, TNDLEN, TCU[20][5];
static __thread u64 c_leaves, c_rows, c_hits, c_s1rows;
static u64 G_leaves, G_rows, G_hits, G_s1rows, G_cnt[8];
static pthread_mutex_t OUTMU = PTHREAD_MUTEX_INITIALIZER;

static void rd(FILE *f, const char *name, int n, u64 *out) {
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, name)) { fprintf(stderr, "param: expected %s got %s\n", name, key); exit(2); }
    for (int i = 0; i < n; i++) if (fscanf(f, "%llu", (unsigned long long *)&out[i]) != 1) { fprintf(stderr, "param: %s[%d]\n", name, i); exit(2); }
}
static void rdi(FILE *f, const char *name, int n, i64 *out) {
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, name)) { fprintf(stderr, "param: expected %s got %s\n", name, key); exit(2); }
    for (int i = 0; i < n; i++) if (fscanf(f, "%lld", (long long *)&out[i]) != 1) { fprintf(stderr, "param: %s[%d]\n", name, i); exit(2); }
}
static inline void mul2(u64 alo, u64 ahi, u64 blo, u64 bhi, u64 *rlo, u64 *rhi) {   /* mod TM1^2, base-TM1 digits */
    u128 p = (u128)alo * blo; u64 q = (u64)(p / TM1); *rlo = (u64)(p - (u128)q * TM1);
    u128 s = (u128)q + (u128)alo * bhi + (u128)ahi * blo;
    *rhi = (u64)(s % TM1);
}
static inline u64 mulm1(u64 a, u64 b) {                 /* a b mod 5^K, a, b < 5^K < 2^46 (q off by <= 1, fixed) */
    u64 q = (u64)((long double)a * (long double)b * TINVM1);
    i64 r = (i64)(a * b - q * TM1);
    if (r < 0) r += (i64)TM1; else if (r >= (i64)TM1) r -= (i64)TM1;
    return (u64)r;
}
static inline u64 mods(i64 x, u64 m) {
    if (x >= 0) return (u64)x < m ? (u64)x : (u64)x % m;
    if ((u64)(-x) < m) return (u64)((i64)m + x);
    i64 r = x % (i64)m; return r ? (u64)(r + (i64)m) : 0;
}
static inline u64 modMm(u64 x) { u64 q = (u64)(((u128)x * TBRM) >> 64), r = x - q * TMm; while (r >= TMm) r -= TMm; return r; }
static u64 invmod(u64 a, u64 m) {
    i128 t = 0, nt = 1, r = (i128)m, nr = (i128)(a % m);
    while (nr) { i128 q = r / nr, x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    return (u64)(t < 0 ? t + (i128)m : t);
}
static inline i128 tfdiv(i128 a, i128 b) { i128 q = a / b; if ((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }
static inline i128 tcdiv(i128 a, i128 b) { return -tfdiv(-a, b); }
static inline i128 tfdiv_f(i128 a, i128 b) {            /* floor(a / b), b > 0: long double estimate, exact fix */
    i128 q = (i128)floorl((long double)a / (long double)b);
    while (b * q > a) q--;
    while (b * (q + 1) <= a) q++;
    return q;
}
static inline i64 floordiv_d(i64 x, i64 c, double ic) {  /* floor(x / c) for c > 0, |x| < 2^53 */
    i64 q = (i64)__builtin_floor((double)x * ic);
    while (c * q > x) q--;
    while (c * (q + 1) <= x) q++;
    return q;
}
static int eliminate(const ineq *src, int ns, ineq *dst, int k) {    /* Fourier-Motzkin: drop coordinate k */
    int nd = 0;
    for (int x = 0; x < ns; x++) if (src[x].a[k] == 0) dst[nd++] = src[x];
    for (int x = 0; x < ns; x++) if (src[x].a[k] > 0) for (int y = 0; y < ns; y++) if (src[y].a[k] < 0) {
        i128 cp = src[x].a[k], cn = -src[y].a[k];
        ineq d; memset(&d, 0, sizeof d);
        for (int j = 0; j < k; j++) d.a[j] = cn * src[x].a[j] + cp * src[y].a[j];
        d.b0 = cn * src[x].b0 + cp * src[y].b0; d.br = cn * src[x].br + cp * src[y].br;
        dst[nd++] = d;
    }
    return nd;
}
static inline int bounds(const ineq *q, int nq, int k, const i64 *x, i128 rr, i128 *lo, i128 *hi) {
    *lo = -((i128)1 << 100); *hi = (i128)1 << 100;
    for (int t = 0; t < nq; t++) {
        i128 rhs = q[t].b0 + q[t].br * rr;
        for (int j = 0; j < k; j++) rhs -= q[t].a[j] * x[j];
        i128 c = q[t].a[k];
        if (c > 0) { i128 v = tfdiv_f(rhs, c); if (v < *hi) *hi = v; }
        else if (c < 0) { i128 v = -tfdiv_f(rhs, -c); if (v > *lo) *lo = v; }
        else if (rhs < 0) return 0;
    }
    return *lo <= *hi;
}
static i64 R2A[32][2], R2B0[32], R2BR[32]; static double R2IC[32]; static int R2OK;   /* RQ2 in i64 if it fits */
static inline int bounds2(i64 i0, i64 r, i64 *lo, i64 *hi) {   /* i1 range given i0 (exact) */
    *lo = -(1LL << 62); *hi = 1LL << 62;
    for (int t = 0; t < n2; t++) {
        i64 rhs = R2B0[t] + R2BR[t] * r - R2A[t][0] * i0, c = R2A[t][1];
        if (c > 0) { i64 v = floordiv_d(rhs, c, R2IC[t]); if (v < *hi) *hi = v; }
        else if (c < 0) { i64 v = -floordiv_d(rhs, -c, R2IC[t]); if (v > *lo) *lo = v; }
        else if (rhs < 0) return 0;
    }
    return *lo <= *hi;
}
static inline int okblock_from(u64 v) {                 /* digits 4 TSKIP + 1 .. K of the K-digit block v */
    for (int c = 0; c < TSKIP; c++) v /= 10000;
    int ok = 1;
    for (int c = TSKIP; c < TNFULL; c++) { ok &= TOK4[v % 10000]; v /= 10000; }
    return TNREM ? ok & TOKR[v] : ok;
}
static int okdigits(u64 v, int n) { for (int i = 0; i < n; i++, v /= 10) if (v % 10 == 0 || (v % 10 == 5 && DM == 0)) return 0; return 1; }

/* T(h0): y mod 5^M with the lowest M digits of h0 + 2^K y admissible; u = y with TCU[i][t] = t 5^i, or u = delta2^-1 y */
static void tnode(int i, u64 u, u64 q) {
    if (i == TM) { TBUF[TN++] = (u32)u; return; }
    for (int t = 0; t < 5; t++) {
        u64 qq = q + (u64)t * TP2K[i], d = qq % 10;
        if (d == 0 || (d == 5 && DM == 0)) continue;
        u64 uu = u + TCU[i][t]; if (uu >= TMm) uu -= TMm;
        tnode(i + 1, uu, qq / 10);
    }
}
static u64 *TTOFF;                                      /* T(h0) = TT[TTOFF[h0] .. TTOFF[h0 + 1]) */
static void build_tt(size_t cap) {
    size_t maxper = (size_t)TMm, nh = (size_t)1 << TK, tot = 0;
    TBUF = malloc(maxper * 4); TTOFF = malloc((nh + 1) * 8);
    for (int i = 0; i < TM; i++) for (int t = 0; t < 5; t++) TCU[i][t] = (u64)t * TL5[i];
    for (size_t h = 0; h < nh; h++) { TTOFF[h] = tot; TN = 0; tnode(0, 0, h); tot += (size_t)TN;
        if (DM == 0 && TN != 1 << (2 * TM)) { fprintf(stderr, "T(%zu) has %d elements, expected 4^%d\n", h, TN, TM); exit(4); } }
    TTOFF[nh] = tot;
    if (tot * 4 > cap || !(TT = malloc(tot * 4))) { TTOK = 0; free(TBUF); TBUF = NULL; return; }
    for (size_t h = 0; h < nh; h++) { TN = 0; tnode(0, 0, h); memcpy(TT + TTOFF[h], TBUF, (size_t)TN * 4); }
    TTOK = 1; free(TBUF); TBUF = NULL;
}
static inline u64 bits64(u64 pos) { u64 q = pos >> 6, r = pos & 63; return r ? (TBMP[q] >> r) | (TBMP[q + 1] << (64 - r)) : TBMP[q]; }
static void clearbm(void) { for (int x = 0; x < TN; x++) { TBMP[TBUF[x] >> 6] = 0; TBMP[(TBUF[x] + TMm) >> 6] = 0; } }

/* stages 3-5 for the points whose last 26 digits are admissible: gsieve's later_stages, batched per thread so that
 * the table rows can be prefetched (survivors arrive in coset order, so a, b, c jump around) */
#define SBMAX 4096
#define PFD 12
typedef struct { u32 a, b, c; } spt;
static __thread spt *SB, *SB3; static __thread int nsb;
static u32 (*T3L)[4], (*T7L)[4];                       /* low 4 limbs of 3^b, 7^c (compact, for the prefilter) */
static void flush_survivors(void) {
    int n3 = 0, n4 = 0;
    for (int k = 0; k < nsb; k++) {                     /* region start x >= 10^XMIN (gsieve main_c), stage 3 */
        u64 a = SB[k].a, b = SB[k].b, c = SB[k].c;
        if (b * LG3 + c * LG7 < XMIN + 1) {
            double need = (XMIN - b * LG3 - c * LG7) / LG2;
            i64 alo = need > AMIN ? (i64)ceil(need - 1e-6) : AMIN;
            if ((i64)a < alo) continue;
        }
        cnt[2]++;
        if (DUMP) emit("S2 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        u64 bin = (a * L2 + b * L3 + c * L7) >> (64 - TB);
        if ((TOP[bin >> 6] >> (bin & 63)) & 1) { if (DUMP) emit("R3 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        cnt[3]++;
        SB3[n3++] = SB[k];
    }
    for (int k = 0; k < n3; k++) {                      /* limb_check's first test from low limbs (same decision) */
        if (k + PFD < n3) { __builtin_prefetch(T3L[SB3[k + PFD].b]); __builtin_prefetch(T7L[SB3[k + PFD].c]); __builtin_prefetch(T2H[SB3[k + PFD].a]); }
        u64 a = SB3[k].a, b = SB3[k].b, c = SB3[k].c;
        u64 m0, m1, m2, m3, acc, cy;
        const u32 *B3 = T3L[b], *B7 = T7L[c];
        acc = (u64)B3[0] * B7[0]; cy = acc / E9; m0 = acc - cy * E9;
        acc = (u64)B3[0] * B7[1] + (u64)B3[1] * B7[0] + cy; cy = acc / E9; m1 = acc - cy * E9;
        acc = (u64)B3[0] * B7[2] + (u64)B3[1] * B7[1] + (u64)B3[2] * B7[0] + cy; cy = acc / E9; m2 = acc - cy * E9;
        acc = (u64)B3[0] * B7[3] + (u64)B3[1] * B7[2] + (u64)B3[2] * B7[1] + (u64)B3[3] * B7[0] + cy; cy = acc / E9; m3 = acc - cy * E9;
        const u32 *AH = T2H[a]; u64 a0 = AH[0], a1 = AH[1], a2 = AH[2], a3 = AH[3];
        acc = a0 * m0; cy = acc / E9;
        acc = a0 * m1 + a1 * m0 + cy; cy = acc / E9;
        acc = a0 * m2 + a1 * m1 + a2 * m0 + cy; cy = acc / E9; u64 l2 = acc - cy * E9;
        acc = a0 * m3 + a1 * m2 + a2 * m1 + a3 * m0 + cy; cy = acc / E9; u64 l3 = acc - cy * E9;
#if KK == 13
        uint8_t f = FL1[l2 / 100000000ULL] | fl9(l3);
#else
#error "tsieve6 needs KK = 13"
#endif
        if (bad_of(f, C_EVEN)) { if (DUMP) emit("R4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        SB3[n4++] = SB3[k];
    }
    for (int k = 0; k < n4; k++) {                      /* the rest of limb_check (full M), then GMP */
        u64 a = SB3[k].a, b = SB3[k].b, c = SB3[k].c;
        u32 Mf[NL]; u64 cy = 0;
        for (int q = 0; q < NL; q++) {
            u64 acc = cy;
            for (int u = 0; u <= q; u++) acc += (u64)T3[b][u] * T7[c][q - u];
            cy = acc / E9; Mf[q] = (u32)(acc - cy * E9);
        }
        if (!limb_check(T2[a], Mf, T2H[a])) { if (DUMP) emit("R4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        cnt[4]++;
        if (DUMP) emit("S4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        full_check(0, a, b, c, C_EVEN);
    }
    nsb = 0;
}
static inline void survivor(i64 a, i64 b, i64 c) {
    SB[nsb++] = (spt){(u32)a, (u32)b, (u32)c};
    if (nsb == SBMAX) flush_survivors();
}

static void build_u(u64 h0, u64 inv2) {                /* U = delta2^-1 T(h0) as a doubled bitmap (positions in TBUF) */
    if (TTOK) {
        const u32 *T = TT + TTOFF[h0]; TN = (int)(TTOFF[h0 + 1] - TTOFF[h0]);
        for (int x = 0; x < TN; x++) TBUF[x] = (u32)modMm((u64)T[x] * inv2);
    } else {
        for (int i = 0; i < TM; i++) { u64 cc = inv2 * TL5[i] % TMm; for (int t = 0; t < 5; t++) TCU[i][t] = cc * t % TMm; }
        TN = 0; tnode(0, 0, h0);
    }
    for (int x = 0; x < TN; x++) { u64 u = TBUF[x], u2 = u + TMm; TBMP[u >> 6] |= 1ULL << (u & 63); TBMP[u2 >> 6] |= 1ULL << (u2 & 63); }
}
static void leaf(u64 r, u64 alo, u64 ahi) {
    c_leaves++;
    u64 h0 = (alo << TK) / TM1;
    u64 d[3]; for (int m = 0; m < 3; m++) d[m] = mulm1(alo, TWEPS[m]);
    u64 inv2 = invmod(d[2] % TMm, TMm);
    int built = 0;
#define BUILD_U() if (!built) { build_u(h0, inv2); built = 1; }
    TNDLEN = 1; TND[0] = 0;
    u64 sd1 = inv2 * (d[1] % TMm) % TMm;
    i64 xi[3]; i128 lo, hi, rr = (i128)r;
    if (bounds(RQ1, n1, 0, xi, rr, &lo, &hi)) for (i64 i0 = (i64)lo; i0 <= (i64)hi; i0++) {
        xi[0] = i0; i128 lo1, hi1;
        if (R2OK) { i64 l, h; if (!bounds2(i0, (i64)r, &l, &h)) continue; lo1 = l; hi1 = h; }
        else if (!bounds(RQ2, n2, 1, xi, rr, &lo1, &hi1)) continue;
        i64 rhs[8];
        for (int q = 0; q < n3; q++) rhs[q] = (i64)(RQ3[q].b0 + RQ3[q].br * rr) - A3[q][0] * i0 - A3[q][1] * (i64)lo1;
        u64 Yb = (ahi + mulm1(mods(i0, TM1), d[0]) + mulm1(mods((i64)lo1, TM1), d[1])) % TM1;
        u64 sb = Yb % TMm * inv2 % TMm;
        i64 plo = 0; u64 Ylo = 0; int have = 0;
        for (i64 i1 = (i64)lo1; i1 <= (i64)hi1; i1++) {
            i64 lo2 = -(1LL << 62), hi2 = 1LL << 62; int empty = 0;
            for (int q = 0; q < n3; q++) {
                i64 cc = A3[q][2];
                if (cc > 0) { i64 v = floordiv_d(rhs[q], cc, IC3[q]); if (v < hi2) hi2 = v; }
                else if (cc < 0) { i64 v = -floordiv_d(rhs[q], -cc, IC3[q]); if (v > lo2) lo2 = v; }
                else if (rhs[q] < 0) empty = 1;
            }
            if (!empty && lo2 <= hi2) {
                BUILD_U();
                c_rows++;
                u64 len = (u64)(hi2 - lo2 + 1);
                c_s1rows += len;
                i64 dl = lo2 - plo; u64 adl = have ? (u64)(dl < 0 ? -dl : dl) + 1 : 0, need = len > adl ? len : adl;
                while (TNDLEN < need && TNDLEN < NDMAX) { u64 t = TND[TNDLEN - 1] + d[2]; TND[TNDLEN++] = t >= TM1 ? t - TM1 : t; }
                if (have && dl >= 0 && dl < NDMAX) { Ylo += TND[dl]; if (Ylo >= TM1) Ylo -= TM1; }
                else if (have && dl < 0 && -dl < NDMAX) { Ylo += TM1 - TND[-dl]; if (Ylo >= TM1) Ylo -= TM1; }
                else Ylo = mulm1(mods(lo2, TM1), d[2]);
                plo = lo2; have = 1;
                u64 Ys = Yb + Ylo; if (Ys >= TM1) Ys -= TM1;
                u64 s = sb + mods(lo2, TMm); if (s >= TMm) s -= TMm;
                for (u64 base = 0; base < len; base += TMm) {
                    u64 wlen = len - base < TMm ? len - base : TMm;
                    for (u64 off = 0; off < wlen; off += 64) {
                        u64 word = bits64(s + off), rem = wlen - off;
                        if (rem < 64) word &= (1ULL << rem) - 1;
                        while (word) {
                            u64 n = base + off + (u64)__builtin_ctzll(word); word &= word - 1;
                            c_hits++;
                            u64 Y = Ys + (n < NDMAX ? TND[n] : mulm1(n % TM1, d[2])); if (Y >= TM1) Y -= TM1;
                            if (okblock_from(h0 + (Y << TK))) {
                                i64 i2 = lo2 + (i64)n;
                                survivor((i64)r + i0 * TV[0][0] + i1 * TV[1][0] + i2 * TV[2][0],
                                         i0 * TV[0][1] + i1 * TV[1][1] + i2 * TV[2][1], i0 * TV[0][2] + i1 * TV[1][2] + i2 * TV[2][2]);
                            }
                        }
                    }
                }
            }
            for (int q = 0; q < n3; q++) rhs[q] -= A3[q][1];
            Yb += d[1]; if (Yb >= TM1) Yb -= TM1;
            sb += sd1; if (sb >= TMm) sb -= TMm;
        }
    }
    if (built) clearbm();
}

static const int INV5[5] = {0, 1, 3, 2, 4};
typedef struct { int j; u64 r, ylo, yhi; } troot;
static troot *ROOTS; static u64 NROOTS, NEXTROOT, RFROM, RTO; static int DR;
static void node(int j, u64 r, u64 ylo, u64 yhi) {
    if (j == DR && ROOTS && NROOTS != (u64)-1 && !TBMP) { ROOTS[NROOTS++] = (troot){j, r, ylo, yhi}; return; }
    if (j == TK) { leaf(r, ylo, yhi); return; }
    u64 m = TL5[j] * 5;
    u64 Z = (u64)((u128)(ylo % m) * TLC[j] % m);
    u64 zlo = Z % TL5[j], q = Z / TL5[j];
    int p = 2 * zlo >= TL5[j];
    int v = (int)(zlo % 5 * TLU[j] % 5);
    int tbad = (DM == 1 && p) ? -1 : (int)(((p ? 2 : 0) + 5 - q) % 5) * INV5[v] % 5;   /* E: odd digits all admissible */
    for (int t = 0; t < 5; t++) {
        if (t == tbad) continue;
        u64 cl = ylo, ch = yhi;
        if (t) mul2(ylo, yhi, TLG[j][t][0], TLG[j][t][1], &cl, &ch);
        node(j + 1, r + (u64)t * TLP[j], cl, ch);
    }
}
static void *worker6(void *arg) {
    (void)arg;
    mpz_init(gx); mpz_init(gm); mpz_init(gt);
    OBC = 1 << 16; OB = malloc(OBC); OBN = 0;
    SB = malloc(sizeof(spt) * SBMAX); SB3 = malloc(sizeof(spt) * SBMAX); nsb = 0;
    TBUF = malloc(sizeof(u32) * TMm); TBMP = calloc(((2 * TMm) >> 6) + 3, 8); TND = malloc(8 * NDMAX);
    for (;;) {
        u64 k = __atomic_fetch_add(&NEXTROOT, 1, __ATOMIC_RELAXED);
        if (k >= RTO) break;
        node(ROOTS[k].j, ROOTS[k].r, ROOTS[k].ylo, ROOTS[k].yhi);
        flush_survivors();
        if (OBN) { pthread_mutex_lock(&OUTMU); fwrite(OB, 1, OBN, stdout); pthread_mutex_unlock(&OUTMU); OBN = 0; }
    }
    pthread_mutex_lock(&OUTMU);
    G_leaves += c_leaves; G_rows += c_rows; G_hits += c_hits; G_s1rows += c_s1rows;
    for (int k = 0; k < 8; k++) G_cnt[k] += cnt[k];
    pthread_mutex_unlock(&OUTMU);
    free(TBUF); free(TBMP); free(TND); free(OB);
    mpz_clear(gx); mpz_clear(gm); mpz_clear(gt);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: tsieve6 PARAMFILE G|E M N [n0=N0] [threads=T] [dump]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    u64 t[16];
    rd(f, "K", 1, t); TK = (int)t[0]; rd(f, "M1", 1, &TM1); rd(f, "PK", 1, &TPK);
    rdi(f, "V0", 3, TV[0]); rdi(f, "V1", 3, TV[1]); rdi(f, "V2", 3, TV[2]); rd(f, "WEPS", 3, TWEPS);
    for (int r = 0; r < 4; r++) { char nm[16]; sprintf(nm, "ROOT%d", r); rd(f, nm, 2, TROOT[r]); }
    for (int j = 1; j < TK; j++) {
        char nm[16]; sprintf(nm, "LEV%d", j); rd(f, nm, 14, t);
        TLP[j] = t[0]; TLC[j] = t[2]; TLU[j] = t[3];
        for (int u = 0; u < 5; u++) { TLG[j][u][0] = t[4 + 2 * u]; TLG[j][u][1] = t[5 + 2 * u]; }
    }
    fclose(f);
    if (TK != KK) { fprintf(stderr, "PARAMFILE has K=%d, this build needs K=%d (gsieve's later stages start at digit 2K+1)\n", TK, KK); return 1; }
    TL5[0] = 1; for (int j = 1; j <= TK; j++) TL5[j] = TL5[j - 1] * 5;
    DM = argv[2][0] == 'E' ? 1 : 0; TM = atoi(argv[3]); NN = atoll(argv[4]);
    int NTH = 1;
    for (int i = 5; i < argc; i++) {
        if (!strncmp(argv[i], "n0=", 3)) N0 = atoll(argv[i] + 3);
        else if (!strncmp(argv[i], "threads=", 8)) NTH = atoi(argv[i] + 8);
        else if (!strcmp(argv[i], "dump")) DUMP = 1;
        else if (!strncmp(argv[i], "roots=", 6)) { RFROM = strtoull(argv[i] + 6, 0, 10); RTO = strtoull(strchr(argv[i], ':') + 1, 0, 10); }
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 1; }
    }
    if (TM < 1 || TM > TK || TM > 13 || NTH < 1) { fprintf(stderr, "need 1 <= M <= 13, threads >= 1\n"); return 1; }
#ifdef PROFALARM
    signal(SIGALRM, (void (*)(int))exit); alarm(PROFALARM);
#endif
    BADEVEN = DM == 0 ? (F_ZERO | F_FIVE) : F_ZERO;
    build_fl(FL1, 1); build_fl(FL2, 2); build_fl(FL3, 3); build_fl(FL4, 4); build_fl(FL5, 5);
    build_top(); build_t2(3 * (u64)NN + 2);
    T3 = malloc((2 * (u64)NN + 2) * sizeof *T3); T7 = malloc(((u64)NN + 2) * sizeof *T7);
    memset(T3[0], 0, sizeof T3[0]); T3[0][0] = 1; memset(T7[0], 0, sizeof T7[0]); T7[0][0] = 1;
    for (u64 b = 1; b <= 2 * (u64)NN + 1; b++) { memcpy(T3[b], T3[b - 1], sizeof T3[b]); mulsmall(T3[b], 3); }
    for (u64 c = 1; c <= (u64)NN + 1; c++) { memcpy(T7[c], T7[c - 1], sizeof T7[c]); mulsmall(T7[c], 7); }
    T3L = malloc((2 * (u64)NN + 2) * sizeof *T3L); T7L = malloc(((u64)NN + 2) * sizeof *T7L);
    for (u64 b = 0; b <= 2 * (u64)NN + 1; b++) memcpy(T3L[b], T3[b], sizeof T3L[b]);
    for (u64 c = 0; c <= (u64)NN + 1; c++) memcpy(T7L[c], T7[c], sizeof T7L[c]);
    TMm = TL5[TM]; TSKIP = TM / 4; TBRM = (u64)(((u128)1 << 64) / TMm);
    for (int i = 0; i < TM; i++) TP2K[i] = 1ULL << (TK - i);
    for (int x = 0; x < 10000; x++) TOK4[x] = okdigits(x, 4);
    TNFULL = TK / 4; TNREM = TK % 4;
    u64 pr = 1; for (int i = 0; i < TNREM; i++) pr *= 10;
    for (int x = 0; x < 10000; x++) TOKR[x] = TNREM ? ((u64)x < pr && okdigits(x, TNREM)) : 0;
    TINVM1 = 1.0L / (long double)TM1;
    /* region a >= AMIN, b, c >= 0, 2a + 3b + 6c <= 6N (> 6 N0) in lattice coordinates: p = (r,0,0) + V^T i */
    i64 nv[5][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {2, 3, 6}, {-2, -3, -6}};
    i64 bv[5] = {-AMIN, 0, 0, 6 * NN, -6 * N0 - 1};
    n3 = N0 >= 0 ? 5 : 4;
    for (int q = 0; q < n3; q++) {
        for (int m = 0; m < 3; m++) RQ3[q].a[m] = nv[q][0] * TV[m][0] + nv[q][1] * TV[m][1] + nv[q][2] * TV[m][2];
        RQ3[q].b0 = bv[q]; RQ3[q].br = -nv[q][0];
    }
    n2 = eliminate(RQ3, n3, RQ2, 2); n1 = eliminate(RQ2, n2, RQ1, 1);
    for (int q = 0; q < n3; q++) { for (int m = 0; m < 3; m++) A3[q][m] = (i64)RQ3[q].a[m]; IC3[q] = A3[q][2] ? 1.0 / (double)(A3[q][2] > 0 ? A3[q][2] : -A3[q][2]) : 0; }
    {   /* RQ2 in i64: every |rhs| stays below 2^52 for coordinates up to 2^26 and r < P_K */
        R2OK = 1;
        for (int t = 0; t < n2; t++) {
            i128 bound = (RQ2[t].b0 < 0 ? -RQ2[t].b0 : RQ2[t].b0) + (RQ2[t].br < 0 ? -RQ2[t].br : RQ2[t].br) * (i128)TPK;
            for (int j = 0; j < 2; j++) bound += (RQ2[t].a[j] < 0 ? -RQ2[t].a[j] : RQ2[t].a[j]) << 26;
            if (bound >= ((i128)1 << 52)) R2OK = 0;
            R2A[t][0] = (i64)RQ2[t].a[0]; R2A[t][1] = (i64)RQ2[t].a[1]; R2B0[t] = (i64)RQ2[t].b0; R2BR[t] = (i64)RQ2[t].br;
            R2IC[t] = R2A[t][1] ? 1.0 / (double)(R2A[t][1] > 0 ? R2A[t][1] : -R2A[t][1]) : 0;
        }
    }
    build_tt((size_t)1 << 30);
    /* points (gsieve's region loop) and the s1 points below 10^XMIN that the cosets contain (subtracted) */
    u64 points = 0, corner = 0;
    for (i64 c = 0; c <= NN && !RFROM && !RTO; c++) {      /* (skipped for a subset of roots: benchmarks) */
        if (6 * NN - 6 * c - 2 * AMIN < 0) break;
        i64 bmax = (6 * NN - 6 * c - 2 * AMIN) / 3;
        for (i64 b = 0; b <= bmax; b++) {
            i64 ahi = (6 * NN - 3 * b - 6 * c) / 2;
            double need = (XMIN - b * LG3 - c * LG7) / LG2;
            i64 alo = need > AMIN ? (i64)ceil(need - 1e-6) : AMIN, a0 = AMIN;
            if (N0 >= 0 && 6 * N0 - 3 * b - 6 * c >= 0) { i64 t0 = (6 * N0 - 3 * b - 6 * c) / 2 + 1; if (t0 > alo) alo = t0; if (t0 > a0) a0 = t0; }
            if (alo <= ahi) points += (u64)(ahi - alo + 1);
            i64 ae = (alo - 1 < ahi ? alo - 1 : ahi);
            if (a0 <= ae) {                                 /* polytope points with a < alo: count those with 13 good digits */
                const u64 m13 = 10000000000000ULL; u128 x = 1;
                for (i64 i = 0; i < a0; i++) x = x * 2 % m13;
                for (i64 i = 0; i < b; i++) x = x * 3 % m13;
                for (i64 i = 0; i < c; i++) x = x * 7 % m13;
                for (i64 a = a0; a <= ae; a++) { if (okdigits((u64)x, 13)) corner++; x = x * 2 % m13; }
            }
        }
    }
    /* subtree roots at depth DR, then the threads */
    DR = TK < 7 ? TK : 7;
    ROOTS = malloc(sizeof(troot) * 4 * TL5[DR - 1]); NROOTS = 0;      /* at most 4 5^(DR-1) roots (E mode) */
    for (int r = 0; r < 4; r++) node(1, (u64)r, TROOT[r][0], TROOT[r][1]);
    if (DM == 0 && NROOTS != 1ULL << (2 * DR)) { fprintf(stderr, "roots %llu, expected 4^%d\n", (unsigned long long)NROOTS, DR); return 4; }
    if (RTO == 0 || RTO > NROOTS) RTO = NROOTS;
    NEXTROOT = RFROM;
    struct timespec w0, w1; clock_gettime(CLOCK_MONOTONIC, &w0); clock_t c0 = clock();
    pthread_t *tid = malloc(NTH * sizeof *tid);
    for (int i = 0; i < NTH; i++) pthread_create(&tid[i], NULL, worker6, NULL);
    for (int i = 0; i < NTH; i++) pthread_join(tid[i], NULL);
    clock_gettime(CLOCK_MONOTONIC, &w1);
    double wall = (double)(w1.tv_sec - w0.tv_sec) + 1e-9 * (double)(w1.tv_nsec - w0.tv_nsec), cpu = (double)(clock() - c0) / CLOCKS_PER_SEC;
    if (DM == 0 && RFROM == 0 && RTO == NROOTS && G_leaves != 1ULL << (2 * TK)) { fprintf(stderr, "leaves %llu, expected 4^%d\n", (unsigned long long)G_leaves, TK); return 4; }
    printf("DONE mode=tree dm=%s N=%lld N0=%lld K=%d M=%d NL=%d points=%llu s1=%llu s2=%llu s3=%llu s4=%llu gmp=%llu full=%llu members=%llu "
           "leaves=%llu rows=%llu hits=%llu tt=%d roots=%llu:%llu/%llu time=%.2f threads=%d wall=%.2f\n",
           DM ? "E" : "G", (long long)NN, (long long)N0, TK, TM, NL, (unsigned long long)points, (unsigned long long)(G_s1rows - corner),
           (unsigned long long)G_cnt[2], (unsigned long long)G_cnt[3], (unsigned long long)G_cnt[4], (unsigned long long)G_cnt[5],
           (unsigned long long)G_cnt[6], (unsigned long long)G_cnt[7], (unsigned long long)G_leaves, (unsigned long long)G_rows,
           (unsigned long long)G_hits, TTOK, (unsigned long long)RFROM, (unsigned long long)RTO, (unsigned long long)NROOTS, cpu, NTH, wall);
    return 0;
}
