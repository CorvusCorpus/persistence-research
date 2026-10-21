/* Tree sieve, version 4 (prototype): cosets of depth K, bitmap hit lists for digits K+1..K+m.
 *
 * Region and output as tsieve.c: every point (a,b,c) with a >= AMIN (>= 2K), b, c >= 0, 2a + 3b + 6c <= 6N (and
 * > 6 N0 if N0 >= 0) whose last 2K digits avoid {0, 5}, without visiting points that fail the first K + m of them.
 *   1. Tree of admissible residues r mod P_K (4^K leaves, lemma L2), as in tsieve.c.
 *   2. Coset of r: points (r,0,0) + i0 V0 + i1 V1 + i2 V2, rows along V2. Digits K+1..2K are the K-digit block
 *      h0 + 2^K Y, Y = Y0 + i.delta (mod 5^K) (lemma L3). Its lowest m digits depend only on Y mod 5^m, and
 *      T = {y mod 5^m : those m digits admissible} has exactly 4^m elements (each 5-adic lift of y moves digit i+1 by
 *      2^(K-i) t mod 10: five digits of one parity, exactly one of them 0 or 5). T is built as a tree from h0.
 *   3. Along a row, Y = Ys + n delta2; the admissible n are n = (u - s) mod 5^m with u in U = delta2^-1 T, s =
 *      delta2^-1 Ys. U is a bitmap of 5^m bits (stored twice, so windows never wrap); a row reads its window
 *      [s, s + len) 64 positions per word, and each set bit is a hit.
 *   4. Each hit gets the rest of its K-digit block checked (branch-free table lookups).
 * Timing builds: -DNOSCAN (rows only), -DNOHIT (rows and bitmap windows), -DNOCHECK (no digit check).
 * Row bounds are exact: the region is a polytope in (i0, i1, i2); i2 bounds come from its inequalities, (i0, i1)
 * bounds from Fourier-Motzkin projections (rows with an empty i2 range are skipped).
 *
 * usage: tsieve4 PARAMFILE M N N0 AMIN [dump]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64; typedef long long i64; typedef unsigned __int128 u128; typedef __int128 i128;
typedef unsigned int u32;

static int K, MM; static u64 M1, PK, Mm; static i64 V[3][3]; static u64 WEPS[3];
static u64 ROOT[4][2];
static u64 LP[20], L5[20], LC[20], LU[20], LG[20][5][2];
static i64 N, N0, AMIN; static int DUMP;
static unsigned char OK4[10000], OKR[10000]; static int NFULL, NREM;
static u64 c_nodes, c_leaves, c_rows, c_hits, c_surv, c_tnodes;

/* inequalities in lattice coordinates: A.i <= B0 + Br * r  (3D: i0, i1, i2; 2D: i0, i1; 1D: i0) */
typedef struct { i128 a[3], b0, br; } ineq;
static ineq Q3[8], Q2[32], Q1[256]; static int n3, n2, n1;

static void rd(FILE *f, const char *name, int n, u64 *out) {
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, name)) { fprintf(stderr, "param: expected %s got %s\n", name, key); exit(2); }
    for (int i = 0; i < n; i++) if (fscanf(f, "%llu", &out[i]) != 1) { fprintf(stderr, "param: %s[%d]\n", name, i); exit(2); }
}
static void rdi(FILE *f, const char *name, int n, i64 *out) {
    char key[64];
    if (fscanf(f, "%63s", key) != 1 || strcmp(key, name)) { fprintf(stderr, "param: expected %s got %s\n", name, key); exit(2); }
    for (int i = 0; i < n; i++) if (fscanf(f, "%lld", &out[i]) != 1) { fprintf(stderr, "param: %s[%d]\n", name, i); exit(2); }
}
static inline void mul2(u64 alo, u64 ahi, u64 blo, u64 bhi, u64 *rlo, u64 *rhi) {
    u128 p = (u128)alo * blo; u64 q = (u64)(p / M1); *rlo = (u64)(p - (u128)q * M1);
    u128 s = (u128)q + (u128)alo * bhi + (u128)ahi * blo;
    *rhi = (u64)(s % M1);
}
static inline i128 fdiv(i128 a, i128 b) { i128 q = a / b; if ((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }
static inline i128 cdiv(i128 a, i128 b) { return -fdiv(-a, b); }
static inline u64 modp(i128 x, u64 m) { i128 r = x % (i128)m; return (u64)(r < 0 ? r + (i128)m : r); }
static inline u64 mods(i64 x, u64 m) {                  /* x mod m in [0, m), fast when |x| < m */
    if (x >= 0) return (u64)x < m ? (u64)x : (u64)x % m;
    if ((u64)(-x) <= m) return (u64)((i64)m + x) == m ? 0 : (u64)((i64)m + x);
    i64 r = x % (i64)m; return r ? (u64)(r + (i64)m) : 0;
}
static inline int okblock_nb(u64 v) {                   /* branch-free: all K digits outside {0,5} */
    int ok = 1;
    for (int c = 0; c < NFULL; c++) { ok &= OK4[v % 10000]; v /= 10000; }
    return NREM ? ok & OKR[v] : ok;
}
static inline int okblock(u64 v) {
    for (int c = 0; c < NFULL; c++) { if (!OK4[v % 10000]) return 0; v /= 10000; }
    return NREM ? OKR[v] : 1;
}
static int okdigits(u64 v, int n) { for (int i = 0; i < n; i++, v /= 10) if (v % 10 == 0 || v % 10 == 5) return 0; return 1; }
static u64 invmod(u64 a, u64 m) {   /* a^-1 mod m, gcd(a, m) = 1 */
    i128 t = 0, nt = 1, r = (i128)m, nr = (i128)(a % m);
    while (nr) { i128 q = r / nr, x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    return (u64)(t < 0 ? t + (i128)m : t);
}

/* Fourier-Motzkin: eliminate coordinate k from inequalities src (dim k+1) into dst (dim k) */
static int eliminate(const ineq *src, int ns, ineq *dst, int k) {
    int nd = 0;
    for (int x = 0; x < ns; x++) if (src[x].a[k] == 0) dst[nd++] = src[x];
    for (int x = 0; x < ns; x++) if (src[x].a[k] > 0) for (int y = 0; y < ns; y++) if (src[y].a[k] < 0) {
        i128 cp = src[x].a[k], cn = -src[y].a[k];       /* cn * (x) + cp * (y) cancels coordinate k */
        ineq d; memset(&d, 0, sizeof d);
        for (int j = 0; j < k; j++) d.a[j] = cn * src[x].a[j] + cp * src[y].a[j];
        d.b0 = cn * src[x].b0 + cp * src[y].b0; d.br = cn * src[x].br + cp * src[y].br;
        dst[nd++] = d;
    }
    return nd;
}
/* bounds on coordinate k given the earlier coordinates (exact; returns 0 if infeasible) */
static inline int bounds(const ineq *q, int nq, int k, const i64 *x, i128 rr, i128 *lo, i128 *hi) {
    *lo = -((i128)1 << 100); *hi = (i128)1 << 100;
    for (int t = 0; t < nq; t++) {
        i128 rhs = q[t].b0 + q[t].br * rr;
        for (int j = 0; j < k; j++) rhs -= q[t].a[j] * x[j];
        i128 c = q[t].a[k];
        if (c > 0) { i128 v = fdiv(rhs, c); if (v < *hi) *hi = v; }
        else if (c < 0) { i128 v = cdiv(rhs, c); if (v > *lo) *lo = v; }
        else if (rhs < 0) return 0;
    }
    return *lo <= *hi;
}

/* U = delta2^-1 T mod 5^m, built directly: T is the tree of y mod 5^m with the lowest m digits of h0 + 2^K y
 * admissible (from y = 0, Q = h0; digit i+1 of a lift y + t 5^i is (Q + t 2^(K-i)) mod 10), and u = delta2^-1 y is
 * carried along (u += t * (delta2^-1 5^i mod 5^m)). */
static u32 *TBUF; static int TN;
static u64 P2K_I[20], CU[20][5];
static void tnode(int i, u64 u, u64 q) {
    c_tnodes++;
    if (i == MM) { TBUF[TN++] = (u32)u; return; }
    for (int t = 0; t < 5; t++) {
        u64 qq = q + (u64)t * P2K_I[i];
        u64 d = qq % 10;
        if (d == 0 || d == 5) continue;
        u64 uu = u + CU[i][t]; if (uu >= Mm) uu -= Mm;
        tnode(i + 1, uu, qq / 10);
    }
}
static long double INVM1;
static inline u64 mulm1(u64 a, u64 b) {                 /* a b mod 5^K for a, b < 5^K < 2^46 (exact: q is off by <= 1) */
    u64 q = (u64)((long double)a * (long double)b * INVM1);
    i64 r = (i64)(a * b - q * M1);
    if (r < 0) r += (i64)M1; else if (r >= (i64)M1) r -= (i64)M1;
    return (u64)r;
}
static inline i64 floordiv_d(i64 x, i64 c, double ic) {  /* floor(x / c), c != 0, |x| < 2^53: float estimate, exact fix */
    i64 q = (i64)__builtin_floor((double)x * ic);
    if (c > 0) { while (c * q > x) q--; while (c * (q + 1) <= x) q++; }
    else { while (c * q < x) q--; while (c * (q + 1) >= x) q++; }
    return q;
}
static u64 *BM; static u64 *ND; static u64 ndlen;
#define NDMAX 4096
static void clearbm(void) { for (int x = 0; x < TN; x++) { BM[TBUF[x] >> 6] = 0; BM[(TBUF[x] + Mm) >> 6] = 0; } }
static inline u64 bits64(u64 pos) { u64 q = pos >> 6, r = pos & 63; return r ? (BM[q] >> r) | (BM[q + 1] << (64 - r)) : BM[q]; }
static int SKIP;                                        /* 4-digit chunks of the block known admissible (M >= 4 SKIP) */
static inline int okblock_from(u64 v) {
    for (int c = 0; c < SKIP; c++) v /= 10000;
    int ok = 1;
    for (int c = SKIP; c < NFULL; c++) { ok &= OK4[v % 10000]; v /= 10000; }
    return NREM ? ok & OKR[v] : ok;
}

static i64 A3[8][3]; static double IC3[8];              /* 3D inequalities as i64, with 1/a2 */
static void leaf(u64 r, u64 alo, u64 ahi) {
    c_leaves++;
    u64 h0 = (alo << K) / M1;
    u64 d[3]; for (int m = 0; m < 3; m++) d[m] = mulm1(alo, WEPS[m]);
    u64 inv2 = invmod(d[2] % Mm, Mm);
    for (int i = 0; i < MM; i++) { u64 c = inv2 * L5[i] % Mm; for (int t = 0; t < 5; t++) CU[i][t] = c * t % Mm; }
    TN = 0; tnode(0, 0, h0);
    if (TN != 1 << (2 * MM)) { fprintf(stderr, "T has %d elements, expected 4^%d\n", TN, MM); exit(4); }
    for (int x = 0; x < TN; x++) { u64 u = TBUF[x], u2 = u + Mm; BM[u >> 6] |= 1ULL << (u & 63); BM[u2 >> 6] |= 1ULL << (u2 & 63); }
    ndlen = 1; ND[0] = 0;
    u64 sd1 = inv2 * (d[1] % Mm) % Mm;
    i64 xi[3]; i128 lo, hi; i128 rr = (i128)r;
    if (!bounds(Q1, n1, 0, xi, rr, &lo, &hi)) { clearbm(); return; }
    for (i64 i0 = (i64)lo; i0 <= (i64)hi; i0++) {
        xi[0] = i0; i128 lo1, hi1;
        if (!bounds(Q2, n2, 1, xi, rr, &lo1, &hi1)) continue;
        i64 rhs[8];
        for (int q = 0; q < n3; q++) rhs[q] = (i64)(Q3[q].b0 + Q3[q].br * rr) - A3[q][0] * i0 - A3[q][1] * (i64)lo1;
        u64 Yb = (ahi + mulm1(mods(i0, M1), d[0]) + mulm1(mods((i64)lo1, M1), d[1])) % M1;
        u64 sb = Yb % Mm * inv2 % Mm;
        i64 plo = 0; u64 Ylo = 0; int have = 0;                 /* Ylo = plo * d2 mod 5^K, updated by small steps */
        for (i64 i1 = (i64)lo1; i1 <= (i64)hi1; i1++) {
            i64 lo2 = -(1LL << 62), hi2 = 1LL << 62; int empty = 0;
            for (int q = 0; q < n3; q++) {
                i64 c = A3[q][2];
                if (c > 0) { i64 v = floordiv_d(rhs[q], c, IC3[q]); if (v < hi2) hi2 = v; }
                else if (c < 0) { i64 v = -floordiv_d(rhs[q], -c, -IC3[q]); if (v > lo2) lo2 = v; }   /* ceil(rhs/c) */
                else if (rhs[q] < 0) empty = 1;
            }
            if (!empty && lo2 <= hi2) {
                c_rows++;
                u64 len = (u64)(hi2 - lo2 + 1);
#ifdef NOSCAN
                c_hits += len; if (1) goto nextrow;
#endif
                i64 dl = lo2 - plo; u64 adl = have ? (u64)(dl < 0 ? -dl : dl) + 1 : 0, need = len > adl ? len : adl;
                while (ndlen < need && ndlen < NDMAX) { u64 t = ND[ndlen - 1] + d[2]; ND[ndlen++] = t >= M1 ? t - M1 : t; }
                if (have && dl >= 0 && dl < NDMAX) { Ylo += ND[dl]; if (Ylo >= M1) Ylo -= M1; }
                else if (have && dl < 0 && -dl < NDMAX) { Ylo += M1 - ND[-dl]; if (Ylo >= M1) Ylo -= M1; }
                else Ylo = mulm1(mods(lo2, M1), d[2]);
                plo = lo2; have = 1;
                u64 Ys = Yb + Ylo; if (Ys >= M1) Ys -= M1;
                u64 s = sb + mods(lo2, Mm); if (s >= Mm) s -= Mm;
                for (u64 base = 0; base < len; base += Mm) {
                    u64 wlen = len - base < Mm ? len - base : Mm;      /* window [s, s + wlen) of the doubled bitmap */
                    for (u64 off = 0; off < wlen; off += 64) {
                        u64 word = bits64(s + off), rem = wlen - off;
                        if (rem < 64) word &= (1ULL << rem) - 1;
#ifdef NOHIT
                        c_hits += (u64)__builtin_popcountll(word); word = 0;
#endif
                        while (word) {
                            u64 n = base + off + (u64)__builtin_ctzll(word); word &= word - 1;
                            c_hits++;
                            u64 Y = Ys + (n < NDMAX ? ND[n] : mulm1(n % M1, d[2])); if (Y >= M1) Y -= M1;
#ifdef NOCHECK
                            int ok = (int)(Y & 1); (void)h0;
#else
                            int ok = okblock_from(h0 + (Y << K));
#endif
                            c_surv += (u64)ok;
                            if (DUMP && ok) {
                                i64 i2 = lo2 + (i64)n;
                                printf("S %lld %lld %lld\n", (i64)r + i0 * V[0][0] + i1 * V[1][0] + i2 * V[2][0],
                                       i0 * V[0][1] + i1 * V[1][1] + i2 * V[2][1], i0 * V[0][2] + i1 * V[1][2] + i2 * V[2][2]);
                            }
                        }
                    }
                }
            }
#ifdef NOSCAN
            nextrow:
#endif
            for (int q = 0; q < n3; q++) rhs[q] -= A3[q][1];
            Yb += d[1]; if (Yb >= M1) Yb -= M1;
            sb += sd1; if (sb >= Mm) sb -= Mm;
        }
    }
    clearbm();
}

static const int INV5[5] = {0, 1, 3, 2, 4};
static void node(int j, u64 r, u64 ylo, u64 yhi) {
    c_nodes++;
    if (j == K) { leaf(r, ylo, yhi); return; }
    u64 m = L5[j] * 5;
    u64 Z = (u64)((u128)(ylo % m) * LC[j] % m);
    u64 zlo = Z % L5[j], q = Z / L5[j];
    int p = 2 * zlo >= L5[j];
    int v = (int)(zlo % 5 * LU[j] % 5);
    int tbad = (int)(((p ? 2 : 0) + 5 - q) % 5) * INV5[v] % 5;
    for (int t = 0; t < 5; t++) {
        if (t == tbad) continue;
        u64 cl = ylo, ch = yhi;
        if (t) mul2(ylo, yhi, LG[j][t][0], LG[j][t][1], &cl, &ch);
        node(j + 1, r + (u64)t * LP[j], cl, ch);
    }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: tsieve4 PARAMFILE M N N0 AMIN [dump]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    u64 t[16];
    rd(f, "K", 1, t); K = (int)t[0]; rd(f, "M1", 1, &M1); rd(f, "PK", 1, &PK);
    rdi(f, "V0", 3, V[0]); rdi(f, "V1", 3, V[1]); rdi(f, "V2", 3, V[2]); rd(f, "WEPS", 3, WEPS);
    for (int r = 0; r < 4; r++) { char nm[16]; sprintf(nm, "ROOT%d", r); rd(f, nm, 2, ROOT[r]); }
    for (int j = 1; j < K; j++) {
        char nm[16]; sprintf(nm, "LEV%d", j); rd(f, nm, 14, t);
        LP[j] = t[0]; L5[j] = t[1]; LC[j] = t[2]; LU[j] = t[3];
        for (int u = 0; u < 5; u++) { LG[j][u][0] = t[4 + 2 * u]; LG[j][u][1] = t[5 + 2 * u]; }
    }
    fclose(f);
    L5[0] = 1; for (int j = 1; j <= K; j++) L5[j] = L5[j - 1] * 5;
    MM = atoi(argv[2]); N = atoll(argv[3]); N0 = atoll(argv[4]); AMIN = atoll(argv[5]);
    DUMP = argc > 6 && !strcmp(argv[6], "dump");
    if (AMIN < 2 * K || MM < 1 || MM > K || MM > 13) { fprintf(stderr, "need AMIN >= 2K and 1 <= M <= min(K, 13)\n"); return 1; }
    Mm = L5[MM]; SKIP = MM / 4;
    for (int i = 0; i < MM; i++) P2K_I[i] = 1ULL << (K - i);
    for (int x = 0; x < 10000; x++) OK4[x] = okdigits(x, 4);
    NFULL = K / 4; NREM = K % 4;
    u64 pr = 1; for (int i = 0; i < NREM; i++) pr *= 10;
    for (int x = 0; x < 10000; x++) OKR[x] = NREM ? ((u64)x < pr && okdigits(x, NREM)) : 0;
    /* region in lattice coordinates: n.p <= beta with p = (r,0,0) + V^T i */
    i64 nv[5][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {2, 3, 6}, {-2, -3, -6}};
    i64 bv[5] = {-AMIN, 0, 0, 6 * N, -6 * N0 - 1};
    n3 = N0 >= 0 ? 5 : 4;
    for (int q = 0; q < n3; q++) {
        for (int m = 0; m < 3; m++) Q3[q].a[m] = nv[q][0] * V[m][0] + nv[q][1] * V[m][1] + nv[q][2] * V[m][2];
        Q3[q].b0 = bv[q]; Q3[q].br = -nv[q][0];
    }
    n2 = eliminate(Q3, n3, Q2, 2); n1 = eliminate(Q2, n2, Q1, 1);
    for (int q = 0; q < n3; q++) { for (int m = 0; m < 3; m++) A3[q][m] = (i64)Q3[q].a[m]; IC3[q] = A3[q][2] ? 1.0 / (double)A3[q][2] : 0; }
    INVM1 = 1.0L / (long double)M1;
    TBUF = malloc(sizeof(u32) << (2 * MM)); BM = calloc(((2 * Mm) >> 6) + 3, 8); ND = malloc(8 * NDMAX);
    clock_t t0 = clock();
    for (int r = 0; r < 4; r++) node(1, (u64)r, ROOT[r][0], ROOT[r][1]);
    double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("DONE K=%d M=%d N=%lld N0=%lld AMIN=%lld nodes=%llu leaves=%llu tnodes=%llu rows=%llu hits=%llu surv=%llu "
           "ineq=%d/%d/%d cpu=%.3f\n", K, MM, N, N0, AMIN, c_nodes, c_leaves, c_tnodes, c_rows, c_hits, c_surv,
           n3, n2, n1, sec);
    return 0;
}
