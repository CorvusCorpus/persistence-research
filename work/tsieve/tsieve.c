/* Tree sieve for chain numbers of form 2^a 3^b 7^c (prototype).
 *
 * Decides, for every point (a,b,c) of the region  a >= AMIN, b >= 0, c >= 0, 2a + 3b + 6c <= 6N  (and
 * 2a + 3b + 6c > 6 N0 if N0 >= 0), whether the last 2K digits of x = 2^a 3^b 7^c avoid {0, 5}, without looking at
 * the points whose last K digits fail:
 *   1. Walk the tree of residues r mod P_k (P_k = 4*5^(k-1)) whose power 2^r has k admissible last digits, k = 1..K.
 *      Every node has exactly 4 admissible children (lemma L2); the excluded child is found from one digit computed
 *      at t = 0, since the five lifts r + t P_k have digits p + 2((q + t v) mod 5).
 *   2. At a leaf r (depth K), the region points with e(a,b,c) = a + b mu3 + c mu7 = r (mod P_K) are the points
 *      (r,0,0) + i0 V0 + i1 V1 + i2 V2 (reduced basis of the lattice L_K). Digits K+1..2K of x are h0 + 2^K Y with
 *      Y = A_hi + delta . i (mod 5^K) (lemma L3 with j = J(p) affine in i), so each step along i2 is one addition.
 * Requires AMIN >= 2K (suffix collapse at level 2K). The digit block is exactly K digits, leading zeros included,
 * which is correct whenever x has more than 2K digits (the production sieve requires x >= 10^130).
 *
 * usage: tsieve PARAMFILE N N0 AMIN [dump]       dump: print "S a b c" for every point that passes 2K digits
 * build: cc -O3 -march=native -o tsieve tsieve.c   (-DCHECK: verify every leaf and every digit rule directly)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned long long u64; typedef long long i64; typedef unsigned __int128 u128;

static int K; static u64 M1, PK; static i64 V[3][3]; static u64 WEPS[3];
static u64 ROOT[4][2];
static u64 LP[20], L5[20], LC[20], LU[20], LG[20][5][2];
static i64 N, N0, AMIN; static int DUMP;
static unsigned char OK4[10000], OKR[10000]; static int NFULL, NREM;
static double MINV[3][3];
static i64 NQ, CM[5][3], BETA[5], NA[5];            /* inequalities CM.i <= BETA - NA*r */
static u64 c_nodes, c_leaves, c_rows, c_boxpairs, c_points, c_surv;

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

/* (a_lo + M1 a_hi)(b_lo + M1 b_hi) mod M1^2, in the same representation */
static inline void mul2(u64 alo, u64 ahi, u64 blo, u64 bhi, u64 *rlo, u64 *rhi) {
    u128 p = (u128)alo * blo; u64 q = (u64)(p / M1); *rlo = (u64)(p - (u128)q * M1);
    u128 s = (u128)q + (u128)alo * bhi + (u128)ahi * blo;
    *rhi = (u64)(s % M1);
}
static inline i64 fdiv(i64 a, i64 b) { i64 q = a / b; if ((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }
static inline i64 cdiv(i64 a, i64 b) { return -fdiv(-a, b); }
static inline u64 modm(i64 x) { i64 r = x % (i64)M1; return (u64)(r < 0 ? r + (i64)M1 : r); }
static inline int okblock(u64 v) {
    for (int c = 0; c < NFULL; c++) { if (!OK4[v % 10000]) return 0; v /= 10000; }
    return NREM ? OKR[v] : 1;
}
static int okdigits(u64 v, int n) { for (int i = 0; i < n; i++, v /= 10) if (v % 10 == 0 || v % 10 == 5) return 0; return 1; }

#ifdef CHECK
static u64 M2lo_pow(u64 base, u64 e, u64 *hi) {   /* base^e mod M1^2 (base < M1) in lo/hi form */
    u64 rl = 1, rh = 0, bl = base % M1, bh = base / M1;
    while (e) { if (e & 1) mul2(rl, rh, bl, bh, &rl, &rh); mul2(bl, bh, bl, bh, &bl, &bh); e >>= 1; }
    *hi = rh; return rl;
}
static u64 powmod64(u64 b, u64 e, u64 m) { u128 r = 1 % m, x = b % m; while (e) { if (e & 1) r = r * x % m; x = x * x % m; e >>= 1; } return (u64)r; }
#endif

static void leaf(u64 r, u64 alo, u64 ahi) {
    c_leaves++;
#ifdef CHECK
    {   /* y = 2^(r-2K): 2^(2K) * y == 2^r (mod 5^(2K)); 2^r' (r' = r mod P_K, r' >= K) has K admissible digits */
        u64 thi, tlo = M2lo_pow(2, 2 * K, &thi), plo, phi, qlo, qhi;
        mul2(alo, ahi, tlo, thi, &plo, &phi); qlo = M2lo_pow(2, r, &qhi);
        if (plo != qlo || phi != qhi) { fprintf(stderr, "CHECK: leaf state wrong at r=%llu\n", r); exit(4); }
        u64 m = 1; for (int i = 0; i < K; i++) m *= 10;
        u64 rr = r; while (rr < (u64)K) rr += PK;
        if (!okdigits(powmod64(2, rr, m), K)) { fprintf(stderr, "CHECK: leaf r=%llu not admissible\n", r); exit(4); }
    }
#endif
    u64 h0 = (alo << K) / M1;
    u64 d[3]; for (int m = 0; m < 3; m++) d[m] = (u64)((u128)alo * WEPS[m] % M1);
    /* bounding box of the region in lattice coordinates (outer simplex vertices; exact tests follow) */
    double vx[4][3] = {{(double)AMIN, 0, 0}, {3.0 * N, 0, 0}, {(double)AMIN, (6.0 * N - 2.0 * AMIN) / 3, 0},
                       {(double)AMIN, 0, (6.0 * N - 2.0 * AMIN) / 6}};
    double lo0 = 1e300, hi0 = -1e300, lo1 = 1e300, hi1 = -1e300;
    for (int q = 0; q < 4; q++) {
        double p[3] = {vx[q][0] - (double)r, vx[q][1], vx[q][2]};
        double i0 = MINV[0][0] * p[0] + MINV[0][1] * p[1] + MINV[0][2] * p[2];
        double i1 = MINV[1][0] * p[0] + MINV[1][1] * p[1] + MINV[1][2] * p[2];
        lo0 = fmin(lo0, i0); hi0 = fmax(hi0, i0); lo1 = fmin(lo1, i1); hi1 = fmax(hi1, i1);
    }
    i64 a0 = (i64)floor(lo0) - 1, b0 = (i64)ceil(hi0) + 1, a1 = (i64)floor(lo1) - 1, b1 = (i64)ceil(hi1) + 1;
    i64 rhs0[5]; for (int q = 0; q < NQ; q++) rhs0[q] = BETA[q] - NA[q] * (i64)r;
    for (i64 i0 = a0; i0 <= b0; i0++) for (i64 i1 = a1; i1 <= b1; i1++) {
        c_boxpairs++;
        i64 lo = -(1LL << 62), hi = 1LL << 62; int empty = 0;
        for (int q = 0; q < NQ; q++) {
            i64 rhs = rhs0[q] - CM[q][0] * i0 - CM[q][1] * i1, c3 = CM[q][2];
            if (c3 > 0) { i64 t = fdiv(rhs, c3); if (t < hi) hi = t; }
            else if (c3 < 0) { i64 t = cdiv(rhs, c3); if (t > lo) lo = t; }
            else if (rhs < 0) { empty = 1; break; }
        }
        if (empty || lo > hi) continue;
        c_rows++; c_points += (u64)(hi - lo + 1);
        u64 Y = (u64)(((u128)ahi + (u128)modm(i0) * d[0] + (u128)modm(i1) * d[1] + (u128)modm(lo) * d[2]) % M1);
        for (i64 i2 = lo; i2 <= hi; i2++) {
            u64 v = h0 + (Y << K);
            if (okblock(v)) {
                c_surv++;
                if (DUMP) {
                    i64 a = (i64)r + i0 * V[0][0] + i1 * V[1][0] + i2 * V[2][0];
                    i64 b = i0 * V[0][1] + i1 * V[1][1] + i2 * V[2][1];
                    i64 c = i0 * V[0][2] + i1 * V[1][2] + i2 * V[2][2];
                    printf("S %lld %lld %lld\n", a, b, c);
                }
            }
            Y += d[2]; if (Y >= M1) Y -= M1;
        }
    }
}

static const int INV5[5] = {0, 1, 3, 2, 4};
static void node(int j, u64 r, u64 ylo, u64 yhi) {
    c_nodes++;
    if (j == K) { leaf(r, ylo, yhi); return; }
    u64 m = L5[j] * 5;                                   /* 5^(j+1) */
    u64 Z = (u64)((u128)(ylo % m) * LC[j] % m);           /* 2^(r - j - 1) mod 5^(j+1): digit j+1 = floor(2Z/5^j) */
    u64 zlo = Z % L5[j], q = Z / L5[j];
    int p = 2 * zlo >= L5[j];
    int v = (int)(zlo % 5 * LU[j] % 5);
    int tbad = (int)(((p ? 2 : 0) + 5 - q) % 5) * INV5[v] % 5;
    for (int t = 0; t < 5; t++) {
        if (t == tbad) continue;
        u64 cl = ylo, ch = yhi;
        if (t) mul2(ylo, yhi, LG[j][t][0], LG[j][t][1], &cl, &ch);
#ifdef CHECK
        {   u64 Zt = (u64)((u128)(cl % m) * LC[j] % m), dig = 2 * Zt / L5[j];
            if (dig != (u64)p + 2 * ((q + (u64)t * v) % 5) || dig == 0 || dig == 5) { fprintf(stderr, "CHECK: digit rule\n"); exit(4); }
            u64 Zb; mul2(ylo, yhi, LG[j][tbad][0], LG[j][tbad][1], &Zb, &(u64){0});
            u64 db = 2 * ((u64)((u128)(Zb % m) * LC[j] % m)) / L5[j];
            if (db != 0 && db != 5) { fprintf(stderr, "CHECK: excluded child admissible\n"); exit(4); } }
#endif
        node(j + 1, r + (u64)t * LP[j], cl, ch);
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: tsieve PARAMFILE N N0 AMIN [dump]\n"); return 1; }
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
    N = atoll(argv[2]); N0 = atoll(argv[3]); AMIN = atoll(argv[4]); DUMP = argc > 5 && !strcmp(argv[5], "dump");
    if (AMIN < 2 * K) { fprintf(stderr, "AMIN must be >= 2K\n"); return 1; }
    for (int x = 0; x < 10000; x++) OK4[x] = okdigits(x, 4);
    NFULL = K / 4; NREM = K % 4;
    for (int x = 0; x < 10000; x++) OKR[x] = NREM ? (x < (int)pow(10, NREM) && okdigits(x, NREM)) : 0;
    /* lattice coordinates: p - (r,0,0) = i0 V0 + i1 V1 + i2 V2 */
    double A[3][3]; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) A[i][j] = (double)V[j][i];
    double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
               + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
        int i1 = (j + 1) % 3, i2 = (j + 2) % 3, j1 = (i + 1) % 3, j2 = (i + 2) % 3;
        MINV[i][j] = (A[i1][j1] * A[i2][j2] - A[i1][j2] * A[i2][j1]) / det;
    }
    i64 nv[5][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {2, 3, 6}, {-2, -3, -6}};
    i64 bv[5] = {-AMIN, 0, 0, 6 * N, -6 * N0 - 1};
    NQ = N0 >= 0 ? 5 : 4;
    for (int q = 0; q < NQ; q++) {
        for (int m = 0; m < 3; m++) CM[q][m] = nv[q][0] * V[m][0] + nv[q][1] * V[m][1] + nv[q][2] * V[m][2];
        BETA[q] = bv[q]; NA[q] = nv[q][0];
    }
    clock_t t0 = clock();
    for (int r = 0; r < 4; r++) node(1, (u64)r, ROOT[r][0], ROOT[r][1]);
    double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("DONE K=%d N=%lld N0=%lld AMIN=%lld nodes=%llu leaves=%llu boxpairs=%llu rows=%llu points=%llu surv=%llu cpu=%.3f\n",
           K, N, N0, AMIN, c_nodes, c_leaves, c_boxpairs, c_rows, c_points, c_surv, sec);
    return 0;
}
