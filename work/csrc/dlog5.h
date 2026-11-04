/* 5-adic discrete logarithm base 3 (3 is a primitive root mod 5^n, n <= 27 so that 5^n < 2^63), shared by csrc/mitm_bottom.c
 * (the bottom-window table) and csrc/mitm_exact.h (exact decision of candidate pairs). Call dlog_init(n), then
 * dlog_fast_init(); dlog3(x) returns e in [0, 4 5^(n-1)) with 3^e = x mod 5^n. */
#ifndef DLOG5_H
#define DLOG5_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static u64 invmod(u64 a, u64 m) {                              /* m = 5^n, a coprime to 5 */
    i128 t = 0, nt = 1, r = m, nr = a % m;
    while (nr) { i128 q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    if (t < 0) t += m;
    return (u64)t;
}

/* Pohlig-Hellman: e in [0, 4 5^(n-1)) with 3^e = x mod 5^n */
static u64 P5[40], G4inv[40][5], GM[5], T4[4];
static int NN;
static void dlog_init(int n) {
    NN = n; P5[0] = 1; for (int i = 1; i < 40; i++) P5[i] = P5[i - 1] * 5;
    u64 M = P5[n];
    for (int e = 0; e < 4; e++) T4[e] = powmod(3, (u64)e * P5[n - 1], M);
    u64 g = powmod(3, 4, M), gi = invmod(g, M); int m = n - 1;
    if (m > 0) { u64 gm = powmod(g, P5[m - 1], M); for (int d = 0; d < 5; d++) GM[d] = powmod(gm, d, M); }
    for (int i = 0; i < m; i++) for (int d = 0; d < 5; d++) G4inv[i][d] = powmod(gi, (u64)d * P5[i], M);
}
static u64 dlog3_ph(u64 x) {
    int n = NN, m = n - 1; u64 M = P5[n];
    u64 xp = powmod(x, P5[n - 1], M); int e4 = -1;
    for (int e = 0; e < 4; e++) if (T4[e] == xp) e4 = e;
    if (e4 < 0) { fprintf(stderr, "dlog: no order-4 part\n"); exit(2); }
    u64 h = powmod(x, 4, M), e5 = 0;
    for (int i = 0; i < m; i++) {
        u64 hi = h; for (int s = 0; s < m - 1 - i; s++) { u64 x2 = mulmod(hi, hi, M), x4 = mulmod(x2, x2, M); hi = mulmod(x4, hi, M); }
        int d = -1; for (int dd = 0; dd < 5; dd++) if (GM[dd] == hi) d = dd;
        if (d < 0) { fprintf(stderr, "dlog: digit not found\n"); exit(2); }
        e5 += (u64)d * P5[i]; h = mulmod(h, G4inv[i][d], M);
    }
    /* e = e5 mod 5^m, e = e4 mod 4 */
    u64 e = e5; while (e % 4 != (u64)e4) e += P5[m];
    return e;
}
/* The same discrete log in about n products: e mod 4 from x mod 5 (3^e mod 5 = 1, 3, 4, 2 for e = 0 .. 3); e mod 5^(n-1) from
 * h = x^4 = 81^e in 1 + 5Z: for s = 1 .. n-1, h = 1 + c 5^s (mod 5^(s+1)) with c in 0..4; dividing h by 1 + c 5^s keeps
 * h = 1 mod 5^(s+1) and adds LG[s][c] = log_81(1 + c 5^s) (mod 5^(n-1), tabulated by dlog3_ph at start). Env MB_CHECK=1
 * compares every value with dlog3_ph. */
static u64 LG[40][5], LGINV[40][5]; static int E4OF[5] = {-1, 0, 3, 1, 2}, MBCHECK;
static void dlog_fast_init(void) {
    int n = NN, m = n - 1; u64 M = P5[n];
    for (int s = 1; s < n; s++) for (int c = 1; c < 5; c++) {
        u64 y = 1 + (u64)c * P5[s];                              /* < 5^(s+1) <= M */
        u64 d = dlog3_ph(y);                                     /* y = 3^d, y = 1 mod 5: 4 | d, y = 81^(d/4) */
        if (d % 4) { fprintf(stderr, "dlog table: 4 does not divide\n"); exit(2); }
        LG[s][c] = (d / 4) % P5[m]; LGINV[s][c] = invmod(y, M);
    }
}
static u64 dlog3(u64 x) {
    int n = NN, m = n - 1; u64 M = P5[n];
    int e4 = E4OF[x % 5];
    if (e4 < 0) { fprintf(stderr, "dlog: 5 divides x\n"); exit(2); }
    u64 x2 = mulmod(x, x, M), h = mulmod(x2, x2, M), e5 = 0;
    for (int s = 1; s < n; s++) {
        const int c = (int)(((h - 1) / P5[s]) % 5);
        if (c) { h = mulmod(h, LGINV[s][c], M); e5 += LG[s][c]; if (e5 >= P5[m]) e5 -= P5[m]; }
    }
    if (h != 1 % M) { fprintf(stderr, "dlog: reduction failed\n"); exit(2); }
    u64 e = e5; while (e % 4 != (u64)e4) e += P5[m];
    if (MBCHECK && e != dlog3_ph(x)) { fprintf(stderr, "dlog: fast %lu != Pohlig-Hellman\n", e); exit(2); }
    return e;
}

static long long rnd_div(i128 n, i128 d) {                     /* nearest integer to n/d (d != 0), ties up */
    if (d < 0) { n = -n; d = -d; }
    i128 q = n >= 0 ? (2 * n + d) / (2 * d) : -((-2 * n + d - 1) / (2 * d));
    return (long long)q;
}

#endif
