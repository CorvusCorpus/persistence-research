#define _GNU_SOURCE
// gsieve.c -- chain sieve for multiplicative persistence (final digit 0 and all others).
//
// G = 7-smooth numbers x >= 10 with no digit 0 and never both a digit 5 and an even digit.
// Every intermediate number c0..c8 of a persistence-12 orbit lies in G.
//
// Forms (x not divisible by 10):  (i)  x = 2^a 3^b 7^c        region 2a + 3b + 6c <= 6N
//                                 (iii) x = 3^a 5^b 7^c, b>=1  region a + 2b + 2c <= 2N
// (region = "digit length <= N": a/3 + b/2 + c (+ #5s) <= N is necessary for x = p(n), len(n) <= N)
//
// Paths: main   form (i), a >= 2K  : bitmap of last K digits (suffix collapse), affine lift for
//                                    digits K+1..2K, leading-digit table, base-1e9 limb engine, GMP.
//        direct form (i), a < 2K   : 3 x 18-digit limbs stepped along b, GMP.
//        odd    form (iii)         : b <= 39 stepped along a; b >= 40 via the unique 2-adic odd thread.
// Only x >= 10^XMIN are handled here; smaller x are enumerated exhaustively elsewhere.
//
// Digit modes: G (default), E (zero-free only; validation), P (permissive; counts only).
// Threads: threads=T (or GSIEVE_THREADS=T) workers share the read-only tables and take c values in increasing
// order; per-c output (dump/MEMBER lines, CDONE) is written in increasing c, so logs do not depend on T.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <immintrin.h>
#include <sys/mman.h>
#include <pthread.h>
#include <stdarg.h>
#include "gconst.h"

// Stage-2 variant of the main path (bitmap scan + affine lift + digit checks of digits K+1..2K):
//   S2_IFMA    AVX-512 F/DQ/BW/VL + IFMA: the production build of 19-20 Oct 2026
//   S2_AVX512  AVX-512 F/DQ/BW/VL without IFMA: floor(v/10^4) by exact 32x32-bit multiply-shift
//   S2_AVX2    AVX2: the same exact arithmetic on 4 lanes; survivors extracted with movemask
//   SCALAR     portable reference: per-candidate % and /
// Default: the first one the compiler target supports. All later stages are shared.
#if !defined(SCALAR) && !defined(S2_IFMA) && !defined(S2_AVX512) && !defined(S2_AVX2)
#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512BW__) && defined(__AVX512VL__) && defined(__AVX512IFMA__)
#define S2_IFMA
#elif defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512BW__) && defined(__AVX512VL__)
#define S2_AVX512
#elif defined(__AVX2__)
#define S2_AVX2
#else
#define SCALAR
#endif
#endif
#if defined(S2_IFMA)
#define VARIANT "ifma"
#elif defined(S2_AVX512)
#define VARIANT "avx512"
#elif defined(S2_AVX2)
#define VARIANT "avx2"
#elif defined(SCALAR)
#define VARIANT "scalar"
#endif
#if defined(SCALAR) + defined(S2_IFMA) + defined(S2_AVX512) + defined(S2_AVX2) != 1
#error choose exactly one of SCALAR, S2_IFMA, S2_AVX512, S2_AVX2
#endif

#ifndef KK
#define KK 13
#endif
#if KK == 13
#define P1 K13_P1
#define QQ K13_Q
#define P2 K13_P2
#define Q2 K13_Q2
#define MU3 K13_MU3
#define MU7 K13_MU7
#define WW K13_W
#define Y0 K13_Y0
#define X1 K13_X1
#define TENK K13_TEN
#define TWOK K13_TWO
#elif KK == 12
#define P1 K12_P1
#define QQ K12_Q
#define P2 K12_P2
#define Q2 K12_Q2
#define MU3 K12_MU3
#define MU7 K12_MU7
#define WW K12_W
#define Y0 K12_Y0
#define X1 K12_X1
#define TENK K12_TEN
#define TWOK K12_TWO
#else
#error KK must be 12 or 13
#endif
#ifndef NL
#define NL 14
#endif
#define AMIN (2 * KK)
#define XMIN (9 * NL + 4)
#define E9 1000000000ULL
#define E18 1000000000000000000ULL

typedef uint64_t u64; typedef int64_t i64; typedef uint32_t u32; typedef unsigned __int128 u128;

enum { F_ZERO = 1, F_FIVE = 2, F_EVEN = 4 };
enum { C_EVEN = 0, C_ODD5 = 1, C_ODDN5 = 2 };
static int DM;                 // 0 = G, 1 = E, 2 = P
static int BADEVEN;            // bad-flag mask for even numbers
static uint8_t FL1[10], FL2[100], FL3[1000], FL4[10000], FL5[100000];
static int NOGMP = 0, DUMP = 0, PROG = 0;
static __thread u64 cnt[8];        // per-thread counters (merged per c by the reporter)
static i64 NN, N0 = -1;
static const double LG2 = 0.30102999566398119521, LG3 = 0.47712125471966243730,
                    LG5 = 0.69897000433601880479, LG7 = 0.84509804001425683071;

static uint8_t dflags(int d) { return d == 0 ? F_ZERO : d == 5 ? F_FIVE : (d % 2 == 0) ? F_EVEN : 0; }
static void build_fl(uint8_t *T, int w) {
    int n = 1; for (int i = 0; i < w; i++) n *= 10;
    for (int v = 0; v < n; v++) { int x = v; uint8_t f = 0; for (int i = 0; i < w; i++) { f |= dflags(x % 10); x /= 10; } T[v] = f; }
}
static inline uint8_t flw(u64 v, int w) {     // flags of v written with exactly w digits (zero padded), w in 1..9
    switch (w) {
        case 1: return FL1[v]; case 2: return FL2[v]; case 3: return FL3[v]; case 4: return FL4[v]; case 5: return FL5[v];
        default: return FL4[v % 10000] | flw(v / 10000, w - 4);
    }
}
static inline uint8_t fl9(u64 v) { return FL4[v % 10000] | FL5[v / 10000]; }
static inline uint8_t fl18(u64 v) { return fl9(v % E9) | fl9(v / E9); }
static inline int bad_of(uint8_t f, int cond) {
    if (DM == 2) return 0;
    if (DM == 1) return f & F_ZERO;
    if (cond == C_EVEN) return f & (F_ZERO | F_FIVE);
    if (cond == C_ODD5) return f & (F_ZERO | F_EVEN);
    return (f & F_ZERO) || ((f & (F_FIVE | F_EVEN)) == (F_FIVE | F_EVEN));
}

// ---------------- exact (GMP) check ----------------
static __thread mpz_t gx, gm, gt;
// per-thread output buffer for the current c (dump and MEMBER lines), handed to the reporter when c is done
static __thread char *OB; static __thread size_t OBN, OBC;
static void emit(const char *fmt, ...) {
    for (;;) {
        va_list ap; va_start(ap, fmt);
        int n = vsnprintf(OB + OBN, OBC - OBN, fmt, ap);
        va_end(ap);
        if (n < 0) { fprintf(stderr, "emit failed\n"); exit(1); }
        if ((size_t)n < OBC - OBN) { OBN += (size_t)n; return; }
        OBC = 2 * OBC + (size_t)n + 256;
        if (!(OB = realloc(OB, OBC))) { fprintf(stderr, "alloc failed\n"); exit(1); }
    }
}
static int full_check(int form, u64 a, u64 b, u64 c, int cond) {
    // form 0: 2^a 3^b 7^c ; form 1: 3^a 5^b 7^c
    u64 pr[3] = { form ? 3 : 2, form ? 5 : 3, 7 }, ex[3] = { a, b, c };
    cnt[5]++;
    if (DM == 2 && NOGMP) return 0;
    // partial: last 600 digits
    mpz_ui_pow_ui(gm, 10, 600); mpz_set_ui(gx, 1);
    for (int i = 0; i < 3; i++) { mpz_set_ui(gt, pr[i]); mpz_powm_ui(gt, gt, ex[i], gm); mpz_mul(gx, gx, gt); mpz_mod(gx, gx, gm); }
    // the number has >= XMIN digits; only check as many low digits as it has
    double lg = (form ? a * LG3 + b * LG5 : a * LG2 + b * LG3) + c * LG7;
    int nd = lg > 650 ? 600 : -1;
    if (nd > 0) {
        char buf[700]; mpz_get_str(buf, 10, gx); int L = strlen(buf); uint8_t f = 0;
        if (L < nd) f |= F_ZERO;     // leading zeros inside the window are genuine zero digits
        for (int i = 0; i < L; i++) f |= dflags(buf[i] - '0');
        if (bad_of(f, cond)) return 0;
    }
    cnt[6]++;
    mpz_set_ui(gx, 1);
    for (int i = 0; i < 3; i++) { mpz_ui_pow_ui(gt, pr[i], ex[i]); mpz_mul(gx, gx, gt); }
    char *s = mpz_get_str(NULL, 10, gx); uint8_t f = 0;
    for (char *p = s; *p; p++) f |= dflags(*p - '0');
    int ok = !bad_of(f, cond);
    if (ok) { cnt[7]++; emit("MEMBER form=%d a=%llu b=%llu c=%llu digits=%zu\n", form, (unsigned long long)a, (unsigned long long)b, (unsigned long long)c, strlen(s)); }
    free(s);
    return ok;
}

// ---------------- tables for the main path ----------------
static void *hmalloc(size_t n) {                 // 2 MB aligned, transparent huge pages
    void *p; size_t al = 1 << 21; n = (n + al - 1) / al * al;
    if (posix_memalign(&p, al, n)) { fprintf(stderr, "alloc failed\n"); exit(1); }
    madvise(p, n, MADV_HUGEPAGE); memset(p, 0, n); return p;
}
static u64 *BM; static u32 *RK; static u64 *LO; static u64 *AB; static uint16_t *H0; static u32 *POS; static u64 NGOOD;
static uint8_t OK4B[10000 + 8], OK5B[100000 + 8];
static u64 *TOP;
static u32 (*T2)[NL]; static u32 (*T2H)[4];

static inline uint8_t flK(u64 x) {               // flags of the K-digit window
#if KK == 13
    return FL4[x % 10000] | FL4[(x / 10000) % 10000] | FL5[x / 100000000];
#else
    return FL4[x % 10000] | FL4[(x / 10000) % 10000] | FL4[x / 100000000];
#endif
}
// Bitmap and stage-2 tables, built in parallel over word-aligned chunks of exponents e in [0, P1):
// x(e) = 2^(P1+e) mod 10^K = X1 2^e and y(e) = 2^(e-2K) mod 5^(2K) = Y0 2^e, both by doubling from the
// chunk start (computed by modular exponentiation). Pass 1 sets the bits of the chunk's own words; after the
// rank prefix sums, pass 2 fills AB/H0/POS (LO for SCALAR) from rank RK[w0] on. Same tables for every T.
static u64 pow2mod(u64 base, u64 e, u64 m) {     // base * 2^e mod m, m < 2^62
    u64 r = base % m, p = 2 % m;
    for (; e; e >>= 1) { if (e & 1) r = (u64)((u128)r * p % m); p = (u64)((u128)p * p % m); }
    return r;
}
typedef struct { u64 w0, w1, n; int bad; } bchunk_t;
static void *bitmap_pass1(void *arg) {
    bchunk_t *ch = arg;
    u64 e0 = ch->w0 * 64, e1 = ch->w1 * 64 < P1 ? ch->w1 * 64 : P1, n = 0;
    u64 x = pow2mod(X1, e0, TENK);
    for (u64 e = e0; e < e1; e++) {
        if (!bad_of(flK(x), C_EVEN)) { BM[e >> 6] |= 1ULL << (e & 63); n++; }
        x <<= 1; if (x >= TENK) x -= TENK;
    }
    ch->n = n; return NULL;
}
static void *bitmap_pass2(void *arg) {
    bchunk_t *ch = arg;
    u64 e0 = ch->w0 * 64, e1 = ch->w1 * 64 < P1 ? ch->w1 * 64 : P1, i = RK[ch->w0];
    u64 y = pow2mod(Y0, e0, Q2);
    for (u64 e = e0; e < e1; e++) {
        if ((BM[e >> 6] >> (e & 63)) & 1) {
#ifdef SCALAR
            LO[i] = y;
#else
            u64 alo = y % QQ, ahi = y / QQ;
            AB[i] = ahi | (((u64)(((u128)((alo * WW) % QQ) << 32) % QQ)) << 32); H0[i] = (uint16_t)((TWOK * alo) / QQ); POS[i] = (u32)e;
#endif
            i++;
        }
        y <<= 1; if (y >= Q2) y -= Q2;
    }
    ch->bad = i != RK[ch->w1]; return NULL;
}
static void run_chunks(void *(*fn)(void *), bchunk_t *ch, int nth) {
    pthread_t *tid = malloc(nth * sizeof *tid);
    for (int t = 0; t < nth; t++) if (pthread_create(&tid[t], NULL, fn, &ch[t])) { fprintf(stderr, "pthread_create failed\n"); exit(1); }
    for (int t = 0; t < nth; t++) pthread_join(tid[t], NULL);
    free(tid);
}
static u64 fnv(u64 h, const void *p, size_t n) {
    const uint8_t *q = p;
    for (size_t i = 0; i < n; i++) { h ^= q[i]; h *= 1099511628211ULL; }
    return h;
}
static void build_bitmap(int nth) {
    u64 nw = (P1 + 63) / 64;
    BM = hmalloc((nw + 1) * 8); RK = hmalloc((nw + 1) * 4);
    bchunk_t *ch = calloc(nth, sizeof *ch);
    for (int t = 0; t < nth; t++) { ch[t].w0 = nw * t / nth; ch[t].w1 = nw * (t + 1) / nth; }
    run_chunks(bitmap_pass1, ch, nth);
    u64 n = 0;
    for (int t = 0; t < nth; t++) n += ch[t].n;
    NGOOD = n;
#ifdef SCALAR
    LO = malloc((n + 1) * 8);
#else
    AB = hmalloc((n + 8) * 8); H0 = hmalloc((n + 8) * 2); POS = hmalloc((n + 8) * 4);
#endif
    u64 run = 0;
    for (u64 w = 0; w < nw; w++) { RK[w] = (u32)run; run += __builtin_popcountll(BM[w]); }
    RK[nw] = (u32)run;
    if (run != n) { fprintf(stderr, "bitmap count mismatch\n"); exit(1); }
    run_chunks(bitmap_pass2, ch, nth);
    for (int t = 0; t < nth; t++) if (ch[t].bad) { fprintf(stderr, "bitmap fill mismatch\n"); exit(1); }
    free(ch);
    for (int v = 0; v < 10000; v++) OK4B[v] = !bad_of(FL4[v], C_EVEN);
    for (int v = 0; v < 100000; v++) OK5B[v] = !bad_of(FL5[v], C_EVEN);
    if (getenv("GSIEVE_TABLE_HASH")) {          // fingerprint of the tables (for checking the parallel build)
        u64 h = fnv(1469598103934665603ULL, BM, (nw + 1) * 8); h = fnv(h, RK, (nw + 1) * 4);
#ifdef SCALAR
        h = fnv(h, LO, n * 8);
#else
        h = fnv(h, AB, n * 8); h = fnv(h, H0, n * 2); h = fnv(h, POS, n * 4);
#endif
        fprintf(stderr, "table hash %016llx (n=%llu)\n", (unsigned long long)h, (unsigned long long)n);
    }
}
#ifndef TB
#define TB 18
#endif
#ifndef TD
#define TD 5
#endif
static void build_top(void) {
    TOP = hmalloc(((1ULL << TB) / 64 + 1) * 8);
    if (DM == 2) return;
    long double d = ldexpl(1.0L, -40);
    u64 lo10 = 1; for (int i = 1; i < TD; i++) lo10 *= 10;
    for (u64 k = 1; k + 1 < (1ULL << TB); k++) {
        long double flo = ldexpl((long double)k, -TB) - d, fhi = ldexpl((long double)(k + 1), -TB) + d;
        long double vlo = exp10l(flo + (TD - 1)) * (1.0L - 1e-15L), vhi = exp10l(fhi + (TD - 1)) * (1.0L + 1e-15L);
        u64 mlo = (u64)floorl(vlo), mhi = (u64)floorl(vhi);
        if (mlo < lo10 || mhi >= 10 * lo10) continue;
        int bad = 1;
        for (u64 m = mlo; m <= mhi; m++) if (!bad_of(FL4[m % 10000] | flw(m / 10000, TD - 4), C_EVEN)) { bad = 0; break; }
        if (bad) TOP[k >> 6] |= 1ULL << (k & 63);
    }
}
static void build_t2(u64 amax) {
    T2 = hmalloc((amax + 1) * sizeof *T2); T2H = hmalloc((amax + 1) * sizeof *T2H);
    memset(T2[0], 0, sizeof T2[0]); T2[0][0] = 1;
    for (u64 a = 1; a <= amax; a++) {
        u64 cy = 0;
        for (int i = 0; i < NL; i++) { u64 v = 2ULL * T2[a - 1][i] + cy; cy = v >= E9; T2[a][i] = (u32)(v - cy * E9); }
    }
    for (u64 a = 0; a <= amax; a++) for (int i = 0; i < 4; i++) T2H[a][i] = T2[a][i];
}
static inline void mulsmall(u32 *M, u64 k) {
    u64 cy = 0;
    for (int i = 0; i < NL; i++) { u64 v = k * M[i] + cy; cy = v / E9; M[i] = (u32)(v - cy * E9); }
}
// digits [2K, 9NL) of (T2[a] * M) mod 10^(9NL); returns 1 if all allowed
static inline int limb_check(const u32 *A, const u32 *B, const u32 *AH) {
    // digits [2K, 9NL) of (2^a * M) mod 10^(9NL); limbs base 1e9; first 4 limbs straight-line
    u64 a0 = AH[0], a1 = AH[1], a2 = AH[2], a3 = AH[3], b0 = B[0], b1 = B[1], b2 = B[2], b3 = B[3];
    u64 acc = a0 * b0, cy = acc / E9;
    acc = a0 * b1 + a1 * b0 + cy; cy = acc / E9;
    acc = a0 * b2 + a1 * b1 + a2 * b0 + cy; cy = acc / E9; u64 l2 = acc - cy * E9;
    acc = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0 + cy; cy = acc / E9; u64 l3 = acc - cy * E9;
#if KK == 13
    uint8_t f = FL1[l2 / 100000000ULL] | fl9(l3);
#else
    uint8_t f = FL3[l2 / 1000000ULL] | fl9(l3);
#endif
    if (bad_of(f, C_EVEN)) return 0;
    for (int k = 4; k < NL; k++) {
        acc = cy;
        for (int u = 0; u <= k; u++) acc += (u64)A[u] * B[k - u];
        cy = acc / E9; u64 limb = acc - cy * E9;
        if (bad_of(fl9(limb), C_EVEN)) return 0;
    }
    return 1;
}

// ---------------- main path ----------------
#ifdef SCALAR
static inline void scan(u64 s, u64 e, u64 j, i64 delta, u64 b, u64 c, u64 tb, const u32 *M) {
    u64 w0 = s >> 6, w1 = (e - 1) >> 6;
    u64 idx = RK[w0] + __builtin_popcountll(BM[w0] & ((1ULL << (s & 63)) - 1));
    for (u64 w = w0; w <= w1; w++) {
        u64 bits = BM[w];
        if (w == w0) bits &= ~0ULL << (s & 63);
        if (w == w1 && (e & 63)) bits &= (1ULL << (e & 63)) - 1;
#ifdef BENCH0
        cnt[1] += __builtin_popcountll(bits); continue;
#endif
#ifdef BENCH1
        while (bits) { bits &= bits - 1; cnt[2] += LO[idx++]; } continue;
#endif
        while (bits) {
            int tz = __builtin_ctzll(bits); bits &= bits - 1;
            u64 p = (w << 6) | (u64)tz;
            u64 y = LO[idx++];
            cnt[1]++;
            u64 alo = y % QQ, ahi = y / QQ;
            u64 h0 = (TWOK * alo) / QQ, beta = (alo * WW) % QQ;
            u64 Y = (ahi + beta * j) % QQ;
            u64 xh = h0 + TWOK * Y;
            if (bad_of(flK(xh), C_EVEN)) continue;
            cnt[2]++;
            u64 a = (u64)((i64)p + delta);
            if (DUMP) emit("S2 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
            u64 t = a * L2 + tb;
            u64 bin = t >> (64 - TB);
            if ((TOP[bin >> 6] >> (bin & 63)) & 1) { if (DUMP) emit("R3 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
            cnt[3]++;
            if (!limb_check(T2[a], M, T2H[a])) { if (DUMP) emit("R4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
            cnt[4]++;
            if (DUMP) emit("S4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
            full_check(0, a, b, c, C_EVEN);
        }
    }
}
#endif
#ifndef SCALAR
static __thread u64 *PBUF, *PA, *PBIN, *HBUF;      // per-thread survivor buffers
#if KK == 13
#define QINV 163364931ULL
#else
#define QINV 816824655ULL
#endif
// Exact floor(v / 10^4) without IFMA: M4 = ceil(2^45 / 10^4) = 3518437209 < 2^32 and M4 * 10^4 - 2^45 = 1168,
// so floor(v * M4 / 2^45) = floor(v / 10^4) for 0 <= v < 2^45 / 1168 (> 3.0e10). Every v divided this way
// is < 2^31, so v * M4 < 2^63. The 13-digit block xh = h0 + 2^K Y (< 2^44) is split as Y = 10^4 sq + r:
// xh mod 10^4 = w mod 10^4 and floor(xh / 10^4) = 2^K sq + floor(w / 10^4), with w = h0 + 2^K r < 2^K 10^4.
#define M4 3518437209ULL
static inline void later_stages(u64 np, i64 delta, u64 b, u64 c, u64 tb, const u32 *M);
#if defined(S2_IFMA) || defined(S2_AVX512)
#ifdef S2_IFMA
#define DIV4(v) _mm512_srli_epi64(_mm512_madd52hi_epu64(z, v, mag), 8)
#else
#define DIV4(v) _mm512_srli_epi64(_mm512_mul_epu32(v, vm4), 45)
#endif
static inline void scan(u64 s, u64 e, u64 j, i64 delta, u64 b, u64 c, u64 tb, const u32 *M) {
    u64 w0 = s >> 6, w1 = (e - 1) >> 6;
    u64 r0 = RK[w0] + __builtin_popcountll(BM[w0] & ((1ULL << (s & 63)) - 1));
    u64 r1 = (e & 63) ? RK[w1] + __builtin_popcountll(BM[w1] & ((1ULL << (e & 63)) - 1)) : RK[w1 + 1];
    if (r1 <= r0) return;
    cnt[1] += r1 - r0;
    const __m512i vq = _mm512_set1_epi64((long long)QQ), vj = _mm512_set1_epi64((long long)j);
    const __m512i vqinv = _mm512_set1_epi64((long long)QINV), m32 = _mm512_set1_epi64(0xFFFFFFFFLL);
    const __m512i v4 = _mm512_set1_epi64(10000), z = _mm512_setzero_si512();
#ifdef S2_IFMA
    const __m512i mag = _mm512_set1_epi64(115292150460685LL);
#else
    const __m512i vm4 = _mm512_set1_epi64((long long)M4);
#endif
    u64 np = 0, nh = 0;
    // phase A: affine lift, low 4 digits of the 13-digit block
    for (u64 i = r0; i < r1; i += 8) {
        __mmask8 lm = (r1 - i >= 8) ? 0xFF : (__mmask8)((1u << (r1 - i)) - 1);
        __m512i ab = _mm512_maskz_loadu_epi64(lm, AB + i);
        __m512i h0 = _mm512_cvtepu16_epi64(_mm_maskz_loadu_epi16(lm, H0 + i));
        __m512i t = _mm512_mul_epu32(_mm512_srli_epi64(ab, 32), vj);
        __m512i mm = _mm512_and_si512(_mm512_mul_epu32(t, vqinv), m32);
        __m512i u = _mm512_srli_epi64(_mm512_add_epi64(t, _mm512_mul_epu32(mm, vq)), 32);
        u = _mm512_mask_sub_epi64(u, _mm512_cmpge_epu64_mask(u, vq), u, vq);
        __m512i Y = _mm512_add_epi64(u, _mm512_and_si512(ab, m32));
        Y = _mm512_mask_sub_epi64(Y, _mm512_cmpge_epu64_mask(Y, vq), Y, vq);
#ifdef S2_IFMA
        __m512i xh = _mm512_add_epi64(h0, _mm512_slli_epi64(Y, KK));
        __m512i q1 = DIV4(xh);
        __m512i d1 = _mm512_sub_epi64(xh, _mm512_mul_epu32(q1, v4));
#else
        __m512i sq = DIV4(Y);
        __m512i w = _mm512_add_epi64(h0, _mm512_slli_epi64(_mm512_sub_epi64(Y, _mm512_mul_epu32(sq, v4)), KK));
        __m512i t1 = DIV4(w);
        __m512i d1 = _mm512_sub_epi64(w, _mm512_mul_epu32(t1, v4));
        __m512i q1 = _mm512_add_epi64(_mm512_slli_epi64(sq, KK), t1);
#endif
        __m256i g1 = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lm, d1, OK4B, 1);
        __mmask8 pass = _mm256_test_epi32_mask(g1, _mm256_set1_epi32(0xFF)) & lm;
        // keep (q1, position) of survivors: q1 < 2^30, position < 2^32 -> pack into one word
        __m512i pk = _mm512_or_si512(_mm512_slli_epi64(q1, 32), _mm512_cvtepu32_epi64(_mm256_maskz_loadu_epi32(pass, POS + i)));
        _mm512_mask_compressstoreu_epi64(HBUF + nh, pass, pk);
        nh += __builtin_popcount(pass);
    }
    // phase B: remaining 9 digits (KK=13) / 8 digits (KK=12)
    for (u64 i = 0; i < nh; i += 8) {
        __mmask8 lm = (nh - i >= 8) ? 0xFF : (__mmask8)((1u << (nh - i)) - 1);
        __m512i pk = _mm512_maskz_loadu_epi64(lm, HBUF + i);
        __m512i q1 = _mm512_srli_epi64(pk, 32);
        __m512i q2 = DIV4(q1);
        __m512i d2 = _mm512_sub_epi64(q1, _mm512_mul_epu32(q2, v4));
        __m256i g2 = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lm, d2, OK4B, 1);
#if KK == 13
        __m512i q3 = DIV4(q2);
        __m512i d3 = _mm512_sub_epi64(q2, _mm512_mul_epu32(q3, v4));
        __m256i g3 = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lm, d3, OK4B, 1);
        __mmask8 top = _mm512_cmpneq_epu64_mask(q3, z);
        if (BADEVEN & F_FIVE) top &= _mm512_cmpneq_epu64_mask(q3, _mm512_set1_epi64(5));
        if (DM == 2) top = 0xFF;
#else
        __m256i g3 = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lm, q2, OK4B, 1);
        __mmask8 top = 0xFF;
#endif
        __m256i g = _mm256_and_si256(_mm256_and_si256(g2, g3), _mm256_set1_epi32(0xFF));
        __mmask8 pass = _mm256_test_epi32_mask(g, g) & lm & top;
        _mm512_mask_compressstoreu_epi64(PBUF + np, pass, _mm512_and_si512(pk, m32));
        np += __builtin_popcount(pass);
    }
#ifdef BENCH2
    cnt[2] += np; return;
#endif
    later_stages(np, delta, b, c, tb, M);
}
#undef DIV4
#elif defined(S2_AVX2)
// left-pack of the selected 64-bit lanes of a 4-lane vector (AVX2 has no compress): PACK4[m] lists, for each
// lane l set in m in increasing order, the 32-bit source indices 2l, 2l+1
static uint32_t PACK4[16][8] __attribute__((aligned(32)));
static void build_pack4(void) {
    for (int m = 0; m < 16; m++) {
        int k = 0;
        for (int l = 0; l < 4; l++) if (m >> l & 1) { PACK4[m][k++] = 2 * l; PACK4[m][k++] = 2 * l + 1; }
        while (k < 8) { PACK4[m][k] = k; k++; }
    }
}
static inline void scan(u64 s, u64 e, u64 j, i64 delta, u64 b, u64 c, u64 tb, const u32 *M) {
    u64 w0 = s >> 6, w1 = (e - 1) >> 6;
    u64 r0 = RK[w0] + __builtin_popcountll(BM[w0] & ((1ULL << (s & 63)) - 1));
    u64 r1 = (e & 63) ? RK[w1] + __builtin_popcountll(BM[w1] & ((1ULL << (e & 63)) - 1)) : RK[w1 + 1];
    if (r1 <= r0) return;
    cnt[1] += r1 - r0;
    const __m256i vq = _mm256_set1_epi64x((long long)QQ), vqm1 = _mm256_set1_epi64x((long long)QQ - 1);
    const __m256i vj = _mm256_set1_epi64x((long long)j), vqinv = _mm256_set1_epi64x((long long)QINV);
    const __m256i m32 = _mm256_set1_epi64x(0xFFFFFFFFLL), vm4 = _mm256_set1_epi64x((long long)M4);
    const __m256i v4 = _mm256_set1_epi64x(10000), z = _mm256_setzero_si256();
    const __m128i lowbyte = _mm_set1_epi32(0xFF), z128 = _mm_setzero_si128();
    const __m256i lanes = _mm256_set_epi64x(3, 2, 1, 0);
#if KK == 13
    const __m256i five = _mm256_set1_epi64x(5);
    const int fivebad = (BADEVEN & F_FIVE) != 0, nocheck = DM == 2;
#endif
    u64 np = 0, nh = 0;
    // phase A: affine lift, low 4 digits of the block. AB, H0 and POS have >= 8 zero-filled padding entries
    // beyond NGOOD, so the unmasked loads stay in bounds and every gather index is a remainder mod 10^4.
    for (u64 i = r0; i < r1; i += 4) {
        unsigned lm = (r1 - i >= 4) ? 0xFu : (1u << (r1 - i)) - 1;
        __m256i ab = _mm256_loadu_si256((const __m256i *)(AB + i));
        __m256i h0 = _mm256_cvtepu16_epi64(_mm_loadl_epi64((const __m128i *)(H0 + i)));
        __m256i t = _mm256_mul_epu32(_mm256_srli_epi64(ab, 32), vj);
        __m256i mm = _mm256_mul_epu32(t, vqinv);                  // low 32 bits: t * QINV mod 2^32
        __m256i u = _mm256_srli_epi64(_mm256_add_epi64(t, _mm256_mul_epu32(mm, vq)), 32);
        u = _mm256_sub_epi64(u, _mm256_and_si256(_mm256_cmpgt_epi64(u, vqm1), vq));     // u < 2Q < 2^63: signed compare exact
        __m256i Y = _mm256_add_epi64(u, _mm256_and_si256(ab, m32));
        Y = _mm256_sub_epi64(Y, _mm256_and_si256(_mm256_cmpgt_epi64(Y, vqm1), vq));
        __m256i sq = _mm256_srli_epi64(_mm256_mul_epu32(Y, vm4), 45);
        __m256i w = _mm256_add_epi64(h0, _mm256_slli_epi64(_mm256_sub_epi64(Y, _mm256_mul_epu32(sq, v4)), KK));
        __m256i t1 = _mm256_srli_epi64(_mm256_mul_epu32(w, vm4), 45);
        __m256i d1 = _mm256_sub_epi64(w, _mm256_mul_epu32(t1, v4));
        __m256i q1 = _mm256_add_epi64(_mm256_slli_epi64(sq, KK), t1);
        __m128i g1 = _mm256_i64gather_epi32((const int *)OK4B, d1, 1);
        unsigned pass = ~(unsigned)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_and_si128(g1, lowbyte), z128))) & lm;
        // keep (q1, position) of survivors: q1 < 2^30, position < 2^32 -> one word; left-pack into HBUF
        __m256i pk = _mm256_or_si256(_mm256_slli_epi64(q1, 32), _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *)(POS + i))));
        _mm256_storeu_si256((__m256i *)(HBUF + nh), _mm256_permutevar8x32_epi32(pk, _mm256_load_si256((const __m256i *)PACK4[pass])));
        nh += __builtin_popcount(pass);
    }
    // phase B: remaining 9 digits (KK=13) / 8 digits (KK=12); lanes beyond nh are loaded as 0
    for (u64 i = 0; i < nh; i += 4) {
        __m256i valid = _mm256_cmpgt_epi64(_mm256_set1_epi64x((long long)(nh - i)), lanes);
        __m256i pk = _mm256_maskload_epi64((const long long *)(HBUF + i), valid);
        unsigned lm = (unsigned)_mm256_movemask_pd(_mm256_castsi256_pd(valid));
        __m256i q1 = _mm256_srli_epi64(pk, 32);
        __m256i q2 = _mm256_srli_epi64(_mm256_mul_epu32(q1, vm4), 45);
        __m256i d2 = _mm256_sub_epi64(q1, _mm256_mul_epu32(q2, v4));
#if KK == 13
        __m256i q3 = _mm256_srli_epi64(_mm256_mul_epu32(q2, vm4), 45);
        __m256i d3 = _mm256_sub_epi64(q2, _mm256_mul_epu32(q3, v4));
        __m256i badtop = _mm256_cmpeq_epi64(q3, z);
        if (fivebad) badtop = _mm256_or_si256(badtop, _mm256_cmpeq_epi64(q3, five));
        unsigned top = nocheck ? 0xFu : ~(unsigned)_mm256_movemask_pd(_mm256_castsi256_pd(badtop)) & 0xFu;
#else
        __m256i d3 = q2;
        unsigned top = 0xFu;
#endif
        __m128i g2 = _mm256_i64gather_epi32((const int *)OK4B, d2, 1);
        __m128i g3 = _mm256_i64gather_epi32((const int *)OK4B, d3, 1);
        __m128i g = _mm_and_si128(_mm_and_si128(g2, g3), lowbyte);
        unsigned pass = ~(unsigned)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(g, z128))) & lm & top;
        __m256i pos = _mm256_and_si256(pk, m32);
        _mm256_storeu_si256((__m256i *)(PBUF + np), _mm256_permutevar8x32_epi32(pos, _mm256_load_si256((const __m256i *)PACK4[pass])));
        np += __builtin_popcount(pass);
    }
#ifdef BENCH2
    cnt[2] += np; return;
#endif
    later_stages(np, delta, b, c, tb, M);
}
#endif
// stages 3-5 for the stage-2 survivors in PBUF (bit positions, ascending)
static inline void later_stages(u64 np, i64 delta, u64 b, u64 c, u64 tb, const u32 *M) {
    for (u64 k = 0; k < np; k++) {
        u64 p = PBUF[k];
        u64 a = (u64)((i64)p + delta);
        u64 bin = (a * L2 + tb) >> (64 - TB);
        PA[k] = a; PBIN[k] = bin;
        // (no prefetch: the small table stays cache resident)
    }
    u64 n3 = 0;
#ifdef BENCH4
    cnt[2] += np; return;
#endif
    for (u64 k = 0; k < np; k++) {
        u64 a = PA[k], bin = PBIN[k];
        if (DUMP) emit("S2 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        if ((TOP[bin >> 6] >> (bin & 63)) & 1) { if (DUMP) emit("R3 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        PA[n3++] = a; __builtin_prefetch(T2H[a]);
    }
    cnt[2] += np; cnt[3] += n3;
#ifdef BENCH3
    return;
#endif
    for (u64 k = 0; k < n3; k++) {
        u64 a = PA[k];
        if (!limb_check(T2[a], M, T2H[a])) { if (DUMP) emit("R4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        cnt[4]++;
        if (DUMP) emit("S4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        full_check(0, a, b, c, C_EVEN);
    }
}
#endif
static void pow7_limbs(u32 *M, u64 c) {        // M = 7^c mod 10^(9 NL) in base-1e9 limbs
    mpz_ui_pow_ui(gm, 10, 9 * NL); mpz_set_ui(gt, 7); mpz_powm_ui(gx, gt, c, gm);
    for (int i = 0; i < NL; i++) M[i] = (u32)mpz_fdiv_q_ui(gx, gx, E9);
}
static void main_c(u64 c) {
    {
        u64 lam = (u64)(((u128)c * MU7) % P2);
        u64 tb = c * L7;
        u32 M[NL]; pow7_limbs(M, c);
        if (6 * NN - 6 * (i64)c - 2 * AMIN < 0) return;
        i64 bmax = (6 * NN - 6 * (i64)c - 2 * AMIN) / 3;
        for (i64 b = 0; b <= bmax; b++) {
            i64 ahi = (6 * NN - 3 * b - 6 * (i64)c) / 2;
            double need = (XMIN - b * LG3 - c * LG7) / LG2;
            i64 alo = need > AMIN ? (i64)ceil(need - 1e-6) : AMIN;
            if (N0 >= 0 && 6 * N0 - 3 * b - 6 * (i64)c >= 0) { i64 t0 = (6 * N0 - 3 * b - 6 * (i64)c) / 2 + 1; if (t0 > alo) alo = t0; }
            if (alo <= ahi) {
                u64 len = (u64)(ahi - alo + 1);
                cnt[0] += len;
                u64 e2 = lam + (u64)alo; if (e2 >= P2) e2 -= P2;
                u64 j = e2 / P1, r = e2 % P1;
                if (r + len <= P1) scan(r, r + len, j, alo - (i64)r, b, c, tb, M);
                else {
                    scan(r, P1, j, alo - (i64)r, b, c, tb, M);
                    u64 j2 = j + 1; if (j2 == QQ) j2 = 0;
                    scan(0, r + len - P1, j2, alo + (i64)(P1 - r), b, c, tb, M);
                }
            }
            lam += MU3; if (lam >= P2) lam -= P2;
            tb += L3;
            mulsmall(M, 3);
        }
    }
}

// ---------------- 3-limb (base 1e18) stepping for direct / odd paths ----------------
static void set3(u64 *L, int form, u64 a, u64 b, u64 c) {
    mpz_ui_pow_ui(gm, 10, 54); mpz_set_ui(gx, 1);
    u64 pr[3] = { form ? 3 : 2, form ? 5 : 3, 7 }, ex[3] = { a, b, c };
    for (int i = 0; i < 3; i++) { mpz_set_ui(gt, pr[i]); mpz_powm_ui(gt, gt, ex[i], gm); mpz_mul(gx, gx, gt); mpz_mod(gx, gx, gm); }
    mpz_t q, e18; mpz_init(q); mpz_init_set_str(e18, "1000000000000000000", 10);
    for (int i = 0; i < 3; i++) { mpz_fdiv_qr(gx, q, gx, e18); L[i] = mpz_get_ui(q); }
    mpz_clear(q); mpz_clear(e18);
}
static inline void mul3(u64 *L) {
    u64 v0 = 3 * L[0], c0 = v0 >= 2 * E18 ? 2 : v0 >= E18 ? 1 : 0; L[0] = v0 - c0 * E18;
    u64 v1 = 3 * L[1] + c0, c1 = v1 >= 2 * E18 ? 2 : v1 >= E18 ? 1 : 0; L[1] = v1 - c1 * E18;
    u64 v2 = 3 * L[2] + c1; v2 -= (v2 >= 2 * E18 ? 2 : v2 >= E18 ? 1 : 0) * E18; L[2] = v2;
}
static inline int pass54(const u64 *L, int cond) {
    if (DM == 2) return 1;
    if (bad_of(FL4[L[0] % 10000], cond == C_ODDN5 ? C_EVEN : cond) && (cond != C_ODDN5 || (FL4[L[0] % 10000] & F_ZERO))) return 0;
    if (cond == C_ODDN5) {          // accumulate flags (bad = zero, or five together with even)
        uint8_t f = 0;
        for (int i = 0; i < 3; i++) { f |= fl18(L[i]); if (bad_of(f, cond)) return 0; }
        return 1;
    }
    for (int i = 0; i < 3; i++) if (bad_of(fl18(L[i]), cond)) return 0;
    return 1;
}
static void direct_c(u64 c) {
        for (u64 a = 0; a < AMIN; a++) {
            if (6 * NN - 2 * (i64)a - 6 * (i64)c < 0) continue;
            i64 bmax = (6 * NN - 2 * (i64)a - 6 * (i64)c) / 3;
            double need = (XMIN - a * LG2 - c * LG7) / LG3;
            i64 blo = need > 0 ? (i64)ceil(need - 1e-6) : 0;
            if (N0 >= 0 && 6 * N0 - 2 * (i64)a - 6 * (i64)c >= 0) { i64 t0 = (6 * N0 - 2 * (i64)a - 6 * (i64)c) / 3 + 1; if (t0 > blo) blo = t0; }
            if (blo > bmax) continue;
            int cond = a == 0 ? C_ODDN5 : C_EVEN;
            u64 L[3]; set3(L, 0, a, (u64)blo, c);
            for (i64 b = blo; b <= bmax; b++) {
                cnt[0]++;
                if (pass54(L, cond)) { cnt[4]++; if (DUMP) emit("D4 0 %llu %lld %llu\n", (unsigned long long)a, (long long)b, (unsigned long long)c); full_check(0, a, (u64)b, c, cond); }
                mul3(L);
            }
        }
}
static void odd_c(u64 c) {
    const u64 MSK = (1ULL << 38) - 1;
    {
        u64 bdirect = (DM == 0) ? 39 : (u64)NN;      // E / P modes: enumerate every b directly
        for (u64 b = 1; b <= bdirect; b++) {
            i64 amax = 2 * NN - 2 * (i64)b - 2 * (i64)c;
            if (amax < 0) break;
            double need = (XMIN - b * LG5 - c * LG7) / LG3;
            i64 alo = need > 0 ? (i64)ceil(need - 1e-6) : 0;
            if (N0 >= 0 && 2 * N0 - 2 * (i64)b - 2 * (i64)c >= 0) { i64 t0 = 2 * N0 - 2 * (i64)b - 2 * (i64)c + 1; if (t0 > alo) alo = t0; }
            if (alo > amax) continue;
            u64 L[3]; set3(L, 1, (u64)alo, b, c);
            u64 x0 = L[0];
            for (i64 a = alo; a <= amax; a++) {
                cnt[0]++;
                if (DM == 2 || (!bad_of(FL4[x0 % 10000], C_ODD5) && !bad_of(fl18(x0), C_ODD5))) {      // low 18 digits pass: check 54 digits exactly
                    set3(L, 1, (u64)a, b, c);
                    if (pass54(L, C_ODD5)) { cnt[4]++; if (DUMP) emit("D4 1 %lld %llu %llu\n", (long long)a, (unsigned long long)b, (unsigned long long)c); full_check(1, (u64)a, b, c, C_ODD5); }
                }
                x0 = 3 * x0; x0 -= (x0 >= 2 * E18 ? 2 : x0 >= E18 ? 1 : 0) * E18;
            }
        }
        if (DM != 0) return;
        // b >= 40: x mod 2^40 = (-1)^(a+c) 5^(b + a f3 + c f7); all-odd forces (a+c) odd and
        // b + a f3 + c f7 == OTHREAD (mod 2^38)  (the unique odd-digit 2-adic thread)
        for (i64 a = 0; 2 * NN - a - 2 * (i64)c >= 80; a++) {
            if (((u64)a + c) % 2 == 0) continue;
            u64 b0 = (OTHREAD - (u64)a * OF3 - c * OF7) & MSK;
            cnt[1]++;
            if (b0 < 40 || (i64)(2 * b0) > 2 * NN - a - 2 * (i64)c) continue;
            if (N0 >= 0 && a + 2 * (i64)b0 + 2 * (i64)c <= 2 * N0) continue;
            if (a * LG3 + b0 * LG5 + c * LG7 < XMIN) continue;     // tiny range handled elsewhere
            cnt[4]++;
            if (DUMP) emit("T4 1 %lld %llu %llu\n", (long long)a, (unsigned long long)b0, (unsigned long long)c);
            full_check(1, (u64)a, b0, c, C_ODD5);
        }
    }
}

// ---------------- threads: dynamic c scheduling, per-c results reported in increasing c ----------------
typedef struct { char *buf; size_t len; u64 d[8]; int ready; } crec_t;
static void (*CFN)(u64);
static u64 C0, C1, NEXTC, NEXTREP, TOT[8];
static crec_t *REC;
static pthread_mutex_t REPMU = PTHREAD_MUTEX_INITIALIZER;
static void report(u64 c, char *buf, size_t len, const u64 *d) {
    pthread_mutex_lock(&REPMU);
    crec_t *r = &REC[c - C0]; r->buf = buf; r->len = len; memcpy(r->d, d, sizeof r->d); r->ready = 1;
    while (NEXTREP < C1 && REC[NEXTREP - C0].ready) {
        crec_t *q = &REC[NEXTREP - C0];
        if (q->len) fwrite(q->buf, 1, q->len, stdout);
        free(q->buf); q->buf = NULL;
        for (int k = 0; k < 8; k++) TOT[k] += q->d[k];
        if (PROG) printf("CDONE %llu points=%llu s4=%llu members=%llu\n", (unsigned long long)NEXTREP, (unsigned long long)TOT[0], (unsigned long long)TOT[4], (unsigned long long)TOT[7]);
        if (PROG || q->len) fflush(stdout);
        NEXTREP++;
    }
    pthread_mutex_unlock(&REPMU);
}
static void *worker(void *arg) {
    (void)arg;
    mpz_init(gx); mpz_init(gm); mpz_init(gt);
#ifndef SCALAR
    if (CFN == main_c) {
        size_t n = (3 * (size_t)NN + 64) * 8;
        PBUF = malloc(n); HBUF = malloc(n); PA = malloc(n); PBIN = malloc(n);
        if (!PBUF || !HBUF || !PA || !PBIN) { fprintf(stderr, "alloc failed\n"); exit(1); }
    }
#endif
    OBC = 1 << 16; OB = malloc(OBC);
    for (;;) {
        u64 c = __atomic_fetch_add(&NEXTC, 1, __ATOMIC_RELAXED);
        if (c >= C1) break;
        u64 snap[8], d[8]; memcpy(snap, cnt, sizeof snap);
        OBN = 0;
        CFN(c);
        for (int k = 0; k < 8; k++) d[k] = cnt[k] - snap[k];
        if (OBN) { report(c, OB, OBN, d); OBC = 1 << 16; OB = malloc(OBC); }
        else report(c, NULL, 0, d);
        if (!OB) { fprintf(stderr, "alloc failed\n"); exit(1); }
    }
    free(OB);
#ifndef SCALAR
    free(PBUF); free(HBUF); free(PA); free(PBIN);
#endif
    mpz_clear(gx); mpz_clear(gm); mpz_clear(gt);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: gsieve mode(main|direct|odd) G|E|P N c0 c1 [dump] [nogmp] [prog] [n0=M] [threads=T]\n"); return 1; }
    const char *mode = argv[1];
    DM = argv[2][0] == 'G' ? 0 : argv[2][0] == 'E' ? 1 : 2;
    BADEVEN = DM == 0 ? (F_ZERO | F_FIVE) : DM == 1 ? F_ZERO : 0;
    NN = atoll(argv[3]); u64 c0 = strtoull(argv[4], 0, 10), c1 = strtoull(argv[5], 0, 10);
    int NTH = getenv("GSIEVE_THREADS") ? atoi(getenv("GSIEVE_THREADS")) : 1;
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "dump")) DUMP = 1; else if (!strcmp(argv[i], "nogmp")) NOGMP = 1; else if (!strcmp(argv[i], "prog")) PROG = 1;
        else if (!strncmp(argv[i], "n0=", 3)) N0 = atoll(argv[i] + 3);
        else if (!strncmp(argv[i], "threads=", 8)) NTH = atoi(argv[i] + 8);
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 1; }
    }
    if (NTH < 1 || NTH > 1024) { fprintf(stderr, "bad thread count %d\n", NTH); return 1; }
    build_fl(FL1, 1); build_fl(FL2, 2); build_fl(FL3, 3); build_fl(FL4, 4); build_fl(FL5, 5);
    clock_t t0 = clock();
    struct timespec w0, w1; clock_gettime(CLOCK_MONOTONIC, &w0);
    if (!strcmp(mode, "main")) {
        build_bitmap(NTH); build_top(); build_t2(3 * (u64)NN + 2);
#ifdef S2_AVX2
        build_pack4();
#endif
        clock_gettime(CLOCK_MONOTONIC, &w1);
        fprintf(stderr, "tables: K=%d variant=%s good=%llu (density %.5f) build %.1fs (wall, %d threads)\n", KK, VARIANT, (unsigned long long)NGOOD,
                (double)NGOOD / P1, (double)(w1.tv_sec - w0.tv_sec) + 1e-9 * (double)(w1.tv_nsec - w0.tv_nsec), NTH);
        CFN = main_c;
    } else if (!strcmp(mode, "direct")) CFN = direct_c;
    else if (!strcmp(mode, "odd")) CFN = odd_c;
    else { fprintf(stderr, "bad mode\n"); return 1; }
    t0 = clock();
    clock_gettime(CLOCK_MONOTONIC, &w0);
    C0 = c0; C1 = c1; NEXTC = c0; NEXTREP = c0;
    REC = calloc(c1 > c0 ? c1 - c0 : 1, sizeof *REC);
    pthread_t *tid = malloc(NTH * sizeof *tid);
    if (!REC || !tid) { fprintf(stderr, "alloc failed\n"); return 1; }
    for (int i = 0; i < NTH; i++) if (pthread_create(&tid[i], NULL, worker, NULL)) { fprintf(stderr, "pthread_create failed\n"); return 1; }
    for (int i = 0; i < NTH; i++) pthread_join(tid[i], NULL);
    if (NEXTREP != c1 && c1 > c0) { fprintf(stderr, "reporter incomplete\n"); return 1; }
    if (PROG && c1 > c0) { printf("CDONE %llu points=%llu s4=%llu members=%llu\n", (unsigned long long)(c1 - 1), (unsigned long long)TOT[0], (unsigned long long)TOT[4], (unsigned long long)TOT[7]); fflush(stdout); }
    double el = (double)(clock() - t0) / CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC, &w1);
    double wall = (double)(w1.tv_sec - w0.tv_sec) + 1e-9 * (double)(w1.tv_nsec - w0.tv_nsec);
    // time = CPU seconds of all threads (core-seconds, as in single-threaded runs); wall = elapsed seconds
    printf("DONE mode=%s dm=%s N=%lld N0=%lld c=[%llu,%llu) K=%d NL=%d variant=%s points=%llu s1=%llu s2=%llu s3=%llu s4=%llu gmp=%llu full=%llu members=%llu time=%.2f threads=%d wall=%.2f\n",
           mode, argv[2], (long long)NN, (long long)N0, (unsigned long long)c0, (unsigned long long)c1, KK, NL, VARIANT,
           (unsigned long long)TOT[0], (unsigned long long)TOT[1], (unsigned long long)TOT[2], (unsigned long long)TOT[3],
           (unsigned long long)TOT[4], (unsigned long long)TOT[5], (unsigned long long)TOT[6], (unsigned long long)TOT[7], el, NTH, wall);
    return 0;
}
