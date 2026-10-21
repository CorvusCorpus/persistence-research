#define _GNU_SOURCE
// z2sieve: decimal zeros of powers of two (OEIS A007377: 2^n has no digit 0 only for the 36 known n <= 86).
// For every n in [N_LO, N_HI) it finds a digit 0 of 2^n, and it reports every n for which it finds none.
// All rejections are exact integer arithmetic, except the leading-digit table, which is built with explicit margins.
//
// Write n = r + j P with P = 4 5^(K-1) (the order of 2 modulo 5^K) and 0 <= r < P (docs/ZEROLESS_POW2.md).
//   stage 1: the last K digits of 2^n are w(r) = 2^(r+P) mod 10^K (n >= K); w is doubled along r.
//   stage 2: digits K+1..2K are U(j) = a0 + 2^K t(j), t(j) = (z_hi + j d) mod 5^K (n >= 2K): along j, U steps by the
//            constant 2^K d modulo 10^K.
//   stage 3: the first TD digits, from frac(n log10 2) in 128-bit fixed point and a table of 2^TB bins; a bin rejects
//            only if every TD-digit prefix that its interval (widened by margins) can hold contains a 0.
//   stage 4: digits 2K+1..9 NL of (2^r mod 10^(9 NL)) (2^(jP) mod 10^(9 NL)), base-10^9 limbs, one column at a time.
//   stage 5: the last 1000 digits (GMP) and the first 40 (MPFR, rigorous enclosure); 2^n in full if n <= FULLMAX.
// Work: chunks of JCH consecutive j (JCH P consecutive n), each split into blocks of RB residues r over the threads.
// A chunk's CDONE line is written only after its integrity checks pass (tables re-hashed, every block's running state
// re-derived by modular exponentiation), so a log can be resumed from its last CDONE line.
// usage: z2sieve N_LO N_HI [threads=T] [dump] [checkall]      (N_LO >= NMIN; output does not depend on T)
//        z2sieve kerneltest SAMPLES SEED                        (every stage against GMP/MPFR at random n)
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef Z2_CUDA
#include "z2sieve_cuda.h"           // stages 2-4 on the GPU (z2sieve_cuda.cu); stage 1 is a host bitmap, stage 5 the host
#define VARIANT "cuda"
#else
#define VARIANT "cpu"
#endif

typedef uint64_t u64; typedef uint32_t u32; typedef unsigned __int128 u128;

#ifndef KK
#define KK 13                      // stage 1 and 2: K digits each
#endif
#ifndef NL
#define NL 18                      // stage 4: the last 9 NL digits
#endif
#ifndef TB
#define TB 24                      // stage 3: 2^TB bins ...
#endif
#ifndef TD
#define TD 7                       // ... testing the first TD digits
#endif
#ifndef NMIN
#define NMIN 1000                  // smallest n accepted: 2^n must be longer than the 9 NL-digit window
#endif
#ifndef JCH
#define JCH 4096                   // j values per chunk (a chunk is JCH P = 4 10^12 consecutive n)
#endif
#define RB (1ULL << 20)            // residues r per block
#define E9 1000000000ULL
#define LDIG 1000                  // stage 5: last LDIG digits ...
#define LTOP 40                    // ... and first LTOP digits
#define FULLMAX 10000000ULL        // 2^n is computed in full when n <= FULLMAX
#if KK < 4 || KK > 13
#error "KK must be in 4..13 (5^(2K) must fit in 63 bits)"
#endif
#if NL > 18 || 9 * NL <= 2 * KK
#error "need 2K < 9 NL and NL <= 18 (column sums of the limb engine stay below 2^64)"
#endif
#if TD < 2 || TD > 9 || TB < 12 || TB > 32
#error "need TD in 2..9 and TB in 12..32"
#endif

static const u64 P10T[19] = { 1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL, 100000000ULL,
    1000000000ULL, 10000000000ULL, 100000000000ULL, 1000000000000ULL, 10000000000000ULL, 100000000000000ULL,
    1000000000000000ULL, 10000000000000000ULL, 100000000000000000ULL, 1000000000000000000ULL };
static const u64 F5T[27] = { 1ULL, 5ULL, 25ULL, 125ULL, 625ULL, 3125ULL, 15625ULL, 78125ULL, 390625ULL, 1953125ULL,
    9765625ULL, 48828125ULL, 244140625ULL, 1220703125ULL, 6103515625ULL, 30517578125ULL, 152587890625ULL,
    762939453125ULL, 3814697265625ULL, 19073486328125ULL, 95367431640625ULL, 476837158203125ULL, 2384185791015625ULL,
    11920928955078125ULL, 59604644775390625ULL, 298023223876953125ULL, 1490116119384765625ULL };
#define TENK P10T[KK]
#define M5 F5T[KK]                 // 5^K
#define M2 F5T[2 * KK]             // 5^(2K)
#define PP (4 * F5T[KK - 1])       // P
#define TWOK (1ULL << KK)
// a K-digit window is kept as three parts a|b|c of WA, WB, WC digits, so that stepping it needs no division
#define WC (KK < 5 ? KK : 5)
#define WB (KK - WC < 4 ? KK - WC : 4)
#define WA (KK - WC - WB)
#define BC P10T[WC]
#define BB P10T[WB]
#define BA P10T[WA]
#define LQ (2 * KK / 9)            // the limb that holds digit 2K+1 (digits counted from 1 at the right) ...
#define LO (2 * KK % 9)            // ... has LO digits below it

static u64 C2K, UU, LLO, LHI;      // 2^(-2K) mod 5^(2K); u with 2^P = 1 + 5^K u (mod 5^(2K)); floor(log10(2) 2^128)
static u64 *ZT[6];                 // ZT[w]: bit v set iff v written with w digits (zero padded) has no 0
static u64 *TOP;                   // stage 3: bit k set iff bin k rejects
static u32 (*Q)[NL];               // Q[i] = 2^((CH_JA + i) P) mod 10^(9 NL) for the current chunk
static int NTH = 1, DUMP = 0, CHKALL = 0;
static u64 CH_JA, CH_JB, CH_N0, CH_N1, NBLK, NEXTB __attribute__((unused));
typedef struct { u64 cnt[8]; char *out; size_t on, oc; int bad; } bres_t;    // one block's counters and output
static bres_t *BRES;
static __thread bres_t *CUR;
static __thread u32 *JBUF __attribute__((unused));
static __thread mpz_t gx, gm, gt; static __thread mpfr_t fa, fb, fc, fd;
// counters: 1 = passed stage 1, 2 = stage 2, 3 = stage 3, 4 = stage 4, 5 = no 0 in the stage-5 windows (or n small
// enough to compute 2^n in full), 6 = CAND (no 0 found, 2^n too large to compute), 7 = MEMBER (2^n has no 0)

static void die(const char *msg) { fprintf(stderr, "%s\n", msg); exit(1); }
static void fail(const char *msg) {        // integrity failure: stop before anything of this chunk is written
    fprintf(stderr, "INTEGRITY CHECK FAILED: %s (memory fault?); nothing of this chunk was written\n", msg);
    exit(3);
}
static void emit(const char *fmt, ...) {
    bres_t *b = CUR;
    for (;;) {
        va_list ap; va_start(ap, fmt);
        int n = vsnprintf(b->out ? b->out + b->on : NULL, b->out ? b->oc - b->on : 0, fmt, ap);
        va_end(ap);
        if (n < 0) die("emit failed");
        if (b->out && (size_t)n < b->oc - b->on) { b->on += (size_t)n; return; }
        b->oc = 2 * b->oc + (size_t)n + 256;
        if (!(b->out = realloc(b->out, b->oc))) die("alloc failed");
    }
}
static u64 fnv(u64 h, const void *p, size_t n) {
    const uint8_t *q = p;
    for (size_t i = 0; i < n; i++) { h ^= q[i]; h *= 1099511628211ULL; }
    return h;
}

// ---------------- digit windows ----------------
#define ZBIT(T, v) ((int)(((T)[(v) >> 6] >> ((v) & 63)) & 1))
static void build_zt(void) {
    for (int w = 0; w <= 5; w++) {
        u64 n = P10T[w];
        ZT[w] = calloc(n / 64 + 1, 8);
        for (u64 v = 0; v < n; v++) {
            u64 x = v; int ok = 1;
            for (int i = 0; i < w; i++) { if (x % 10 == 0) ok = 0; x /= 10; }
            if (ok) ZT[w][v >> 6] |= 1ULL << (v & 63);
        }
    }
}
static inline int zfw(u64 v, int w) {           // v written with w <= 9 digits has no 0
    return w <= 5 ? ZBIT(ZT[w], v) : ZBIT(ZT[4], v % 10000) & zfw(v / 10000, w - 4);
}
static inline int zf9(u64 v) { return ZBIT(ZT[4], v % 10000) & ZBIT(ZT[5], v / 10000); }
typedef struct { u64 a, b, c; } kd;            // a K-digit window as parts a|b|c
static inline kd kd_split(u64 x) { kd k; k.c = x % BC; x /= BC; k.b = x % BB; k.a = x / BB; return k; }
static inline u64 kd_join(kd k) { return (k.a * BB + k.b) * BC + k.c; }
static inline int kd_ok(kd k) { return ZBIT(ZT[WA], k.a) & ZBIT(ZT[WB], k.b) & ZBIT(ZT[WC], k.c); }
static inline void kd_add(kd *k, kd d) {        // k = (k + d) mod 10^K
    u64 c = k->c + d.c, cy = c >= BC; c -= cy ? BC : 0;
    u64 b = k->b + d.b + cy; cy = b >= BB; b -= cy ? BB : 0;
    u64 a = k->a + d.a + cy; a -= a >= BA ? BA : 0;
    k->a = a; k->b = b; k->c = c;
}

// ---------------- arithmetic ----------------
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1 % m; b %= m;
    for (; e; e >>= 1) { if (e & 1) r = mulmod(r, b, m); b = mulmod(b, b, m); }
    return r;
}
static void init_consts(void) {
    // C2K = (2^(2K))^(-1) mod 5^(2K): 2^(-1) = (5^(2K) + 1) / 2
    C2K = powmod((M2 + 1) / 2, 2 * KK, M2);
    if (mulmod(C2K, powmod(2, 2 * KK, M2), M2) != 1) die("C2K");
    u64 g = powmod(2, PP, M2);
    if ((g - 1) % M5) die("2^P is not 1 mod 5^K");
    UU = (g - 1) / M5 % M5;
    if (UU % 5 == 0) die("2 is not a primitive root mod 5^(K+1)");
    mpfr_t x; mpz_t z; mpfr_init2(x, 320); mpz_init(z);
    mpfr_set_ui(x, 2, MPFR_RNDN); mpfr_log10(x, x, MPFR_RNDD); mpfr_mul_2ui(x, x, 128, MPFR_RNDD); mpfr_get_z(z, x, MPFR_RNDD);
    LLO = mpz_getlimbn(z, 0); LHI = mpz_getlimbn(z, 1);
    if (mpz_sizeinbase(z, 2) > 128 || sizeof(mp_limb_t) != 8) die("log10(2) constant");
    mpfr_clear(x); mpz_clear(z);
}
static u64 const_hash(void) { u64 v[4] = { C2K, UU, LLO, LHI }; return fnv(1469598103934665603ULL, v, sizeof v); }
static inline u64 frac_hi(u64 n) {              // top 64 bits of frac(n log10 2), error < 2^-63
    u128 p = (u128)n * LLO;
    return (u64)(p >> 64) + n * LHI;
}
static void to_limbs(mpz_t x, u32 *L) {         // x mod 10^(9 NL) as NL base-10^9 limbs (x is destroyed)
    for (int i = 0; i < NL; i++) L[i] = (u32)mpz_fdiv_q_ui(x, x, E9);
}
static void pow2_limbs(u64 e, u32 *L) {         // 2^e mod 10^(9 NL)
    mpz_set_ui(gm, 10); mpz_pow_ui(gm, gm, 9 * NL); mpz_set_ui(gt, 2); mpz_powm_ui(gx, gt, e, gm); to_limbs(gx, L);
}
static inline void limbs_mul_small(u32 *L, u64 k) {     // L = L k mod 10^(9 NL), k < 2^30
    u64 cy = 0;
    for (int i = 0; i < NL; i++) { u64 v = (u64)L[i] * k + cy; cy = v / E9; L[i] = (u32)(v - cy * E9); }
}

// ---------------- stage 3: leading-digit table ----------------
typedef struct { u64 w0, w1; } trange_t;
static void *top_part(void *arg) {
    trange_t *t = arg;
    long double d = ldexpl(1.0L, -40);
    u64 lo10 = P10T[TD - 1];
    for (u64 k = t->w0 * 64; k < t->w1 * 64 && k < (1ULL << TB); k++) {
        if (k == 0 || k + 1 == (1ULL << TB)) continue;          // bins next to the wrap-around never reject
        long double flo = ldexpl((long double)k, -TB) - d, fhi = ldexpl((long double)(k + 1), -TB) + d;
        long double vlo = exp10l(flo + (TD - 1)) * (1.0L - 1e-15L), vhi = exp10l(fhi + (TD - 1)) * (1.0L + 1e-15L);
        u64 mlo = (u64)floorl(vlo), mhi = (u64)floorl(vhi);
        if (mlo < lo10 || mhi >= 10 * lo10) continue;
        int bad = 1;
        for (u64 m = mlo; m <= mhi && bad; m++) if (zfw(m, TD)) bad = 0;
        if (bad) TOP[k >> 6] |= 1ULL << (k & 63);
    }
    return NULL;
}
static void build_top(void) {
    u64 nw = (1ULL << TB) / 64 + 1;
    TOP = calloc(nw, 8);
    pthread_t tid[64]; trange_t tr[64]; int nt = NTH < 64 ? NTH : 64;
    for (int i = 0; i < nt; i++) { tr[i].w0 = nw * i / nt; tr[i].w1 = nw * (i + 1) / nt; pthread_create(&tid[i], NULL, top_part, &tr[i]); }
    for (int i = 0; i < nt; i++) pthread_join(tid[i], NULL);
}
static u64 tables_hash(void) {
    u64 h = fnv(1469598103934665603ULL, TOP, ((1ULL << TB) / 64 + 1) * 8);
    for (int w = 0; w <= 5; w++) h = fnv(h, ZT[w], (P10T[w] / 64 + 1) * 8);
    return fnv(h, &(u64){ const_hash() }, 8);
}

// ---------------- stage 4: limb engine ----------------
// digits 2K+1..9 NL of (A B) mod 10^(9 NL), with A = 2^r and B = 2^(jP) modulo 10^(9 NL): returns 1 if none is 0.
// With chk, the limbs below digit 2K+1 must equal those of w(r) + 10^K U(j), already known from stages 1 and 2.
static inline int limb_engine(const u32 *A, const u32 *B, u64 w, u64 U, int chk) {
    u64 e18 = 0, ehi = 0;
    if (chk) { e18 = w + TENK * (U % P10T[18 - KK]); ehi = U / P10T[18 - KK]; }
    u64 cy = 0;
    for (int k = 0; k < NL; k++) {
        u64 acc = cy;
        for (int u = 0; u <= k; u++) acc += (u64)A[u] * B[k - u];
        cy = acc / E9; u64 limb = acc - cy * E9;
        if (k < LQ) {
            if (chk && limb != (k == 0 ? e18 % E9 : k == 1 ? e18 / E9 : ehi)) fail("limb engine disagrees with stages 1-2");
        } else if (k == LQ) {
            if (chk && limb % P10T[LO] != (k == 0 ? e18 % E9 : k == 1 ? e18 / E9 : ehi)) fail("limb engine disagrees with stages 1-2");
            if (!zfw(limb / P10T[LO], 9 - LO)) return 0;
        } else if (!zf9(limb)) return 0;
    }
    return 1;
}

// ---------------- stage 5 ----------------
static int lead_digits(u64 n, char *out) {      // the leading digits of 2^n that are certain (at most LTOP)
    out[0] = 0;
    mpfr_set_ui(fa, 2, MPFR_RNDN); mpfr_log10(fa, fa, MPFR_RNDD); mpfr_mul_ui(fa, fa, n, MPFR_RNDD);   // fa <= n log10 2
    mpfr_set_ui(fb, 2, MPFR_RNDN); mpfr_log10(fb, fb, MPFR_RNDU); mpfr_mul_ui(fb, fb, n, MPFR_RNDU);   // fb >= n log10 2
    mpfr_floor(fc, fa); mpfr_floor(fd, fb);
    if (!mpfr_equal_p(fc, fd)) return 0;
    mpfr_sub(fa, fa, fc, MPFR_RNDD); mpfr_sub(fb, fb, fc, MPFR_RNDU);
    mpfr_add_ui(fa, fa, LTOP - 1, MPFR_RNDD); mpfr_add_ui(fb, fb, LTOP - 1, MPFR_RNDU);
    mpfr_exp10(fa, fa, MPFR_RNDD); mpfr_exp10(fb, fb, MPFR_RNDU);
    mpfr_get_z(gx, fa, MPFR_RNDD); mpfr_get_z(gt, fb, MPFR_RNDD);
    char sa[LTOP + 8], sb[LTOP + 8];
    if (mpz_sizeinbase(gx, 10) > LTOP + 1 || mpz_sizeinbase(gt, 10) > LTOP + 1) return 0;
    mpz_get_str(sa, 10, gx); mpz_get_str(sb, 10, gt);
    if (strlen(sa) != LTOP || strlen(sb) != LTOP) return 0;
    int k = 0;
    while (k < LTOP && sa[k] == sb[k]) k++;
    memcpy(out, sa, k); out[k] = 0;
    return k;
}
static void full_check(u64 n) {
    mpz_ui_pow_ui(gx, 2, n);
    char *s = mpz_get_str(NULL, 10, gx);
    int zero = strchr(s, '0') != NULL;
    free(s);
    if (zero) { if (DUMP) emit("R5 %llu\n", (unsigned long long)n); return; }
    CUR->cnt[7]++; emit("MEMBER %llu\n", (unsigned long long)n);
}
static void final_check(u64 n) {
    if (n <= 4000) { CUR->cnt[5]++; full_check(n); return; }     // at most 1205 digits: compute 2^n in full
    mpz_set_ui(gm, 10); mpz_pow_ui(gm, gm, LDIG); mpz_set_ui(gt, 2); mpz_powm_ui(gx, gt, n, gm);
    char s[LDIG + 8]; mpz_get_str(s, 10, gx);
    if (strlen(s) < LDIG || strchr(s, '0')) { if (DUMP) emit("R5 %llu\n", (unsigned long long)n); return; }   // 2^n > 10^1205
    char t[LTOP + 8]; lead_digits(n, t);
    if (strchr(t, '0')) { if (DUMP) emit("R5 %llu\n", (unsigned long long)n); return; }
    CUR->cnt[5]++;
    if (n <= FULLMAX) { full_check(n); return; }
    CUR->cnt[6]++; emit("CAND %llu\n", (unsigned long long)n);
}

// ---------------- the sieve ----------------
#ifndef Z2_CUDA                                 // the CPU variant's stages 2-4
static void do_r(u64 r, u64 w, u64 s, const u32 *R, u64 j0, u64 j1) {
    bres_t *B = CUR;
    u64 z0 = mulmod(s, C2K, M2), zlo = z0 % M5, zhi = z0 / M5;
    u64 a0 = TWOK * zlo / M5;                            // 2^K zlo < 10^K and 2^(2K) zlo < 2^K 10^K < 2^64
    if (TWOK * TWOK * zlo % TENK != w) fail("stage-2 setup disagrees with stage 1");
    u64 d = mulmod(zlo, UU, M5), t0 = (zhi + mulmod(j0 % M5, d, M5)) % M5;
    kd u = kd_split(a0 + TWOK * t0), dd = kd_split(TWOK * d);
    u32 nc = 0;
    for (u64 j = j0; j < j1; j++) { JBUF[nc] = (u32)(j - j0); nc += kd_ok(u); kd_add(&u, dd); }
    B->cnt[1] += j1 - j0; B->cnt[2] += nc;
    for (u32 i = 0; i < nc; i++) {
        u64 off = JBUF[i], j = j0 + off, n = r + j * PP;
        if (DUMP) emit("S2 %llu\n", (unsigned long long)n);
        u64 bin = frac_hi(n) >> (64 - TB);
        if (ZBIT(TOP, bin)) { if (DUMP) emit("R3 %llu\n", (unsigned long long)n); continue; }
        B->cnt[3]++;
        int chk = CHKALL || DUMP || !(n & 63);
        u64 U = chk ? a0 + TWOK * ((t0 + off * d) % M5) : 0;
        if (!limb_engine(R, Q[j - CH_JA], w, U, chk)) { if (DUMP) emit("R4 %llu\n", (unsigned long long)n); continue; }
        B->cnt[4]++;
        if (DUMP) emit("S4 %llu\n", (unsigned long long)n);
        final_check(n);
    }
}
static void do_block(u64 b) {
    u64 ra = b * RB, rb = ra + RB < PP ? ra + RB : PP;
    int any = 0;
    for (u64 j = CH_JA; j < CH_JB && !any; j++) if (ra + j * PP < CH_N1 && rb + j * PP > CH_N0) any = 1;
    if (!any) return;
    kd W = kd_split(powmod(2, ra + PP, TENK));          // w(r) = 2^(r+P) mod 10^K
    u64 rl = ra, s = powmod(2, ra, M2);                  // s = 2^rl mod 5^(2K) and R = 2^rl mod 10^(9 NL), advanced lazily
    u32 R[NL]; pow2_limbs(ra, R);
    for (u64 r = ra; r < rb; r++, kd_add(&W, W)) {
        if (!kd_ok(W)) continue;
        u64 jlo = CH_N0 > r ? (CH_N0 - r + PP - 1) / PP : 0, jhi = CH_N1 > r ? (CH_N1 - r + PP - 1) / PP : 0;
        if (jlo < CH_JA) jlo = CH_JA;
        if (jhi > CH_JB) jhi = CH_JB;
        if (jlo >= jhi) continue;
        while (rl < r) { u64 g = r - rl < 29 ? r - rl : 29; s = (u64)(((u128)s << g) % M2); limbs_mul_small(R, 1ULL << g); rl += g; }
        do_r(r, kd_join(W), s, R, jlo, jhi);
    }
    // the running state must agree with modular exponentiation at the end of the block
    u32 R2[NL]; pow2_limbs(rl, R2);
    if (kd_join(W) != powmod(2, rb + PP, TENK) || s != powmod(2, rl, M2) || memcmp(R, R2, sizeof R)) CUR->bad = 1;
}
static void *worker(void *arg) {
    (void)arg;
    mpz_inits(gx, gm, gt, NULL); mpfr_inits2(256, fa, fb, fc, fd, (mpfr_ptr)0);
    JBUF = malloc(JCH * sizeof *JBUF);
    for (;;) {
        u64 b = __atomic_fetch_add(&NEXTB, 1, __ATOMIC_RELAXED);
        if (b >= NBLK) break;
        CUR = &BRES[b]; do_block(b);
    }
    free(JBUF); mpz_clears(gx, gm, gt, NULL); mpfr_clears(fa, fb, fc, fd, (mpfr_ptr)0); mpfr_free_cache();
    return NULL;
}
#endif
static void build_q(void) {                     // Q for the chunk [CH_JA, CH_JB), by modular exponentiation
    mpz_t m, x, two, e; mpz_inits(m, x, two, e, NULL);
    mpz_set_ui(m, 10); mpz_pow_ui(m, m, 9 * NL); mpz_set_ui(two, 2);
    for (u64 i = 0; i < JCH; i++) { mpz_set_ui(e, CH_JA + i); mpz_mul_ui(e, e, PP); mpz_powm(x, two, e, m); to_limbs(x, Q[i]); }
    mpz_clears(m, x, two, e, NULL);
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static u64 TH0;                                 // hash of the tables, taken after they were built

// ---------------- CUDA variant: host tables, one chunk on the GPU, stage 5 on host threads ----------------
#ifdef Z2_CUDA
static u64 *B1, NW1, S1TOT, NRS, *RSW, *RSS, PW13[1024], PSV[1024]; static u32 *RK1, *RSR, PWL[1024 * NL];
typedef struct { u64 a, b; } span_t;
static void run_spans(void *(*fn)(void *), u64 n) {     // fn over [0, n) split into NTH spans
    pthread_t tid[64]; span_t sp[64]; int nt = NTH < 64 ? NTH : 64;
    for (int t = 0; t < nt; t++) { sp[t].a = n * t / nt; sp[t].b = n * (t + 1) / nt; if (pthread_create(&tid[t], NULL, fn, &sp[t])) die("pthread_create failed"); }
    for (int t = 0; t < nt; t++) pthread_join(tid[t], NULL);
}
static void *b1_part(void *arg) {               // stage 1 for the words [a, b) of the bitmap
    span_t *s = arg;
    kd W = kd_split(powmod(2, s->a * 64 + PP, TENK));
    for (u64 w = s->a; w < s->b; w++) {
        u64 bits = 0;
        for (int k = 0; k < 64; k++, kd_add(&W, W)) if (w * 64 + k < PP && kd_ok(W)) bits |= 1ULL << k;
        B1[w] = bits;
    }
    return NULL;
}
static void *rs_part(void *arg) {               // the state at r = 1024 k for k in [a, b)
    span_t *s = arg;
    mpz_inits(gx, gm, gt, NULL);
    for (u64 k = s->a; k < s->b; k++) { RSW[k] = powmod(2, 1024 * k + PP, TENK); RSS[k] = powmod(2, 1024 * k, M2); pow2_limbs(1024 * k, RSR + k * NL); }
    mpz_clears(gx, gm, gt, NULL);
    return NULL;
}
static void gpu_setup(void) {
    double t0 = now();
    NW1 = (PP + 63) / 64; B1 = calloc(NW1 + 1, 8); RK1 = malloc((NW1 + 1) * 4);
    run_spans(b1_part, NW1);
    u64 run = 0;
    for (u64 w = 0; w < NW1; w++) { RK1[w] = (u32)run; run += __builtin_popcountll(B1[w]); }
    RK1[NW1] = (u32)run; S1TOT = run;
    if (run >> 32) die("too many stage-1 survivors for 32-bit ranks");
    NRS = (PP + 1023) / 1024; RSW = malloc(NRS * 8); RSS = malloc(NRS * 8); RSR = malloc(NRS * NL * 4);
    run_spans(rs_part, NRS);
    mpz_inits(gx, gm, gt, NULL); mpfr_inits2(256, fa, fb, fc, fd, (mpfr_ptr)0);    // stage 5 in the main thread
    for (int g = 0; g < 1024; g++) { PW13[g] = powmod(2, g, TENK); PSV[g] = powmod(2, g, M2); pow2_limbs(g, PWL + g * NL); }
    static u64 zt[1740]; u64 o = 0;                   // static: kept as the host copy of the device table
    for (int w = 0; w <= 5; w++) for (u64 i = 0; i < P10T[w] / 64 + 1; i++) zt[o++] = ZT[w][i];
    if (o != 1740) die("digit table size");
    z2c_t c = { C2K, UU, LLO, LHI };
    fprintf(stderr, "host tables: stage-1 bitmap %llu survivors of %llu (%.5f), run starts %llu, built in %.1fs\n",
            (unsigned long long)S1TOT, (unsigned long long)PP, (double)S1TOT / PP, (unsigned long long)NRS, now() - t0);
    z2g_init(&c, zt, TOP, (1ULL << TB) / 64 + 1, B1, RK1, NW1, RSW, RSS, RSR, NRS, PW13, PSV, PWL);
}
static u64 rank1(u64 r) {                       // stage-1 survivors below r
    u64 w = r >> 6;
    return RK1[w] + ((r & 63) ? (u64)__builtin_popcountll(B1[w] & ((1ULL << (r & 63)) - 1)) : 0);
}
static int cmpu64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
typedef struct { const u64 *v; u64 n; bres_t b; } s5_t;
static void *s5_worker(void *arg) {
    s5_t *s = arg; CUR = &s->b;
    mpz_inits(gx, gm, gt, NULL); mpfr_inits2(256, fa, fb, fc, fd, (mpfr_ptr)0);
    for (u64 i = 0; i < s->n; i++) final_check(s->v[i]);
    mpz_clears(gx, gm, gt, NULL); mpfr_clears(fa, fb, fc, fd, (mpfr_ptr)0); mpfr_free_cache();
    return NULL;
}
static void gpu_chunk(u64 *cc, u64 qh) {
    u64 iv[4]; int niv;                         // survivor indices whose r has some n of [CH_N0, CH_N1)
    if (CH_N1 - CH_N0 >= PP) { iv[0] = 0; iv[1] = S1TOT; niv = 1; }
    else {
        u64 ja = CH_N0 / PP, jb = (CH_N1 - 1) / PP;
        if (ja == jb) { iv[0] = rank1(CH_N0 - ja * PP); iv[1] = rank1(CH_N1 - ja * PP); niv = 1; }
        else { iv[0] = rank1(CH_N0 - ja * PP); iv[1] = S1TOT; iv[2] = 0; iv[3] = rank1(CH_N1 - jb * PP); niv = 2; }
    }
    z2r_t res;
    z2g_chunk(CH_JA, CH_JB, CH_N0, CH_N1, &Q[0][0], iv, niv, DUMP, &res);
    if (res.flags) { fprintf(stderr, "gpu flags %u\n", res.flags); fail("the GPU's stages disagree with the host tables, or a buffer overflowed"); }
    for (int k = 0; k < 4; k++) cc[k + 1] = res.ctr[k];
    if (res.ns4 != cc[4]) fail("the GPU's stage-4 list disagrees with its counter");
    bres_t B = { 0 };
    if (DUMP) {                                 // all lines in increasing n (stage 5 in this thread)
        CUR = &B; qsort(res.dmp, res.ndmp, 8, cmpu64);
        u64 n4 = 0;
        for (u64 i = 0; i < res.ndmp; i++) {
            u64 n = res.dmp[i] >> 2; int cls = (int)(res.dmp[i] & 3);
            emit("S2 %llu\n", (unsigned long long)n);
            if (cls == 1) emit("R3 %llu\n", (unsigned long long)n);
            else if (cls == 2) emit("R4 %llu\n", (unsigned long long)n);
            else { emit("S4 %llu\n", (unsigned long long)n); final_check(n); n4++; }
        }
        if (res.ndmp != cc[2] || n4 != cc[4]) fail("the GPU dump disagrees with its counters");
    } else {                                    // stage 5 on NTH host threads, output in increasing n
        qsort(res.s4, res.ns4, 8, cmpu64);
        s5_t sp[64]; pthread_t tid[64]; int nt = NTH < 64 ? NTH : 64;
        for (int t = 0; t < nt; t++) { memset(&sp[t], 0, sizeof sp[t]); sp[t].v = res.s4 + res.ns4 * t / nt; sp[t].n = res.ns4 * (t + 1) / nt - res.ns4 * t / nt; if (pthread_create(&tid[t], NULL, s5_worker, &sp[t])) die("pthread_create failed"); }
        for (int t = 0; t < nt; t++) {
            pthread_join(tid[t], NULL);
            for (int k = 5; k < 8; k++) B.cnt[k] += sp[t].b.cnt[k];
            if (sp[t].b.on) { CUR = &B; emit("%.*s", (int)sp[t].b.on, sp[t].b.out); }
            free(sp[t].b.out);
        }
    }
    for (int k = 5; k < 8; k++) cc[k] = B.cnt[k];
    if (tables_hash() != TH0) fail("the tables changed");
    if (fnv(1469598103934665603ULL, Q, JCH * sizeof *Q) != qh) fail("the Q table changed");
    if (z2g_verify()) fail("a device table no longer equals its host copy");
    if (B.on) fwrite(B.out, 1, B.on, stdout);
    free(B.out);
}
#endif

// ---------------- kernel test ----------------
static int kerneltest(u64 samples, u64 seed) {
    mpz_inits(gx, gm, gt, NULL); mpfr_inits2(256, fa, fb, fc, fd, (mpfr_ptr)0);
    bres_t dummy = { 0 }; CUR = &dummy;
    mpz_t x, m, two; mpz_inits(x, m, two, NULL); mpz_set_ui(m, 10); mpz_pow_ui(m, m, 9 * NL); mpz_set_ui(two, 2);
    u64 bad[6] = { 0 }, rej3 = 0, rej4 = 0, st = seed * 0x9E3779B97F4A7C15ULL + 1;
    for (u64 i = 0; i < samples; i++) {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        u64 n = NMIN + st % ((1ULL << (20 + i % 36)) - NMIN), r = n % PP, j = n / PP;
        mpz_powm_ui(x, two, n, m);
        u32 X[NL]; to_limbs(x, X);                       // 2^n mod 10^(9 NL), exactly
        u128 low = 0;                                     // 2^n mod 10^(2K) from the limbs (2K <= 9 min(NL, 3))
        for (int k = (NL < 3 ? NL : 3) - 1; k >= 0; k--) low = low * E9 + X[k];
        u64 w = powmod(2, r + PP, TENK), lowK = (u64)(low % TENK), U = (u64)((low / TENK) % TENK);
        bad[0] += w != lowK;                              // stage 1
        u64 s = powmod(2, r, M2), z0 = mulmod(s, C2K, M2), zlo = z0 % M5, zhi = z0 / M5;
        u64 a0 = (u64)(((u128)TWOK * zlo) / M5), d = mulmod(zlo, UU, M5), t = (zhi + mulmod(j % M5, d, M5)) % M5;
        bad[1] += a0 + TWOK * t != U;                     // stage 2
        u32 A[NL], B[NL], C[NL]; pow2_limbs(r, A);
        mpz_set_ui(gx, j); mpz_mul_ui(gx, gx, PP); mpz_powm(gx, two, gx, m); to_limbs(gx, B);
        u64 cy = 0;                                       // the limb product, all limbs
        for (int k = 0; k < NL; k++) { u64 acc = cy; for (int u = 0; u <= k; u++) acc += (u64)A[u] * B[k - u]; cy = acc / E9; C[k] = (u32)(acc - cy * E9); }
        bad[2] += memcmp(C, X, sizeof C) != 0;
        int exact4 = 1;                                   // digits 2K+1..9 NL of 2^n have no 0
        for (int k = LQ; k < NL; k++) if (!(k == LQ ? zfw(X[k] / P10T[LO], 9 - LO) : zf9(X[k]))) exact4 = 0;
        int eng = limb_engine(A, B, w, U, 1);
        bad[3] += eng != exact4; rej4 += !eng;
        mpfr_set_ui(fa, 2, MPFR_RNDN); mpfr_log10(fa, fa, MPFR_RNDN); mpfr_mul_ui(fa, fa, n, MPFR_RNDN);
        mpfr_frac(fa, fa, MPFR_RNDN); mpfr_mul_2ui(fa, fa, 64, MPFR_RNDN);
        u64 fh = frac_hi(n), fx = mpfr_get_uj(fa, MPFR_RNDD), dl = fh - fx;
        bad[4] += dl > 2 && 0 - dl > 2;                    // stage-3 fixed point within 2^-63 of the exact fraction (mod 1)
        if (ZBIT(TOP, fh >> (64 - TB))) {                  // a stage-3 rejection must be right
            rej3++; char t2[LTOP + 8]; int k = lead_digits(n, t2);
            bad[5] += !(k >= TD && memchr(t2, '0', TD));
        }
    }
    printf("kerneltest K=%d NL=%d TB=%d TD=%d samples=%llu seed=%llu: mismatches stage1 %llu stage2 %llu limbs %llu engine %llu frac %llu top %llu"
           " (stage-3 rejections %llu, stage-4 rejections %llu)\n", KK, NL, TB, TD, (unsigned long long)samples, (unsigned long long)seed,
           (unsigned long long)bad[0], (unsigned long long)bad[1], (unsigned long long)bad[2], (unsigned long long)bad[3],
           (unsigned long long)bad[4], (unsigned long long)bad[5], (unsigned long long)rej3, (unsigned long long)rej4);
    return bad[0] + bad[1] + bad[2] + bad[3] + bad[4] + bad[5] ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 3) die("usage: z2sieve N_LO N_HI [threads=T] [dump] [checkall] | z2sieve kerneltest SAMPLES SEED");
    for (int i = 3; i < argc; i++) {
        if (!strncmp(argv[i], "threads=", 8)) NTH = atoi(argv[i] + 8);
        else if (!strcmp(argv[i], "dump")) DUMP = 1;
        else if (!strcmp(argv[i], "checkall")) CHKALL = 1;
        else if (strcmp(argv[1], "kerneltest")) die("unknown option");
    }
    if (NTH < 1) NTH = 1;
    init_consts(); build_zt();
    double t0 = now();
    build_top();
    u64 nrej = 0; for (u64 w = 0; w < (1ULL << TB) / 64 + 1; w++) nrej += __builtin_popcountll(TOP[w]);
    fprintf(stderr, "z2sieve K=%d NL=%d TB=%d TD=%d NMIN=%d JCH=%d variant=%s: P=%llu, log10(2) 2^128 = %016llx%016llx; table: %llu of %llu bins reject (%.4f), built in %.1fs\n",
            KK, NL, TB, TD, NMIN, JCH, VARIANT, (unsigned long long)PP, (unsigned long long)LHI, (unsigned long long)LLO,
            (unsigned long long)nrej, 1ULL << TB, (double)nrej / (double)(1ULL << TB), now() - t0);
    if (!strcmp(argv[1], "kerneltest")) return kerneltest(strtoull(argv[2], NULL, 10), strtoull(argv[3 < argc ? 3 : 2], NULL, 10));
    u64 nlo = strtoull(argv[1], NULL, 10), nhi = strtoull(argv[2], NULL, 10);
    if (nlo < NMIN || nhi <= nlo || nhi > (1ULL << 62)) die("need NMIN <= N_LO < N_HI <= 2^62");
    TH0 = tables_hash();
    Q = malloc(JCH * sizeof *Q);
    NBLK = (PP + RB - 1) / RB; BRES = calloc(NBLK, sizeof *BRES);
#ifdef Z2_CUDA
    gpu_setup();
#endif
    u64 tot[8] = { 0 }; double w0 = now(); clock_t c0 = clock(); u64 nch = 0;
    for (u64 c = nlo / PP / JCH; c * JCH * PP < nhi; c++) {
        CH_JA = c * JCH; CH_JB = CH_JA + JCH;
        CH_N0 = CH_JA * PP > nlo ? CH_JA * PP : nlo; CH_N1 = CH_JB * PP < nhi ? CH_JB * PP : nhi;
        if (CH_N0 >= CH_N1) continue;
        build_q();
        u64 qh = fnv(1469598103934665603ULL, Q, JCH * sizeof *Q);
        u64 cc[8] = { 0 };
#ifdef Z2_CUDA
        gpu_chunk(cc, qh);
#else
        memset(BRES, 0, NBLK * sizeof *BRES); NEXTB = 0;
        pthread_t *tid = malloc(NTH * sizeof *tid);
        for (int t = 0; t < NTH; t++) if (pthread_create(&tid[t], NULL, worker, NULL)) die("pthread_create failed");
        for (int t = 0; t < NTH; t++) pthread_join(tid[t], NULL);
        free(tid);
        // integrity: tables and Q unchanged, every block's running state re-derived
        if (tables_hash() != TH0) fail("the tables changed");
        if (fnv(1469598103934665603ULL, Q, JCH * sizeof *Q) != qh) fail("the Q table changed");
        for (u64 b = 0; b < NBLK; b++) { if (BRES[b].bad) fail("a block's running state disagrees with modular exponentiation"); for (int k = 0; k < 8; k++) cc[k] += BRES[b].cnt[k]; }
        for (u64 b = 0; b < NBLK; b++) if (BRES[b].on) { fwrite(BRES[b].out, 1, BRES[b].on, stdout); }
        for (u64 b = 0; b < NBLK; b++) free(BRES[b].out);
#endif
        printf("CDONE j=[%llu,%llu) n=[%llu,%llu) s1=%llu s2=%llu s3=%llu s4=%llu s5=%llu cand=%llu members=%llu\n",
               (unsigned long long)CH_JA, (unsigned long long)CH_JB, (unsigned long long)CH_N0, (unsigned long long)CH_N1,
               (unsigned long long)cc[1], (unsigned long long)cc[2], (unsigned long long)cc[3], (unsigned long long)cc[4],
               (unsigned long long)cc[5], (unsigned long long)cc[6], (unsigned long long)cc[7]);
        fflush(stdout);
        for (int k = 0; k < 8; k++) tot[k] += cc[k];
        nch++;
    }
    double wall = now() - w0, cpu = (double)(clock() - c0) / CLOCKS_PER_SEC;
    printf("DONE n=[%llu,%llu) K=%d NL=%d TB=%d TD=%d variant=%s chunks=%llu s1=%llu s2=%llu s3=%llu s4=%llu s5=%llu cand=%llu members=%llu time=%.2f threads=%d wall=%.2f\n",
           (unsigned long long)nlo, (unsigned long long)nhi, KK, NL, TB, TD, VARIANT, (unsigned long long)nch,
           (unsigned long long)tot[1], (unsigned long long)tot[2], (unsigned long long)tot[3], (unsigned long long)tot[4],
           (unsigned long long)tot[5], (unsigned long long)tot[6], (unsigned long long)tot[7], cpu, NTH, wall);
    return 0;
}
