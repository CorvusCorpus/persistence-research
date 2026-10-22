/* Fixed-width unsigned multiprecision arithmetic for the GPU bottom solver (bottomgpu.cu): L limbs of 32 bits,
 * little endian, every loop bound a compile-time constant so that the limbs stay in registers. Host and device.
 * Exact integer arithmetic only:
 *   - Montgomery multiplication modulo an odd m < R = 2^(32 L) (CIOS; any m < R, so one L serves 5^N, 5^G and 5^(G-1));
 *     mont(a, b) = a b R^-1 mod m for a, b < m. With one operand stored as B R mod m ("Montgomery form") the result is the
 *     plain residue a B mod m, which is how the tables are used.
 *   - truncated products modulo 2^(32 L) (the 2-adic side), masks to 2^k.
 *   - small helpers: add/sub (mod m), shifts, exact division by a small integer, trailing zeros. */
#ifndef BIGFIX_H
#define BIGFIX_H
#include <stdint.h>
#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD static inline
#endif
typedef uint32_t u32; typedef uint64_t u64;

template <int L> struct Big { u32 w[L]; };

template <int L> HD void bzero_(Big<L> &a) {
#pragma unroll
    for (int i = 0; i < L; i++) a.w[i] = 0;
}
template <int L> HD void bset(Big<L> &a, u64 v) { bzero_(a); a.w[0] = (u32)v; if (L > 1) a.w[1] = (u32)(v >> 32); }
template <int L> HD int bcmp_(const Big<L> &a, const Big<L> &b) {   /* -1, 0, 1 */
    int r = 0;
#pragma unroll
    for (int i = 0; i < L; i++) r = a.w[i] < b.w[i] ? -1 : (a.w[i] > b.w[i] ? 1 : r);
    return r;
}
template <int L> HD int bzero_p(const Big<L> &a) {
    u32 o = 0;
#pragma unroll
    for (int i = 0; i < L; i++) o |= a.w[i];
    return o == 0;
}
template <int L> HD u32 badd(Big<L> &r, const Big<L> &a, const Big<L> &b) {   /* returns the carry */
    u64 c = 0;
#pragma unroll
    for (int i = 0; i < L; i++) { c += (u64)a.w[i] + b.w[i]; r.w[i] = (u32)c; c >>= 32; }
    return (u32)c;
}
template <int L> HD u32 bsub(Big<L> &r, const Big<L> &a, const Big<L> &b) {   /* returns the borrow */
    u64 br = 0;
#pragma unroll
    for (int i = 0; i < L; i++) { u64 d = (u64)a.w[i] - b.w[i] - br; r.w[i] = (u32)d; br = (d >> 32) & 1; }
    return (u32)br;
}
template <int L> HD void baddmod(Big<L> &r, const Big<L> &a, const Big<L> &b, const Big<L> &m) {   /* a, b < m */
    Big<L> t; u32 c = badd(t, a, b); Big<L> s; u32 br = bsub(s, t, m);
    r = (c || !br) ? s : t;
}
template <int L> HD void bsubmod(Big<L> &r, const Big<L> &a, const Big<L> &b, const Big<L> &m) {   /* a, b < m */
    Big<L> t; u32 br = bsub(t, a, b); Big<L> s; badd(s, t, m);
    r = br ? s : t;
}
/* r = a * s (+ carry out), s < 2^32 */
template <int L> HD u32 bmul_u32(Big<L> &r, const Big<L> &a, u32 s) {
    u64 c = 0;
#pragma unroll
    for (int i = 0; i < L; i++) { c += (u64)a.w[i] * s; r.w[i] = (u32)c; c >>= 32; }
    return (u32)c;
}
/* Montgomery product a b R^-1 mod m, R = 2^(32 L); a, b < m < R, m odd, mp = -m^-1 mod 2^32 (CIOS) */
template <int L> HD void bmont(Big<L> &r, const Big<L> &a, const Big<L> &b, const Big<L> &m, u32 mp) {
    u32 t[L + 2];
#pragma unroll
    for (int i = 0; i < L + 2; i++) t[i] = 0;
#pragma unroll
    for (int i = 0; i < L; i++) {
        u64 c = 0;
#pragma unroll
        for (int j = 0; j < L; j++) { c += (u64)a.w[j] * b.w[i] + t[j]; t[j] = (u32)c; c >>= 32; }
        c += t[L]; t[L] = (u32)c; t[L + 1] = (u32)(c >> 32);
        const u32 q = t[0] * mp;
        c = ((u64)q * m.w[0] + t[0]) >> 32;
#pragma unroll
        for (int j = 1; j < L; j++) { c += (u64)q * m.w[j] + t[j]; t[j - 1] = (u32)c; c >>= 32; }
        c += t[L]; t[L - 1] = (u32)c; t[L] = t[L + 1] + (u32)(c >> 32); t[L + 1] = 0;
    }
    /* t < 2m: subtract m once if t >= m */
    Big<L> x, y;
#pragma unroll
    for (int i = 0; i < L; i++) x.w[i] = t[i];
    u32 br = bsub(y, x, m);
    r = (t[L] || !br) ? y : x;
}
/* truncated product a b mod 2^(32 L) */
template <int L> HD void bmul_lo(Big<L> &r, const Big<L> &a, const Big<L> &b) {
    u32 t[L];
#pragma unroll
    for (int i = 0; i < L; i++) t[i] = 0;
#pragma unroll
    for (int i = 0; i < L; i++) {
        u64 c = 0;
#pragma unroll
        for (int j = 0; j + i < L; j++) { c += (u64)a.w[i] * b.w[j] + t[i + j]; t[i + j] = (u32)c; c >>= 32; }
    }
#pragma unroll
    for (int i = 0; i < L; i++) r.w[i] = t[i];
}
/* full product into 2L limbs */
template <int L> HD void bmul_full(u32 *r, const Big<L> &a, const Big<L> &b) {
#pragma unroll
    for (int i = 0; i < 2 * L; i++) r[i] = 0;
#pragma unroll
    for (int i = 0; i < L; i++) {
        u64 c = 0;
#pragma unroll
        for (int j = 0; j < L; j++) { c += (u64)a.w[i] * b.w[j] + r[i + j]; r[i + j] = (u32)c; c >>= 32; }
        r[i + L] = (u32)c;
    }
}
/* keep the low k bits (0 <= k <= 32 L) */
template <int L> HD void bmask(Big<L> &a, int k) {
#pragma unroll
    for (int i = 0; i < L; i++) {
        const int lo = 32 * i;
        a.w[i] = k <= lo ? 0u : (k >= lo + 32 ? a.w[i] : a.w[i] & ((1u << (k - lo)) - 1u));
    }
}
/* logical shifts by s bits, 0 <= s < 32 L: whole limbs by a barrel of compile-time offsets (no dynamic indexing, so the
 * limbs stay in registers on the GPU), then the bit part */
template <int L> HD void bshr(Big<L> &r, const Big<L> &a, int s) {
    Big<L> t = a;
    const int q = s >> 5, b = s & 31;
#pragma unroll
    for (int k = 1; k < L; k <<= 1) {
        if (q & k) {
#pragma unroll
            for (int i = 0; i < L; i++) t.w[i] = i + k < L ? t.w[i + k] : 0u;
        }
    }
    if (b) {
#pragma unroll
        for (int i = 0; i < L; i++) t.w[i] = (t.w[i] >> b) | (i + 1 < L ? t.w[i + 1] << (32 - b) : 0u);
    }
    r = t;
}
template <int L> HD void bshl(Big<L> &r, const Big<L> &a, int s) {
    Big<L> t = a;
    const int q = s >> 5, b = s & 31;
#pragma unroll
    for (int k = 1; k < L; k <<= 1) {
        if (q & k) {
#pragma unroll
            for (int i = L - 1; i >= 0; i--) t.w[i] = i - k >= 0 ? t.w[i - k] : 0u;
        }
    }
    if (b) {
#pragma unroll
        for (int i = L - 1; i >= 0; i--) t.w[i] = (t.w[i] << b) | (i >= 1 ? t.w[i - 1] >> (32 - b) : 0u);
    }
    r = t;
}
HD int ctz32_(u32 x) {
#ifdef __CUDA_ARCH__
    return __ffs((int)x) - 1;
#else
    return __builtin_ctz(x);
#endif
}
/* number of trailing zero bits (32 L if a = 0) */
template <int L> HD int bctz(const Big<L> &a) {
    int r = 32 * L;
#pragma unroll
    for (int i = L - 1; i >= 0; i--) if (a.w[i]) r = 32 * i + ctz32_(a.w[i]);
    return r;
}
/* a mod 2^32 * k ... small residues: a mod D for a compile-time D (the division by D compiles to multiplications) */
template <u32 D, int L> HD u32 bmod_small(const Big<L> &a) {
    u64 rem = 0;
#pragma unroll
    for (int i = L - 1; i >= 0; i--) rem = ((rem << 32) | a.w[i]) % D;
    return (u32)rem;
}
#endif
