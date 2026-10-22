/* Kernel of the GPU top filter: classifies the cases of csrc/top3.c's loop, so that the host runs top3.c's exact code only on
 * the cases that are not certainly killed (docs/GPU_EVEN_TARGETS.md section 3).
 *
 * top3.c kills a case at level 1 when dist(m theta, Z) > rho and at level 2 by the next Fincke-Pohst level, with
 * m = round(-C lambda), lambda = log(1 + z). Since m = -C lambda + delta (|delta| <= 1/2 + 2^-50) and theta, theta2 are tiny
 * (the host checks |theta|, |theta2| < 2^-60), m theta = -Theta lambda + delta theta with Theta = C theta: both levels are
 * decided by psi1 = frac(-Theta lambda) and psi2 = frac(-Theta2 lambda) up to eta = 0.6 max(|theta|, |theta2|). The kernel
 * computes psi1, psi2 in fixed point with a rigorous error bound and reports a case as killed only when the exact test is
 * certain to kill it; it also counts the certain level-1 passes among the cases it kills at level 2 (top3.c's f1).
 *
 * Fixed point: values in [0, 2^32) with FL = 32 LF fraction bits in L = LF + 1 limbs; products truncated (error < 1 ulp).
 *   prefix (digits c_1 .. c_{j-1} at heights S_i = g_0 + .. + g_{i-1}): level 1 from host tables log(1 + c 10^-S),
 *   1/(1 + c 10^-S); level i >= 2: w = c_i 10^-S_i / A_{i-1}, lambda_i = lambda_{i-1} + log1p(w), 1/A_i = (1/A_{i-1}) / (1 + w)
 *   by Horner: log1p(w) = w (1 - w (1/2 - w (1/3 - ...))), 1/(1 + w) = 1 - w (1 - w (1 - ...)), every intermediate in (0, 1].
 *   case (last digit c at height e): lambda = lambda_A + log1p(c 10^-e / A).
 *   Digits at heights below Smin = 5 (|w| could exceed 2^-10; at most four digits, heights 1 .. 4) are not summed: the leading
 *   digits below Smin take lambda_A and 1/A from a host table over every such configuration (8^4 entries), and so does a case
 *   whose last digit is below Smin. nterm[S] series terms: |w|^(n+1) / (1 - |w|) < 2^-FL.
 *   frac(|Theta| lambda): the low LF + 2 limbs of TH Lambda (TH = round(|Theta| 2^64), Lambda = lambda 2^FL), whose top two limbs
 *   are the top 64 fraction bits (only TH mod 2^(FL + 64) matters).
 * Error (ulps of 2^-FL): tables 1 each; per level w 3, log1p and 1/(1+w) 2 nterm + 2 each, 1/A_i 3; per case 3 + 2 nterm.
 *   The host takes FL >= bits(Theta) + 80, so |Theta| err(lambda) < 2^-66 for any j <= 40, and adds lambda 2^-65 (rounding of
 *   TH), 2^-64 (truncated product), eta and a factor 4 to every margin.
 * Level 1 along a prefix (recurrence; the cases of one prefix have consecutive heights e): frac(|Theta| lambda) =
 *   frac(phi_A + T_1 - T_2 + T_3 ...) with phi_A = frac(|Theta| lambda_A) (96 bits of TH Lambda_A) and T_n = |Theta| u^n / n,
 *   u = c 10^-e / A, in fixed point with 96 fraction bits. At the prefix's first height with nrec[e] > 0: T_1 = TH U / 2^(FL - 32)
 *   (U = u 2^FL), T_n = T_(n-1) U (n - 1) / n; at each next height T_n = floor(T_n / 10^n), as u(e + 1) = u(e) / 10 exactly.
 *   nrec[e] (host): the fewest N <= TP_NREC terms with omitted tail |Theta| (72 10^-e)^(N+1) / (N + 1) / (1 - 2^-10) < 2^-90, else
 *   0 (direct series). Errors: T_1 at the start |Theta| err(U) 2^-FL + u 2^-65 + 2^-96 < 2^-70 (err(U) a few ulps plus 2^-6 err(1/A),
 *   |Theta| < 2^(FL - 80)), T_n alike but smaller; each step adds at most 2^-96 and divides the earlier error by 10^n; tail < 2^-90;
 *   the 96-bit sum is truncated to 64 bits. So psi1 stays within the budget of the direct path. Level 2 (after a certain level-1
 *   pass, 1-2% of cases) and the dump mode's psi2 use the direct lambda. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bigfix.h"
#include "topgpu.h"

#ifdef TP_EMU
#define DF static inline
#define LDG(p) (*(p))
#else
#define DF __device__ __forceinline__
#define LDG(p) __ldg(p)
#endif

template <int L> DF void ldb(Big<L> &x, const u32 *p) {
#pragma unroll
    for (int i = 0; i < L; i++) x.w[i] = LDG(p + i);
}
/* fixed point product: limbs [L - 1, 2L - 1) of a b (truncation); FL = 32 (L - 1) */
template <int L> DF void fmul_body(Big<L> &r, const Big<L> &a, const Big<L> &b) {
    u32 t[2 * L]; bmul_full<L>(t, a, b);
#pragma unroll
    for (int i = 0; i < L; i++) r.w[i] = t[L - 1 + i];
}
#ifndef TP_EMU
/* for L >= 24 the product stays out of line: nvcc 13.3 crashes (cicc, exit 139) on the fully inlined kernel */
template <int L> __device__ __noinline__ void fmul_ni(Big<L> &r, const Big<L> &a, const Big<L> &b) { fmul_body<L>(r, a, b); }
#endif
template <int L> DF void fmul(Big<L> &r, const Big<L> &a, const Big<L> &b) {
#ifndef TP_EMU
    if constexpr (L >= 24) { fmul_ni<L>(r, a, b); return; }
#endif
    fmul_body<L>(r, a, b);
}
/* log1p(w), 0 <= w < 2^-10, N >= 1 terms: w (1 - w (1/2 - w (1/3 - ... - w / N))); inv + n L = 1/n */
template <int L> DF void flog1p(Big<L> &r, const Big<L> &w, int N, const u32 *inv) {
    Big<L> Pp, c, t; ldb(Pp, inv + (size_t)N * L);
    for (int n = N - 1; n >= 1; n--) { fmul(t, w, Pp); ldb(c, inv + (size_t)n * L); bsub(Pp, c, t); }
    fmul(r, w, Pp);
}
/* 1/(1 + w), 0 <= w < 2^-10, N terms: 1 - w (1 - w (1 - ...)) */
template <int L> DF void finv1p(Big<L> &r, const Big<L> &w, int N) {
    Big<L> one, Q, t; bzero_(one); one.w[L - 1] = 1u; Q = one;
    for (int n = 0; n < N; n++) { fmul(t, w, Q); bsub(Q, one, t); }
    r = Q;
}

struct TCtx { const tp_params *P; const u32 *T; };

/* top 64 bits of frac(-Theta lambda); th[0] = 1 if Theta < 0, th[1 ..] = low L + 1 limbs of round(|Theta| 2^64) */
template <int L> DF u64 frac_theta(const Big<L> &lam, const u32 *th) {
    Big<L + 1> a, b, p;
#pragma unroll
    for (int i = 0; i < L + 1; i++) { a.w[i] = LDG(th + 1 + i); b.w[i] = i < L ? lam.w[i] : 0u; }
    bmul_lo<L + 1>(p, a, b);                             /* (TH Lambda) mod 2^(32 (L + 1)) */
    const u64 f = ((u64)p.w[L] << 32) | p.w[L - 1];      /* top 64 bits of frac(|Theta| lambda) */
    return LDG(th) ? f : (u64)0 - f;                     /* Theta < 0: frac(-Theta lambda) = frac(|Theta| lambda) */
}

/* the top 96 bits of frac(|Theta| lambda), low limb first (as frac_theta, before the sign) */
template <int L> DF void frac_theta96(u32 f[3], const Big<L> &lam, const u32 *th) {
    Big<L + 1> a, b, p;
#pragma unroll
    for (int i = 0; i < L + 1; i++) { a.w[i] = LDG(th + 1 + i); b.w[i] = i < L ? lam.w[i] : 0u; }
    bmul_lo<L + 1>(p, a, b);
    f[0] = p.w[L - 2]; f[1] = p.w[L - 1]; f[2] = p.w[L];
}
/* recurrence start: T = floor(|Theta| u 2^96) from TH = round(|Theta| 2^64) and U = u 2^FL: limbs [L - 2, 2L - 2) of TH U
 * (|Theta| < 2^(FL - 80) and u < 2^-10, so the value is below 2^(32 L - 26) and the limbs above are zero) */
template <int L> DF void th_mul_body(Big<L> &r, const u32 *th, const Big<L> &U) {
    Big<L + 1> a, b; u32 t[2 * L + 2];
#pragma unroll
    for (int i = 0; i < L + 1; i++) { a.w[i] = LDG(th + 1 + i); b.w[i] = i < L ? U.w[i] : 0u; }
    bmul_full<L + 1>(t, a, b);
#pragma unroll
    for (int i = 0; i < L; i++) r.w[i] = t[L - 2 + i];
}
#ifndef TP_EMU
template <int L> __device__ __noinline__ void th_mul_ni(Big<L> &r, const u32 *th, const Big<L> &U) { th_mul_body<L>(r, th, U); }
#endif
template <int L> DF void th_mul(Big<L> &r, const u32 *th, const Big<L> &U) {
#ifndef TP_EMU
    if constexpr (L >= 24) { th_mul_ni<L>(r, th, U); return; }
#endif
    th_mul_body<L>(r, th, U);
}
/* a = floor(a / D), D < 2^32 a compile-time constant */
template <u32 D, int L> DF void bdivc(Big<L> &a) {
    u64 r = 0;
#pragma unroll
    for (int i = L - 1; i >= 0; i--) { const u64 cur = (r << 32) | a.w[i]; const u64 q = cur / D; a.w[i] = (u32)q; r = cur - q * D; }
}
/* f +/- the low 3 limbs of t, mod 2^96 */
template <int L> DF void add96(u32 f[3], const Big<L> &t) {
    u64 c = (u64)f[0] + t.w[0]; f[0] = (u32)c; c >>= 32; c += (u64)f[1] + t.w[1]; f[1] = (u32)c; c >>= 32; f[2] += t.w[2] + (u32)c;
}
template <int L> DF void sub96(u32 f[3], const Big<L> &t) {
    u64 d = (u64)f[0] - t.w[0]; f[0] = (u32)d; u32 br = (u32)(d >> 32) & 1u;
    d = (u64)f[1] - t.w[1] - br; f[1] = (u32)d; br = (u32)(d >> 32) & 1u; f[2] = f[2] - t.w[2] - br;
}

/* level 1 from psi1: 0 certainly killed; 2 the host decides (within the margins of rho); 3 certain pass (level 2 next) */
DF int level1(u64 psi1, const tp_params &P) {
    const u64 d1 = psi1 < (1ULL << 63) ? psi1 : (u64)0 - psi1;
    return d1 > P.rho_hi ? 0 : (d1 >= P.rho_lo ? 2 : 3);
}
/* level 2 after a certain level-1 pass: 1 certainly killed (counted in f1); 2 the host decides */
DF int level2(u64 psi1, u64 psi2, const tp_params &P) {
    const long long s = (long long)psi1;                 /* signed distance of psi1 to Z (|s| < rho < 2^62), units 2^-64 */
#ifdef TP_EMU
    const long long ms = (long long)(((__int128)P.mu63 * s) >> 63);
#else
    const long long hi = __mul64hi(P.mu63, s); const unsigned long long lo = (unsigned long long)P.mu63 * (unsigned long long)s;
    const long long ms = (long long)(((unsigned long long)hi << 1) | (lo >> 63));
#endif
    const u64 y = (u64)0 - psi2 - (u64)ms;               /* frac(-m theta2 - mu s) */
    const u64 dy = y < (1ULL << 63) ? y : (u64)0 - y;
    const u64 as = s < 0 ? (u64)(-s) : (u64)s;
    const u64 a = dy > P.dy ? dy - P.dy : 0, bq = as > P.ds ? as - P.ds : 0;   /* lower bounds of dist(y), |s| */
    const u64 a32 = a >> 32, b32 = bq >> 32;             /* truncated: still lower bounds; a32, b32 < 2^31 */
    const u64 x1 = a32 * a32, x2 = b32 * b32;            /* < 2^62 */
#ifdef TP_EMU
    const unsigned __int128 lhs = (unsigned __int128)x1 * P.g2l + (unsigned __int128)x2 * P.g3l;
    return lhs > ((unsigned __int128)1 << 96) ? 1 : 2;
#else
    const u64 l1 = x1 * P.g2l, h1 = __umul64hi(x1, P.g2l), l2 = x2 * P.g3l, h2 = __umul64hi(x2, P.g3l);
    const u64 lo2 = l1 + l2, hi2 = h1 + h2 + (lo2 < l1 ? 1 : 0);
    return (hi2 > (1ULL << 32) || (hi2 == (1ULL << 32) && lo2 > 0)) ? 1 : 2;
#endif
}

DF unsigned long long bump(unsigned long long *n) {
#ifdef TP_EMU
    return (*n)++;
#else
    return atomicAdd(n, 1ULL);
#endif
}
DF void putr(tp_out *out, long long outcap, unsigned long long *nout, unsigned prefix, unsigned gl, unsigned kind) {
    const unsigned long long k = bump(nout);
    if ((long long)k < outcap) { out[k].prefix = prefix; out[k].gl = gl; out[k].kind = kind; out[k].pad = 0; out[k].psi1 = out[k].psi2 = 0; }
}
/* lambda of the case with last digit c (index cl) at height e >= Smin: lambda_A + log1p(c 10^-e / A) */
template <int L> DF void case_lambda(Big<L> &la, const Big<L> &lam, const Big<L> &ia, int cl, int e, const TCtx &C) {
    const tp_params &P = *C.P;
    Big<L> cw, u, lu; ldb(cw, C.T + P.o_c10 + ((size_t)cl * (P.smax + 1) + e) * L);
    fmul(u, cw, ia);                                     /* u = c 10^-e / A */
    flog1p(lu, u, P.nterm[e], C.T + P.o_inv); badd(la, lam, lu);
}

template <int L>
DF void do_prefix(const TCtx &C, const tp_prefix &pf, unsigned idx, tp_out *out, long long outcap, unsigned long long *nout,
                  unsigned long long cnt[2]) {
    const tp_params &P = *C.P; const u32 *T = C.T;
    const int j = P.j, SS = P.smax + 1;
    const u32 *csw = T + P.o_cs + (size_t)pf.o * ((j + 3) / 4);
    Big<L> lam, ia;                                      /* lambda_A = log A, 1/A */
    int acc = 0, i = 0;
    unsigned sidx = 0;                                   /* the digits at heights below Smin: one entry of the small table */
    cnt[0] += (unsigned long long)(pf.ghi - pf.glo + 1);
    for (; i < j - 1 && acc + pf.g[i] < P.Smin; i++) {
        acc += pf.g[i];
        sidx |= (((LDG(csw + (i >> 2)) >> (8 * (i & 3))) & 0xFFu) + 1u) << (3 * (acc - 1));
    }
    if (i == 0 && j >= 2) {                              /* first digit at height >= Smin: level-1 tables */
        acc = pf.g[0];
        const int c1 = (int)((LDG(csw) >> 0) & 0xFF);
        ldb(lam, T + P.o_lam1 + ((size_t)c1 * SS + acc) * L);
        ldb(ia, T + P.o_ia1 + ((size_t)c1 * SS + acc) * L);
        i = 1;
    } else {                                             /* entry 0 (no digit, j = 1): lambda 0, 1/A = 1 */
        ldb(lam, T + P.o_sl + (size_t)sidx * L);
        ldb(ia, T + P.o_sia + (size_t)sidx * L);
    }
    for (; i < j - 1; i++) {                             /* heights >= Smin: |w| <= 72 10^-5 < 2^-10 */
        acc += pf.g[i];
        const int ci = (int)((LDG(csw + (i >> 2)) >> (8 * (i & 3))) & 0xFF);
        Big<L> cw, w, lw, q, t; ldb(cw, T + P.o_c10 + ((size_t)ci * SS + acc) * L);
        fmul(w, cw, ia);                                 /* w = c_i 10^-S_i / A_{i-1} */
        const int N = P.nterm[acc];
        flog1p(lw, w, N, T + P.o_inv); badd(lam, lam, lw);
        finv1p(q, w, N); fmul(t, ia, q); ia = t;
    }
    const int cl = (int)((LDG(csw + ((j - 1) >> 2)) >> (8 * ((j - 1) & 3))) & 0xFF);
    const u32 *th = T + P.o_th;
    /* recurrence over the consecutive heights e with nrec[e] > 0 (all heights above the first such one): T_n = |Theta| u^n / n,
     * u = c 10^-e / A, 96 fraction bits; T_n(e + 1) = T_n(e) / 10^n; frac(|Theta| lambda) = frac(phi_A + T_1 - T_2 + T_3) */
    Big<L> T1, T2, T3; u32 phA[3] = {0u, 0u, 0u}; int nr = 0;
    for (int gl = pf.glo; gl <= pf.ghi; gl++) {
        const int e = acc + gl;
        Big<L> la; int have = 0, path; u64 psi1;
        if (e < P.Smin) {                                /* every digit below Smin (then i = j - 1): the small table */
            ldb(la, T + P.o_sl + (size_t)(sidx | ((unsigned)(cl + 1) << (3 * (e - 1)))) * L);
            have = 1; path = 0; psi1 = frac_theta<L>(la, th);
        } else if (nr == 0 && P.nrec[e] == 0) {          /* more series terms than the recurrence carries: direct */
            case_lambda<L>(la, lam, ia, cl, e, C); have = 1; path = 1; psi1 = frac_theta<L>(la, th);
        } else {
            if (nr == 0) {                               /* start: u from the tables, the T_n by products */
                nr = P.nrec[e]; path = 2;
                frac_theta96<L>(phA, lam, th);
                Big<L> cw, U; ldb(cw, T + P.o_c10 + ((size_t)cl * SS + e) * L);
                fmul(U, cw, ia);                         /* U = u 2^FL */
                th_mul<L>(T1, th, U);
                if (nr >= 2) { fmul(T2, T1, U); bshr(T2, T2, 1); }
                if (TP_NREC >= 3 && nr >= 3) { fmul(T3, T2, U); bmul_u32(T3, T3, 2u); bdivc<3u>(T3); }
            } else {                                     /* step e - 1 -> e: u / 10 */
                path = 3;
                bdivc<10u>(T1);
                if (nr >= 2) bdivc<100u>(T2);
                if (TP_NREC >= 3 && nr >= 3) bdivc<1000u>(T3);
            }
            u32 F[3] = {phA[0], phA[1], phA[2]};
            add96(F, T1);
            if (nr >= 2) sub96(F, T2);
            if (TP_NREC >= 3 && nr >= 3) add96(F, T3);
            const u64 f = ((u64)F[2] << 32) | F[1];      /* top 64 bits of frac(|Theta| lambda) */
            psi1 = LDG(th) ? f : (u64)0 - f;
        }
        if (P.dump) {                                    /* validation only: psi1 (and its path), psi2 of every case */
            if (!have) case_lambda<L>(la, lam, ia, cl, e, C);
            const unsigned long long k = bump(nout);
            if ((long long)k < outcap) {
                out[k].prefix = idx; out[k].gl = (unsigned)gl; out[k].kind = 3; out[k].pad = (unsigned)path;
                out[k].psi1 = psi1; out[k].psi2 = frac_theta<L>(la, T + P.o_th2);
            }
            continue;
        }
        const int k1 = level1(psi1, P);
        if (k1 == 0) continue;                           /* certainly killed at level 1 */
        if (k1 == 2) { putr(out, outcap, nout, idx, (unsigned)gl, 1); continue; }
        if (!have) case_lambda<L>(la, lam, ia, cl, e, C);   /* certain level-1 pass: level 2 needs psi2 */
        if (level2(psi1, frac_theta<L>(la, T + P.o_th2), P) == 1) cnt[1]++;
        else putr(out, outcap, nout, idx, (unsigned)gl, 1);
    }
}

/* ---------------- size classes: L = LF + 1 limbs ---------------- */
#define TP_CLASSES(X) X(0, 6) X(1, 8) X(2, 10) X(3, 12) X(4, 14) X(5, 16) X(6, 20) X(7, 22) X(8, 24) X(9, 28) X(10, 32)

#ifdef TP_EMU
static tp_params HP; static const u32 *HT;
extern "C" int tp_init(const tp_params *p, const unsigned int *tables) { HP = *p; HT = tables; return 0; }
extern "C" long long tp_run(const tp_prefix *pre, int n, tp_out *out, long long outcap, unsigned long long counts[2]) {
    unsigned long long nout = 0, cnt[2] = {0, 0}; TCtx C = {&HP, HT};
    for (int i = 0; i < n; i++) {
        switch (HP.cls) {
#define CASE(k, l) case k: do_prefix<l>(C, pre[i], (unsigned)i, out, outcap, &nout, cnt); break;
            TP_CLASSES(CASE)
#undef CASE
        }
    }
    if ((long long)nout > outcap) return -1;
    counts[0] += cnt[0]; counts[1] += cnt[1];
    return (long long)nout;
}
extern "C" const char *tp_device(void) { return "emulation"; }
#else
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(5); } } while (0)
__constant__ tp_params DP;
template <int L>
__global__ void k_prefixes(const u32 *T, const tp_prefix *pre, int n, tp_out *out, long long outcap, unsigned long long *nout,
                           unsigned long long *gcnt) {
    unsigned long long cnt[2] = {0, 0}; TCtx C = {&DP, T};
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x)
        do_prefix<L>(C, pre[i], (unsigned)i, out, outcap, nout, cnt);
    if (cnt[0]) atomicAdd(&gcnt[0], cnt[0]);
    if (cnt[1]) atomicAdd(&gcnt[1], cnt[1]);
}
static tp_params HP; static u32 *DT; static tp_prefix *DPF; static tp_out *DOUT; static unsigned long long *DCNT;
static int DPFCAP; static long long DOUTCAP; static char DEVNAME[256]; static int NSM;
extern "C" int tp_init(const tp_params *p, const unsigned int *tables) {
    HP = *p;
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr, 0)); NSM = pr.multiProcessorCount; snprintf(DEVNAME, sizeof DEVNAME, "%s", pr.name);
    CK(cudaMemcpyToSymbol(DP, p, sizeof *p));
    CK(cudaMalloc(&DT, (size_t)p->words * 4)); CK(cudaMemcpy(DT, tables, (size_t)p->words * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DCNT, 3 * 8));
    return 0;
}
extern "C" long long tp_run(const tp_prefix *pre, int n, tp_out *out, long long outcap, unsigned long long counts[2]) {
    if (n > DPFCAP) { if (DPF) cudaFree(DPF); DPFCAP = n; CK(cudaMalloc(&DPF, (size_t)n * sizeof(tp_prefix))); }
    if (outcap > DOUTCAP) { if (DOUT) cudaFree(DOUT); DOUTCAP = outcap; CK(cudaMalloc(&DOUT, (size_t)outcap * sizeof(tp_out))); }
    CK(cudaMemcpy(DPF, pre, (size_t)n * sizeof(tp_prefix), cudaMemcpyHostToDevice));
    CK(cudaMemset(DCNT, 0, 3 * 8));
    const int nt = 128; int grid = (n + nt - 1) / nt; if (grid > 64 * NSM) grid = 64 * NSM;
    switch (HP.cls) {
#define CASE(k, l) case k: k_prefixes<l><<<grid, nt>>>(DT, DPF, n, DOUT, outcap, DCNT + 2, DCNT); break;
        TP_CLASSES(CASE)
#undef CASE
    }
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    unsigned long long h[3]; CK(cudaMemcpy(h, DCNT, 3 * 8, cudaMemcpyDeviceToHost));
    long long no = (long long)h[2];
    if (no > outcap) return -1;
    counts[0] += h[0]; counts[1] += h[1];
    if (no) CK(cudaMemcpy(out, DOUT, (size_t)no * sizeof(tp_out), cudaMemcpyDeviceToHost));
    return no;
}
extern "C" const char *tp_device(void) { return DEVNAME; }
#endif
