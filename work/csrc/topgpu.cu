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
 *   Heights below Smin (|w| > 2^-10) are left to the host. nterm[S] series terms: |w|^(n+1) / (1 - |w|) < 2^-FL.
 *   frac(|Theta| lambda): the low LF + 2 limbs of TH Lambda (TH = round(|Theta| 2^64), Lambda = lambda 2^FL), whose top two limbs
 *   are the top 64 fraction bits (only TH mod 2^(FL + 64) matters).
 * Error (ulps of 2^-FL): tables 1 each; per level w 3, log1p and 1/(1+w) 2 nterm + 2 each, 1/A_i 3; per case 3 + 2 nterm.
 *   The host takes FL >= bits(Theta) + 80, so |Theta| err(lambda) < 2^-66 for any j <= 40, and adds lambda 2^-65 (rounding of
 *   TH), 2^-64 (truncated product), eta and a factor 4 to every margin. */
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

/* 0: certainly killed at level 1; 1: certain level-1 pass, certainly killed at level 2; 2: the host decides */
template <int L> DF int classify(const Big<L> &lam, const TCtx &C) {
    const tp_params &P = *C.P;
    const u64 psi1 = frac_theta<L>(lam, C.T + P.o_th);
    const u64 d1 = psi1 < (1ULL << 63) ? psi1 : (u64)0 - psi1;
    if (d1 > P.rho_hi) return 0;
    if (d1 >= P.rho_lo) return 2;                        /* within the margins of rho */
    const u64 psi2 = frac_theta<L>(lam, C.T + P.o_th2);
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

template <int L>
DF void do_prefix(const TCtx &C, const tp_prefix &pf, unsigned idx, tp_out *out, long long outcap, unsigned long long *nout,
                  unsigned long long cnt[2]) {
    const tp_params &P = *C.P; const u32 *T = C.T;
    const int j = P.j, SS = P.smax + 1;
    const u32 *csw = T + P.o_cs + (size_t)pf.o * ((j + 3) / 4);
    Big<L> lam, ia;                                      /* lambda_A = log A, 1/A */
    int acc = 0;
    cnt[0] += (unsigned long long)(pf.ghi - pf.glo + 1);
    if (j >= 2) {
        acc = pf.g[0];
        const int c1 = (int)((LDG(csw) >> 0) & 0xFF);
        ldb(lam, T + P.o_lam1 + ((size_t)c1 * SS + acc) * L);
        ldb(ia, T + P.o_ia1 + ((size_t)c1 * SS + acc) * L);
        for (int i = 1; i < j - 1; i++) {
            acc += pf.g[i];
            if (acc < P.Smin) { putr(out, outcap, nout, idx, 0, 2); return; }   /* the host does the whole prefix */
            const int ci = (int)((LDG(csw + (i >> 2)) >> (8 * (i & 3))) & 0xFF);
            Big<L> cw, w, lw, q, t; ldb(cw, T + P.o_c10 + ((size_t)ci * SS + acc) * L);
            fmul(w, cw, ia);                             /* w = c_i 10^-S_i / A_{i-1} */
            const int N = P.nterm[acc];
            flog1p(lw, w, N, T + P.o_inv); badd(lam, lam, lw);
            finv1p(q, w, N); fmul(t, ia, q); ia = t;
        }
    } else {
        bzero_(lam); bzero_(ia); ia.w[L - 1] = 1u;
    }
    const int cl = (int)((LDG(csw + ((j - 1) >> 2)) >> (8 * ((j - 1) & 3))) & 0xFF);
    for (int gl = pf.glo; gl <= pf.ghi; gl++) {
        const int e = acc + gl;
        if (e < P.Smin) { putr(out, outcap, nout, idx, (unsigned)gl, 1); continue; }
        Big<L> cw, u, lu, la; ldb(cw, T + P.o_c10 + ((size_t)cl * SS + e) * L);
        fmul(u, cw, ia);                                 /* u = c 10^-e / A */
        flog1p(lu, u, P.nterm[e], T + P.o_inv); badd(la, lam, lu);
        if (P.dump) {                                    /* validation only: psi1, psi2 of every case */
            const unsigned long long k = bump(nout);
            if ((long long)k < outcap) {
                out[k].prefix = idx; out[k].gl = (unsigned)gl; out[k].kind = 3; out[k].pad = 0;
                out[k].psi1 = frac_theta<L>(la, T + P.o_th); out[k].psi2 = frac_theta<L>(la, T + P.o_th2);
            }
            continue;
        }
        const int k = classify<L>(la, C);
        if (k == 1) cnt[1]++;
        else if (k == 2) putr(out, outcap, nout, idx, (unsigned)gl, 1);
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
