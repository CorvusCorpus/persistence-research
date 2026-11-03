/* Kernel of the GPU MITM top probe: csrc/mitm_top's join on the GPU (docs/MITM_STEPS.md section 2).
 *
 * Per top configuration (a prefix of j - 1 digits and the last digit at height e): lambda = log(1 + z) in fixed point exactly as
 * topgpu.cu's direct path computes it (prefix: small-height table or level-1 tables, then log1p / 1/(1 + w) series per digit;
 * last digit: case_lambda, or the small table when e < Smin), then P_i = frac(Theta_i lambda), i = 0, 1, 2 (top 64 bits),
 * by topgpu's frac_theta with the host's sign flag of -Theta_i. The box [P_i - hw_i, P_i + hw_i] covers the cells
 * cell(P_i - hw_i) .. cell(P_i + hw_i) of each axis (cyclically; cell(x) = floor(x ncell_i / 2^64)), which is every cell it
 * meets when its width is below (ncell_i - 1)/ncell_i of the axis; wider boxes (host flag full) take every cell of the axis.
 * Each of these cells is looked up in an open-addressing hash (linear probing) of the bottom table's cells;
 * every entry of a cell found with |P_i - Q_i| <= hw_i (signed 64-bit differences) for i = 0, 1, 2 and digits compatible with
 * the configuration's (per value: configuration + entry <= bag) is a candidate, decided exactly by the host side (Python).
 * Error of P_i: topgpu.cu's analysis of the direct path (FL >= bits(Theta) + 80): |Theta| err(lambda) < 2^-66, rounding of
 * TH lambda 2^-65, truncation 2^-64, in all below 2 units of 2^-64; the host adds 64 units to hw_i, as topgpu's margin. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bigfix.h"
#include "mitmgpu.h"

#ifdef MG_EMU
#define DF static inline
#define LDG(p) (*(p))
#else
#define DF __device__ __forceinline__
#define LDG(p) __ldg(p)
#endif

/* ---------------- fixed point (topgpu.cu, unchanged) ---------------- */
template <int L> DF void ldb(Big<L> &x, const u32 *p) {
#pragma unroll
    for (int i = 0; i < L; i++) x.w[i] = LDG(p + i);
}
template <int L> DF void fmul_body(Big<L> &r, const Big<L> &a, const Big<L> &b) {
    u32 t[2 * L]; bmul_full<L>(t, a, b);
#pragma unroll
    for (int i = 0; i < L; i++) r.w[i] = t[L - 1 + i];
}
#ifndef MG_EMU
template <int L> __device__ __noinline__ void fmul_ni(Big<L> &r, const Big<L> &a, const Big<L> &b) { fmul_body<L>(r, a, b); }
#endif
template <int L> DF void fmul(Big<L> &r, const Big<L> &a, const Big<L> &b) {
#ifndef MG_EMU
    if constexpr (L >= 24) { fmul_ni<L>(r, a, b); return; }
#endif
    fmul_body<L>(r, a, b);
}
template <int L> DF void flog1p(Big<L> &r, const Big<L> &w, int N, const u32 *inv) {
    Big<L> Pp, c, t; ldb(Pp, inv + (size_t)N * L);
    for (int n = N - 1; n >= 1; n--) { fmul(t, w, Pp); ldb(c, inv + (size_t)n * L); bsub(Pp, c, t); }
    fmul(r, w, Pp);
}
template <int L> DF void finv1p(Big<L> &r, const Big<L> &w, int N) {
    Big<L> one, Q, t; bzero_(one); one.w[L - 1] = 1u; Q = one;
    for (int n = 0; n < N; n++) { fmul(t, w, Q); bsub(Q, one, t); }
    r = Q;
}
/* top 64 bits of frac(-Theta lambda); th[0] = 1 if Theta < 0, th[1 ..] = low L + 1 limbs of round(|Theta| 2^64) */
template <int L> DF u64 frac_theta(const Big<L> &lam, const u32 *th) {
    Big<L + 1> a, b, p;
#pragma unroll
    for (int i = 0; i < L + 1; i++) { a.w[i] = LDG(th + 1 + i); b.w[i] = i < L ? lam.w[i] : 0u; }
    bmul_lo<L + 1>(p, a, b);
    const u64 f = ((u64)p.w[L] << 32) | p.w[L - 1];
    return LDG(th) ? f : (u64)0 - f;
}

struct MCtx { const mg_params *P; const u32 *T; const mg_entry *E; const mg_slot *H; };

template <int L> DF void case_lambda(Big<L> &la, const Big<L> &lam, const Big<L> &ia, int cl, int e, const MCtx &C) {
    const mg_params &P = *C.P;
    Big<L> cw, u, lu; ldb(cw, C.T + P.o_c10 + ((size_t)cl * (P.smax + 1) + e) * L);
    fmul(u, cw, ia);                                     /* u = c 10^-e / A */
    flog1p(lu, u, P.nterm[e], C.T + P.o_inv); badd(la, lam, lu);
}

DF unsigned long long bump(unsigned long long *n) {
#ifdef MG_EMU
    return (*n)++;
#else
    return atomicAdd(n, 1ULL);
#endif
}
DF u64 cellof(u64 x, u64 n) {                            /* floor(x n / 2^64) */
#ifdef MG_EMU
    return (u64)(((unsigned __int128)x * n) >> 64);
#else
    return __umul64hi(x, n);
#endif
}
DF u64 hash_slot(u64 key, int hbits) { return (key * 0x9E3779B97F4A7C15ULL) >> (64 - hbits); }

/* the entries of the cells of P's box: candidates to out */
DF void probe(const MCtx &C, const u64 Pt[3], u64 sigT, unsigned idx, unsigned gl, mg_out *out, long long outcap,
              unsigned long long *nout, unsigned long long cnt[2]) {
    const mg_params &P = *C.P;
    u64 lo[3], span[3];
#pragma unroll
    for (int a = 0; a < 3; a++) {                        /* box width < (n - 1)/n of the axis: the cyclic run lo .. hi */
        const u64 n = P.ncell[a];
        if (P.full[a]) { lo[a] = 0; span[a] = n; continue; }
        lo[a] = cellof(Pt[a] - P.hw[a], n);
        const u64 hi = cellof(Pt[a] + P.hw[a], n);
        span[a] = (hi + n - lo[a]) % n + 1;
    }
    const u64 mask = (1ULL << P.hbits) - 1;
    for (u64 s0 = 0; s0 < span[0]; s0++) {
        const u64 c0 = (lo[0] + s0) % P.ncell[0];
        for (u64 s1 = 0; s1 < span[1]; s1++) {
            const u64 c1 = (lo[1] + s1) % P.ncell[1];
            for (u64 s2 = 0; s2 < span[2]; s2++) {
                const u64 c2 = (lo[2] + s2) % P.ncell[2];
                const u64 key = (c0 * P.ncell[1] + c1) * P.ncell[2] + c2;
                for (u64 h = hash_slot(key, P.hbits);; h = (h + 1) & mask) {
                    const u64 k = LDG(&C.H[h].key);
                    if (k == ~0ULL) break;
                    if (k != key) continue;
                    const unsigned first = LDG(&C.H[h].first), count = LDG(&C.H[h].count);
                    for (unsigned x = first; x < first + count; x++) {
                        cnt[1]++;
                        int in = 1;
#pragma unroll
                        for (int a = 0; a < 3; a++) {
                            const u64 dd = Pt[a] - LDG(&C.E[x].q[a]);
                            const u64 ad = (long long)dd < 0 ? (u64)0 - dd : dd;
                            if (ad > P.hw[a]) in = 0;
                        }
                        if (!in) continue;
                        const u64 sigE = LDG(&C.E[x].sig);
                        for (int c = 0; c < MG_NC && in; c++)
                            if (((sigT >> (6 * c)) & 63) + ((sigE >> (6 * c)) & 63) > ((P.cntsig >> (6 * c)) & 63)) in = 0;
                        if (!in) continue;
                        const unsigned long long r = bump(nout);
                        if ((long long)r < outcap) { out[r].prefix = idx; out[r].gl = gl; out[r].idx = x; out[r].p[0] = Pt[0]; out[r].p[1] = Pt[1]; out[r].p[2] = Pt[2]; }
                    }
                    break;
                }
            }
        }
    }
}

template <int L>
DF void do_prefix(const MCtx &C, const mg_prefix &pf, unsigned idx, mg_out *out, long long outcap, unsigned long long *nout,
                  unsigned long long cnt[2]) {
    const mg_params &P = *C.P; const u32 *T = C.T;
    const int j = P.j, SS = P.smax + 1;
    const u32 *csw = T + P.o_cs + (size_t)pf.o * ((j + 3) / 4);
    const u64 sigT = ((u64)LDG(T + P.o_sig + 2 * (size_t)pf.o + 1) << 32) | LDG(T + P.o_sig + 2 * (size_t)pf.o);
    Big<L> lam, ia;                                      /* lambda_A = log A, 1/A (topgpu.cu's prefix code) */
    int acc = 0, i = 0;
    unsigned sidx = 0;
    cnt[0] += (unsigned long long)(pf.ghi - pf.glo + 1);
    for (; i < j - 1 && acc + pf.g[i] < P.Smin; i++) {
        acc += pf.g[i];
        sidx |= (((LDG(csw + (i >> 2)) >> (8 * (i & 3))) & 0xFFu) + 1u) << (3 * (acc - 1));
    }
    if (i == 0 && j >= 2) {
        acc = pf.g[0];
        const int c1 = (int)((LDG(csw) >> 0) & 0xFF);
        ldb(lam, T + P.o_lam1 + ((size_t)c1 * SS + acc) * L);
        ldb(ia, T + P.o_ia1 + ((size_t)c1 * SS + acc) * L);
        i = 1;
    } else {
        ldb(lam, T + P.o_sl + (size_t)sidx * L);
        ldb(ia, T + P.o_sia + (size_t)sidx * L);
    }
    for (; i < j - 1; i++) {
        acc += pf.g[i];
        const int ci = (int)((LDG(csw + (i >> 2)) >> (8 * (i & 3))) & 0xFF);
        Big<L> cw, w, lw, q, t; ldb(cw, T + P.o_c10 + ((size_t)ci * SS + acc) * L);
        fmul(w, cw, ia);
        const int N = P.nterm[acc];
        flog1p(lw, w, N, T + P.o_inv); badd(lam, lam, lw);
        finv1p(q, w, N); fmul(t, ia, q); ia = t;
    }
    const int cl = (int)((LDG(csw + ((j - 1) >> 2)) >> (8 * ((j - 1) & 3))) & 0xFF);
    for (int gl = pf.glo; gl <= pf.ghi; gl++) {
        const int e = acc + gl;
        Big<L> la;
        if (e < P.Smin) ldb(la, T + P.o_sl + (size_t)(sidx | ((unsigned)(cl + 1) << (3 * (e - 1)))) * L);
        else case_lambda<L>(la, lam, ia, cl, e, C);
        u64 Pt[3];
#pragma unroll
        for (int a = 0; a < 3; a++) Pt[a] = frac_theta<L>(la, T + P.o_th[a]);
        if (P.dump) {
            const unsigned long long r = bump(nout);
            if ((long long)r < outcap) { out[r].prefix = idx; out[r].gl = (unsigned)gl; out[r].idx = 0; out[r].p[0] = Pt[0]; out[r].p[1] = Pt[1]; out[r].p[2] = Pt[2]; }
            continue;
        }
        probe(C, Pt, sigT, idx, (unsigned)gl, out, outcap, nout, cnt);
    }
}

/* ---------------- size classes: L = LF + 1 limbs (topgpu's) ---------------- */
#define MG_CLASSES(X) X(0, 6) X(1, 8) X(2, 10) X(3, 12) X(4, 14) X(5, 16) X(6, 20) X(7, 22) X(8, 24) X(9, 28) X(10, 32)

#ifdef MG_EMU
static mg_params HP; static const u32 *HT; static const mg_entry *HE; static const mg_slot *HH;
extern "C" int mg_init(const mg_params *p, const unsigned int *tables, const mg_entry *E, long long nE, const mg_slot *H) {
    (void)nE; HP = *p; HT = tables; HE = E; HH = H; return 0;
}
extern "C" long long mg_run(const mg_prefix *pre, int n, mg_out *out, long long outcap, unsigned long long counts[2]) {
    unsigned long long nout = 0, cnt[2] = {0, 0}; MCtx C = {&HP, HT, HE, HH};
    for (int i = 0; i < n; i++) {
        switch (HP.cls) {
#define CASE(k, l) case k: do_prefix<l>(C, pre[i], (unsigned)i, out, outcap, &nout, cnt); break;
            MG_CLASSES(CASE)
#undef CASE
        }
    }
    if ((long long)nout > outcap) return -1;
    counts[0] += cnt[0]; counts[1] += cnt[1];
    return (long long)nout;
}
extern "C" const char *mg_device(void) { return "emulation"; }
#else
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(5); } } while (0)
__constant__ mg_params DP;
template <int L>
__global__ void k_prefixes(const u32 *T, const mg_entry *E, const mg_slot *H, const mg_prefix *pre, int n, mg_out *out,
                           long long outcap, unsigned long long *nout, unsigned long long *gcnt) {
    unsigned long long cnt[2] = {0, 0}; MCtx C = {&DP, T, E, H};
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x)
        do_prefix<L>(C, pre[i], (unsigned)i, out, outcap, nout, cnt);
    if (cnt[0]) atomicAdd(&gcnt[0], cnt[0]);
    if (cnt[1]) atomicAdd(&gcnt[1], cnt[1]);
}
static mg_params HP; static u32 *DT; static mg_entry *DE; static mg_slot *DH; static mg_prefix *DPF; static mg_out *DOUT;
static unsigned long long *DCNT; static int DPFCAP; static long long DOUTCAP; static char DEVNAME[256]; static int NSM;
extern "C" int mg_init(const mg_params *p, const unsigned int *tables, const mg_entry *E, long long nE, const mg_slot *H) {
    HP = *p;
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr, 0)); NSM = pr.multiProcessorCount; snprintf(DEVNAME, sizeof DEVNAME, "%s", pr.name);
    CK(cudaMemcpyToSymbol(DP, p, sizeof *p));
    CK(cudaMalloc(&DT, (size_t)p->words * 4)); CK(cudaMemcpy(DT, tables, (size_t)p->words * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DE, (size_t)(nE ? nE : 1) * sizeof(mg_entry)));
    if (nE) CK(cudaMemcpy(DE, E, (size_t)nE * sizeof(mg_entry), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DH, ((size_t)1 << p->hbits) * sizeof(mg_slot)));
    CK(cudaMemcpy(DH, H, ((size_t)1 << p->hbits) * sizeof(mg_slot), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DCNT, 3 * 8));
    CK(cudaDeviceSynchronize());                         /* pageable copies: complete before any launch reads them */
    return 0;
}
extern "C" long long mg_run(const mg_prefix *pre, int n, mg_out *out, long long outcap, unsigned long long counts[2]) {
    if (n > DPFCAP) { if (DPF) cudaFree(DPF); DPFCAP = n; CK(cudaMalloc(&DPF, (size_t)n * sizeof(mg_prefix))); }
    if (outcap > DOUTCAP) { if (DOUT) cudaFree(DOUT); DOUTCAP = outcap; CK(cudaMalloc(&DOUT, (size_t)outcap * sizeof(mg_out))); }
    CK(cudaMemcpy(DPF, pre, (size_t)n * sizeof(mg_prefix), cudaMemcpyHostToDevice));
    CK(cudaMemset(DCNT, 0, 3 * 8));
    const int nt = 128; int grid = (n + nt - 1) / nt; if (grid > 64 * NSM) grid = 64 * NSM;
    switch (HP.cls) {
#define CASE(k, l) case k: k_prefixes<l><<<grid, nt>>>(DT, DE, DH, DPF, n, DOUT, outcap, DCNT + 2, DCNT); break;
        MG_CLASSES(CASE)
#undef CASE
    }
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    unsigned long long h[3]; CK(cudaMemcpy(h, DCNT, 3 * 8, cudaMemcpyDeviceToHost));
    long long no = (long long)h[2];
    if (no > outcap) return -1;
    counts[0] += h[0]; counts[1] += h[1];
    if (no) CK(cudaMemcpy(out, DOUT, (size_t)no * sizeof(mg_out), cudaMemcpyDeviceToHost));
    return no;
}
extern "C" const char *mg_device(void) { return DEVNAME; }
#endif
