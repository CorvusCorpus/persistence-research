/* Kernels of the GPU tree sieve: one thread block per coset of L_K (see tsieve6.c for the method).
 *
 * Per coset, in phases separated by barriers:
 *   1. U = delta2^-1 T(h0) as a doubled bitmap of 2 5^M bits in shared memory. T(h0) is enumerated leaf by leaf: in G mode
 *      leaf l picks, at level i, the (l >> 2i & 3)-th of the 4 admissible lifts; in E mode l runs over all 5^M paths and
 *      the paths through a digit 0 are dropped.
 *   2. For a chunk of i0 values, each thread finds the i1 range of its i0 (exact projection), then thread 0 forms
 *      prefix sums, so that rows (i0, i1) can be handed out by index.
 *   3. Rows: exact i2 range, then the window [s, s + len) of U, 32 positions per word; each set bit is a hit, whose
 *      K-digit block (digits K+1..2K) is checked; points that pass go through gsieve's later stages inline: the
 *      leading-digit table (stage 3), then limb_check's digits 27..9NL (low limbs first, then the full product).
 *      Points close to 10^XMIN are handed to the host (TG_CORNER), which applies gsieve's region rule.
 * The same source compiles as C++ for the CPU with -DTG_EMU: blocks run one after another and every phase is a loop
 * over the thread ids, so the emulation follows the kernel's exact control flow (used to validate it without a GPU).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tsieve_gpu.h"
typedef unsigned long long u64; typedef long long i64; typedef unsigned int u32;

#define IMAX 1024
#define E9 1000000000ULL
#define LG3 0.47712125471966243730
#define LG7 0.84509804001425683071

typedef struct { const u32 *T2, *T3, *T7; const u64 *TOP; } tg_tabs;

#ifdef TG_EMU
#define DEV static inline
#define MULHI64(x, y) ((u64)(((unsigned __int128)(x) * (y)) >> 64))
#define CTZ32(w) __builtin_ctz(w)
#define ATOMIC_OR32(p, v) (*(p) |= (v))
#define ATOMIC_ADD64(p, v) ([&]() { u64 o_ = *(p); *(p) += (v); return o_; }())
static tg_params P;
static int NTE = 32;                                   /* emulated threads per block */
#define PHASE_BEGIN for (int tid = 0; tid < NTE; tid++) { const int nth = NTE; (void)nth;
#define PHASE_END }
#else
#define DEV __device__ __forceinline__
#define MULHI64(x, y) __umul64hi((x), (y))
#define CTZ32(w) (__ffs(w) - 1)
#define ATOMIC_OR32(p, v) atomicOr((p), (v))
#define ATOMIC_ADD64(p, v) atomicAdd((unsigned long long *)(p), (unsigned long long)(v))
__constant__ tg_params P;
#define PHASE_BEGIN { const int tid = threadIdx.x; const int nth = blockDim.x; (void)nth;
#define PHASE_END } __syncthreads();
#endif

DEV u64 mod_br(u64 x, u64 m, u64 br) { u64 q = MULHI64(x, br), r = x - q * m; while (r >= m) r -= m; return r; }   /* x mod m */
DEV u64 mods_br(i64 x, u64 m, u64 br) { return x >= 0 ? mod_br((u64)x, m, br) : (m - mod_br((u64)(-x), m, br)) % m; }
DEV i64 floordiv(i64 x, i64 c) {                       /* floor(x / c), c > 0, |x| < 2^52 (float estimate, exact fix) */
    i64 q = (i64)floor((double)x / (double)c);
    while (c * q > x) q--;
    while (c * (q + 1) <= x) q++;
    return q;
}
DEV int ok9(u64 v, const unsigned char *ok4) { return ok4[v % 10000] & ok4[(v / 10000) % 10000] & P.ok1[v / 100000000ULL]; }

/* gsieve's stages 3-4 for a point whose last 2K digits are admissible; returns its decision (0 = S4 and further) */
DEV int later(u64 a, u64 b, u64 c, tg_tabs tb, const unsigned char *ok4, u64 *l2c, u64 *l3c, u64 *l4c) {
    if (b * LG3 + c * LG7 < P.cornerlg) return TG_CORNER;
    (*l2c)++;
    u64 bin = (a * P.l2 + b * P.l3 + c * P.l7) >> (64 - P.tb);
    if ((tb.TOP[bin >> 6] >> (bin & 63)) & 1) return TG_R3;
    (*l3c)++;
    const int NL = P.nl;
    const u32 *B3 = tb.T3 + b * NL, *B7 = tb.T7 + c * NL, *A = tb.T2 + a * NL;
    u64 m[16], acc, cy;
    acc = (u64)B3[0] * B7[0]; cy = acc / E9; m[0] = acc - cy * E9;
    acc = (u64)B3[0] * B7[1] + (u64)B3[1] * B7[0] + cy; cy = acc / E9; m[1] = acc - cy * E9;
    acc = (u64)B3[0] * B7[2] + (u64)B3[1] * B7[1] + (u64)B3[2] * B7[0] + cy; cy = acc / E9; m[2] = acc - cy * E9;
    acc = (u64)B3[0] * B7[3] + (u64)B3[1] * B7[2] + (u64)B3[2] * B7[1] + (u64)B3[3] * B7[0] + cy; cy = acc / E9; m[3] = acc - cy * E9;
    u64 a0 = A[0], a1 = A[1], a2 = A[2], a3 = A[3];
    acc = a0 * m[0]; cy = acc / E9;
    acc = a0 * m[1] + a1 * m[0] + cy; cy = acc / E9;
    acc = a0 * m[2] + a1 * m[1] + a2 * m[0] + cy; cy = acc / E9; u64 l2 = acc - cy * E9;
    acc = a0 * m[3] + a1 * m[2] + a2 * m[1] + a3 * m[0] + cy; cy = acc / E9; u64 l3 = acc - cy * E9;
    if (!(P.ok1[l2 / 100000000ULL] & ok9(l3, ok4))) return TG_R4;              /* digits 27..36 (limb_check, first test) */
    cy = 0;
    for (int k = 0; k < NL; k++) {                                              /* M = 3^b 7^c mod 10^(9 NL) in full */
        u64 ac = cy;
        for (int u = 0; u <= k; u++) ac += (u64)B3[u] * B7[k - u];
        cy = ac / E9; m[k] = ac - cy * E9;
    }
    cy = 0;
    for (int k = 0; k < NL; k++) {                                              /* limbs of 2^a M; digits 27.. checked */
        u64 ac = cy;
        for (int u = 0; u <= k; u++) ac += (u64)A[u] * m[k - u];
        cy = ac / E9; u64 limb = ac - cy * E9;
        if (k == 2 && !P.ok1[limb / 100000000ULL]) return TG_R4;
        if (k >= 3 && !ok9(limb, ok4)) return TG_R4;
    }
    (*l4c)++;
    return TG_S4;
}

#ifdef TG_EMU
static void k_coset(const tg_coset *cos, int ncos, tg_tabs tb, int dump, tg_point *out, u64 outcap, u64 *nout, u64 *gcnt, char *smem) {
    const int bid0 = 0, bstep = 1;
#else
__global__ void k_coset(const tg_coset *cos, int ncos, tg_tabs tb, int dump, tg_point *out, u64 outcap, u64 *nout, u64 *gcnt) {
    extern __shared__ __align__(16) char smem[];
    const int bid0 = blockIdx.x, bstep = gridDim.x;
#endif
    const u64 Mm = P.mm, M1 = P.m1;
    const int nwords = (int)((2 * Mm) / 32 + 3);
    u32 *bm = (u32 *)smem;
    i64 *i1lo = (i64 *)(smem + ((nwords * 4 + 15) & ~15));
    u32 *rcnt = (u32 *)(i1lo + IMAX), *rstart = rcnt + IMAX;
    u64 *bc = (u64 *)(rstart + IMAX + 2);                  /* block counters: s1, rows, hits, s2, s3, s4 */
    unsigned char *ok4 = (unsigned char *)(bc + 8);
    PHASE_BEGIN
        for (int x = tid; x < 10000; x += nth) ok4[x] = P.ok4[x];
        for (int w = tid; w < nwords; w += nth) bm[w] = 0;
        if (tid < 8) bc[tid] = 0;
    PHASE_END
    for (int ci = bid0; ci < ncos; ci += bstep) {
        const tg_coset C = cos[ci];
        PHASE_BEGIN                                            /* 1. U */
            for (u64 l = tid; l < P.nleaf; l += nth) {
                u64 Q = C.h0, y = 0, ll = l; int ok = 1;
                for (int i = 0; i < P.m; i++) {
                    int t;
                    if (P.dm == 0) {
                        int want = (int)(ll & 3), seen = -1; ll >>= 2;
                        for (t = 0; t < 5; t++) { u64 d = (Q + (u64)t * P.p2k[i]) % 10; if (d != 0 && d != 5 && ++seen == want) break; }
                    } else {
                        t = (int)(ll % 5); ll /= 5;
                        if ((Q + (u64)t * P.p2k[i]) % 10 == 0) { ok = 0; break; }
                    }
                    Q = (Q + (u64)t * P.p2k[i]) / 10; y += (u64)t * P.p5[i];
                }
                if (!ok) continue;
                u64 u = mod_br(y * (u64)C.inv2, Mm, P.brm), u2 = u + Mm;
                ATOMIC_OR32(&bm[u >> 5], 1u << (u & 31)); ATOMIC_OR32(&bm[u2 >> 5], 1u << (u2 & 31));
            }
        PHASE_END
        for (i64 i0c = C.i0lo; i0c <= C.i0hi; i0c += IMAX) {
            const int nk = (int)(C.i0hi - i0c + 1 < IMAX ? C.i0hi - i0c + 1 : IMAX);
            PHASE_BEGIN                                        /* 2. i1 ranges */
                for (int k = tid; k < nk; k += nth) {
                    i64 i0 = i0c + k, lo = -(1LL << 62), hi = 1LL << 62; int empty = 0;
                    for (int t = 0; t < P.n2; t++) {
                        i64 rhs = P.b2_0[t] + P.b2_r[t] * (i64)C.r - P.a2[t][0] * i0, cc = P.a2[t][1];
                        if (cc > 0) { i64 v = floordiv(rhs, cc); if (v < hi) hi = v; }
                        else if (cc < 0) { i64 v = -floordiv(rhs, -cc); if (v > lo) lo = v; }
                        else if (rhs < 0) empty = 1;
                    }
                    rcnt[k] = (empty || lo > hi) ? 0 : (u32)(hi - lo + 1); i1lo[k] = lo;
                }
            PHASE_END
            PHASE_BEGIN
                if (tid == 0) { rstart[0] = 0; for (int k = 0; k < nk; k++) rstart[k + 1] = rstart[k] + rcnt[k]; }
            PHASE_END
            PHASE_BEGIN                                        /* 3. rows */
                u64 ls1 = 0, lrows = 0, lhits = 0, l2c = 0, l3c = 0, l4c = 0;
                const u32 total = rstart[nk];
                for (u32 row = tid; row < total; row += nth) {
                    int lo = 0, hi = nk;
                    while (hi - lo > 1) { int mid = (lo + hi) >> 1; if (rstart[mid] <= row) lo = mid; else hi = mid; }
                    const i64 i0 = i0c + lo, i1 = i1lo[lo] + (i64)(row - rstart[lo]);
                    i64 lo2 = -(1LL << 62), hi2 = 1LL << 62; int empty = 0;
                    for (int q = 0; q < P.n3; q++) {
                        i64 rhs = P.b3_0[q] + P.b3_r[q] * (i64)C.r - P.a3[q][0] * i0 - P.a3[q][1] * i1, cc = P.a3[q][2];
                        if (cc > 0) { i64 v = floordiv(rhs, cc); if (v < hi2) hi2 = v; }
                        else if (cc < 0) { i64 v = -floordiv(rhs, -cc); if (v > lo2) lo2 = v; }
                        else if (rhs < 0) empty = 1;
                    }
                    if (empty || lo2 > hi2) continue;
                    const u64 len = (u64)(hi2 - lo2 + 1);
                    ls1 += len; lrows++;
                    u64 Ys = C.ybase + mod_br(mods_br(i0, M1, P.br1) * C.d0, M1, P.br1) + mod_br(mods_br(i1, M1, P.br1) * C.d1, M1, P.br1)
                           + mod_br(mods_br(lo2, M1, P.br1) * C.d2, M1, P.br1);
                    Ys = mod_br(Ys, M1, P.br1);
                    u64 s = C.s0 + mod_br(mods_br(i0, Mm, P.brm) * C.si0, Mm, P.brm) + mod_br(mods_br(i1, Mm, P.brm) * C.si1, Mm, P.brm)
                          + mods_br(lo2, Mm, P.brm);
                    s = mod_br(s, Mm, P.brm);
                    for (u64 base = 0; base < len; base += Mm) {
                        const u64 wlen = len - base < Mm ? len - base : Mm;
                        for (u64 off = 0; off < wlen; off += 32) {
                            const u64 pos = s + off, q32 = pos >> 5, sh = pos & 31;
                            u32 w = sh ? (bm[q32] >> sh) | (bm[q32 + 1] << (32 - sh)) : bm[q32];
                            const u64 rem = wlen - off;
                            if (rem < 32) w &= (1u << rem) - 1;
                            while (w) {
                                const u64 n = base + off + (u64)CTZ32(w); w &= w - 1;
                                lhits++;
                                const u64 Y = mod_br(Ys + n * C.d2, M1, P.br1);
                                u64 v = C.h0 + (Y << P.k);
                                for (int cc = 0; cc < P.skip; cc++) v /= 10000;
                                int ok = 1;
                                for (int cc = P.skip; cc < P.nfull; cc++) { ok &= ok4[v % 10000]; v /= 10000; }
                                if (P.nrem) ok &= P.okr[v];
                                if (!ok) continue;
                                const i64 i2 = lo2 + (i64)n;
                                const u64 a = (u64)((i64)C.r + i0 * P.v[0][0] + i1 * P.v[1][0] + i2 * P.v[2][0]);
                                const u64 b = (u64)(i0 * P.v[0][1] + i1 * P.v[1][1] + i2 * P.v[2][1]);
                                const u64 c = (u64)(i0 * P.v[0][2] + i1 * P.v[1][2] + i2 * P.v[2][2]);
                                const int code = later(a, b, c, tb, ok4, &l2c, &l3c, &l4c);
                                if (dump || code == TG_S4 || code == TG_CORNER) {
                                    u64 idx = ATOMIC_ADD64(nout, 1);
                                    if (idx < outcap) { out[idx].a = (u32)a; out[idx].b = (u32)b; out[idx].c = (u32)c; out[idx].code = (u32)code; }
                                }
                            }
                        }
                    }
                }
                ATOMIC_ADD64(&bc[0], ls1); ATOMIC_ADD64(&bc[1], lrows); ATOMIC_ADD64(&bc[2], lhits);
                ATOMIC_ADD64(&bc[3], l2c); ATOMIC_ADD64(&bc[4], l3c); ATOMIC_ADD64(&bc[5], l4c);
            PHASE_END
        }
        PHASE_BEGIN                                            /* clear U for the next coset */
            for (int w = tid; w < nwords; w += nth) bm[w] = 0;
        PHASE_END
    }
    PHASE_BEGIN
        if (tid < 6) ATOMIC_ADD64(&gcnt[tid], bc[tid]);
    PHASE_END
}

static size_t smem_bytes(const tg_params *p) {
    size_t nwords = (size_t)((2 * p->mm) / 32 + 3);
    return ((nwords * 4 + 15) & ~(size_t)15) + IMAX * 8 + (2 * IMAX + 2) * 4 + 8 * 8 + 10000;
}

#ifdef TG_EMU
/* ---------------- CPU emulation ---------------- */
struct tg_ctx { tg_point *out; u64 outcap; char *smem; };
static tg_tabs TABS;
extern "C" int tg_init(const tg_params *p, const u32 *T2, const u32 *T3, const u32 *T7, u64 n2, u64 n3, u64 n7, const u64 *TOP, u64 topwords) {
    (void)n2; (void)n3; (void)n7; (void)topwords;
    P = *p; TABS.T2 = T2; TABS.T3 = T3; TABS.T7 = T7; TABS.TOP = TOP;
    if (getenv("TG_EMU_THREADS")) NTE = atoi(getenv("TG_EMU_THREADS"));
    return 0;
}
extern "C" tg_ctx *tg_ctx_new(int maxcos, u64 outcap) {
    (void)maxcos;
    tg_ctx *x = (tg_ctx *)calloc(1, sizeof *x);
    x->out = (tg_point *)malloc(outcap * sizeof(tg_point)); x->outcap = outcap; x->smem = (char *)aligned_alloc(64, (smem_bytes(&P) + 63) & ~(size_t)63);
    return x;
}
extern "C" int tg_run(tg_ctx *x, const tg_coset *cos, int n, int dump, tg_point **out, u64 *nout, u64 cnt[6]) {
    u64 no = 0, gc[6] = {0, 0, 0, 0, 0, 0};
    k_coset(cos, n, TABS, dump, x->out, x->outcap, &no, gc, x->smem);
    for (int k = 0; k < 6; k++) cnt[k] += gc[k];
    *out = x->out; *nout = no < x->outcap ? no : x->outcap;
    return no > x->outcap;
}
extern "C" const char *tg_device(void) { return "emulation"; }
#else
/* ---------------- CUDA ---------------- */
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(5); } } while (0)
struct tg_ctx { cudaStream_t st; tg_coset *dcos, *hcos; tg_point *dout, *hout; u64 *dno, *dcnt, *hno; u64 outcap; int maxcos; };
static tg_tabs DTABS; static tg_params HP; static int NSM; static char DEVNAME[256];
extern "C" int tg_init(const tg_params *p, const u32 *T2, const u32 *T3, const u32 *T7, u64 n2, u64 n3, u64 n7, const u64 *TOP, u64 topwords) {
    HP = *p;
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr, 0)); NSM = pr.multiProcessorCount; snprintf(DEVNAME, sizeof DEVNAME, "%s", pr.name);
    CK(cudaMemcpyToSymbol(P, p, sizeof *p));
    u32 *d2, *d3, *d7; u64 *dtop;
    CK(cudaMalloc(&d2, n2 * 4)); CK(cudaMalloc(&d3, n3 * 4)); CK(cudaMalloc(&d7, n7 * 4)); CK(cudaMalloc(&dtop, topwords * 8));
    CK(cudaMemcpy(d2, T2, n2 * 4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d3, T3, n3 * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d7, T7, n7 * 4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(dtop, TOP, topwords * 8, cudaMemcpyHostToDevice));
    DTABS.T2 = d2; DTABS.T3 = d3; DTABS.T7 = d7; DTABS.TOP = dtop;
    size_t sm = smem_bytes(p);
    if (sm > 48 * 1024) CK(cudaFuncSetAttribute(k_coset, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)sm));
    return 0;
}
extern "C" tg_ctx *tg_ctx_new(int maxcos, u64 outcap) {
    tg_ctx *x = (tg_ctx *)calloc(1, sizeof *x);
    x->maxcos = maxcos; x->outcap = outcap;
    CK(cudaStreamCreateWithFlags(&x->st, cudaStreamNonBlocking));
    CK(cudaMalloc(&x->dcos, maxcos * sizeof(tg_coset))); CK(cudaMallocHost(&x->hcos, maxcos * sizeof(tg_coset)));
    CK(cudaMalloc(&x->dout, outcap * sizeof(tg_point))); CK(cudaMallocHost(&x->hout, outcap * sizeof(tg_point)));
    CK(cudaMalloc(&x->dno, 8 * 8)); CK(cudaMalloc(&x->dcnt, 8 * 8)); CK(cudaMallocHost(&x->hno, 8 * 8));
    return x;
}
extern "C" int tg_run(tg_ctx *x, const tg_coset *cos, int n, int dump, tg_point **out, u64 *nout, u64 cnt[6]) {
    if (n > x->maxcos) { fprintf(stderr, "tg_run: %d cosets > %d\n", n, x->maxcos); exit(5); }
    memcpy(x->hcos, cos, n * sizeof(tg_coset));
    CK(cudaMemcpyAsync(x->dcos, x->hcos, n * sizeof(tg_coset), cudaMemcpyHostToDevice, x->st));
    CK(cudaMemsetAsync(x->dno, 0, 8 * 8, x->st)); CK(cudaMemsetAsync(x->dcnt, 0, 8 * 8, x->st));
    int grid = n < 8 * NSM ? n : 8 * NSM;
    if (grid > 0) k_coset<<<grid, HP.nt, smem_bytes(&HP), x->st>>>(x->dcos, n, DTABS, dump, x->dout, x->outcap, x->dno, x->dcnt);
    CK(cudaGetLastError());
    CK(cudaMemcpyAsync(x->hno, x->dno, 8, cudaMemcpyDeviceToHost, x->st));
    CK(cudaMemcpyAsync(x->hno + 1, x->dcnt, 6 * 8, cudaMemcpyDeviceToHost, x->st));
    CK(cudaStreamSynchronize(x->st));                      /* pinned buffers, synchronized before reading (CLAUDE.md) */
    u64 no = x->hno[0];
    for (int k = 0; k < 6; k++) cnt[k] += x->hno[1 + k];
    u64 m = no < x->outcap ? no : x->outcap;
    if (m) { CK(cudaMemcpyAsync(x->hout, x->dout, m * sizeof(tg_point), cudaMemcpyDeviceToHost, x->st)); CK(cudaStreamSynchronize(x->st)); }
    *out = x->hout; *nout = m;
    return no > x->outcap;
}
extern "C" const char *tg_device(void) { return DEVNAME; }
#endif
