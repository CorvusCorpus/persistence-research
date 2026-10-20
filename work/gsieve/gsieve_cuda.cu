// gsieve_cuda.cu -- stages 1-2 of the chain sieve's main path on an NVIDIA GPU (the S2_CUDA variant of gsieve.c).
//
// Stage 1 (last K digits, lemma L1): the stage-1 survivors of a segment [s, e) of exponent residues are the set
// bits of the bitmap BM in that range, i.e. the dense index range [r0, r1) of the compacted tables AB, H0, POS,
// with r0 = rank(s), r1 = rank(e) from the prefix counts RK and a popcount -- as in the CPU variants.
// Stage 2 (digits K+1..2K, lemma L3): for table index i, AB[i] = A_hi | betaM << 32 with betaM = beta 2^32 mod 5^K
// (Montgomery form) and H0[i] = h0. Y = (A_hi + beta j) mod 5^K by one Montgomery multiply, and the K-digit block
// xh = h0 + 2^K Y is checked 4 digits at a time against the allowed-window table. The integer arithmetic is that of
// the S2_AVX512 / S2_AVX2 paths: floor(v / 10^4) = (v * M4) >> 45, M4 = ceil(2^45 / 10^4), exact for v < 2^31
// (proof in docs/CHAIN_SIEVE.md, exhaustive check in val_kernel.c), applied to the split Y = 10^4 sq + r,
// w = h0 + 2^K r, xh = 2^K 10^4 sq + w, so that every divided value is < 2^31. No floating point.
//
// One thread block per segment. Each tile of TILE consecutive indices is decided, and its survivors' bit positions
// POS[i] are appended in increasing i to a shared-memory buffer (warp ballots + one warp-level scan per tile), which
// is written to the launch's output buffer when the segment ends (or earlier, as an extra chunk, if it fills up).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cuda_runtime.h>
#include "gconst.h"
#include "gsieve_cuda.h"

#ifndef KK
#define KK 13
#endif
#if KK == 13
#define QQ 1220703125u                 // 5^13
#define QINV 163364931u                // -5^-13 mod 2^32 (as in gsieve.c)
static_assert(QQ == K13_Q, "Q");
#elif KK == 12
#define QQ 244140625u                  // 5^12
#define QINV 816824655u
static_assert(QQ == K12_Q, "Q");
#else
#error KK must be 12 or 13
#endif
#define M4 3518437209ull               // ceil(2^45 / 10^4)
static_assert((unsigned)(QQ * QINV) == 0xFFFFFFFFu, "QINV = -Q^-1 mod 2^32");
static_assert(M4 * 10000ull - (1ull << 45) == 1168ull, "M4");
static_assert(2ull * QQ < (1ull << 32) && QQ < (1u << 31), "Q < 2^31");

#define BLK 256                        // threads per block (one block per segment)
#define NWARP (BLK / 32)
#define ITEMS 4                        // indices per thread per tile
#define TILE (BLK * ITEMS)
#define SCAP (3 * TILE)                // shared survivor buffer; written out early only past 2 TILE survivors
#define OKW ((10000 + 31) / 32)        // allowed 4-digit windows as a bitmask
static_assert(ITEMS * NWARP == 32, "the per-tile scan uses one warp");

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error: %s (%s:%d: %s)\n", cudaGetErrorString(e_), __FILE__, __LINE__, #call); exit(1); } } while (0)

typedef unsigned long long ull;
struct tabs_t { const ull *BM; const unsigned *RK; const ull *AB; const unsigned short *H0; const unsigned *POS;
                const unsigned *OKM; int top; };
static tabs_t T;                       // device tables (read only after gpu_init)

__device__ __forceinline__ unsigned div4(unsigned v) { return (unsigned)(((ull)v * M4) >> 45); }    // v < 2^31
__device__ __forceinline__ unsigned okw(const unsigned *okm, unsigned d) { return (okm[d >> 5] >> (d & 31)) & 1u; }

// stage-2 decision for table entry (ab, h0) at lift index j: 1 if digits K+1..2K are all allowed
__device__ __forceinline__ unsigned stage2(ull ab, unsigned h0, unsigned j, const unsigned *okm, int top) {
    unsigned ahi = (unsigned)ab, bm = (unsigned)(ab >> 32);
    ull t = (ull)bm * j;                                     // < 2^62
    unsigned mm = (unsigned)t * QINV;                        // t * (-Q^-1) mod 2^32
    unsigned u = (unsigned)((t + (ull)mm * QQ) >> 32);       // = beta j mod Q, in [0, 2Q); t + mm Q < 2^64
    if (u >= QQ) u -= QQ;
    unsigned Y = u + ahi; if (Y >= QQ) Y -= QQ;               // (A_hi + beta j) mod Q < 2^31
    unsigned sq = div4(Y);
    unsigned w = h0 + ((Y - 10000u * sq) << KK);              // < 2^K 10^4 < 2^27; xh = 2^K 10^4 sq + w
    unsigned t1 = div4(w);
    if (!okw(okm, w - 10000u * t1)) return 0;                 // digits K+1..K+4 (xh mod 10^4 = w mod 10^4)
    unsigned q1 = (sq << KK) + t1;                            // floor(xh / 10^4) < 10^9
    unsigned q2 = div4(q1);
    if (!okw(okm, q1 - 10000u * q2)) return 0;                // digits K+5..K+8
#if KK == 13
    unsigned q3 = div4(q2);
    if (!okw(okm, q2 - 10000u * q3)) return 0;                // digits 22..25
    return top == 2 || (q3 != 0 && (top != 0 || q3 != 5));    // digit 26
#else
    return okw(okm, q2);                                      // digits 21..24
#endif
}

__global__ void __launch_bounds__(BLK) s12_kernel(tabs_t t, const gseg_t *__restrict__ seg, unsigned *__restrict__ s1o,
        unsigned *__restrict__ offo, unsigned *__restrict__ cnto, unsigned *__restrict__ pos, ull poscap,
        ull *__restrict__ ctr, govf_t *__restrict__ ovf, ull ovfcap) {
    __shared__ unsigned okm[OKW], sbuf[SCAP], wc[32];
    __shared__ unsigned sh_r0, sh_r1, sh_tot;
    __shared__ ull sh_off;
    const unsigned tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, k = blockIdx.x;
    for (unsigned v = tid; v < OKW; v += BLK) okm[v] = t.OKM[v];
    const gseg_t sg = seg[k];
    if (tid == 0) {                                           // stage 1: rank range of the set bits in [s, e)
        unsigned w0 = sg.s >> 6, w1 = (sg.e - 1) >> 6;
        unsigned r0 = t.RK[w0] + __popcll(t.BM[w0] & ((1ull << (sg.s & 63)) - 1));
        unsigned r1 = (sg.e & 63) ? t.RK[w1] + __popcll(t.BM[w1] & ((1ull << (sg.e & 63)) - 1)) : t.RK[w1 + 1];
        if (r1 < r0) r1 = r0;
        sh_r0 = r0; sh_r1 = r1; s1o[k] = r1 - r0;
    }
    __syncthreads();
    const unsigned r0 = sh_r0, r1 = sh_r1, j = sg.j, below = (1u << lane) - 1;
    unsigned nbuf = 0, seq = 0;                               // block-uniform
    // write sbuf[0, n) to the output: the segment's final chunk, or an early chunk recorded in ovf
#define FLUSH(n, fin) do {                                                                                         \
        if (tid == 0) sh_off = (n) ? atomicAdd(&ctr[0], (ull)(n)) : 0;                                           \
        __syncthreads();                                                                                         \
        const ull off_ = sh_off;                                                                                 \
        if (off_ + (n) <= poscap) for (unsigned v = tid; v < (n); v += BLK) pos[off_ + v] = sbuf[v];             \
        if (tid == 0) {                                                                                          \
            if (fin) { offo[k] = (unsigned)off_; cnto[k] = (n); }                                                \
            else { ull q_ = atomicAdd(&ctr[1], 1ull); if (q_ < ovfcap) ovf[q_] = govf_t{ k, seq, (unsigned)off_, (n) }; } \
        }                                                                                                        \
        __syncthreads();                                                                                         \
    } while (0)
    for (unsigned base = r0; base < r1; base += TILE) {
        if (nbuf + TILE > SCAP) { FLUSH(nbuf, 0); nbuf = 0; seq++; }
        ull ab[ITEMS]; unsigned h0[ITEMS], bal[ITEMS];
#pragma unroll
        for (int q = 0; q < ITEMS; q++) {
            unsigned i = base + q * BLK + tid;
            ab[q] = i < r1 ? t.AB[i] : 0; h0[q] = i < r1 ? t.H0[i] : 0;
        }
#pragma unroll
        for (int q = 0; q < ITEMS; q++) {
            unsigned i = base + q * BLK + tid;
            bal[q] = __ballot_sync(0xFFFFFFFFu, i < r1 && stage2(ab[q], h0[q], j, okm, t.top));
            if (lane == 0) wc[q * NWARP + warp] = __popc(bal[q]);
        }
        __syncthreads();
        if (warp == 0) {                                      // exclusive scan of the 32 counts in (q, warp) order
            unsigned v = wc[lane], x = v;
#pragma unroll
            for (int d = 1; d < 32; d <<= 1) { unsigned y = __shfl_up_sync(0xFFFFFFFFu, x, d); if (lane >= (unsigned)d) x += y; }
            wc[lane] = x - v;
            if (lane == 31) sh_tot = x;
        }
        __syncthreads();
#pragma unroll
        for (int q = 0; q < ITEMS; q++)                       // index base + q BLK + tid has rank (q, warp, lane)
            if ((bal[q] >> lane) & 1) sbuf[nbuf + wc[q * NWARP + warp] + __popc(bal[q] & below)] = t.POS[base + q * BLK + tid];
        nbuf += sh_tot;
        __syncthreads();                                      // wc, sh_tot are reused by the next tile
    }
    FLUSH(nbuf, 1);
#undef FLUSH
}

// ---------------------------------------------------------------- host side
// Every device table has a host copy that stays allocated, so that the device copy can be compared with it again
// at any time (gpu_verify_tables): GeForce memory has no ECC, and the tables are never written after the upload,
// so a memory fault that changes them persists until such a check finds it.
struct tbl_t { const void *h; const void *d; size_t bytes; const char *name; };
static tbl_t TBL[8]; static int NTBL;
static char *VSTAGE; static cudaStream_t VSTREAM;            // pinned staging buffer and stream of the checks
#define VCH ((size_t)64 << 20)
static int cmp_device(const tbl_t &t) {                    // 0 if the device copy equals the host copy
    for (size_t o = 0; o < t.bytes; o += VCH) {
        size_t m = t.bytes - o < VCH ? t.bytes - o : VCH;
        CK(cudaMemcpyAsync(VSTAGE, (const char *)t.d + o, m, cudaMemcpyDeviceToHost, VSTREAM));
        CK(cudaStreamSynchronize(VSTREAM));
        if (memcmp(VSTAGE, (const char *)t.h + o, m)) return 1;
    }
    return 0;
}
template <class X> static X *dupload(const X *h, size_t n, const char *name, size_t *total) {
    X *d; size_t bytes = n * sizeof(X);
    CK(cudaMalloc(&d, bytes));
    CK(cudaMemcpy(d, h, bytes, cudaMemcpyHostToDevice));
    CK(cudaDeviceSynchronize());     // a pageable upload may return before its DMA lands; the checks use another stream
    TBL[NTBL] = tbl_t{ h, d, bytes, name };
    if (cmp_device(TBL[NTBL])) { fprintf(stderr, "gpu: device copy of %s differs from the host table\n", name); exit(1); }
    NTBL++; *total += bytes;
    return d;
}
// compare every device table with its host copy; returns 0 if all are equal, else prints which one differs
extern "C" int gpu_verify_tables(void) {
    for (int i = 0; i < NTBL; i++) if (cmp_device(TBL[i])) { fprintf(stderr, "gpu: device table %s no longer equals the host copy\n", TBL[i].name); return 1; }
    return 0;
}

extern "C" void gpu_init(const uint64_t *BM, const uint32_t *RK, uint64_t nw, const uint64_t *AB, const uint16_t *H0,
                         const uint32_t *POS, uint64_t ngood, const uint8_t *ok4b, int topmode) {
    // workers waiting for their launch block instead of spinning: the CPU time goes to stages 3+
    CK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    cudaDeviceProp p; CK(cudaSetDevice(0)); CK(cudaGetDeviceProperties(&p, 0));
    CK(cudaHostAlloc((void **)&VSTAGE, VCH, cudaHostAllocDefault));
    CK(cudaStreamCreateWithFlags(&VSTREAM, cudaStreamNonBlocking));
    static unsigned okm[OKW];                                  // static: kept as the host copy of OKM
    memset(okm, 0, sizeof okm);
    for (unsigned v = 0; v < 10000; v++) if (ok4b[v]) okm[v >> 5] |= 1u << (v & 31);
    size_t tot = 0;
    T.BM = dupload((const ull *)BM, nw + 1, "BM", &tot);
    T.RK = dupload((const unsigned *)RK, nw + 1, "RK", &tot);
    T.AB = dupload((const ull *)AB, ngood + 8, "AB", &tot);
    T.H0 = dupload((const unsigned short *)H0, ngood + 8, "H0", &tot);
    T.POS = dupload((const unsigned *)POS, ngood + 8, "POS", &tot);
    T.OKM = dupload(okm, OKW, "OKM", &tot);
    T.top = topmode;
    size_t fr, all; CK(cudaMemGetInfo(&fr, &all));
    fprintf(stderr, "gpu: %s (cc %d.%d, %d SMs), K=%d tables %.3f GB in device memory, verified; %.2f GB free\n",
            p.name, p.major, p.minor, p.multiProcessorCount, KK, tot / 1e9, fr / 1e9);
}

static void ctx_alloc_seg(gctx_t *g, uint64_t n) {
    CK(cudaHostAlloc((void **)&g->h_seg, n * sizeof(gseg_t), cudaHostAllocDefault));
    CK(cudaHostAlloc((void **)&g->h_s1, n * 4, cudaHostAllocDefault));
    CK(cudaHostAlloc((void **)&g->h_off, n * 4, cudaHostAllocDefault));
    CK(cudaHostAlloc((void **)&g->h_cnt, n * 4, cudaHostAllocDefault));
    CK(cudaMalloc((void **)&g->d_seg, n * sizeof(gseg_t)));
    CK(cudaMalloc((void **)&g->d_s1, n * 4)); CK(cudaMalloc((void **)&g->d_off, n * 4)); CK(cudaMalloc((void **)&g->d_cnt, n * 4));
    g->segcap = n;
}
static void ctx_alloc_pos(gctx_t *g, uint64_t n) {
    if (n >= (1ull << 32)) { fprintf(stderr, "gpu: survivor buffer would exceed 2^32 entries\n"); exit(1); }
    CK(cudaHostAlloc((void **)&g->h_pos, n * 4, cudaHostAllocDefault)); CK(cudaMalloc((void **)&g->d_pos, n * 4));
    g->poscap = n;
}
static void ctx_alloc_ovf(gctx_t *g, uint64_t n) {
    CK(cudaHostAlloc((void **)&g->h_ovf, n * sizeof(govf_t), cudaHostAllocDefault)); CK(cudaMalloc((void **)&g->d_ovf, n * sizeof(govf_t)));
    g->ovfcap = n;
}
extern "C" gctx_t *gpu_ctx_new(uint64_t segcap, uint64_t poscap) {
    gctx_t *g = (gctx_t *)calloc(1, sizeof *g);
    cudaStream_t s; CK(cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking)); g->stream = s;
    ctx_alloc_seg(g, segcap); ctx_alloc_pos(g, poscap); ctx_alloc_ovf(g, 1024);
    CK(cudaHostAlloc((void **)&g->h_ctr, 2 * sizeof(ull), cudaHostAllocDefault)); CK(cudaMalloc((void **)&g->d_ctr, 2 * sizeof(ull)));
    return g;
}
extern "C" void gpu_run(gctx_t *g, uint64_t nseg) {
    if (nseg == 0) { g->h_ctr[0] = g->h_ctr[1] = 0; return; }
    if (nseg > g->segcap || nseg >= (1ull << 31)) { fprintf(stderr, "gpu_run: %llu segments exceed the buffers\n", (ull)nseg); exit(1); }
    cudaStream_t s = (cudaStream_t)g->stream;
    CK(cudaMemcpyAsync(g->d_seg, g->h_seg, nseg * sizeof(gseg_t), cudaMemcpyHostToDevice, s));
    for (;;) {
        CK(cudaMemsetAsync(g->d_ctr, 0, 2 * sizeof(ull), s));
        s12_kernel<<<(unsigned)nseg, BLK, 0, s>>>(T, g->d_seg, g->d_s1, g->d_off, g->d_cnt, g->d_pos, g->poscap, g->d_ctr, g->d_ovf, g->ovfcap);
        CK(cudaGetLastError());
        CK(cudaMemcpyAsync(g->h_ctr, g->d_ctr, 2 * sizeof(ull), cudaMemcpyDeviceToHost, s));
        CK(cudaStreamSynchronize(s));
        ull npos = g->h_ctr[0], novf = g->h_ctr[1];
        if (npos <= g->poscap && novf <= g->ovfcap) break;
        // more survivors than the buffers hold: enlarge them and run the launch again (nothing was dropped silently)
        if (npos > g->poscap) {
            CK(cudaFreeHost(g->h_pos)); CK(cudaFree(g->d_pos));
            ctx_alloc_pos(g, npos + npos / 2 + 4096);
        }
        if (novf > g->ovfcap) {
            CK(cudaFreeHost(g->h_ovf)); CK(cudaFree(g->d_ovf));
            ctx_alloc_ovf(g, 2 * novf + 1024);
        }
    }
    CK(cudaMemcpyAsync(g->h_s1, g->d_s1, nseg * 4, cudaMemcpyDeviceToHost, s));
    CK(cudaMemcpyAsync(g->h_off, g->d_off, nseg * 4, cudaMemcpyDeviceToHost, s));
    CK(cudaMemcpyAsync(g->h_cnt, g->d_cnt, nseg * 4, cudaMemcpyDeviceToHost, s));
    if (g->h_ctr[0]) CK(cudaMemcpyAsync(g->h_pos, g->d_pos, g->h_ctr[0] * 4, cudaMemcpyDeviceToHost, s));
    if (g->h_ctr[1]) CK(cudaMemcpyAsync(g->h_ovf, g->d_ovf, g->h_ctr[1] * sizeof(govf_t), cudaMemcpyDeviceToHost, s));
    CK(cudaStreamSynchronize(s));
}
extern "C" void gpu_ctx_free(gctx_t *g) {
    CK(cudaStreamSynchronize((cudaStream_t)g->stream));
    CK(cudaFreeHost(g->h_seg)); CK(cudaFreeHost(g->h_s1)); CK(cudaFreeHost(g->h_off)); CK(cudaFreeHost(g->h_cnt));
    CK(cudaFreeHost(g->h_pos)); CK(cudaFreeHost(g->h_ovf)); CK(cudaFreeHost(g->h_ctr));
    CK(cudaFree(g->d_seg)); CK(cudaFree(g->d_s1)); CK(cudaFree(g->d_off)); CK(cudaFree(g->d_cnt));
    CK(cudaFree(g->d_pos)); CK(cudaFree(g->d_ovf)); CK(cudaFree(g->d_ctr));
    CK(cudaStreamDestroy((cudaStream_t)g->stream));
    free(g);
}
