// z2sieve_cuda.cu: stages 2-4 of z2sieve on an NVIDIA GPU (docs/ZEROLESS_POW2.md, "CUDA variant").
// Stage 1 is the host's bitmap over r; stage 5 stays on the host. One thread per surviving residue r: it derives the
// state of r from the run-start tables (one multiplication), then walks j over the chunk in lockstep with its block.
// Stage-2 survivors go to a shared queue that the whole block drains through stages 3-4 every JSTEP steps of j.
// Same decisions as the CPU code: the dumps are compared line by line (val_z2.sh).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cuda_runtime.h>
#include "z2sieve_cuda.h"

#ifndef KK
#define KK 13
#endif
#ifndef NL
#define NL 18
#endif
#ifndef TB
#define TB 24
#endif
#ifndef TD
#define TD 7
#endif
#ifndef JCH
#define JCH 4096
#endif
typedef unsigned long long ull; typedef unsigned __int128 u128;
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error: %s (%s:%d: %s)\n", cudaGetErrorString(e_), __FILE__, __LINE__, #call); exit(1); } } while (0)

__host__ __device__ constexpr ull p10(int k) { return k ? 10ull * p10(k - 1) : 1ull; }
__host__ __device__ constexpr ull f5(int k) { return k ? 5ull * f5(k - 1) : 1ull; }
constexpr ull TENK = p10(KK), M5 = f5(KK), M2 = f5(2 * KK), PP = 4 * f5(KK - 1), TWOK = 1ull << KK, E9 = 1000000000ull;
constexpr int WC = KK < 5 ? KK : 5, WB = KK - WC < 4 ? KK - WC : 4, WA = KK - WC - WB;
constexpr unsigned BC = (unsigned)p10(WC), BB = (unsigned)p10(WB), BA = (unsigned)p10(WA);
constexpr int LQ = 2 * KK / 9, LO = 2 * KK % 9;
// digit bit tables ZT[0..5] concatenated: words per width 1, 1, 2, 16, 157, 1563
__host__ __device__ constexpr int zof(int w) { return w == 0 ? 0 : w == 1 ? 1 : w == 2 ? 2 : w == 3 ? 4 : w == 4 ? 20 : 177; }
#define ZTOT 1740
#define BLK 128
#define JSB 3                           // JSTEP = 2^JSB steps of j between queue drains
#define JSTEP (1 << JSB)
#define QCAP (BLK * JSTEP)
#define LCH (1ull << 20)                // survivors per launch (a few ms each)
#define S4CAP (1u << 22)
#define DCAP (1u << 25)
static_assert(BLK << JSB <= 65536, "queue entries pack the thread and j - jb0 into 16 bits");

struct gdev_t {
    const ull *zt, *top, *b1; const unsigned *rk; ull nw;
    const ull *rsw, *rss; const unsigned *rsr; const ull *pw13, *ps; const unsigned *pw;
    const unsigned *q; ull *ctr; ull *s4; unsigned *ns4; ull *dmp; unsigned *ndmp; unsigned *flags;
};
struct prm_t { ull ja, jb, n0, n1, C2K, UU, LLO, LHI; int dump, probe; };
// Timing probes exist only in builds with -DZ2PROBE (never production: they skip stages): env Z2_PROBE = 1 stops after
// the setup, 2 after stage 2, 3 after stage 3.
#ifdef Z2PROBE
#define PROBE(k) (p.probe == (k))
#else
#define PROBE(k) 0
#endif

__device__ __forceinline__ ull mulmod(ull a, ull b, ull m) { return (ull)((u128)a * b % m); }
// 32-bit digit windows: 64-bit integer division is emulated on the GPU, 32-bit division by a constant is a multiply
__device__ __forceinline__ int zbit(const ull *zs, int w, unsigned v) { return (int)((zs[zof(w) + (v >> 6)] >> (v & 63)) & 1); }
__device__ __forceinline__ int zfw(const ull *zs, unsigned v, int w) {  // v written with w <= 9 digits has no 0
    int ok = 1;
    while (w > 5) { ok &= zbit(zs, 4, v % 10000u); v /= 10000u; w -= 4; }
    return ok & zbit(zs, w, v);
}
__device__ __forceinline__ int zf9(const ull *zs, unsigned v) { return zbit(zs, 4, v % 10000u) & zbit(zs, 5, v / 10000u); }
// acc mod 10^9, and q = acc div 10^9, for every acc < 2^64: q = (acc >> 9) div 5^9 by the multiplier m = ceil(2^76/5^9).
// With e = m 5^9 - 2^76 = 799614 the error of (acc >> 9) m / 2^76 is below 2^55 e / (5^9 2^76) < 1/5^9, so the floor is
// exact (also checked at 1.4 10^6 values with edge cases, NOTEBOOK P.5)
__device__ __forceinline__ unsigned divE9(ull acc, ull &q) {
    q = __umul64hi(acc >> 9, 38685626227668134ull) >> 12;
    return (unsigned)(acc - q * E9);
}
// a K-digit window as parts a|b|c
__device__ __forceinline__ void kd_split(ull x, unsigned &a, unsigned &b, unsigned &c) {
    c = (unsigned)(x % BC); x /= BC; b = (unsigned)(x % BB); a = (unsigned)(x / BB);
}
__device__ __forceinline__ int kd_ok(const ull *zs, unsigned a, unsigned b, unsigned c) {
    return zbit(zs, WA, a) & zbit(zs, WB, b) & zbit(zs, WC, c);
}
__device__ __forceinline__ void kd_add(unsigned &a, unsigned &b, unsigned &c, unsigned da, unsigned db, unsigned dc) {
    unsigned x = c + dc, cy = x >= BC; c = cy ? x - BC : x;
    unsigned y = b + db + cy; cy = y >= BB; b = cy ? y - BB : y;
    unsigned z = a + da + cy; a = z >= BA ? z - BA : z;
}
// stage 4: the CPU limb_engine on the first KEND limbs (A = 2^r limbs at stride BLK, B = a Q row, both in shared
// memory). With chk (n a multiple of 64, or dump mode, as on the CPU) the limbs below digit 2K+1 must equal
// w(r) + 10^K U(j) from stages 1-2. limb_cols<PRE> is an exact pre-filter (digits 2K+1..9 PRE, no early exit, so
// whole warps run it); limb_cols<NL> is the full engine.
constexpr int PRE = LQ + 2 < NL ? LQ + 2 : NL;
template <int KEND> __device__ __forceinline__ int limb_cols(const ull *zs, const unsigned *A, const unsigned *B, ull w, ull U, int chk, unsigned *flags) {
    ull cy = 0, e18 = 0, ehi = 0;
    if (chk) { e18 = w + TENK * (U % p10(18 - KK)); ehi = U / p10(18 - KK); }
#pragma unroll
    for (int k = 0; k < KEND; k++) {
        ull acc = cy;
#pragma unroll
        for (int u = 0; u <= k; u++) acc += (ull)A[u * BLK] * B[k - u];
        unsigned limb = divE9(acc, cy);
        if (k < LQ) {
            if (chk && limb != (k == 0 ? e18 % E9 : k == 1 ? e18 / E9 : ehi)) atomicOr(flags, 4u);
        } else if (k == LQ) {
            if (chk && limb % (unsigned)p10(LO) != (k == 0 ? e18 % E9 : k == 1 ? e18 / E9 : ehi)) atomicOr(flags, 4u);
            if (!zfw(zs, limb / (unsigned)p10(LO), 9 - LO)) return 0;
        } else if (!zf9(zs, limb)) return 0;
    }
    return 1;
}
__device__ __forceinline__ void dput(const gdev_t &D, ull n, unsigned cls) {
    unsigned pos = atomicAdd(D.ndmp, 1u);
    if (pos < DCAP) D.dmp[pos] = n << 2 | cls; else atomicOr(D.flags, 16u);
}

__global__ void __launch_bounds__(BLK, 3) ksieve(gdev_t D, prm_t p, ull i0, ull i1) {
    __shared__ ull zs[ZTOT];
    __shared__ unsigned short q[QCAP], q2[QCAP];               // queue entries: thread << JSB | (j - jb0)
    __shared__ unsigned qn, qn2, qn3;
    __shared__ unsigned sr[BLK], szhi[BLK], sd[BLK], sa0[BLK], sR[NL * BLK], sQ[JSTEP * NL];   // sR: limb-major, sR[u BLK + t]
    __shared__ ull sw[BLK], jmn, jmx, bc[4];
    for (int t = threadIdx.x; t < ZTOT; t += BLK) zs[t] = D.zt[t];
    if (threadIdx.x == 0) { qn = qn2 = qn3 = 0; jmn = ~0ull; jmx = 0; bc[0] = bc[1] = bc[2] = bc[3] = 0; }
    __syncthreads();
    ull i = i0 + (ull)blockIdx.x * BLK + threadIdx.x, jlo = 0, jhi = 0;
    int act = i < i1;
    unsigned ua = 0, ub = 0, uc = 0, da = 0, db = 0, dc = 0;
    ull r = 0;
    if (act) {
        ull lo = 0, hi = D.nw;                                  // the word holding survivor i: rk[lo] <= i < rk[lo+1]
        while (hi - lo > 1) { ull mid = (lo + hi) / 2; if (D.rk[mid] <= i) lo = mid; else hi = mid; }
        ull wv = D.b1[lo]; unsigned k = (unsigned)(i - D.rk[lo]); int bit = -1;
        for (int b = 0; b < 64; b++) if ((wv >> b) & 1) { if (!k) { bit = b; break; } k--; }
        if (bit < 0) { atomicOr(D.flags, 1u); act = 0; }
        r = lo * 64 + (unsigned)(bit < 0 ? 0 : bit);
    }
    if (act) {
        ull kr = r >> 10; unsigned g = (unsigned)(r & 1023);
        ull w = mulmod(D.rsw[kr], D.pw13[g], TENK), s = mulmod(D.rss[kr], D.ps[g], M2);
        const unsigned *A = D.rsr + kr * NL, *B = D.pw + (ull)g * NL;    // R = 2^r mod 10^(9 NL)
        ull cy = 0;
        for (int c = 0; c < NL; c++) { ull acc = cy; for (int u = 0; u <= c; u++) acc += (ull)A[u] * B[c - u]; sR[c * BLK + threadIdx.x] = divE9(acc, cy); }
        unsigned a, b2, c2; kd_split(w, a, b2, c2);
        if (!kd_ok(zs, a, b2, c2)) atomicOr(D.flags, 1u);      // the bitmap must agree with w(r)
        jlo = p.n0 > r ? (p.n0 - r + PP - 1) / PP : 0; jhi = p.n1 > r ? (p.n1 - r + PP - 1) / PP : 0;
        if (jlo < p.ja) jlo = p.ja;
        if (jhi > p.jb) jhi = p.jb;
        ull z0 = mulmod(s, p.C2K, M2), zlo = z0 % M5, zhi = z0 / M5, a0 = TWOK * zlo / M5;
        if (TWOK * TWOK * zlo % TENK != w) atomicOr(D.flags, 2u);
        ull d = mulmod(zlo, p.UU, M5), t0 = (zhi + mulmod(jlo % M5, d, M5)) % M5;
        kd_split(a0 + TWOK * t0, ua, ub, uc); kd_split(TWOK * d, da, db, dc);
        sr[threadIdx.x] = (unsigned)r; sw[threadIdx.x] = w; szhi[threadIdx.x] = (unsigned)zhi; sd[threadIdx.x] = (unsigned)d; sa0[threadIdx.x] = (unsigned)a0;
        if (jlo < jhi) { atomicMin(&jmn, jlo); atomicMax(&jmx, jhi); atomicAdd(&bc[0], jhi - jlo); }
        else act = 0;
    }
    __syncthreads();
    ull j0 = jmn, j1 = jmx, n3 = 0, n4 = 0;                    // this thread's stage-3 and stage-4 passes
    if (PROBE(1)) j1 = j0;
    for (ull jb0 = j0; jb0 < j1; jb0 += JSTEP) {
        for (int jj = 0; jj < JSTEP; jj++) {                   // stage 2, lockstep over the block
            ull j = jb0 + jj;
            int in = act && j >= jlo && j < jhi, pred = in && kd_ok(zs, ua, ub, uc);
            unsigned m = __ballot_sync(0xffffffffu, pred);
            if (m) {
                int lane = threadIdx.x & 31, ldr = __ffs(m) - 1; unsigned base = 0;
                if (lane == ldr) base = atomicAdd(&qn, (unsigned)__popc(m));
                base = __shfl_sync(0xffffffffu, base, ldr);
                if (pred) q[base + __popc(m & ((1u << lane) - 1))] = (unsigned short)(threadIdx.x << JSB | jj);
            }
            if (in) kd_add(ua, ub, uc, da, db, dc);
        }
        for (int t = threadIdx.x; t < JSTEP * NL; t += BLK) {  // the Q rows of this window, for stage 4
            ull row = jb0 - p.ja + t / NL;
            sQ[t] = row < JCH ? D.q[row * NL + t % NL] : 0;
        }
        __syncthreads();
        unsigned nq = qn;
        if (threadIdx.x == 0) bc[1] += nq;
        if (PROBE(2)) nq = 0;
        for (unsigned b0 = 0; b0 < nq; b0 += BLK) {            // stage 3; its survivors go to q2 (whole warps take part)
            unsigned e = b0 + threadIdx.x, it = e < nq ? q[e] : 0; int pass = 0;
            if (e < nq) {
                ull n = sr[it >> JSB] + (jb0 + (it & (JSTEP - 1))) * PP;
                ull fh = __umul64hi(n, p.LLO) + n * p.LHI, bin = fh >> (64 - TB);
                if ((__ldg(D.top + (bin >> 6)) >> (bin & 63)) & 1) { if (p.dump) dput(D, n, 1); } else { pass = 1; n3++; }
            }
            unsigned m = __ballot_sync(0xffffffffu, pass);
            if (m) {
                int lane = threadIdx.x & 31, ldr = __ffs(m) - 1; unsigned base = 0;
                if (lane == ldr) base = atomicAdd(&qn2, (unsigned)__popc(m));
                base = __shfl_sync(0xffffffffu, base, ldr);
                if (pass) q2[base + __popc(m & ((1u << lane) - 1))] = (unsigned short)it;
            }
        }
        __syncthreads();
        unsigned nq2 = PROBE(3) ? 0 : qn2;
        for (unsigned b0 = 0; b0 < nq2; b0 += BLK) {           // stage 4a: the first PRE limbs; survivors back to q
            unsigned e = b0 + threadIdx.x, it = e < nq2 ? q2[e] : 0; int pass = 0;
            if (e < nq2) {
                unsigned lane = it >> JSB, jj = it & (JSTEP - 1);
                ull j = jb0 + jj, n = sr[lane] + j * PP;
                int chk = p.dump || !(n & 63);
                ull U = chk ? sa0[lane] + TWOK * ((szhi[lane] + (j % M5) * sd[lane] % M5) % M5) : 0;
                if (limb_cols<PRE>(zs, sR + lane, sQ + jj * NL, sw[lane], U, chk, D.flags)) pass = 1;
                else if (p.dump) dput(D, n, 2);
            }
            unsigned m = __ballot_sync(0xffffffffu, pass);
            if (m) {
                int lane = threadIdx.x & 31, ldr = __ffs(m) - 1; unsigned base = 0;
                if (lane == ldr) base = atomicAdd(&qn3, (unsigned)__popc(m));
                base = __shfl_sync(0xffffffffu, base, ldr);
                if (pass) q[base + __popc(m & ((1u << lane) - 1))] = (unsigned short)it;    // q was drained by stage 3
            }
        }
        __syncthreads();
        unsigned nq3 = qn3;
        for (unsigned e = threadIdx.x; e < nq3; e += BLK) {    // stage 4b: all NL limbs, early exit
            unsigned it = q[e], lane = it >> JSB, jj = it & (JSTEP - 1);
            ull n = sr[lane] + (jb0 + jj) * PP;
            if (!limb_cols<NL>(zs, sR + lane, sQ + jj * NL, 0, 0, 0, D.flags)) { if (p.dump) dput(D, n, 2); continue; }
            n4++;
            unsigned pos = atomicAdd(D.ns4, 1u);
            if (pos < S4CAP) D.s4[pos] = n; else atomicOr(D.flags, 8u);
            if (p.dump) dput(D, n, 3);
        }
        __syncthreads();
        if (threadIdx.x == 0) qn = qn2 = qn3 = 0;
        __syncthreads();
    }
    for (int o = 16; o; o >>= 1) { n3 += __shfl_down_sync(0xffffffffu, n3, o); n4 += __shfl_down_sync(0xffffffffu, n4, o); }
    if ((threadIdx.x & 31) == 0) { atomicAdd(&bc[2], n3); atomicAdd(&bc[3], n4); }
    __syncthreads();
    if (threadIdx.x == 0) for (int k = 0; k < 4; k++) atomicAdd(D.ctr + k, bc[k]);
}

// ---------------------------------------------------------------- host side
// Every device table keeps its host copy, so that z2g_verify can compare them again (GeForce memory has no ECC).
struct tbl_t { const void *h; const void *d; size_t bytes; const char *name; };
static tbl_t TBL[16]; static int NTBL;
static char *VSTAGE; static cudaStream_t VSTREAM;
#define VCH ((size_t)64 << 20)
static int cmp_device(const tbl_t &t) {
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
    CK(cudaDeviceSynchronize());                     // a pageable upload may return before its DMA lands
    TBL[NTBL] = tbl_t{ h, d, bytes, name };
    if (cmp_device(TBL[NTBL])) { fprintf(stderr, "gpu: device copy of %s differs from the host table\n", name); exit(1); }
    NTBL++; *total += bytes;
    return d;
}
static gdev_t D; static z2c_t C;
static unsigned *QD; static const uint32_t *QH;      // Q of the current chunk (device, host)
static ull *HS4, *HDMP;
extern "C" int z2g_verify(void) {
    for (int i = 0; i < NTBL; i++) if (cmp_device(TBL[i])) { fprintf(stderr, "gpu: device table %s no longer equals the host copy\n", TBL[i].name); return 1; }
    if (QH && cmp_device(tbl_t{ QH, QD, (size_t)JCH * NL * 4, "Q" })) { fprintf(stderr, "gpu: device Q no longer equals the host copy\n"); return 1; }
    return 0;
}
extern "C" void z2g_init(const z2c_t *c, const uint64_t *zt, const uint64_t *top, uint64_t topw, const uint64_t *b1,
                         const uint32_t *rk, uint64_t nw, const uint64_t *rsw, const uint64_t *rss, const uint32_t *rsr,
                         uint64_t nrs, const uint64_t *pw13, const uint64_t *ps, const uint32_t *pw) {
    CK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    cudaDeviceProp pr; CK(cudaSetDevice(0)); CK(cudaGetDeviceProperties(&pr, 0));
    CK(cudaHostAlloc((void **)&VSTAGE, VCH, cudaHostAllocDefault));
    CK(cudaStreamCreateWithFlags(&VSTREAM, cudaStreamNonBlocking));
    C = *c; size_t tot = 0;
    D.zt = dupload((const ull *)zt, ZTOT, "ZT", &tot);
    D.top = dupload((const ull *)top, topw, "TOP", &tot);
    D.b1 = dupload((const ull *)b1, nw, "B1", &tot);
    D.rk = dupload((const unsigned *)rk, nw + 1, "RK", &tot); D.nw = nw;
    D.rsw = dupload((const ull *)rsw, nrs, "RSW", &tot);
    D.rss = dupload((const ull *)rss, nrs, "RSS", &tot);
    D.rsr = dupload((const unsigned *)rsr, nrs * NL, "RSR", &tot);
    D.pw13 = dupload((const ull *)pw13, 1024, "PW13", &tot);
    D.ps = dupload((const ull *)ps, 1024, "PS", &tot);
    D.pw = dupload((const unsigned *)pw, 1024 * NL, "PW", &tot);
    CK(cudaMalloc(&QD, (size_t)JCH * NL * 4)); D.q = QD;
    CK(cudaMalloc(&D.ctr, 4 * 8)); CK(cudaMalloc(&D.ns4, 4)); CK(cudaMalloc(&D.ndmp, 4)); CK(cudaMalloc(&D.flags, 4));
    CK(cudaMalloc(&D.s4, (size_t)S4CAP * 8)); HS4 = (ull *)malloc((size_t)S4CAP * 8);
    D.dmp = NULL; HDMP = NULL;
    size_t fr, all; CK(cudaMemGetInfo(&fr, &all));
    fprintf(stderr, "gpu: %s (cc %d.%d, %d SMs): tables %.3f GB in device memory, verified; %.2f GB free\n",
            pr.name, pr.major, pr.minor, pr.multiProcessorCount, tot / 1e9, fr / 1e9);
}
extern "C" void z2g_chunk(uint64_t ja, uint64_t jb, uint64_t n0, uint64_t n1, const uint32_t *q, const uint64_t *iv,
                          int niv, int dump, z2r_t *res) {
    if (dump && !D.dmp) { CK(cudaMalloc(&D.dmp, (size_t)DCAP * 8)); HDMP = (ull *)malloc((size_t)DCAP * 8); }
    QH = q;
    CK(cudaMemcpy(QD, q, (size_t)JCH * NL * 4, cudaMemcpyHostToDevice)); CK(cudaDeviceSynchronize());
    if (cmp_device(tbl_t{ QH, QD, (size_t)JCH * NL * 4, "Q" })) { fprintf(stderr, "gpu: device copy of Q differs from the host table\n"); exit(1); }
    CK(cudaMemset(D.ctr, 0, 32)); CK(cudaMemset(D.ns4, 0, 4)); CK(cudaMemset(D.ndmp, 0, 4)); CK(cudaMemset(D.flags, 0, 4));
    prm_t p = { ja, jb, n0, n1, C.C2K, C.UU, C.LLO, C.LHI, dump, 0 };
#ifdef Z2PROBE
    if (getenv("Z2_PROBE")) { p.probe = atoi(getenv("Z2_PROBE")); fprintf(stderr, "Z2PROBE build, probe %d: timing only, the counts are not valid\n", p.probe); }
#endif
    for (int k = 0; k < niv; k++)
        for (ull a = iv[2 * k]; a < iv[2 * k + 1]; a += LCH) {
            ull b = a + LCH < iv[2 * k + 1] ? a + LCH : iv[2 * k + 1];
            ksieve<<<(unsigned)((b - a + BLK - 1) / BLK), BLK>>>(D, p, a, b);
            CK(cudaGetLastError());
        }
    CK(cudaDeviceSynchronize());
    unsigned ns4, ndmp, fl;
    CK(cudaMemcpy(res->ctr, D.ctr, 32, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&ns4, D.ns4, 4, cudaMemcpyDeviceToHost)); CK(cudaMemcpy(&ndmp, D.ndmp, 4, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&fl, D.flags, 4, cudaMemcpyDeviceToHost));
    if (ns4 > S4CAP) { fl |= 8; ns4 = S4CAP; }
    if (ndmp > DCAP) { fl |= 16; ndmp = DCAP; }
    CK(cudaMemcpy(HS4, D.s4, (size_t)ns4 * 8, cudaMemcpyDeviceToHost));
    if (dump) CK(cudaMemcpy(HDMP, D.dmp, (size_t)ndmp * 8, cudaMemcpyDeviceToHost));
    res->s4 = (uint64_t *)HS4; res->ns4 = ns4; res->dmp = (uint64_t *)HDMP; res->ndmp = dump ? ndmp : 0; res->flags = fl;
}
