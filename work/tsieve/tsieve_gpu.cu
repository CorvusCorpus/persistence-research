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
#include <pthread.h>
#include "tsieve_gpu.h"
typedef unsigned long long u64; typedef long long i64; typedef unsigned int u32;

#ifndef IMAX
#define IMAX 1024                                      /* i0 values per chunk (small values in tests force chunking) */
#endif
#ifndef QCAP
#define QCAP 2048                                      /* queued candidates per round (more are handled inline) */
#endif
#define E9 1000000000ULL
#define LG3 0.47712125471966243730
#define LG7 0.84509804001425683071

typedef struct { const u32 *T2, *T3, *T7; const u64 *TOP; } tg_tabs;

#ifdef TG_EMU
#define DEV static inline
#define MULHI64(x, y) ((u64)(((unsigned __int128)(x) * (y)) >> 64))
#define CTZ32(w) __builtin_ctz(w)
#define ATOMIC_OR32(p, v) (*(p) |= (v))
#define ATOMIC_ADD32(p, v) ([&]() { u32 o_ = *(p); *(p) += (v); return o_; }())
#define ATOMIC_ADD64(p, v) ([&]() { u64 o_ = *(p); *(p) += (v); return o_; }())
static tg_params P;
static int NTE = 32;                                   /* emulated threads per block (<= 1024) */
#define PHASE_BEGIN for (int tid = 0; tid < NTE; tid++) { const int nth = NTE; (void)nth;
#define PHASE_END }
#define TLS_DECL(v) u64 v##_t[1024]; for (int t_ = 0; t_ < NTE; t_++) v##_t[t_] = 0;   /* per-thread across phases */
#define TLS(v) v##_t[tid]
#else
#define DEV __device__ __forceinline__
#define MULHI64(x, y) __umul64hi((x), (y))
#define CTZ32(w) (__ffs(w) - 1)
#define ATOMIC_OR32(p, v) atomicOr((p), (v))
#define ATOMIC_ADD32(p, v) atomicAdd((p), (v))
#define ATOMIC_ADD64(p, v) atomicAdd((unsigned long long *)(p), (unsigned long long)(v))
__constant__ tg_params P;
#define PHASE_BEGIN { const int tid = threadIdx.x; const int nth = blockDim.x; (void)nth;
#define PHASE_END } __syncthreads();
#define TLS_DECL(v) u64 v##_t = 0;
#define TLS(v) v##_t
#endif

DEV u64 mod_br(u64 x, u64 m, u64 br) { u64 q = MULHI64(x, br), r = x - q * m; while (r >= m) r -= m; return r; }   /* x mod m */
DEV u64 mods_br(i64 x, u64 m, u64 br) { if (x >= 0) return mod_br((u64)x, m, br); u64 r = mod_br((u64)(-x), m, br); return r ? m - r : 0; }
DEV i64 floordiv(i64 x, i64 c, u64 mc) {               /* floor(x / c), c > 0, mc = floor((2^64 - 1) / c): exact, no FP64 */
    const u64 xa = x < 0 ? (u64)(-x) : (u64)x;
    u64 q = MULHI64(xa, mc), r = xa - q * (u64)c;       /* q <= floor(xa / c) <= q + 1 */
    while (r >= (u64)c) { q++; r -= (u64)c; }
    return x >= 0 ? (i64)q : (r ? -(i64)q - 1 : -(i64)q);
}
DEV int ok9(u64 v, const unsigned char *ok4) { return ok4[v % 10000] & ok4[(v / 10000) % 10000] & P.ok1[v / 100000000ULL]; }

/* gsieve's stages 3-4 for a point whose last 2K digits are admissible; returns its decision (0 = S4 and further) */
DEV int later(u32 a, u32 b, u32 c, tg_tabs tb, const unsigned char *ok4, u64 *l2c, u64 *l3c, u64 *l4c) {
    if (b < P.cornerb && c < P.cornerc) return TG_CORNER;      /* covers b lg3 + c lg7 < cornerlg; no FP64 (slow on GeForce) */
    (*l2c)++;
    u64 bin = ((u64)a * P.l2 + (u64)b * P.l3 + (u64)c * P.l7) >> (64 - P.tb);
    if ((tb.TOP[bin >> 6] >> (bin & 63)) & 1) return TG_R3;
    (*l3c)++;
    const int NL = P.nl;
    const u32 *B3 = tb.T3 + b * NL, *B7 = tb.T7 + c * NL, *A = tb.T2 + a * NL;
    u64 m[16], acc, cy, cm;
    acc = (u64)B3[0] * B7[0]; cy = acc / E9; m[0] = acc - cy * E9;
    acc = (u64)B3[0] * B7[1] + (u64)B3[1] * B7[0] + cy; cy = acc / E9; m[1] = acc - cy * E9;
    acc = (u64)B3[0] * B7[2] + (u64)B3[1] * B7[1] + (u64)B3[2] * B7[0] + cy; cy = acc / E9; m[2] = acc - cy * E9;
    acc = (u64)B3[0] * B7[3] + (u64)B3[1] * B7[2] + (u64)B3[2] * B7[1] + (u64)B3[3] * B7[0] + cy; cm = acc / E9; m[3] = acc - cm * E9;
    u64 a0 = A[0], a1 = A[1], a2 = A[2], a3 = A[3];
    acc = a0 * m[0]; cy = acc / E9;
    acc = a0 * m[1] + a1 * m[0] + cy; cy = acc / E9;
    acc = a0 * m[2] + a1 * m[1] + a2 * m[0] + cy; cy = acc / E9; u64 l2 = acc - cy * E9;
    acc = a0 * m[3] + a1 * m[2] + a2 * m[1] + a3 * m[0] + cy; cy = acc / E9; u64 l3 = acc - cy * E9;
    if (!(P.ok1[l2 / 100000000ULL] & ok9(l3, ok4))) return TG_R4;              /* digits 27..36 (limb_check, first test) */
    /* the remaining limbs of M = 3^b 7^c and of 2^a M, in order and only as far as needed: the same limbs as limb_check
     * (carries cm, cy continue from limb 3), so the decision is the same AND over limbs 4..NL-1 */
    for (int k = 4; k < NL; k++) {
        u64 ac = cm;
        for (int u = 0; u <= k; u++) ac += (u64)B3[u] * B7[k - u];
        cm = ac / E9; m[k] = ac - cm * E9;
        ac = cy;
        for (int u = 0; u <= k; u++) ac += (u64)A[u] * m[k - u];
        cy = ac / E9;
        if (!ok9(ac - cy * E9, ok4)) return TG_R4;
    }
    (*l4c)++;
    return TG_S4;
}

#ifdef TG_EMU
static void k_coset(const tg_coset *cos, int ncos, tg_tabs tb, int dump, tg_point *out, u64 outcap, u64 *nout, u64 *gcnt, char *smem) {
    const int bid0 = 0, bstep = 1, nthr = NTE;
#else
__global__ void k_coset(const tg_coset *cos, int ncos, tg_tabs tb, int dump, tg_point *out, u64 outcap, u64 *nout, u64 *gcnt) {
    extern __shared__ __align__(16) char smem[];
    const int bid0 = blockIdx.x, bstep = gridDim.x, nthr = blockDim.x;
#endif
    const u64 Mm = P.mm, M1 = P.m1;
    const int nwords = (int)((2 * Mm) / 32 + 3);
    u32 *bm = (u32 *)smem;
    i64 *i1lo = (i64 *)(smem + ((nwords * 4 + 15) & ~15));
    u32 *rcnt = (u32 *)(i1lo + IMAX), *rstart = rcnt + IMAX;
    u64 *bc = (u64 *)(rstart + IMAX + 2);                  /* block counters: s1, rows, hits, s2, s3, s4 */
    u64 *dt = bc + 8;                                      /* dt[k] = k d2 mod 5^K (k < 32), dt[32] = 32 d2 mod 5^K */
    unsigned char *ok4 = (unsigned char *)(dt + 34);
    u64 *qd = (u64 *)(ok4 + 10000);                        /* queued candidates: a << 42 | b << 21 | c */
    u32 *qn = (u32 *)(qd + QCAP);                          /* queue lengths of even / odd rounds */
    PHASE_BEGIN
        for (int x = tid; x < 10000; x += nth) ok4[x] = P.ok4[x];
        for (int w = tid; w < nwords; w += nth) bm[w] = 0;
        if (tid < 8) bc[tid] = 0;
    PHASE_END
    for (int ci = bid0; ci < ncos; ci += bstep) {
        const tg_coset C = cos[ci];
        PHASE_BEGIN                                            /* 1. U, and the hit offsets dt */
            for (int k = tid; k <= 32; k += nth) dt[k] = mod_br((u64)k * C.d2, M1, P.br1);
            for (u64 l = tid; l < P.nleaf; l += nth) {
                u32 Q = (u32)C.h0; u64 y = 0, ll = l; int ok = 1;
                for (int i = 0; i < P.m; i++) {
                    int t;
                    if (P.dm == 0) {
                        int want = (int)(ll & 3), seen = -1; ll >>= 2;
                        for (t = 0; t < 5; t++) { u32 d = (Q + (u32)t * (u32)P.p2k[i]) % 10; if (d != 0 && d != 5 && ++seen == want) break; }
                    } else {
                        t = (int)(ll % 5); ll /= 5;
                        if ((Q + (u32)t * (u32)P.p2k[i]) % 10 == 0) { ok = 0; break; }
                    }
                    Q = (Q + (u32)t * (u32)P.p2k[i]) / 10; y += (u64)t * P.p5[i];
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
                        if (cc > 0) { i64 v = floordiv(rhs, cc, P.mc2[t]); if (v < hi) hi = v; }
                        else if (cc < 0) { i64 v = -floordiv(rhs, -cc, P.mc2[t]); if (v > lo) lo = v; }
                        else if (rhs < 0) empty = 1;
                    }
                    rcnt[k] = (empty || lo > hi) ? 0 : (u32)(hi - lo + 1); i1lo[k] = lo;
                }
            PHASE_END
            PHASE_BEGIN
                if (tid == 0) { rstart[0] = 0; for (int k = 0; k < nk; k++) rstart[k + 1] = rstart[k] + rcnt[k]; qn[0] = qn[1] = 0; }
            PHASE_END
            /* 3. rows, in rounds of one row per thread: (a) each thread walks its row and queues the points that pass
             *    digits K+1..2K (inline if the queue is full), (b) all threads take the queued points through later(), so
             *    that its branches run with full warps instead of the few lanes that had a candidate */
            TLS_DECL(ls1) TLS_DECL(lrows) TLS_DECL(lhits) TLS_DECL(l2c) TLS_DECL(l3c) TLS_DECL(l4c)
            const u32 total = rstart[nk];
            for (u32 r0 = 0, rnd = 0; r0 < total; r0 += (u32)nthr, rnd++) {
              PHASE_BEGIN                                      /* 3a */
                if (tid == 0) qn[(rnd + 1) & 1] = 0;               /* next round's counter: last used in round rnd - 1 */
                const u32 row = r0 + (u32)tid;
                do {
                    if (row >= total) break;
                    int lo = 0, hi = nk;
                    while (hi - lo > 1) { int mid = (lo + hi) >> 1; if (rstart[mid] <= row) lo = mid; else hi = mid; }
                    const i64 i0 = i0c + lo, i1 = i1lo[lo] + (i64)(row - rstart[lo]);
                    i64 lo2 = -(1LL << 62), hi2 = 1LL << 62; int empty = 0;
                    for (int q = 0; q < P.n3; q++) {
                        i64 rhs = P.b3_0[q] + P.b3_r[q] * (i64)C.r - P.a3[q][0] * i0 - P.a3[q][1] * i1, cc = P.a3[q][2];
                        if (cc > 0) { i64 v = floordiv(rhs, cc, P.mc3[q]); if (v < hi2) hi2 = v; }
                        else if (cc < 0) { i64 v = -floordiv(rhs, -cc, P.mc3[q]); if (v > lo2) lo2 = v; }
                        else if (rhs < 0) empty = 1;
                    }
                    if (empty || lo2 > hi2) break;
                    const u64 len = (u64)(hi2 - lo2 + 1);
                    TLS(ls1) += len; TLS(lrows)++;
                    u64 Ys = C.ybase + mod_br(mods_br(i0, M1, P.br1) * C.d0, M1, P.br1) + mod_br(mods_br(i1, M1, P.br1) * C.d1, M1, P.br1)
                           + mod_br(mods_br(lo2, M1, P.br1) * C.d2, M1, P.br1);
                    Ys = mod_br(Ys, M1, P.br1);
                    u64 s = C.s0 + mod_br(mods_br(i0, Mm, P.brm) * C.si0, Mm, P.brm) + mod_br(mods_br(i1, Mm, P.brm) * C.si1, Mm, P.brm)
                          + mods_br(lo2, Mm, P.brm);
                    s = mod_br(s, Mm, P.brm);
                    /* (a, b, c) of position n: row base + n v2, mod 2^32 (exact: a <= 3N, b <= 2N, c <= N are < 2^32) */
                    const u32 ar = (u32)((i64)C.r + i0 * P.v[0][0] + i1 * P.v[1][0] + lo2 * P.v[2][0]);
                    const u32 brw = (u32)(i0 * P.v[0][1] + i1 * P.v[1][1] + lo2 * P.v[2][1]);
                    const u32 crw = (u32)(i0 * P.v[0][2] + i1 * P.v[1][2] + lo2 * P.v[2][2]);
                    const u32 va = (u32)P.v[2][0], vb = (u32)P.v[2][1], vc = (u32)P.v[2][2];
                    for (u64 base = 0; base < len; base += Mm) {
                        const u64 wlen = len - base < Mm ? len - base : Mm;
                        u64 Yw = base ? mod_br(Ys + mod_br(base, M1, P.br1) * C.d2, M1, P.br1) : Ys;   /* Y at position base + off */
                        for (u64 off = 0; off < wlen; off += 32, Yw = Yw + dt[32] >= M1 ? Yw + dt[32] - M1 : Yw + dt[32]) {
                            const u64 pos = s + off, q32 = pos >> 5, sh = pos & 31;
                            u32 w = sh ? (bm[q32] >> sh) | (bm[q32 + 1] << (32 - sh)) : bm[q32];
                            const u64 rem = wlen - off;
                            if (rem < 32) w &= (1u << rem) - 1;
                            while (w) {
                                const int bit = CTZ32(w); w &= w - 1;
                                const u64 n = base + off + (u64)bit;
                                TLS(lhits)++;
                                const u64 Y = Yw + dt[bit] >= M1 ? Yw + dt[bit] - M1 : Yw + dt[bit];
                                u64 v = C.h0 + (Y << P.k);
                                for (int cc = 0; cc < P.skip; cc++) v /= 10000;
                                int ok = 1;
                                for (int cc = P.skip; cc < P.nfull; cc++) { ok &= ok4[v % 10000]; v /= 10000; }
                                if (P.nrem) ok &= P.okr[v];
                                if (!ok) continue;
                                const u32 nn = (u32)n, a = ar + nn * va, b = brw + nn * vb, c = crw + nn * vc;
                                const u32 qi = ATOMIC_ADD32(&qn[rnd & 1], 1u);
                                if (qi < QCAP) { qd[qi] = (u64)a << 42 | (u64)b << 21 | c; continue; }
                                const int code = later(a, b, c, tb, ok4, &TLS(l2c), &TLS(l3c), &TLS(l4c));   /* queue full */
                                if (dump || code == TG_S4 || code == TG_CORNER) {
                                    u64 idx = ATOMIC_ADD64(nout, 1);
                                    if (idx < outcap) { out[idx].a = a; out[idx].b = b; out[idx].c = c; out[idx].code = (u32)code; }
                                }
                            }
                        }
                    }
                } while (0);
              PHASE_END
              PHASE_BEGIN                                      /* 3b */
                const u32 nq = qn[rnd & 1] < QCAP ? qn[rnd & 1] : QCAP;
                for (u32 x = (u32)tid; x < nq; x += (u32)nth) {
                    const u64 e = qd[x];
                    const u32 a = (u32)(e >> 42), b = (u32)(e >> 21) & 0x1FFFFF, c = (u32)e & 0x1FFFFF;
                    const int code = later(a, b, c, tb, ok4, &TLS(l2c), &TLS(l3c), &TLS(l4c));
                    if (dump || code == TG_S4 || code == TG_CORNER) {
                        u64 idx = ATOMIC_ADD64(nout, 1);
                        if (idx < outcap) { out[idx].a = a; out[idx].b = b; out[idx].c = c; out[idx].code = (u32)code; }
                    }
                }
              PHASE_END
            }
            PHASE_BEGIN
                ATOMIC_ADD64(&bc[0], TLS(ls1)); ATOMIC_ADD64(&bc[1], TLS(lrows)); ATOMIC_ADD64(&bc[2], TLS(lhits));
                ATOMIC_ADD64(&bc[3], TLS(l2c)); ATOMIC_ADD64(&bc[4], TLS(l3c)); ATOMIC_ADD64(&bc[5], TLS(l4c));
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

static u64 fnv(u64 h, const void *p, size_t n) { const unsigned char *c = (const unsigned char *)p; for (size_t i = 0; i < n; i++) h = (h ^ c[i]) * 1099511628211ULL; return h; }
static const u32 *HT2, *HT3, *HT7; static const u64 *HTOP; static u64 HN2, HN3, HN7, HNTOP, HHASH; static tg_params HPAR;
static u64 host_hash(void) {
    u64 h = fnv(1469598103934665603ULL, &HPAR, sizeof HPAR);
    h = fnv(h, HT2, HN2 * 4); h = fnv(h, HT3, HN3 * 4); h = fnv(h, HT7, HN7 * 4); return fnv(h, HTOP, HNTOP * 8);
}
static void keep_host(const tg_params *p, const u32 *T2, const u32 *T3, const u32 *T7, u64 n2, u64 n3, u64 n7, const u64 *TOP, u64 topwords) {
    HPAR = *p; HT2 = T2; HT3 = T3; HT7 = T7; HTOP = TOP; HN2 = n2; HN3 = n3; HN7 = n7; HNTOP = topwords; HHASH = host_hash();
}
static size_t smem_bytes(const tg_params *p) {
    size_t nwords = (size_t)((2 * p->mm) / 32 + 3);
    return ((nwords * 4 + 15) & ~(size_t)15) + IMAX * 8 + (2 * IMAX + 2) * 4 + 8 * 8 + 34 * 8 + 10000 + QCAP * 8 + 2 * 4;
}

#ifdef TG_EMU
/* ---------------- CPU emulation ---------------- */
struct tg_ctx { tg_point *out; u64 outcap; char *smem; };
static tg_tabs TABS;
extern "C" int tg_init(const tg_params *p, const u32 *T2, const u32 *T3, const u32 *T7, u64 n2, u64 n3, u64 n7, const u64 *TOP, u64 topwords) {
    P = *p; TABS.T2 = T2; TABS.T3 = T3; TABS.T7 = T7; TABS.TOP = TOP;
    keep_host(p, T2, T3, T7, n2, n3, n7, TOP, topwords);
    if (getenv("TG_EMU_THREADS")) NTE = atoi(getenv("TG_EMU_THREADS"));
    if (NTE < 1 || NTE > 1024) { fprintf(stderr, "TG_EMU_THREADS must be 1..1024\n"); return 1; }
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
extern "C" int tg_verify(void) { return host_hash() != HHASH || memcmp(&P, &HPAR, sizeof P) != 0; }
#else
/* ---------------- CUDA ---------------- */
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(5); } } while (0)
struct tg_ctx { cudaStream_t st; tg_coset *dcos, *hcos; tg_point *dout, *hout; u64 *dno, *dcnt, *hno; u64 outcap; int maxcos; };
static tg_tabs DTABS; static tg_params HP; static int NSM; static char DEVNAME[256];
extern "C" int tg_init(const tg_params *p, const u32 *T2, const u32 *T3, const u32 *T7, u64 n2, u64 n3, u64 n7, const u64 *TOP, u64 topwords) {
    HP = *p; keep_host(p, T2, T3, T7, n2, n3, n7, TOP, topwords);
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr, 0)); NSM = pr.multiProcessorCount; snprintf(DEVNAME, sizeof DEVNAME, "%s", pr.name);
    CK(cudaMemcpyToSymbol(P, p, sizeof *p));
    u32 *d2, *d3, *d7; u64 *dtop;
    CK(cudaMalloc(&d2, n2 * 4)); CK(cudaMalloc(&d3, n3 * 4)); CK(cudaMalloc(&d7, n7 * 4)); CK(cudaMalloc(&dtop, topwords * 8));
    CK(cudaMemcpy(d2, T2, n2 * 4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d3, T3, n3 * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d7, T7, n7 * 4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(dtop, TOP, topwords * 8, cudaMemcpyHostToDevice));
    DTABS.T2 = d2; DTABS.T3 = d3; DTABS.T7 = d7; DTABS.TOP = dtop;
    size_t sm = smem_bytes(p);
    if (sm > 48 * 1024) CK(cudaFuncSetAttribute(k_coset, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)sm));
    cudaFuncAttributes fa; CK(cudaFuncGetAttributes(&fa, k_coset));
    if (p->nt < 32 || p->nt > fa.maxThreadsPerBlock || p->nt % 32) {
        fprintf(stderr, "TG_NT=%d: this kernel takes a multiple of 32 up to %d threads per block (%d registers each)\n", p->nt, fa.maxThreadsPerBlock, fa.numRegs);
        return 1;
    }
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
static int dev_equal(const void *dptr, const void *hptr, size_t bytes) {   /* device buffer == host buffer, in 64 MB pieces */
    static char *buf; static size_t cap; int eq = 1;
    if (!buf) { cap = (size_t)64 << 20; CK(cudaMallocHost(&buf, cap)); }
    for (size_t off = 0; off < bytes && eq; off += cap) {
        size_t n = bytes - off < cap ? bytes - off : cap;
        CK(cudaMemcpy(buf, (const char *)dptr + off, n, cudaMemcpyDeviceToHost));        /* pinned, synchronous */
        eq = memcmp(buf, (const char *)hptr + off, n) == 0;
    }
    return eq;
}
extern "C" int tg_verify(void) {
    static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mu);
    tg_params dp; CK(cudaMemcpyFromSymbol(&dp, P, sizeof dp));
    int bad = host_hash() != HHASH || memcmp(&dp, &HPAR, sizeof dp) != 0
           || !dev_equal(DTABS.T2, HT2, HN2 * 4) || !dev_equal(DTABS.T3, HT3, HN3 * 4)
           || !dev_equal(DTABS.T7, HT7, HN7 * 4) || !dev_equal(DTABS.TOP, HTOP, HNTOP * 8);
    pthread_mutex_unlock(&mu);
    return bad;
}
#endif
