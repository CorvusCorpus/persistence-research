/* Kernel of the GPU bottom solver: the per-solve loop of csrc/bottom3.c ('cum' mode) in fixed-width exact arithmetic.
 *
 * One thread per lower configuration ("Flow" = -1 + sum c_q 10^pos_q over the digits below the top one), as in bottom3.c's
 * cum_rec: the thread computes low(Flow) once, then solves for every height p of the top digit. Every quantity is the same
 * residue bottom3.c computes (mod 5^M with M = G, mod 2^n with n = M - t), obtained without exact divisions:
 *   log1p(p^v y) = sum_n (-1)^(n+1) p^(n v - e) y^n / n'   (n = p^e n')      -> Horner in y with tabulated coefficients
 *   plog(x)      = log(x^(p^s)) / p^s = sum_n (-1)^(n+1) p^(n vz - e - s) z'^n / n',  x^(p^s) = 1 + p^vz z'
 *   Flow^-1 mod 5^G, mod 2^K by Newton; D5 = T / 5 and z / 5^vz as products with 5^-1, 5^-vz mod 2^(32 L) (exact quotients).
 *   CRT by Garner: beta = beta2 + 2^(n-2) ((beta5 - beta2) 2^-(n-2) mod 5^(G-1)).
 * 5-adic products are Montgomery products modulo 5^N (plog's powering, N = G + s), 5^G and 5^(G-1), all with the same
 * R = 2^(32 L5); tables hold one operand in Montgomery form so that products come out as plain residues.
 * Lattice step: the coset point nearest the box centre has coordinates floor(K - off0 R) / 2^F in the reduced basis, with
 * R = b/det and K tabulated in fixed point (F = bits of Q + 64, error < 2^-64); when the fractional part is within 2^-60 of an
 * integer the solve is handed to the host (exact GMP, bottom3.c's own solutions()). Candidates are then exact.
 * Compiled for the CPU with -DBG_EMU (same functions, a loop instead of the grid) for validation without a GPU. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bigfix.h"
#include "bottomgpu.h"

#ifdef BG_EMU
#define DEVK static
#define DF static inline
#define LDG(p) (*(p))
#else
#define DEVK __device__
#define DF __device__ __forceinline__
#define LDG(p) __ldg(p)
#endif

template <int L> DF void ldb(Big<L> &x, const u32 *p) {
#pragma unroll
    for (int i = 0; i < L; i++) x.w[i] = LDG(p + i);
}
template <int A, int B> DF void widen(Big<B> &r, const Big<A> &a) {   /* zero-extend or truncate */
#pragma unroll
    for (int i = 0; i < B; i++) r.w[i] = i < A ? a.w[i] : 0u;
}
template <int L> DF void bneg(Big<L> &r, const Big<L> &a) { Big<L> z; bzero_(z); bsub(r, z, a); }
template <int L> DF int bneg_p(const Big<L> &a) { return (int)(a.w[L - 1] >> 31); }   /* two's complement sign */
template <int L> DF void bsub1(Big<L> &a) { Big<L> one; bset(one, 1); bsub(a, a, one); }

struct Ctx { const bg_params *P; const u32 *T; };

DF unsigned long long bump(unsigned long long *nout) {
#ifdef BG_EMU
    return (*nout)++;
#else
    return atomicAdd(nout, 1ULL);
#endif
}
DF void put(bg_out *out, long long outcap, unsigned long long *nout, unsigned flow, int p, int t, int kind, int seq) {
    const unsigned long long k = bump(nout);
    if ((long long)k < outcap) { out[k].flow = flow; out[k].p = (unsigned)p; out[k].t = (unsigned)t; out[k].kind = (unsigned)kind; out[k].seq = (unsigned)seq; }
}

/* sum_{n=1..N} c_n y^n mod m, c_n plain residues at cf[(n-1) L], yh = y R mod m */
template <int L> DF void horner5(Big<L> &r, const u32 *cf, int N, const Big<L> &yh, const Big<L> &m, u32 mp) {
    if (N <= 0) { bzero_(r); return; }
    Big<L> acc, c; ldb(acc, cf + (size_t)(N - 1) * L);
    for (int n = N - 1; n >= 1; n--) { bmont(acc, acc, yh, m, mp); ldb(c, cf + (size_t)(n - 1) * L); baddmod(acc, acc, c, m); }
    bmont(r, acc, yh, m, mp);
}
/* sum_{n=1..N} c_n y^n mod 2^(32 L) */
template <int L> DF void horner2(Big<L> &r, const u32 *cf, int N, const Big<L> &y) {
    if (N <= 0) { bzero_(r); return; }                   /* e.g. log1p(2^v y) with v >= G: 0 mod 2^G */
    Big<L> acc, c; ldb(acc, cf + (size_t)(N - 1) * L);
    for (int n = N - 1; n >= 1; n--) { bmul_lo(acc, acc, y); ldb(c, cf + (size_t)(n - 1) * L); badd(acc, acc, c); }
    bmul_lo(r, acc, y);
}
/* plog(x^4) mod 5^G for a residue x mod 5^N (x not divisible by 5) */
template <int L5> DF void plog5(Big<L5> &r, const Big<L5> &x, const Ctx &C) {
    const bg_params &P = *C.P; const u32 *T = C.T;
    Big<L5> mN, R2N, mG, R2G; ldb(mN, T + P.o_mN); ldb(R2N, T + P.o_R2N); ldb(mG, T + P.o_mG); ldb(R2G, T + P.o_R2G);
    Big<L5> y, y2;
    bmont(y, x, R2N, mN, P.mpN); bmont(y2, y, y, mN, P.mpN); bmont(y, y2, y2, mN, P.mpN);   /* x^4 R */
    for (int i = 0; i < P.s5; i++) { bmont(y2, y, y, mN, P.mpN); bmont(y2, y2, y2, mN, P.mpN); bmont(y, y2, y, mN, P.mpN); }
    Big<L5> one; bset(one, 1);
    bmont(y, y, one, mN, P.mpN);                         /* X = x^(4 5^s) mod 5^N, plain */
    bsub1(y);                                            /* z = X - 1 = 5^vz z' */
    Big<L5> iz; ldb(iz, T + P.o_invz5); bmul_lo(y, y, iz);   /* z' < 5^(N - vz) = 5^(G - 1) */
    bmont(y2, y, R2G, mG, P.mpG);                        /* z' R mod 5^G */
    horner5(r, T + P.o_plog5, P.np5, y2, mG, P.mpG);
}
/* plog(x) mod 2^G for x = 1 mod 4 (x a residue mod 2^(32 L2), meaningful mod 2^(G+2)) */
template <int L2> DF void plog2(Big<L2> &r, const Big<L2> &x, const Ctx &C) {
    const bg_params &P = *C.P; const u32 *T = C.T;
    Big<L2> X = x;
    for (int i = 0; i < P.s2; i++) bmul_lo(X, X, X);
    bmask(X, P.G + P.s2);                                /* x^(2^s) mod 2^N, N = G + s */
    bsub1(X);                                            /* z = 2^vz z' */
    bshr(X, X, P.vz2);
    horner2(r, T + P.o_plog2, P.np2, X);
    bmask(r, P.G);
}

/* every solve (Flow, c, p), p in [lo, phi), of one lower configuration, given Flow mod 5^N, mod 5^G, mod 2^(32 L2), mod 5 */
template <int L5, int L2, int LQ, int LH>
DF void do_core(const Ctx &C, const Big<L5> &F5N, const Big<L5> &F5G, const Big<L2> &F2, const int f5, const int ci0,
                const int lo, const int phi, unsigned flow_idx, bg_out *out, long long outcap, unsigned long long *nout,
                unsigned long long cnt[3]) {
    const bg_params &P = *C.P; const u32 *T = C.T;
    const int G = P.G;
    Big<L5> mN, mG, mQ, t5; ldb(mN, T + P.o_mN); ldb(mG, T + P.o_mG); ldb(mQ, T + P.o_mQ);
    Big<L2> t2;
    const int tl = bctz(F2);
    if (tl > P.tlmax) {                                  /* rare: the host redoes (and counts) every solve of this flow */
        for (int p = lo; p < phi; p++) put(out, outcap, nout, flow_idx, p, 0, 3, 0);
        return;
    }
    /* low(Flow) */
    Big<L2> odd; bshr(odd, F2, tl);
    const int sig = (odd.w[0] & 3u) == 1u ? 1 : -1;
    Big<L5> l5low; plog5(l5low, F5N, C);
    Big<L5> invh;                                        /* Flow^-1 R mod 5^G by Newton, from Flow^-1 mod 5 */
    {
        Big<L5> R2G, twoR, Fh, e, d; ldb(R2G, T + P.o_R2G); ldb(twoR, T + P.o_twoRG);
        bmont(Fh, F5G, R2G, mG, P.mpG);
        Big<L5> y0; bset(y0, (u64)(f5 * f5 * f5 % 5)); bmont(invh, y0, R2G, mG, P.mpG);   /* f5^-1 = f5^3 mod 5 */
        for (int k = 1; k < 2 * G; k <<= 1) { bmont(e, Fh, invh, mG, P.mpG); bsubmod(d, twoR, e, mG); bmont(invh, invh, d, mG, P.mpG); }
    }
    Big<L2> l2low, inv2;
    {
        Big<L2> x = odd; if (sig < 0) bneg(x, odd);
        plog2(l2low, x, C);
        Big<L2> e, d, two; bset(two, 2); inv2 = odd;     /* odd^-1 mod 2^(32 L2): correct mod 8 at the start */
        for (int k = 3; k < 32 * L2; k <<= 1) { bmul_lo(e, odd, inv2); bsub(d, two, e); bmul_lo(inv2, inv2, d); }
    }
    const int c0 = P.cval[ci0];
    for (int p = lo; p < phi; p++) {
        cnt[0]++;
        Big<L2> Fp2; ldb(t2, T + P.o_T10_2 + ((size_t)ci0 * G + p) * L2); badd(Fp2, F2, t2);
        const int t = bctz(Fp2);
        if (t >= G || t > P.tm) continue;
        const int f5F = p >= 1 ? f5 : (f5 + c0) % 5;
        if (f5F == 0) continue;
        const int n = G - t;
        if (!LDG(T + P.o_tok + t)) {                     /* lattice too fine for this t: bottom3.c decides (and stops) */
            if (p == 0) cnt[1]++;
            if (!(tl == t && p + P.v2c[ci0] - tl >= 2)) cnt[2]++;
            put(out, outcap, nout, flow_idx, p, t, 2, 0);
            continue;
        }
        Big<L5> l5F;
        if (p >= 1) {
            Big<L5> yh, Lg; ldb(t5, T + P.o_Tc5 + ((size_t)ci0 * G + p) * L5);
            bmont(yh, t5, invh, mG, P.mpG);              /* c 2^p Flow^-1 R */
            const int N5 = (int)LDG(T + P.o_n5 + p); const u32 off = LDG(T + P.o_off5 + p);
            horner5(Lg, T + P.o_coef5 + (size_t)off * L5, N5, yh, mG, P.mpG);
            baddmod(Lg, Lg, Lg, mG); baddmod(Lg, Lg, Lg, mG); baddmod(l5F, Lg, l5low, mG);
        } else {
            cnt[1]++;
            Big<L5> Fp5; ldb(t5, T + P.o_T10N + ((size_t)ci0 * G + 0) * L5); baddmod(Fp5, F5N, t5, mN);
            plog5(l5F, Fp5, C);
        }
        const int vx2 = p + P.v2c[ci0] - tl;
        Big<L2> l2F; int sigma;
        if (tl == t && vx2 >= 2) {
            Big<L2> y2; ldb(t2, T + P.o_Tc2 + ((size_t)ci0 * G + p) * L2); bmul_lo(y2, t2, inv2); bmask(y2, G);
            const int N2 = (int)LDG(T + P.o_n2 + vx2); const u32 off = LDG(T + P.o_off2 + vx2);
            horner2(l2F, T + P.o_coef2 + (size_t)off * L2, N2, y2);
            badd(l2F, l2F, l2low); bmask(l2F, n);
            sigma = sig;
        } else {
            cnt[2]++;
            Big<L2> Tt; bshr(Tt, Fp2, t); bmask(Tt, n + 2);
            sigma = (Tt.w[0] & 3u) == 1u ? 1 : -1;
            if (sigma < 0) { bneg(Tt, Tt); bmask(Tt, n + 2); }
            plog2(l2F, Tt, C); bmask(l2F, n);
        }
        /* beta5 = ((l5F - t L5_16) / 5) A5inv mod 5^(G-1); beta2 = (l2F >> 2) Ainv mod 2^(n-2) */
        Big<L5> T5, D5, beta5; ldb(t5, T + P.o_L16 + (size_t)t * L5); bsubmod(T5, l5F, t5, mG);
        Big<L5> i5; ldb(i5, T + P.o_inv5); bmul_lo(D5, T5, i5);
        ldb(t5, T + P.o_A5inv + (size_t)t * L5); bmont(beta5, D5, t5, mQ, P.mpQ);
        Big<L2> D2, beta2; bshr(D2, l2F, 2); ldb(t2, T + P.o_Ainv + (size_t)t * L2); bmul_lo(beta2, D2, t2); bmask(beta2, n - 2);
        Big<L5> b25, diff, h; widen<L2, L5>(b25, beta2); bsubmod(diff, beta5, b25, mQ);
        ldb(t5, T + P.o_q2inv5 + (size_t)t * L5); bmont(h, diff, t5, mQ, P.mpQ);
        Big<LQ> beta, tq; widen<L5, LQ>(beta, h); bshl(beta, beta, n - 2); widen<L2, LQ>(tq, beta2); badd(beta, beta, tq);
        /* lattice classes w0 = 0..3 (bottom3.c) */
        int r5 = f5F; for (int i = 0; i < (t & 3); i++) r5 = r5 * 3 % 5;          /* F 2^-t mod 5 (3 = 2^-1, order 4) */
        const int dl5 = r5 == 1 ? 0 : r5 == 2 ? 1 : r5 == 4 ? 2 : 3;               /* dlog base 2 mod 5 */
        const int a2m4 = (int)LDG(T + P.o_a2m4 + t), b2m4 = (int)(beta2.w[0] & 3u);
        Big<LQ> Q; ldb(Q, T + P.o_Q + (size_t)t * LQ);
        for (int w0 = 0; w0 < 4; w0++) {
            const int u0 = (a2m4 * w0 + b2m4) & 3;
            if (((u0 + w0) & 1) != (sigma == 1 ? 0 : 1)) continue;
            if (((3 * u0 + w0) & 3) != dl5) continue;
            Big<LQ> off0, aw; ldb(aw, T + P.o_alphaW + ((size_t)t * 4 + w0) * LQ); baddmod(off0, aw, beta, Q);
            /* ci = floor((K_i - off0 R_i) / 2^F), cj = floor((K_j - off0 R_j) / 2^F); K in 2LQ limbs, R in LQ limbs (signed) */
            Big<LH> ci, cj; int amb = 0;
#pragma unroll
            for (int s = 0; s < 2; s++) {
                Big<LQ> R; ldb(R, T + (s ? P.o_Rj : P.o_Ri) + (size_t)t * LQ);
                const int rneg = bneg_p(R); if (rneg) bneg(R, R);
                Big<2 * LQ> K, prod; ldb(K, T + (s ? P.o_Kj : P.o_Ki) + ((size_t)t * 4 + w0) * 2 * LQ);
                bmul_full<LQ>(prod.w, off0, R);
                if (rneg) badd(K, K, prod); else bsub(K, K, prod);
                /* the 60 fraction bits below 2^F all equal: within 2^-60 of an integer -> the host decides */
                Big<2 * LQ> fr; bshr(fr, K, P.F - 60);
                const u32 lo32 = fr.w[0], hi28 = fr.w[1] & 0x0FFFFFFFu;
                if ((lo32 == 0u && hi28 == 0u) || (lo32 == 0xFFFFFFFFu && hi28 == 0x0FFFFFFFu)) amb = 1;
                Big<2 * LQ> q; bshr(q, K, P.F);
                if (bneg_p(K)) { Big<2 * LQ> ones, hi; bzero_(ones); bsub1(ones); bshl(hi, ones, 64 * LQ - P.F); badd(q, q, hi); }
                if (s) widen<2 * LQ, LH>(cj, q); else widen<2 * LQ, LH>(ci, q);
            }
            if (amb) { put(out, outcap, nout, flow_idx, p, t, 2, 0); break; }   /* host redoes the whole solve */
            /* up = off0 + ci b1x + cj b2x, ww = w0 + ci b1y + cj b2y (two's complement, 2 LH limbs), then the 4 neighbours */
            Big<2 * LH> up, ww, pr, bx[4];
            widen<LQ, 2 * LH>(up, off0); bset(ww, (u64)w0);
#pragma unroll
            for (int i = 0; i < 4; i++) {                /* i = 0: ci b1x, 1: ci b1y, 2: cj b2x, 3: cj b2y */
                Big<LH> bv, a = (i >> 1) ? cj : ci; ldb(bv, T + P.o_b + ((size_t)t * 4 + i) * LH);
                const int ng = bneg_p(bv);
#pragma unroll
                for (int k = 0; k < 2 * LH; k++) bx[i].w[k] = k < LH ? bv.w[k] : (ng ? 0xFFFFFFFFu : 0u);   /* sign-extended */
                const int sn = bneg_p(a) ^ ng;
                if (bneg_p(a)) bneg(a, a);
                if (ng) bneg(bv, bv);
                bmul_full<LH>(pr.w, a, bv); if (sn) bneg(pr, pr);
                if (i & 1) badd(ww, ww, pr); else badd(up, up, pr);
            }
            Big<2 * LH> U2, W2, two;
            { Big<LQ> tU; ldb(tU, T + P.o_U); widen<LQ, 2 * LH>(U2, tU); ldb(tU, T + P.o_W); widen<LQ, 2 * LH>(W2, tU); }
            bset(two, 2);
            for (int di = 0; di <= 1; di++) for (int dj = 0; dj <= 1; dj++) {
                Big<2 * LH> u = up, w = ww;
                if (di) { badd(u, u, bx[0]); badd(w, w, bx[1]); }
                if (dj) { badd(u, u, bx[2]); badd(w, w, bx[3]); }
                if (bneg_p(u) || bneg_p(w)) continue;
                if (bcmp_(u, two) < 0 || bcmp_(u, U2) > 0 || bcmp_(w, W2) > 0) continue;
                const unsigned long long k = bump(nout);
                if ((long long)k < outcap) {
                    out[k].flow = flow_idx; out[k].p = (unsigned)p; out[k].t = (unsigned)t; out[k].kind = 1; out[k].seq = (unsigned)(4 * w0 + 2 * di + dj);
#pragma unroll
                    for (int i = 0; i < 32; i++) { out[k].up[i] = i < 2 * LH ? u.w[i] : 0u; out[k].w[i] = i < 2 * LH ? w.w[i] : 0u; }
                }
            }
        }
    }
}

/* a lower configuration given by its digits and heights: Flow = -1 + sum c_q 10^pos_q from the tables */
template <int L5, int L2, int LQ, int LH>
DF void do_flow(const Ctx &C, const unsigned char *ord, const bg_flow &fl, unsigned flow_idx, bg_out *out, long long outcap,
                unsigned long long *nout, unsigned long long cnt[3]) {
    const bg_params &P = *C.P; const u32 *T = C.T;
    const int b = P.b, G = P.G;
    const unsigned char *oc = ord + (size_t)fl.o * b;
    Big<L5> mN, mG, t5; ldb(mN, T + P.o_mN); ldb(mG, T + P.o_mG);
    Big<L5> F5N = mN, F5G = mG; bsub1(F5N); bsub1(F5G);
    Big<L2> F2, t2; bzero_(F2); bsub1(F2);
    int f5 = 4;
    for (int r = 1; r < b; r++) {
        const int ci = oc[r], pos = fl.pos[r];
        ldb(t5, T + P.o_T10N + ((size_t)ci * G + pos) * L5); baddmod(F5N, F5N, t5, mN);
        ldb(t5, T + P.o_T10G + ((size_t)ci * G + pos) * L5); baddmod(F5G, F5G, t5, mG);
        ldb(t2, T + P.o_T10_2 + ((size_t)ci * G + pos) * L2); badd(F2, F2, t2);
        if (pos == 0) f5 = (f5 + P.cval[ci]) % 5;
    }
    do_core<L5, L2, LQ, LH>(C, F5N, F5G, F2, f5, oc[0], fl.lo, fl.phi, flow_idx, out, outcap, nout, cnt);
}
/* kernel test: an arbitrary Flow */
template <int L5, int L2, int LQ, int LH>
DF void do_query(const Ctx &C, const bg_query &q, unsigned idx, bg_out *out, long long outcap, unsigned long long *nout,
                 unsigned long long cnt[3]) {
    Big<L5> F5N, F5G; Big<L2> F2;
#pragma unroll
    for (int i = 0; i < L5; i++) { F5N.w[i] = q.f5n[i]; F5G.w[i] = q.f5g[i]; }
#pragma unroll
    for (int i = 0; i < L2; i++) F2.w[i] = q.f2[i];
    do_core<L5, L2, LQ, LH>(C, F5N, F5G, F2, q.f5, q.ci, q.p, q.p + 1, idx, out, outcap, nout, cnt);
}

/* ---------------- size classes ---------------- */
#define BG_CLASSES(X) X(0, 6, 5, 10, 7) X(1, 8, 6, 13, 9) X(2, 11, 7, 17, 11) X(3, 13, 9, 21, 13) X(4, 16, 10, 24, 14)

#ifdef BG_EMU
static bg_params HP; static const u32 *HT; static const unsigned char *HORD;
extern "C" int bg_init(const bg_params *p, const unsigned int *tables, const unsigned char *orders_ci, int norders) {
    (void)norders; HP = *p; HT = tables; HORD = orders_ci; return 0;
}
extern "C" long long bg_run(const bg_flow *flows, int n, bg_out *out, long long outcap, unsigned long long counts[3]) {
    unsigned long long nout = 0, cnt[3] = {0, 0, 0}; Ctx C = {&HP, HT};
    for (int i = 0; i < n; i++) {
        switch (HP.cls) {
#define CASE(k, a, b, c, d) case k: do_flow<a, b, c, d>(C, HORD, flows[i], (unsigned)i, out, outcap, &nout, cnt); break;
            BG_CLASSES(CASE)
#undef CASE
        }
    }
    if ((long long)nout > outcap) return -1;             /* counts untouched: the caller reruns the batch */
    for (int k = 0; k < 3; k++) counts[k] += cnt[k];
    return (long long)nout;
}
extern "C" long long bg_runq(const bg_query *q, int n, bg_out *out, long long outcap, unsigned long long counts[3]) {
    unsigned long long nout = 0, cnt[3] = {0, 0, 0}; Ctx C = {&HP, HT};
    for (int i = 0; i < n; i++) {
        switch (HP.cls) {
#define CASE(k, a, b, c, d) case k: do_query<a, b, c, d>(C, q[i], (unsigned)i, out, outcap, &nout, cnt); break;
            BG_CLASSES(CASE)
#undef CASE
        }
    }
    if ((long long)nout > outcap) return -1;
    for (int k = 0; k < 3; k++) counts[k] += cnt[k];
    return (long long)nout;
}
extern "C" const char *bg_device(void) { return "emulation"; }
#else
#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(5); } } while (0)
__constant__ bg_params DP;
template <int L5, int L2, int LQ, int LH>
__global__ void k_flows(const u32 *T, const unsigned char *ord, const bg_flow *flows, int n, bg_out *out, long long outcap,
                        unsigned long long *nout, unsigned long long *gcnt) {
    unsigned long long cnt[3] = {0, 0, 0}; Ctx C = {&DP, T};
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x)
        do_flow<L5, L2, LQ, LH>(C, ord, flows[i], (unsigned)i, out, outcap, nout, cnt);
    for (int k = 0; k < 3; k++) if (cnt[k]) atomicAdd(&gcnt[k], cnt[k]);
}
template <int L5, int L2, int LQ, int LH>
__global__ void k_queries(const u32 *T, const bg_query *q, int n, bg_out *out, long long outcap, unsigned long long *nout,
                          unsigned long long *gcnt) {
    unsigned long long cnt[3] = {0, 0, 0}; Ctx C = {&DP, T};
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x)
        do_query<L5, L2, LQ, LH>(C, q[i], (unsigned)i, out, outcap, nout, cnt);
    for (int k = 0; k < 3; k++) if (cnt[k]) atomicAdd(&gcnt[k], cnt[k]);
}
static bg_params HP; static u32 *DT; static unsigned char *DORD; static bg_flow *DFL; static bg_out *DOUT;
static unsigned long long *DCNT; static int DFLCAP; static long long DOUTCAP; static char DEVNAME[256]; static int NSM;
extern "C" int bg_init(const bg_params *p, const unsigned int *tables, const unsigned char *orders_ci, int norders) {
    HP = *p;
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr, 0)); NSM = pr.multiProcessorCount; snprintf(DEVNAME, sizeof DEVNAME, "%s", pr.name);
    CK(cudaMemcpyToSymbol(DP, p, sizeof *p));
    CK(cudaMalloc(&DT, (size_t)p->words * 4)); CK(cudaMemcpy(DT, tables, (size_t)p->words * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DORD, (size_t)norders * p->b + 1)); CK(cudaMemcpy(DORD, orders_ci, (size_t)norders * p->b, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&DCNT, 4 * 8));
    return 0;
}
extern "C" long long bg_run(const bg_flow *flows, int n, bg_out *out, long long outcap, unsigned long long counts[3]) {
    if (n > DFLCAP) { if (DFL) cudaFree(DFL); DFLCAP = n; CK(cudaMalloc(&DFL, (size_t)n * sizeof(bg_flow))); }
    if (outcap > DOUTCAP) { if (DOUT) cudaFree(DOUT); DOUTCAP = outcap; CK(cudaMalloc(&DOUT, (size_t)outcap * sizeof(bg_out))); }
    CK(cudaMemcpy(DFL, flows, (size_t)n * sizeof(bg_flow), cudaMemcpyHostToDevice));
    CK(cudaMemset(DCNT, 0, 4 * 8));
    const int nt = 128; int grid = (n + nt - 1) / nt; if (grid > 64 * NSM) grid = 64 * NSM;
    switch (HP.cls) {
#define CASE(k, a, b, c, d) case k: k_flows<a, b, c, d><<<grid, nt>>>(DT, DORD, DFL, n, DOUT, outcap, DCNT + 3, DCNT); break;
        BG_CLASSES(CASE)
#undef CASE
    }
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    unsigned long long h[4]; CK(cudaMemcpy(h, DCNT, 4 * 8, cudaMemcpyDeviceToHost));
    long long no = (long long)h[3];
    if (no > outcap) return -1;                          /* counts untouched: the caller reruns the batch */
    for (int k = 0; k < 3; k++) counts[k] += h[k];
    if (no) CK(cudaMemcpy(out, DOUT, (size_t)no * sizeof(bg_out), cudaMemcpyDeviceToHost));
    return no;
}
extern "C" long long bg_runq(const bg_query *q, int n, bg_out *out, long long outcap, unsigned long long counts[3]) {
    bg_query *dq; CK(cudaMalloc(&dq, (size_t)n * sizeof(bg_query))); CK(cudaMemcpy(dq, q, (size_t)n * sizeof(bg_query), cudaMemcpyHostToDevice));
    if (outcap > DOUTCAP) { if (DOUT) cudaFree(DOUT); DOUTCAP = outcap; CK(cudaMalloc(&DOUT, (size_t)outcap * sizeof(bg_out))); }
    CK(cudaMemset(DCNT, 0, 4 * 8));
    const int nt = 128; int grid = (n + nt - 1) / nt; if (grid > 64 * NSM) grid = 64 * NSM;
    switch (HP.cls) {
#define CASE(k, a, b, c, d) case k: k_queries<a, b, c, d><<<grid, nt>>>(DT, dq, n, DOUT, outcap, DCNT + 3, DCNT); break;
        BG_CLASSES(CASE)
#undef CASE
    }
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize()); CK(cudaFree(dq));
    unsigned long long h[4]; CK(cudaMemcpy(h, DCNT, 4 * 8, cudaMemcpyDeviceToHost));
    long long no = (long long)h[3];
    if (no > outcap) return -1;
    for (int k = 0; k < 3; k++) counts[k] += h[k];
    if (no) CK(cudaMemcpy(out, DOUT, (size_t)no * sizeof(bg_out), cudaMemcpyDeviceToHost));
    return no;
}
extern "C" const char *bg_device(void) { return DEVNAME; }
#endif
