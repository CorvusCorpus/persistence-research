/* Interface between the host part of the GPU bottom solver (bottomgpu_host.c, C, includes bottom3.c) and its kernel
 * (bottomgpu.cu: CUDA, or with -DBG_EMU the same source compiled for the CPU). See bottomgpu.cu for the method. */
#ifndef BOTTOMGPU_H
#define BOTTOMGPU_H
#ifdef __cplusplus
extern "C" {
#endif

#define BG_MAXB 16                 /* digits in a bottom configuration (k - j) */
#define BG_NC 7                    /* digit coefficients c = 9 (d - 1), d in 2 3 4 6 7 8 9 */

typedef struct {                   /* one lower configuration (Flow) and the range of heights p of the top digit */
    unsigned int o;                /* ordering index (the digits come from the ordering table) */
    unsigned char pos[BG_MAXB];    /* pos[r], r = 1 .. b-1: heights of digits d_{j+1+r} (pos[0] unused) */
    unsigned char lo, phi;         /* p in [lo, phi) */
    unsigned char pad[2];
} bg_flow;

/* kind 1: a candidate (seq = 4 w0 + 2 di + dj, bottom3.c's print order within a solve); kind 2: the host must redo this
 * solve (counted on the device; drop its kind-1 records); kind 3: the host must redo and count this solve */
typedef struct { unsigned int flow, p, t, kind, seq; unsigned int up[32], w[32]; } bg_out;

/* table offsets (in u32 words) inside one table buffer, and the sizes; filled by the host (bottomgpu_host.c) */
typedef struct {
    int G, tm, b, cls;             /* cls: size class, selects the kernel instantiation (L5, L2, LQ, LH) */
    int L5, L2, LQ, LH, F;         /* limbs; F = fractional bits of the lattice estimates */
    int s5, vz5, s2, vz2;          /* plog parameters at precision G: s = isqrt(G) + 2, vz = s + 1 (p = 5), s + 2 (p = 2) */
    int np5, np2;                  /* terms of the plog series (p = 5, 2) */
    int tlmax;                     /* largest v2(Flow) the kernel handles (larger: the host) */
    unsigned int mpN, mpG, mpQ;    /* -m^-1 mod 2^32 for m = 5^N, 5^G, 5^(G-1) */
    long long o_mN, o_R2N, o_mG, o_R2G, o_twoRG, o_mQ, o_R2Q, o_inv5, o_invz5;
    long long o_T10N, o_T10G, o_T10_2, o_Tc5, o_Tc2;              /* [ci][pos] */
    long long o_coef5, o_n5, o_off5, o_coef2, o_n2, o_off2;       /* series coefficients and their counts / offsets */
    long long o_plog5, o_plog2, o_L16;
    long long o_Ainv, o_A5inv, o_q2inv5, o_a2m4, o_alphaW, o_Q, o_b, o_Ri, o_Rj, o_Ki, o_Kj, o_U, o_W;   /* per t */
    long long o_tok;               /* per t: 1 if the lattice data exist (bottom3.c's mt(G, t) would not stop) */
    long long words;               /* total table size */
    int cval[BG_NC], v2c[BG_NC];   /* c and v2(c) per coefficient index */
} bg_params;

typedef struct {                   /* kernel test: one solve with an arbitrary Flow (residues computed by the host) */
    unsigned int f5n[32], f5g[32], f2[32];          /* Flow mod 5^N, mod 5^G, mod 2^(32 L2) (L5 / L2 limbs used) */
    int f5, ci, p, pad;                             /* Flow mod 5, coefficient index of the top digit, its height */
} bg_query;

int bg_init(const bg_params *p, const unsigned int *tables, const unsigned char *orders_ci, int norders);
/* kernel test: n queries (flow index of a record = query index) */
long long bg_runq(const bg_query *q, int n, bg_out *out, long long outcap, unsigned long long counts[3]);
/* Runs n flows. Returns the number of output records (candidates, and solves the host must redo); counts += nsolve,
 * ndirect5, ndirect2 of the flows' solves done on the device. Records beyond outcap: returns -1 (rerun with fewer flows). */
long long bg_run(const bg_flow *flows, int n, bg_out *out, long long outcap, unsigned long long counts[3]);
const char *bg_device(void);

#ifdef __cplusplus
}
#endif
#endif
