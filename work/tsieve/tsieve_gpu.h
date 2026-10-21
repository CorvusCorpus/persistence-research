/* Interface between the host part of the GPU tree sieve (tsieve_gpu.c, C) and its kernels (tsieve_gpu.cu, CUDA or,
 * with -DTG_EMU, the same source compiled as C++ for the CPU). */
#ifndef TSIEVE_GPU_H
#define TSIEVE_GPU_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {                    /* one coset of L_K: constants computed on the host (see tsieve6.c, leaf()) */
    unsigned long long r, h0, d0, d1, d2, ybase;    /* ybase = A_hi: Y = ybase + i0 d0 + i1 d1 + i2 d2 (mod 5^K) */
    unsigned int s0, si0, si1, inv2;                /* s = s0 + i0 si0 + i1 si1 + i2 (mod 5^M): position in U; inv2 = delta2^-1 */
    int i0lo, i0hi;                                 /* exact range of i0 over the region (empty cosets are skipped) */
} tg_coset;

typedef struct { unsigned int a, b, c, code; } tg_point;
enum { TG_R3 = 3, TG_R4 = 4, TG_S4 = 5, TG_CORNER = 9 };     /* decisions; CORNER: the host decides (x near 10^XMIN) */

typedef struct {
    int k, m, dm, nl, nt, n3, n2, skip, nfull, nrem;
    unsigned long long m1, mm, br1, brm, nleaf;     /* br1 = floor(2^64 / 5^K), brm = floor(2^64 / 5^M) */
    unsigned long long p2k[16], p5[16];             /* 2^(K-i), 5^i */
    long long v[3][3];
    long long a3[8][3], b3_0[8], b3_r[8];           /* region: a3 . (i0, i1, i2) <= b3_0 + b3_r r */
    long long a2[32][2], b2_0[32], b2_r[32];        /* its projection on (i0, i1) */
    double cornerlg;                                /* points with b lg3 + c lg7 below this go to the host */
    unsigned long long l2, l3, l7; int tb;          /* leading-digit table: bin = (a l2 + b l3 + c l7) >> (64 - tb) */
    unsigned char ok4[10000], ok1[10], okr[10000];  /* 4 digits / 1 digit / the top NREM digits admissible */
} tg_params;

typedef struct tg_ctx tg_ctx;
int tg_init(const tg_params *p, const unsigned int *T2, const unsigned int *T3, const unsigned int *T7,
            unsigned long long n2, unsigned long long n3, unsigned long long n7, const unsigned long long *TOP, unsigned long long topwords);
tg_ctx *tg_ctx_new(int maxcos, unsigned long long outcap);
/* Runs n cosets. cnt[0..5] += s1 points, rows, hits, s2, s3, s4 (non-corner points); out gets the S4 and CORNER points,
 * and with dump every s2 point with its decision. Returns 0, or nonzero if out overflowed (rerun with a larger outcap). */
int tg_run(tg_ctx *x, const tg_coset *cos, int n, int dump, tg_point **out, unsigned long long *nout, unsigned long long cnt[6]);
const char *tg_device(void);
/* Fail-stop check for non-ECC memory: the device tables and constants must equal the host copies, and the host copies
 * must still have the hash they had at tg_init. Returns 0 if everything is unchanged. */
int tg_verify(void);

#ifdef __cplusplus
}
#endif
#endif
