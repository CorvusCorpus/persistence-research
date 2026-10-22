/* Interface between the host part of the GPU top filter (topgpu_host.c: csrc/top3.c's logic, exact, MPFR) and its kernel
 * (topgpu.cu: CUDA, or with -DTG_EMU... -DTP_EMU the same source compiled for the CPU). See topgpu.cu for the method. */
#ifndef TOPGPU_H
#define TOPGPU_H
#ifdef __cplusplus
extern "C" {
#endif

#define TP_MAXJ 40
#define TP_NC 7                    /* digit coefficients c = 9 (d - 1), d in 2 3 4 6 7 8 9 */

typedef struct {                   /* one prefix (all digits but the last) and the range of the last gap */
    unsigned int o;                /* ordering index */
    unsigned short glo, ghi;       /* last gap gl in [glo, ghi] */
    unsigned short g[TP_MAXJ];     /* gaps g_0 .. g_{j-2} (16 bits: gaps reach G ~ 411; 8 bits wrapped at 256) */
} tp_prefix;

/* a case the host decides exactly (top3.c's code): kind 1 = not certainly killed (or uncertain), kind 2 = the whole prefix
 * (small heights: the series would converge too slowly) */
typedef struct { unsigned int prefix, gl, kind, pad; unsigned long long psi1, psi2; } tp_out;   /* psi: dump mode (kind 3) */

typedef struct {
    int j, cls, LW, LT;            /* digits; size class (kernel instantiation); limbs of the fixed-point values / of Theta */
    int FL, FT;                    /* fraction bits: lambda, 1/A, c 10^-S (FL); Theta (FT) */
    int Smin;                      /* heights below Smin (levels >= 2 and the last digit) go to the host */
    int smax;                      /* largest height in the tables */
    unsigned long long rho_hi, rho_lo;   /* level 1: certainly killed if dist > rho_hi, certainly passes if dist < rho_lo (2^-64) */
    unsigned long long dy, ds;     /* level 2 margins on dist(y) and |s| (2^-64) */
    long long mu63;                /* mu = P/Q in units 2^-63 (|mu| <= 1/2) */
    unsigned long long g2l, g3l;   /* lower bounds of g2/R2, g3/R2 with 32 fraction bits */
    int dump;                      /* validation: one kind-3 record with psi1, psi2 per case (no classification) */
    int nterm[600];                /* log1p / inverse series terms needed at height S (|w| <= 72 10^-S) */
    long long o_cs, o_lam1, o_ia1, o_c10, o_inv, o_th, o_th2;   /* table offsets (u32 words): c-index per (ordering, digit),
                                   packed 4 per word; level-1 log and inverse per (c, S); c 10^-S per (c, S); 1/n; Theta, Theta2
                                   (word 0: 1 if negative, then the low L + 1 limbs of round(|Theta| 2^64)) */
    long long words;
} tp_params;

int tp_init(const tp_params *p, const unsigned int *tables);
/* Runs n prefixes; counts[0] += cases, counts[1] += certain level-1 passes among the cases not forwarded. Returns the number of
 * records (forwarded cases / prefixes), or -1 if they exceed outcap (counts untouched: rerun with a larger buffer). */
long long tp_run(const tp_prefix *pre, int n, tp_out *out, long long outcap, unsigned long long counts[2]);
const char *tp_device(void);

#ifdef __cplusplus
}
#endif
#endif
