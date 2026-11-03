/* Interface between the host part of the GPU MITM top probe (mitmgpu_host.c: job, tables, enumeration, output as csrc/mitm_top)
 * and its kernel (mitmgpu.cu: CUDA, or with -DMG_EMU the same source compiled for the CPU). See mitmgpu.cu for the method. */
#ifndef MITMGPU_H
#define MITMGPU_H
#ifdef __cplusplus
extern "C" {
#endif

#define MG_MAXJ 40
#define MG_NC 7                    /* digit coefficients c = 9 (d - 1), d in 2 3 4 6 7 8 9 */

typedef struct {                   /* one prefix (all digits but the last) and the range of the last gap, as topgpu.h's tp_prefix */
    unsigned int o;                /* ordering index */
    unsigned short glo, ghi;       /* last gap gl in [glo, ghi] */
    unsigned short g[MG_MAXJ];     /* gaps g_0 .. g_{j-2} */
} mg_prefix;

typedef struct { unsigned int prefix, gl; unsigned long long idx, p[3]; } mg_out;   /* a candidate (prefix, last gap, table
                                   index); dump mode: every case with its P (idx unused) */
typedef struct { unsigned long long q[3], sig; } mg_entry;                       /* bottom-window entry: Q (units 2^-64), digits
                                   (6 bits per value, mitm_bottom's sig) */
typedef struct { unsigned long long key; unsigned int first, count; } mg_slot;   /* hash slot: cell key -> entries
                                   [first, first + count) of the key-sorted table; key ~0 = empty */

typedef struct {
    int j, cls, LW, FL;            /* digits; size class; limbs and fraction bits of the fixed-point values */
    int Smin, smax;                /* = 5 (small-height table, as topgpu); largest height in the tables */
    int dump;                      /* validation: one record per case with P_0, P_1, P_2 (no probe) */
    int hbits;                     /* 2^hbits hash slots */
    int nterm[600];                /* log1p series terms at height S, as topgpu */
    long long o_cs, o_lam1, o_ia1, o_c10, o_inv, o_sl, o_sia;   /* tables as topgpu.h (u32 words) */
    long long o_th[3];             /* word 0: 1 if -Theta_i < 0, then the low L + 1 limbs of round(|Theta_i| 2^64) */
    long long o_sig;               /* per ordering: its digits (6 bits per value), two words, low first */
    unsigned long long hw[3];      /* box half-widths in units 2^-64, margins included (saturated at 2^64 - 1) */
    unsigned long long ncell[3];   /* cells per axis; key (c0 ncell1 + c1) ncell2 + c2 mod 2^64, as mitm_bottom */
    int full[3];                   /* the box covers at least (ncell - 1) / ncell of the axis: probe every cell of it */
    unsigned long long cntsig;     /* the bag's digit counts, 6 bits per value */
    long long words;
} mg_params;

int mg_init(const mg_params *p, const unsigned int *tables, const mg_entry *E, long long nE, const mg_slot *H);
/* Runs n prefixes; counts[0] += cases, counts[1] += table entries examined. Returns the number of records (candidates, or cases
 * in dump mode), or -1 if they exceed outcap (counts untouched: rerun with a larger buffer). */
long long mg_run(const mg_prefix *pre, int n, mg_out *out, long long outcap, unsigned long long counts[2]);
const char *mg_device(void);

#ifdef __cplusplus
}
#endif
#endif
