// z2sieve_cuda.h: interface between z2sieve.c (host) and z2sieve_cuda.cu (stages 2-4 on the GPU)
#ifndef Z2SIEVE_CUDA_H
#define Z2SIEVE_CUDA_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { uint64_t C2K, UU, LLO, LHI; } z2c_t;          // run-time constants (see z2sieve.c)
typedef struct {                                                // results of one chunk
    uint64_t ctr[4];                                            // s1, s2, s3, s4
    uint64_t *s4; uint64_t ns4;                                 // n that passed stage 4 (unordered)
    uint64_t *dmp; uint64_t ndmp;                               // dump: n << 2 | class (1 = R3, 2 = R4, 3 = S4)
    uint32_t flags;                                             // device consistency failures (0 = none)
} z2r_t;
// Upload the tables: zt = the digit bit tables ZT[0..5] concatenated (ZTOT words), top = stage-3 table, b1/rk = the
// stage-1 bitmap over r and its word ranks, rsw/rss/rsr = w, 2^r mod 5^(2K) and 2^r mod 10^(9 NL) at r = 1024 k,
// pw13/ps/pw = 2^g mod 10^K, mod 5^(2K) and mod 10^(9 NL) for g < 1024. Every table keeps its host copy.
void z2g_init(const z2c_t *c, const uint64_t *zt, const uint64_t *top, uint64_t topw, const uint64_t *b1,
              const uint32_t *rk, uint64_t nw, const uint64_t *rsw, const uint64_t *rss, const uint32_t *rsr,
              uint64_t nrs, const uint64_t *pw13, const uint64_t *ps, const uint32_t *pw);
// One chunk: Q = 2^(jP) mod 10^(9 NL) for j in [ja, ja + JCH), survivors [iv[2k], iv[2k+1]) of the bitmap
void z2g_chunk(uint64_t ja, uint64_t jb, uint64_t n0, uint64_t n1, const uint32_t *q, const uint64_t *iv, int niv,
               int dump, z2r_t *res);
int z2g_verify(void);                                           // 0 if every device table equals its host copy
#ifdef __cplusplus
}
#endif
#endif
