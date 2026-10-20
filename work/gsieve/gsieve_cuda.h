// gsieve_cuda.h -- interface between gsieve.c (host side of the S2_CUDA variant) and gsieve_cuda.cu (stages 1-2
// of the main path on the GPU). See docs/CHAIN_SIEVE.md, "CUDA variant".
//
// The host plans each piece of work (rows b0 <= b < b1 of one c) exactly as main_c does and hands the GPU one
// segment per row piece: a bitmap range [s, e) of exponent residues mod P1 and the lift index j (lemma L3). The
// GPU returns, per segment, s1 = the number of stage-1 survivors (set bits of the range) and the bit positions of
// the stage-2 survivors in increasing order -- the contents of PBUF in the CPU variants. Stages 3+ stay on the CPU.
#ifndef GSIEVE_CUDA_H
#define GSIEVE_CUDA_H
#include <stdint.h>

typedef struct { uint32_t s, e, j; } gseg_t;               // segment: bit range [s, e), 0 <= s < e <= P1; j < 5^K
typedef struct { uint32_t seg, seq, off, cnt; } govf_t;    // an early chunk of a segment's survivors (rare)

typedef struct {                                          // one per worker thread: its own stream and buffers
    void *stream;
    uint64_t segcap, poscap, ovfcap;
    gseg_t *h_seg; uint32_t *h_s1, *h_off, *h_cnt, *h_pos; govf_t *h_ovf; unsigned long long *h_ctr;   // pinned
    gseg_t *d_seg; uint32_t *d_s1, *d_off, *d_cnt, *d_pos; govf_t *d_ovf; unsigned long long *d_ctr;   // device
} gctx_t;

#ifdef __cplusplus
extern "C" {
#endif
// Upload the tables (after build_bitmap) and check the device copies byte for byte against the host copies.
// ok4b[v] = 1 if the 4-digit window v (zero padded) is allowed; topmode: 0 = G (top digit not 0 or 5),
// 1 = E (top digit not 0), 2 = P (no check). Prints the device name and memory use to stderr.
void gpu_init(const uint64_t *BM, const uint32_t *RK, uint64_t nw, const uint64_t *AB, const uint16_t *H0,
              const uint32_t *POS, uint64_t ngood, const uint8_t *ok4b, int topmode);
// Compare every device table with its host copy (which must stay allocated); 0 if all are equal. GeForce memory has
// no ECC: gsieve.c runs this periodically and writes results only after a check that started after they were done.
int gpu_verify_tables(void);
gctx_t *gpu_ctx_new(uint64_t segcap, uint64_t poscap);
// Run stages 1-2 on the segments h_seg[0..nseg). On return: h_s1[k], and the survivors of segment k are the
// chunks of h_ovf (entries with seg == k, in increasing seq) followed by h_pos[h_off[k] .. h_off[k] + h_cnt[k]).
// h_ctr[0] = total survivors, h_ctr[1] = number of h_ovf entries. Buffers grow (and the kernel is re-run) when
// a piece has more survivors than they hold, so nothing is ever dropped.
void gpu_run(gctx_t *g, uint64_t nseg);
void gpu_ctx_free(gctx_t *g);
#ifdef __cplusplus
}
#endif
#endif
