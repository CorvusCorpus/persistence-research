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
// Stages 3-4 on the GPU (GSIEVE_GPU34=1): per segment the row data that later_stages() uses. The exponent of 2 is
// a = position + delta, tb = c L7 + b L3 (mod 2^64) and M = 3^b 7^c mod 10^(9 NL) (row mrow of the M array).
typedef struct { long long delta; unsigned long long tb; uint32_t mrow, pad; } gsegx_t;
#define GPU34_CODE(p) ((p) >> 30)                          // stage-3/4 class of an output entry: 0 R3, 1 R4, 2 S4
#define GPU34_POS(p) ((p) & 0x3FFFFFFFu)                   // its bit position (P1 < 2^30)

typedef struct {                                          // one per worker thread: its own stream and buffers
    void *stream;
    uint64_t segcap, poscap, ovfcap;
    gseg_t *h_seg; uint32_t *h_s1, *h_off, *h_cnt, *h_pos; govf_t *h_ovf; unsigned long long *h_ctr;   // pinned
    gseg_t *d_seg; uint32_t *d_s1, *d_off, *d_cnt, *d_pos; govf_t *d_ovf; unsigned long long *d_ctr;   // device
    // stages 3-4 on the GPU: per segment gsegx_t and counts s2, s3, s4; M rows (NL limbs each, segcap rows)
    gsegx_t *h_segx, *d_segx; uint32_t *h_m, *d_m, *h_s2, *d_s2, *h_s3, *d_s3, *h_s4, *d_s4;
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
// Stages 3-4 on the GPU: upload the leading-digit table TOP (ntopw 64-bit words, TB = 18) and the limb table T2
// (nrows rows of nl base-1e9 limbs, T2[a] = 2^a mod 10^(9 nl)); both are included in gpu_verify_tables().
void gpu_init_later(const uint64_t *TOP, uint64_t ntopw, const uint32_t *T2, uint64_t nrows, int nl);
gctx_t *gpu_ctx_new(uint64_t segcap, uint64_t poscap);
void gpu_ctx_later(gctx_t *g);                            // allocate the stages 3-4 buffers (h_segx, h_m, counts)
// Run the segments h_seg[0..nseg). On return: h_s1[k], and the output entries of segment k are the chunks of h_ovf
// (entries with seg == k, in increasing seq) followed by h_pos[h_off[k] .. h_off[k] + h_cnt[k]). h_ctr[0] = total
// entries, h_ctr[1] = number of h_ovf entries. Buffers grow (and the kernel is re-run) when a piece has more
// entries than they hold, so nothing is ever dropped.
//   later = 0: stages 1-2; the entries are the stage-2 survivors' bit positions.
//   later = 1: also stages 3-4 (needs gpu_init_later, h_segx[0..nseg) and nm M rows in h_m). h_s2, h_s3, h_s4 count
//     the stage-2, -3, -4 survivors of each segment; the entries carry their class in the top two bits (GPU34_CODE)
//     and are all stage-2 survivors if dumpall, else only the stage-4 survivors (class 2).
void gpu_run(gctx_t *g, uint64_t nseg, uint64_t nm, int later, int dumpall);
void gpu_ctx_free(gctx_t *g);
#ifdef __cplusplus
}
#endif
#endif
