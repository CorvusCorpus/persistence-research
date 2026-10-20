// Exhaustive check of the exact division used by the IFMA-free stage-2 paths (S2_AVX512, S2_AVX2):
//   (1) floor(v * M4 / 2^45) == floor(v / 10^4) for every 0 <= v < 2^31 (M4 = ceil(2^45 / 10^4));
//   (2) the split of the K-digit block xh = h0 + 2^K Y (K = 12, 13; h0 < 2^K; Y < 5^K):
//       xh mod 10^4 == w mod 10^4 and floor(xh / 10^4) == 2^K sq + floor(w / 10^4), where Y = 10^4 sq + r and
//       w = h0 + 2^K r, for every Y < 5^K and h0 in {0, 1, 2^K - 1} plus 2^26 random (Y, h0) pairs.
// Build and run: gcc -O2 -o val_kernel val_kernel.c && ./val_kernel   (prints OK; exit status 1 on failure)
#include <stdio.h>
#include <stdint.h>
typedef uint64_t u64;
#define M4 3518437209ULL
static inline u64 div4(u64 v) { return (v * M4) >> 45; }
static u64 rng = 88172645463325252ULL;
static inline u64 xorshift(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int check_split(int K, u64 Y, u64 h0) {
    u64 xh = h0 + (Y << K), sq = div4(Y), r = Y - 10000 * sq, w = h0 + (r << K), t1 = div4(w);
    return (w - 10000 * t1) == xh % 10000 && ((sq << K) + t1) == xh / 10000 && w < (1ULL << 31);
}
int main(void) {
    u64 bad = 0;
    for (u64 v = 0; v < (1ULL << 31); v++) bad += div4(v) != v / 10000;
    printf("(1) floor(v*M4/2^45) == v/10^4 for all v < 2^31: %s\n", bad ? "FAIL" : "OK");
    int fail = bad != 0;
    for (int K = 12; K <= 13; K++) {
        u64 Q = 1; for (int i = 0; i < K; i++) Q *= 5;
        u64 nb = 0, H[3] = { 0, 1, (1ULL << K) - 1 };
        for (u64 Y = 0; Y < Q; Y++) for (int h = 0; h < 3; h++) nb += !check_split(K, Y, H[h]);
        for (u64 n = 0; n < (1ULL << 26); n++) nb += !check_split(K, xorshift() % Q, xorshift() % (1ULL << K));
        printf("(2) K=%d split of xh = h0 + 2^K Y: all Y < 5^K (h0 = 0, 1, 2^K-1) + 2^26 random pairs: %s\n", K, nb ? "FAIL" : "OK");
        fail |= nb != 0;
    }
    puts(fail ? "KERNEL CHECK FAILED" : "kernel arithmetic: OK");
    return fail;
}
