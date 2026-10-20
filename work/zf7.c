/*
 * zf7.c -- exhaustive search for zero-free 7-smooth numbers.
 *
 * Every number n with multiplicative persistence >= 2 maps to a digit product
 * P(n) that is 7-smooth.  A 7-smooth number divisible by 10 ends in 0, so the
 * only 7-smooth numbers that can appear inside a long chain are
 *      family A:  2^a 3^b 7^c          (a,b,c >= 0)
 *      family B:  5^a 3^b 7^c          (a >= 1)
 * and they must contain no digit 0.
 *
 * This program enumerates every member of both families with
 * 10^(K+1) <= x < 10^D and prints those whose lowest K decimal digits are all
 * nonzero.  All filtering is exact modular arithmetic on the low digits
 * (x mod 10^36 kept incrementally, x mod 10^72 and 10^144 recomputed from
 * exact tables), so no candidate is ever rejected unless it really contains
 * a 0.  Numbers below 10^(K+1) are handled exactly by verify.py.
 *
 * usage: zf7 D K        (K in {36,72,144})
 * output: lines "A a b c" / "B a b c" for survivors; stats on stderr.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define E9  1000000000ULL
#define E18 1000000000000000000ULL
#define NL  16                    /* 16 limbs of base 1e9 = 144 digits */

static uint8_t zf6[1000000];      /* zf6[v] = 1 iff v, zero-padded to 6 digits, has no 0 */
#define ZF6(v) (zf6[(v)])

static inline int zf9(uint64_t v) {          /* v < 1e9, 9 digits zero-padded */
    return ZF6(v % 1000000) & ZF6(111000 + v / 1000000);
}
static inline int zf18(uint64_t v) {         /* v < 1e18, 18 digits zero-padded (branchless) */
    uint64_t q = v / 1000000;
    return ZF6(v - q * 1000000) & ZF6(q % 1000000) & ZF6(q / 1000000);
}

/* R = A*B mod 10^(9n), limbs base 1e9, little endian. Column sums < 16e18 < 2^64. */
static void mulmod(const uint64_t *A, const uint64_t *B, uint64_t *R, int n) {
    uint64_t acc[NL] = {0};
    for (int i = 0; i < n; i++) {
        uint64_t ai = A[i];
        if (!ai) continue;
        for (int j = 0; j < n - i; j++) acc[i + j] += ai * B[j];
    }
    uint64_t carry = 0;
    for (int k = 0; k < n; k++) { uint64_t v = acc[k] + carry; R[k] = v % E9; carry = v / E9; }
}

static uint64_t *powtab(int base, long cnt) {      /* base^i mod 10^144, i < cnt */
    uint64_t *t = calloc((size_t)cnt * NL, sizeof(uint64_t));
    t[0] = 1;
    for (long i = 1; i < cnt; i++) {
        uint64_t carry = 0;
        for (int k = 0; k < NL; k++) {
            uint64_t v = t[(i - 1) * NL + k] * base + carry;
            t[i * NL + k] = v % E9; carry = v / E9;
        }
    }
    return t;
}

typedef struct { unsigned long long cand, s1, s2, out; } stats;

/* inner multiplier m (2 or 3) over exponent i; rows over (e1 >= e1min, e2 >= 0) */
static void run_family(char fam, int m, double lin, int f1, double lf1, int e1min,
                       double D, int K, stats *st)
{
    const double l7 = log10(7.0);
    const double T = K + 0.5;
    long nin = (long)(D / lin) + 3, n1 = (long)(D / lf1) + 3, n2 = (long)(D / l7) + 3;
    uint64_t *Pin = powtab(m, nin), *P1 = powtab(f1, n1), *P2 = powtab(7, n2);
    unsigned long long cand = 0, s1 = 0, s2 = 0, out = 0;

    #pragma omp parallel for schedule(dynamic, 1) reduction(+:cand, s1, s2, out)
    for (long e2 = 0; e2 < n2 - 1; e2++) {
        if (e2 * l7 >= D + 1) continue;
        uint64_t R[NL], X[NL];
        for (long e1 = e1min; e1 < n1 - 1 && e1 * lf1 + e2 * l7 < D + 1; e1++) {
            double L = e1 * lf1 + e2 * l7;
            long ilo = (long)ceil((T - L) / lin); if (ilo < 0) ilo = 0;
            long ihi = (long)floor((D - L) / lin) + 1; if (ihi > nin - 1) ihi = nin - 1;
            if (ilo > ihi) continue;
            mulmod(P1 + e1 * NL, P2 + e2 * NL, R, NL);
            mulmod(R, Pin + ilo * NL, X, NL);
            uint64_t l0 = X[0] + X[1] * E9, l1 = X[2] + X[3] * E9;   /* x mod 10^36 */
            for (long i = ilo; i <= ihi; i++) {
                cand++;
                if (__builtin_expect(zf18(l0) & zf18(l1), 0)) {
                    s1++;
                    int ok = 1;
                    if (K > 36) {
                        uint64_t Y[NL];
                        mulmod(R, Pin + i * NL, Y, 8);
                        if (Y[0] + Y[1] * E9 != l0 || Y[2] + Y[3] * E9 != l1) {
                            fprintf(stderr, "CONSISTENCY FAILURE %c %ld %ld %ld\n", fam, e1, e2, i);
                            exit(2);
                        }
                        for (int k = 4; k < 8 && ok; k++) ok = zf9(Y[k]);
                        if (ok) {
                            s2++;
                            if (K > 72) {
                                mulmod(R, Pin + i * NL, Y, NL);
                                for (int k = 8; k < NL && ok; k++) ok = zf9(Y[k]);
                            }
                        }
                    }
                    if (ok) {
                        out++;
                        long a, b, c = e2;
                        if (fam == 'A') { a = i; b = e1; } else { a = e1; b = i; }
                        #pragma omp critical
                        printf("%c %ld %ld %ld\n", fam, a, b, c);
                    }
                }
                /* advance: (l1,l0) *= m mod 10^36 */
                if (m == 2) {
                    uint64_t v = l0 << 1; uint64_t c0 = v >= E18; l0 = v - (c0 ? E18 : 0);
                    uint64_t w = (l1 << 1) + c0; l1 = w >= E18 ? w - E18 : w;
                } else {
                    uint64_t v = l0 * 3; uint64_t c0 = (v >= E18) + (v >= 2 * E18); l0 = v - c0 * E18;
                    uint64_t w = l1 * 3 + c0; uint64_t c1 = (w >= E18) + (w >= 2 * E18); l1 = w - c1 * E18;
                }
            }
        }
    }
    st->cand += cand; st->s1 += s1; st->s2 += s2; st->out += out;
    free(Pin); free(P1); free(P2);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: zf7 D K\n"); return 1; }
    double D = atof(argv[1]); int K = atoi(argv[2]);
    if (K != 36 && K != 72 && K != 144) { fprintf(stderr, "K must be 36, 72 or 144\n"); return 1; }
    for (int v = 0; v < 1000000; v++) {
        int ok = 1, x = v;
        for (int d = 0; d < 6; d++) { if (x % 10 == 0) ok = 0; x /= 10; }
        zf6[v] = ok;
    }
    stats st = {0};
    double t0 = omp_get_wtime();
    run_family('A', 2, log10(2.0), 3, log10(3.0), 0, D, K, &st);
    run_family('B', 3, log10(3.0), 5, log10(5.0), 1, D, K, &st);
    fflush(stdout);
    fprintf(stderr, "D=%g K=%d candidates=%llu pass36=%llu pass72=%llu out=%llu time=%.1fs threads=%d\n",
            D, K, st.cand, st.s1, st.s2, st.out, omp_get_wtime() - t0, omp_get_max_threads());
    return 0;
}
