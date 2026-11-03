/* mitm_top: top probe of the MITM top step (mitm_ref.mitm_join in C; docs/MITM_STEPS.md section 2), CPU version.
 * Every top configuration (ordering of j digits in prove_cum.ordered's lexicographic order, depths D_0 < ... < D_(j-1) with
 * D_i < r_i): lambda = log(1 + sum c_i 10^-D_i) (MPFR), P_i = frac(Theta_i lambda) (i = 0, 1, 2; top 64 bits), its cell; the
 * 27 neighbouring cells of the bottom table (mitm_bottom's sorted output) are searched, and every entry that is
 * multiset-compatible and within the box |P_i - Q_i| <= hw_i (units of 2^-64, margins included by mitm_fast.py) is printed:
 *   CAND o D_0 .. D_(j-1) idx            (o: ordering index, idx: table index)
 * then DONE ncases nprobes ncand.        usage: mitm_top JOB [threads]       build: see work/Makefile (OpenMP, MPFR) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#include <omp.h>

typedef uint64_t u64; typedef unsigned __int128 u128;
static const int VALS[7] = {2, 3, 4, 6, 7, 8, 9};
static int prec, j, cnt[7], rng[64];
static char *theta_s[3];
static u64 hw[3], ncell[3];

typedef struct { u64 key; u64 q[3]; u64 hcode; u64 dcode; uint32_t t; uint32_t b; u64 sig; } mb_entry;
static mb_entry *T; static size_t nT;

static int nord; static int (*ords)[64];
static int rem_[7], cur_[64], capo;
static void gen(int r) {
    if (r == j) { if (nord == capo) { capo = capo ? 2 * capo : 1024; ords = realloc(ords, capo * sizeof *ords); }
                  memcpy(ords[nord++], cur_, sizeof cur_); return; }
    for (int c = 0; c < 7; c++) if (rem_[c]) { rem_[c]--; cur_[r] = VALS[c]; gen(r + 1); rem_[c]++; }
}

static size_t lower(u64 key) { size_t lo = 0, hi = nT; while (lo < hi) { size_t m = (lo + hi) / 2; if (T[m].key < key) lo = m + 1; else hi = m; } return lo; }

static u64 frac64(mpfr_t x, mpfr_t tmp) {                     /* top 64 bits of frac(x) in [0, 1) */
    mpfr_frac(tmp, x, MPFR_RNDN);
    if (mpfr_sgn(tmp) < 0) mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_mul_2ui(tmp, tmp, 64, MPFR_RNDN);
    mpfr_floor(tmp, tmp);
    mpz_t z; mpz_init(z); mpfr_get_z(z, tmp, MPFR_RNDN);
    u64 r = 0; if (mpz_sgn(z) > 0) { mpz_t m; mpz_init_set_ui(m, 1); mpz_mul_2exp(m, m, 64); mpz_mod(z, z, m); r = mpz_get_ui(z); mpz_clear(m); }
    mpz_clear(z); return r;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: mitm_top JOB [threads]\n"); return 2; }
    if (argc > 2) omp_set_num_threads(atoi(argv[2]));
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("job"); return 2; }
    char w[64], path[4096] = "";
    static char buf[3][100000];
    while (fscanf(f, "%63s", w) == 1) {
        if (!strcmp(w, "prec")) fscanf(f, "%d", &prec);
        else if (!strcmp(w, "j")) fscanf(f, "%d", &j);
        else if (!strcmp(w, "cnt")) for (int c = 0; c < 7; c++) fscanf(f, "%d", &cnt[c]);
        else if (!strcmp(w, "ranges")) for (int i = 0; i < j; i++) fscanf(f, "%d", &rng[i]);
        else if (!strcmp(w, "theta0") || !strcmp(w, "theta1") || !strcmp(w, "theta2")) { int i = w[5] - '0'; fscanf(f, "%99999s", buf[i]); theta_s[i] = buf[i]; }
        else if (!strcmp(w, "hw")) fscanf(f, "%lu %lu %lu", &hw[0], &hw[1], &hw[2]);
        else if (!strcmp(w, "ncell")) fscanf(f, "%lu %lu %lu", &ncell[0], &ncell[1], &ncell[2]);
        else if (!strcmp(w, "table")) fscanf(f, "%4095s", path);
        else { fprintf(stderr, "unknown key %s\n", w); return 2; }
    }
    fclose(f);
    FILE *tb = fopen(path, "rb"); if (!tb) { perror("table"); return 2; }
    if (fread(&nT, sizeof nT, 1, tb) != 1) return 2;
    T = malloc((nT ? nT : 1) * sizeof *T); if (nT && fread(T, sizeof *T, nT, tb) != nT) { fprintf(stderr, "short table\n"); return 2; }
    fclose(tb);
    for (int c = 0; c < 7; c++) rem_[c] = cnt[c];
    gen(0);
    unsigned long long ncase = 0, nprobe = 0, ncand = 0;
    #pragma omp parallel reduction(+:ncase, nprobe, ncand)
    {
        mpfr_t th[3], z, lam, x, tmp, p10;
        for (int i = 0; i < 3; i++) { mpfr_init2(th[i], prec); mpfr_set_str(th[i], theta_s[i], 10, MPFR_RNDN); }
        mpfr_inits2(prec, z, lam, x, tmp, p10, (mpfr_ptr)0);
        char *out = malloc(1 << 20); size_t olen = 0;
        int D[64];
        #pragma omp for schedule(dynamic, 1)
        for (int o = 0; o < nord; o++) {
            u64 sigT = 0; for (int i = 0; i < j; i++) for (int c = 0; c < 7; c++) if (VALS[c] == ords[o][i]) sigT += (u64)1 << (6 * c);
            /* depth compositions: D_0 >= 1, strictly increasing, D_i < rng[i] */
            int r = 0; D[0] = 0;
            while (r >= 0) {
                D[r]++;
                if (D[r] >= rng[r]) { r--; continue; }
                if (r < j - 1) { r++; D[r] = D[r - 1]; continue; }
                /* a complete configuration */
                ncase++;
                mpfr_set_ui(z, 0, MPFR_RNDN);
                for (int i = 0; i < j; i++) {
                    mpfr_set_ui(p10, 10, MPFR_RNDN); mpfr_pow_ui(p10, p10, (unsigned long)D[i], MPFR_RNDN);
                    mpfr_ui_div(p10, 9 * (ords[o][i] - 1), p10, MPFR_RNDN); mpfr_add(z, z, p10, MPFR_RNDN);
                }
                mpfr_log1p(lam, z, MPFR_RNDN);
                u64 P[3], cell[3];
                for (int i = 0; i < 3; i++) { mpfr_mul(x, lam, th[i], MPFR_RNDN); P[i] = frac64(x, tmp); cell[i] = (u64)(((u128)P[i] * ncell[i]) >> 64); }
                u64 seen[27]; int ns = 0;
                for (int d0 = -1; d0 <= 1; d0++) for (int d1 = -1; d1 <= 1; d1++) for (int d2 = -1; d2 <= 1; d2++) {
                    u64 c0 = (cell[0] + ncell[0] + d0) % ncell[0], c1 = (cell[1] + ncell[1] + d1) % ncell[1], c2 = (cell[2] + ncell[2] + d2) % ncell[2];
                    u64 key = (c0 * ncell[1] + c1) * ncell[2] + c2; int dup = 0;
                    for (int s = 0; s < ns; s++) if (seen[s] == key) dup = 1;
                    if (dup) continue;
                    seen[ns++] = key;
                    for (size_t idx = lower(key); idx < nT && T[idx].key == key; idx++) {
                        nprobe++;
                        int okc = 1; for (int c = 0; c < 7 && okc; c++) if ((int)((sigT >> (6 * c)) & 63) + (int)((T[idx].sig >> (6 * c)) & 63) > cnt[c]) okc = 0;
                        if (!okc) continue;
                        int inbox = 1;
                        for (int i = 0; i < 3 && inbox; i++) { u64 dd = P[i] - T[idx].q[i]; u64 ad = (int64_t)dd < 0 ? (u64)0 - dd : dd; if (ad > hw[i]) inbox = 0; }
                        if (!inbox) continue;
                        ncand++;
                        if (olen > (1 << 20) - 4096) {
                            #pragma omp critical
                            { fwrite(out, 1, olen, stdout); }
                            olen = 0;
                        }
                        olen += sprintf(out + olen, "CAND %d", o);
                        for (int i = 0; i < j; i++) olen += sprintf(out + olen, " %d", D[i]);
                        olen += sprintf(out + olen, " %zu\n", idx);
                    }
                }
            }
        }
        #pragma omp critical
        { fwrite(out, 1, olen, stdout); }
        free(out);
        for (int i = 0; i < 3; i++) mpfr_clear(th[i]);
        mpfr_clears(z, lam, x, tmp, p10, (mpfr_ptr)0);
    }
    printf("DONE %llu %llu %llu\n", ncase, nprobe, ncand);
    return 0;
}
