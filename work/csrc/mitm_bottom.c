/* mitm_bottom: the bottom-window table of the MITM top step (mitm_ref.bottom_entries in C; docs/MITM_STEPS.md section 2).
 * Every window: the lowest b <= bmax digits at heights < Gp (lowest digit at height <= 1, the bottom facts), with their
 * values, b = 0 included (no digit below Gp: only when bot[k] > Gp, i.e. Gp = 1); F = -1 + sum 9(d-1) 10^h; t = v2(F) if
 * v2(F) < Gp, else every t in [Gp, tm]; dl = dlog_3(F 2^-t mod 5^Gp) (digit by digit, dlog5.h; Pohlig-Hellman for its
 * tables and for env MB_CHECK=1); v = (dl, 0) reduced by the basis (b1, b2) of {u + e7 w = 0 mod phi}; the torus point
 * Q = c0 + v_u Au + v_w Aw + t At (mod 1, 128-bit fixed point, constants from mitm_fast.py) and its cell.
 * Output: binary file of entries sorted by cell key (struct mb_entry below); the window as hcode (heights, bit mask) and
 * dcode (3-bit digit indices, at most 21 digits). The earlier code (6 bits per height) overflowed for b > 10 (4 Nov 2026).
 * usage: mitm_bottom JOB OUT            build: see work/Makefile (needs Gp <= 27: 5^Gp < 2^63) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "dlog5.h"

static const int VALS[7] = {2, 3, 4, 6, 7, 8, 9};
static int Gp, tm, k, bmax, cnt[7], botv[128];
static u64 M5, phi, e7;
static long long b1u, b1w, b2u, b2w;
static u64 ncell[3];
static u128 C0[3], AU[3], AW[3], AT[3];

typedef struct { u64 key; u64 q[3]; u64 hcode; u64 dcode; uint32_t t; uint32_t b; u64 sig; } mb_entry;
static mb_entry *E; static size_t nE, capE;

static u128 rd128(const char *s) { u128 x = 0; while (*s >= '0' && *s <= '9') x = x * 10 + (u128)(*s++ - '0'); return x; }

static int hi_[64], vals_[64], hts_[64], b_;
static i128 P10[40]; static u64 INV2T[256];                     /* 10^h; 2^-t mod 5^Gp (t <= tm < 256) */
static void emit(void) {
    i128 F = -1;
    for (int r = 0; r < b_; r++) F += (i128)(9 * (vals_[r] - 1)) * P10[hts_[r]];
    u128 Fu = (u128)F; int v = (u64)Fu ? __builtin_ctzll((u64)Fu) : 64 + __builtin_ctzll((u64)(Fu >> 64));
    int tlo = v < Gp ? v : Gp, thi = v < Gp ? v : tm;
    i128 Fr = F % (i128)M5; u64 Fm = (u64)(Fr < 0 ? Fr + (i128)M5 : Fr);   /* F = -1 for b = 0 */
    u64 sig = 0; for (int r = 0; r < b_; r++) for (int c = 0; c < 7; c++) if (VALS[c] == vals_[r]) sig += (u64)1 << (6 * c);
    /* window codes: hcode = the set of heights as a bit mask (heights < Gp <= 27, read from the top down), dcode = the digits
     * in that order as VALS indices, 3 bits each (at most 21 digits; a wider window stops the program) */
    if (b_ > 21) { fprintf(stderr, "window of %d digits: dcode holds at most 21\n", b_); exit(2); }
    u64 hcode = 0, dcode = 0;
    for (int r = 0; r < b_; r++) { int ci = 0; for (int c = 0; c < 7; c++) if (VALS[c] == vals_[r]) ci = c; hcode |= (u64)1 << hts_[r]; dcode |= (u64)ci << (3 * r); }
    for (int t = tlo; t <= thi; t++) {
        if (t > tm) break;
        u64 x = mulmod(Fm, INV2T[t], M5);
        u64 dl = dlog3(x);
        i128 det = (i128)b1u * b2w - (i128)b1w * b2u;
        long long c1 = rnd_div((i128)dl * b2w, det), c2 = rnd_div(-(i128)b1w * dl, det);
        i128 vu = (i128)dl - (i128)c1 * b1u - (i128)c2 * b2u, vw = -(i128)c1 * b1w - (i128)c2 * b2w;
        if (((vu + (i128)e7 * vw - (i128)dl) % (i128)phi) != 0) { fprintf(stderr, "coset check failed\n"); exit(2); }
        if (nE == capE) { capE = capE ? 2 * capE : 1 << 20; E = realloc(E, capE * sizeof *E); if (!E) { fprintf(stderr, "oom\n"); exit(2); } }
        mb_entry *e = &E[nE++];
        u64 cell[3];
        for (int i = 0; i < 3; i++) {
            u128 q = C0[i] + (u128)(i128)vu * AU[i] + (u128)(i128)vw * AW[i] + (u128)(u64)t * AT[i];
            e->q[i] = (u64)(q >> 64);
            cell[i] = (u64)(((u128)e->q[i] * ncell[i]) >> 64);
        }
        e->key = (cell[0] * ncell[1] + cell[1]) * ncell[2] + cell[2];
        e->hcode = hcode; e->dcode = dcode; e->t = (uint32_t)t; e->b = (uint32_t)b_; e->sig = sig;
    }
}

static int rem_[7];
static void values(int r) {
    if (r == b_) { emit(); return; }
    for (int c = 0; c < 7; c++) if (rem_[c]) { rem_[c]--; vals_[r] = VALS[c]; values(r + 1); rem_[c]++; }
}
static void heights(int r, int top) {
    if (r == b_) { for (int c = 0; c < 7; c++) rem_[c] = cnt[c]; values(0); return; }
    for (int a = top - 1; a >= 0; a--) {
        if (a >= hi_[r]) continue;
        if (r == b_ - 1 && a > 1) continue;                       /* a_k <= 1 */
        hts_[r] = a; heights(r + 1, a);
    }
}

static int cmp_key(const void *a, const void *b) { u64 x = ((const mb_entry *)a)->key, y = ((const mb_entry *)b)->key; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: mitm_bottom JOB OUT\n"); return 2; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("job"); return 2; }
    char w[64], s1[128], s2[128], s3[128];
    for (int i = 0; i < 128; i++) botv[i] = -1;
    while (fscanf(f, "%63s", w) == 1) {
        if (!strcmp(w, "Gp")) fscanf(f, "%d", &Gp);
        else if (!strcmp(w, "tm")) fscanf(f, "%d", &tm);
        else if (!strcmp(w, "k")) fscanf(f, "%d", &k);
        else if (!strcmp(w, "bmax")) fscanf(f, "%d", &bmax);
        else if (!strcmp(w, "cnt")) for (int c = 0; c < 7; c++) fscanf(f, "%d", &cnt[c]);
        else if (!strcmp(w, "bot")) { int i, v; fscanf(f, "%d %d", &i, &v); botv[i] = v; }
        else if (!strcmp(w, "phi")) fscanf(f, "%lu", &phi);
        else if (!strcmp(w, "e7")) fscanf(f, "%lu", &e7);
        else if (!strcmp(w, "b1")) fscanf(f, "%lld %lld", &b1u, &b1w);
        else if (!strcmp(w, "b2")) fscanf(f, "%lld %lld", &b2u, &b2w);
        else if (!strcmp(w, "ncell")) fscanf(f, "%lu %lu %lu", &ncell[0], &ncell[1], &ncell[2]);
        else if (!strcmp(w, "c0") || !strcmp(w, "Au") || !strcmp(w, "Aw") || !strcmp(w, "At")) {
            fscanf(f, "%127s %127s %127s", s1, s2, s3);
            u128 *d = !strcmp(w, "c0") ? C0 : !strcmp(w, "Au") ? AU : !strcmp(w, "Aw") ? AW : AT;
            d[0] = rd128(s1); d[1] = rd128(s2); d[2] = rd128(s3);
        } else { fprintf(stderr, "unknown key %s\n", w); return 2; }
    }
    fclose(f);
    if (Gp < 1 || Gp > 27) { fprintf(stderr, "Gp out of range\n"); return 2; }
    M5 = 1; for (int i = 0; i < Gp; i++) M5 *= 5;
    if (tm < 0 || tm > 255) { fprintf(stderr, "tm out of range\n"); return 2; }
    dlog_init(Gp); dlog_fast_init(); MBCHECK = getenv("MB_CHECK") && atoi(getenv("MB_CHECK"));
    P10[0] = 1; for (int i = 1; i < 40; i++) P10[i] = P10[i - 1] * 10;
    { u64 inv2 = invmod(2, M5); INV2T[0] = 1 % M5; for (int t = 1; t < 256; t++) INV2T[t] = mulmod(INV2T[t - 1], inv2, M5); }
    for (int b = 0; b <= bmax; b++) {                                       /* b = 0: only when bot[k] > Gp (Gp = 1) */
        if (b < k && botv[k - b] >= 0 && botv[k - b] <= Gp) continue;       /* digit k-b must lie inside the window */
        int ok = 1;
        for (int i = k - b + 1; i <= k; i++) {
            int h = Gp - (i - (k - b + 1));
            for (int ip = k - b + 1; ip <= i; ip++) if (botv[ip] >= 0 && botv[ip] - (i - ip) < h) h = botv[ip] - (i - ip);
            hi_[i - (k - b + 1)] = h; if (h <= 0) ok = 0;
        }
        if (!ok) continue;
        b_ = b; heights(0, Gp);
    }
    qsort(E, nE, sizeof *E, cmp_key);
    FILE *o = fopen(argv[2], "wb"); if (!o) { perror("out"); return 2; }
    fwrite(&nE, sizeof nE, 1, o); fwrite(E, sizeof *E, nE, o); fclose(o);
    printf("DONE %zu\n", nE);
    return 0;
}
