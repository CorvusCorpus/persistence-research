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
#include <time.h>
#include <omp.h>
#include "dlog5.h"

static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static const int VALS[7] = {2, 3, 4, 6, 7, 8, 9};
static int Gp, tm, k, bmax, cnt[7], botv[128];
static u64 M5, phi, e7;
static long long b1u, b1w, b2u, b2w;
static u64 ncell[3];
static u128 C0[3], AU[3], AW[3], AT[3];

typedef struct { u64 key; u64 q[3]; u64 hcode; u64 dcode; uint32_t t; uint32_t b; u64 sig; } mb_entry;
static mb_entry *E; static size_t nE;

static u128 rd128(const char *s) { u128 x = 0; while (*s >= '0' && *s <= '9') x = x * 10 + (u128)(*s++ - '0'); return x; }

/* The enumeration runs in OpenMP threads: one task per window size b and heights of its top (up to three) digits; each thread
 * appends to its own buffer (tstate). The table is then sorted by a total order (key, then the rest of the entry), so the
 * file does not depend on the threads. */
static int HI[64][64];                                           /* HI[b][r]: the bound a_r < HI for windows of b digits */
typedef struct { mb_entry *E; size_t n, cap; int vals[64], hts[64], rem[7], b; } tstate;
static i128 P10[40]; static u64 INV2T[256];                     /* 10^h; 2^-t mod 5^Gp (t <= tm < 256) */
static void emit(tstate *s) {
    const int b_ = s->b; const int *vals_ = s->vals, *hts_ = s->hts;
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
        if (s->n == s->cap) { s->cap = s->cap ? 2 * s->cap : 1 << 16; s->E = realloc(s->E, s->cap * sizeof *s->E); if (!s->E) { fprintf(stderr, "oom\n"); exit(2); } }
        mb_entry *e = &s->E[s->n++];
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

static void values(tstate *s, int r) {
    if (r == s->b) { emit(s); return; }
    for (int c = 0; c < 7; c++) if (s->rem[c]) { s->rem[c]--; s->vals[r] = VALS[c]; values(s, r + 1); s->rem[c]++; }
}
static int height_ok(int b, int r, int a) { return a < HI[b][r] && !(r == b - 1 && a > 1); }   /* a_k <= 1 */
static void heights(tstate *s, int r, int top) {
    if (r == s->b) { for (int c = 0; c < 7; c++) s->rem[c] = cnt[c]; values(s, 0); return; }
    for (int a = top - 1; a >= 0; a--) if (height_ok(s->b, r, a)) { s->hts[r] = a; heights(s, r + 1, a); }
}

static int cmp_entry(const void *a_, const void *b_) {           /* total order: key, then the rest of the entry */
    const mb_entry *a = a_, *b = b_;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    for (int i = 0; i < 3; i++) if (a->q[i] != b->q[i]) return a->q[i] < b->q[i] ? -1 : 1;
    if (a->hcode != b->hcode) return a->hcode < b->hcode ? -1 : 1;
    if (a->dcode != b->dcode) return a->dcode < b->dcode ? -1 : 1;
    return a->t < b->t ? -1 : (a->t > b->t);
}
static void merge_runs(const mb_entry *a, size_t na, const mb_entry *b, size_t nb, mb_entry *out) {
    size_t i = 0, j = 0, o = 0;
    while (i < na && j < nb) out[o++] = cmp_entry(&a[i], &b[j]) <= 0 ? a[i++] : b[j++];
    while (i < na) out[o++] = a[i++];
    while (j < nb) out[o++] = b[j++];
}
static void psort(mb_entry *A, size_t n) {                       /* chunks sorted in parallel, then merged pairwise */
    int P = omp_get_max_threads();
    if (n < (1 << 16) || P < 2) { qsort(A, n, sizeof *A, cmp_entry); return; }
    size_t *off = malloc((size_t)(P + 1) * sizeof *off); for (int i = 0; i <= P; i++) off[i] = n * (size_t)i / (size_t)P;
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < P; i++) qsort(A + off[i], off[i + 1] - off[i], sizeof *A, cmp_entry);
    mb_entry *tmp = malloc(n * sizeof *A), *src = A, *dst = tmp;
    if (!tmp) { fprintf(stderr, "oom\n"); exit(2); }
    for (int w = 1; w < P; w *= 2) {
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < P; i += 2 * w) {
            const int m = i + w < P ? i + w : P, e = i + 2 * w < P ? i + 2 * w : P;
            merge_runs(src + off[i], off[m] - off[i], src + off[m], off[e] - off[m], dst + off[i]);
        }
        mb_entry *x = src; src = dst; dst = x;
    }
    if (src != A) memcpy(A, src, n * sizeof *A);
    free(tmp); free(off);
}

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
    const double t0 = now_s();
    M5 = 1; for (int i = 0; i < Gp; i++) M5 *= 5;
    if (tm < 0 || tm > 255) { fprintf(stderr, "tm out of range\n"); return 2; }
    dlog_init(Gp); dlog_fast_init(); MBCHECK = getenv("MB_CHECK") && atoi(getenv("MB_CHECK"));
    P10[0] = 1; for (int i = 1; i < 40; i++) P10[i] = P10[i - 1] * 10;
    { u64 inv2 = invmod(2, M5); INV2T[0] = 1 % M5; for (int t = 1; t < 256; t++) INV2T[t] = mulmod(INV2T[t - 1], inv2, M5); }
    if (bmax > 63) { fprintf(stderr, "bmax out of range\n"); return 2; }
    /* tasks: window size b and the heights of its top min(b, 3) digits */
    int ntask = 0, capt = 1024; int (*task)[4] = malloc(capt * sizeof *task);
    for (int b = 0; b <= bmax; b++) {                                       /* b = 0: only when bot[k] > Gp (Gp = 1) */
        if (b < k && botv[k - b] >= 0 && botv[k - b] <= Gp) continue;       /* digit k-b must lie inside the window */
        int ok = 1;
        for (int i = k - b + 1; i <= k; i++) {
            int h = Gp - (i - (k - b + 1));
            for (int ip = k - b + 1; ip <= i; ip++) if (botv[ip] >= 0 && botv[ip] - (i - ip) < h) h = botv[ip] - (i - ip);
            HI[b][i - (k - b + 1)] = h; if (h <= 0) ok = 0;
        }
        if (!ok) continue;
        const int d = b < 3 ? b : 3;
        int pre[3] = {Gp, Gp, Gp}, r = 0;                                   /* odometer over descending prefixes */
        if (d == 0) { task[ntask][0] = 0; ntask++; continue; }
        pre[0] = Gp;
        while (r >= 0) {
            pre[r]--;
            if (pre[r] < 0) { r--; continue; }
            if (!height_ok(b, r, pre[r])) continue;
            if (r < d - 1) { r++; pre[r] = pre[r - 1]; continue; }
            if (ntask == capt) { capt *= 2; task = realloc(task, capt * sizeof *task); }
            task[ntask][0] = b; task[ntask][1] = pre[0]; task[ntask][2] = pre[1]; task[ntask][3] = pre[2]; ntask++;
        }
    }
    const int T = omp_get_max_threads(); tstate *S = calloc((size_t)T, sizeof *S);
    #pragma omp parallel
    {
        tstate *s = &S[omp_get_thread_num()];
        #pragma omp for schedule(dynamic, 1)
        for (int q = 0; q < ntask; q++) {
            s->b = task[q][0]; const int d = s->b < 3 ? s->b : 3;
            for (int r = 0; r < d; r++) s->hts[r] = task[q][1 + r];
            heights(s, d, d ? s->hts[d - 1] : Gp);
        }
    }
    nE = 0; for (int i = 0; i < T; i++) nE += S[i].n;
    E = malloc((nE ? nE : 1) * sizeof *E); if (!E) { fprintf(stderr, "oom\n"); return 2; }
    { size_t o = 0; for (int i = 0; i < T; i++) { if (S[i].n) memcpy(E + o, S[i].E, S[i].n * sizeof *E); o += S[i].n; free(S[i].E); } }
    free(S); free(task);
    const double t1 = now_s();
    psort(E, nE);
    const double t2 = now_s();
    FILE *o = fopen(argv[2], "wb"); if (!o) { perror("out"); return 2; }
    fwrite(&nE, sizeof nE, 1, o); fwrite(E, sizeof *E, nE, o); fclose(o);
    printf("DONE %zu\n", nE);
    if (getenv("MB_TIME")) fprintf(stderr, "MBTIME %zu entries: setup and enumeration %.2f s, sort %.2f s, write %.2f s\n", nE, t1 - t0, t2 - t1, now_s() - t2);
    return 0;
}
