/* mitm_bottom: the bottom-window table of the MITM top step (mitm_ref.bottom_entries in C; docs/MITM_STEPS.md section 2).
 * Every window: the lowest b <= bmax digits at heights < Gp (lowest digit at height <= 1, the bottom facts), with their
 * values; F = -1 + sum 9(d-1) 10^h; t = v2(F) if v2(F) < Gp, else every t in [Gp, tm]; dl = dlog_3(F 2^-t mod 5^Gp)
 * (Pohlig-Hellman); v = (dl, 0) reduced by the basis (b1, b2) of {u + e7 w = 0 mod phi}; the torus point
 * Q = c0 + v_u Au + v_w Aw + t At (mod 1, 128-bit fixed point, constants from mitm_fast.py) and its cell.
 * Output: binary file of entries sorted by cell key (struct mb_entry below).
 * usage: mitm_bottom JOB OUT            build: see work/Makefile (needs Gp <= 27: 5^Gp < 2^63) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;

static const int VALS[7] = {2, 3, 4, 6, 7, 8, 9};
static int Gp, tm, k, bmax, cnt[7], botv[128];
static u64 M5, phi, e7;
static long long b1u, b1w, b2u, b2w;
static u64 ncell[3];
static u128 C0[3], AU[3], AW[3], AT[3];

typedef struct { u64 key; u64 q[3]; u64 hcode; u64 dcode; uint32_t t; uint32_t b; u64 sig; } mb_entry;
static mb_entry *E; static size_t nE, capE;

static u128 rd128(const char *s) { u128 x = 0; while (*s >= '0' && *s <= '9') x = x * 10 + (u128)(*s++ - '0'); return x; }
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static u64 invmod(u64 a, u64 m) {                              /* m = 5^n, a coprime to 5 */
    i128 t = 0, nt = 1, r = m, nr = a % m;
    while (nr) { i128 q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    if (t < 0) t += m;
    return (u64)t;
}

/* Pohlig-Hellman: e in [0, 4 5^(n-1)) with 3^e = x mod 5^n */
static u64 P5[40], G4inv[40][5], GM[5], T4[4];
static int NN;
static void dlog_init(int n) {
    NN = n; P5[0] = 1; for (int i = 1; i < 40; i++) P5[i] = P5[i - 1] * 5;
    u64 M = P5[n];
    for (int e = 0; e < 4; e++) T4[e] = powmod(3, (u64)e * P5[n - 1], M);
    u64 g = powmod(3, 4, M), gi = invmod(g, M); int m = n - 1;
    if (m > 0) { u64 gm = powmod(g, P5[m - 1], M); for (int d = 0; d < 5; d++) GM[d] = powmod(gm, d, M); }
    for (int i = 0; i < m; i++) for (int d = 0; d < 5; d++) G4inv[i][d] = powmod(gi, (u64)d * P5[i], M);
}
static u64 dlog3(u64 x) {
    int n = NN, m = n - 1; u64 M = P5[n];
    u64 xp = powmod(x, P5[n - 1], M); int e4 = -1;
    for (int e = 0; e < 4; e++) if (T4[e] == xp) e4 = e;
    if (e4 < 0) { fprintf(stderr, "dlog: no order-4 part\n"); exit(2); }
    u64 h = powmod(x, 4, M), e5 = 0;
    for (int i = 0; i < m; i++) {
        u64 hi = h; for (int s = 0; s < m - 1 - i; s++) { u64 x2 = mulmod(hi, hi, M), x4 = mulmod(x2, x2, M); hi = mulmod(x4, hi, M); }
        int d = -1; for (int dd = 0; dd < 5; dd++) if (GM[dd] == hi) d = dd;
        if (d < 0) { fprintf(stderr, "dlog: digit not found\n"); exit(2); }
        e5 += (u64)d * P5[i]; h = mulmod(h, G4inv[i][d], M);
    }
    /* e = e5 mod 5^m, e = e4 mod 4 */
    u64 e = e5; while (e % 4 != (u64)e4) e += P5[m];
    return e;
}

static long long rnd_div(i128 n, i128 d) {                     /* nearest integer to n/d (d != 0), ties up */
    if (d < 0) { n = -n; d = -d; }
    i128 q = n >= 0 ? (2 * n + d) / (2 * d) : -((-2 * n + d - 1) / (2 * d));
    return (long long)q;
}

static int hi_[64], vals_[64], hts_[64], b_;
static void emit(void) {
    i128 F = -1; i128 p10[40]; p10[0] = 1; for (int i = 1; i < 40; i++) p10[i] = p10[i - 1] * 10;
    for (int r = 0; r < b_; r++) F += (i128)(9 * (vals_[r] - 1)) * p10[hts_[r]];
    u128 Fu = (u128)F; int v = (u64)Fu ? __builtin_ctzll((u64)Fu) : 64 + __builtin_ctzll((u64)(Fu >> 64));
    int tlo = v < Gp ? v : Gp, thi = v < Gp ? v : tm;
    u64 Fm = (u64)(F % (i128)M5);
    u64 sig = 0; for (int r = 0; r < b_; r++) for (int c = 0; c < 7; c++) if (VALS[c] == vals_[r]) sig += (u64)1 << (6 * c);
    u64 hcode = 0, dcode = 0; for (int r = 0; r < b_; r++) { hcode |= (u64)hts_[r] << (6 * r); dcode |= (u64)vals_[r] << (4 * r); }
    u64 inv2 = invmod(2, M5);
    for (int t = tlo; t <= thi; t++) {
        if (t > tm) break;
        u64 x = mulmod(Fm, powmod(inv2, t, M5), M5);
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
    dlog_init(Gp);
    for (int b = 1; b <= bmax; b++) {
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
