/* mitmgpu: csrc/mitm_top with the probe on the GPU (kernel csrc/mitmgpu.cu). Same job file, same output:
 *   CAND o D_0 .. D_(j-1) idx     every multiset-compatible table entry within the box (margins: the job's, plus 64 units of
 *                                 2^-64 for the kernel's fixed point), sorted;
 *   DONE ncases nprobes ncand     nprobes: table entries examined (the cells of each box, not mitm_top's 27 neighbours).
 * Tables and the fixed-point lambda are topgpu_host.c's (class with FL >= bits(max |Theta_i|) + 80; beyond the largest class
 * the job goes to mitm_top: env MG_CPU, else mitm_top beside this binary); orderings and the enumeration of depths (D_0 >= 1, D_i < D_(i+1), D_i < ranges[i]) are mitm_top's, as prefixes plus a last gap.
 * usage: mitmgpu JOB      env MG_BATCH (prefixes per launch, default 2^18), MG_DUMP=1 (validation: "PT o D_0 .. P_0 P_1 P_2" for
 *                         every case, no probe), MG_TIME=1 (timing on stderr)
 * build: see work/Makefile (mitmgpu: CUDA; mitmgpu_emu: the same kernel compiled for the CPU) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "mitmgpu.h"

static void die(const char *m) { fprintf(stderr, "%s\n", m); exit(3); }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

#define NCLS 11
static const int CLS[NCLS] = {6, 8, 10, 12, 14, 16, 20, 22, 24, 28, 32};
static const int VALS[MG_NC] = {2, 3, 4, 6, 7, 8, 9};
static const int CVAL[MG_NC] = {9, 18, 27, 45, 54, 63, 72};

static unsigned int *TB; static long long TW, TCAP;
static long long talloc(long long words) { long long o = TW; TW += words; if (TW > TCAP) { TCAP = 2 * TW + 4096; TB = realloc(TB, TCAP * 4); } memset(TB + o, 0, words * 4); return o; }
static void put_mpz(long long off, const mpz_t x, int limbs) {
    size_t cnt = 0; unsigned int tmp[512];
    if (mpz_sgn(x) < 0 || mpz_sizeinbase(x, 2) > (size_t)32 * limbs) die("put_mpz: value out of range");
    memset(tmp, 0, sizeof tmp); mpz_export(tmp, &cnt, -1, 4, 0, 0, x);
    for (int i = 0; i < limbs; i++) TB[off + i] = tmp[i];
}

/* orderings: mitm_top's gen (lexicographic in VALS order) */
static int j_, cnt_[MG_NC], nord, capo, rem_[MG_NC], cur_[MG_MAXJ]; static int (*ords)[MG_MAXJ];
static void gen(int r) {
    if (r == j_) { if (nord == capo) { capo = capo ? 2 * capo : 1024; ords = realloc(ords, capo * sizeof *ords); }
                   memcpy(ords[nord++], cur_, sizeof cur_); return; }
    for (int c = 0; c < MG_NC; c++) if (rem_[c]) { rem_[c]--; cur_[r] = c; gen(r + 1); rem_[c]++; }
}

typedef struct { unsigned long long key, q[3], hcode, dcode; unsigned int t, b; unsigned long long sig; } mb_entry;   /* mitm_bottom */

static mg_out *OUT; static long long OUTCAP = 1 << 16; static int DUMP;
static unsigned long long NCASE, NPROBE, NCAND; static double T_GPU;
static int cmp_out(const void *A_, const void *B_) {
    const mg_out *a = (const mg_out *)A_, *b = (const mg_out *)B_;
    if (a->prefix != b->prefix) return a->prefix < b->prefix ? -1 : 1;
    if (a->gl != b->gl) return a->gl < b->gl ? -1 : 1;
    return a->idx < b->idx ? -1 : (a->idx > b->idx);
}
static void run_batch(const mg_prefix *pf, long n) {
    unsigned long long cn[2] = {0, 0}; long long no;
    if (!OUT) OUT = malloc(sizeof(mg_out) * OUTCAP);
    const double t0 = now_s();
    while ((no = mg_run(pf, (int)n, OUT, OUTCAP, cn)) < 0) { OUTCAP *= 4; OUT = realloc(OUT, sizeof(mg_out) * OUTCAP); }
    T_GPU += now_s() - t0;
    NCASE += cn[0]; NPROBE += cn[1];
    qsort(OUT, (size_t)no, sizeof(mg_out), cmp_out);
    for (long long r = 0; r < no; r++) {
        const mg_prefix *q = &pf[OUT[r].prefix];
        printf(DUMP ? "PT %u" : "CAND %u", q->o);
        long acc = 0;
        for (int s = 0; s < j_ - 1; s++) { acc += q->g[s]; printf(" %ld", acc); }
        printf(" %ld", acc + (long)OUT[r].gl);
        if (DUMP) printf(" %llu %llu %llu\n", OUT[r].p[0], OUT[r].p[1], OUT[r].p[2]);
        else { printf(" %llu\n", OUT[r].idx); NCAND++; }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: mitmgpu JOB\n"); return 2; }
    const double t_start = now_s();
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("job"); return 2; }
    char w[64], path[4096] = ""; static char tbuf[3][100000];
    long prec = 0, rng[MG_MAXJ + 1]; unsigned long long hw[3] = {0, 0, 0}, ncell[3] = {0, 0, 0}; int have_th = 0;
    while (fscanf(f, "%63s", w) == 1) {
        if (!strcmp(w, "prec")) { if (fscanf(f, "%ld", &prec) != 1) die("bad prec"); }
        else if (!strcmp(w, "j")) { if (fscanf(f, "%d", &j_) != 1 || j_ < 1 || j_ > MG_MAXJ) die("bad j"); }
        else if (!strcmp(w, "cnt")) { for (int c = 0; c < MG_NC; c++) if (fscanf(f, "%d", &cnt_[c]) != 1) die("bad cnt"); }
        else if (!strcmp(w, "ranges")) { for (int i = 0; i < j_; i++) if (fscanf(f, "%ld", &rng[i]) != 1) die("bad ranges"); }
        else if (!strcmp(w, "theta0") || !strcmp(w, "theta1") || !strcmp(w, "theta2")) {
            int i = w[5] - '0'; if (fscanf(f, "%99999s", tbuf[i]) != 1) die("bad theta"); have_th |= 1 << i; }
        else if (!strcmp(w, "hw")) { if (fscanf(f, "%llu %llu %llu", &hw[0], &hw[1], &hw[2]) != 3) die("bad hw"); }
        else if (!strcmp(w, "ncell")) { if (fscanf(f, "%llu %llu %llu", &ncell[0], &ncell[1], &ncell[2]) != 3) die("bad ncell"); }
        else if (!strcmp(w, "table")) { if (fscanf(f, "%4095s", path) != 1) die("bad table"); }
        else { fprintf(stderr, "unknown key %s\n", w); return 2; }
    }
    fclose(f);
    if (have_th != 7 || prec < 64 || !path[0]) die("incomplete job");
    long smax = 0; for (int i = 0; i < j_; i++) { if (rng[i] > smax) smax = rng[i]; }
    if (smax >= 590) die("depth range beyond the tables (590)");
    for (int c = 0; c < MG_NC; c++) if (cnt_[c] < 0 || cnt_[c] > 63) die("bad digit count");

    /* bottom table: entries sorted by key (mitm_bottom); the hash maps each key to its run */
    FILE *tb = fopen(path, "rb"); if (!tb) { perror("table"); return 2; }
    unsigned long long nT = 0; if (fread(&nT, sizeof nT, 1, tb) != 1) die("short table");
    if (nT >= (1ULL << 32)) die("table too large");
    mg_entry *E = malloc((nT ? nT : 1) * sizeof *E); unsigned long long *K = malloc((nT ? nT : 1) * sizeof *K);
    unsigned long long ndist = 0;
    {
        const size_t CH = 1 << 16; mb_entry *buf = malloc(CH * sizeof *buf); unsigned long long done = 0;
        while (done < nT) {
            size_t want = nT - done < CH ? (size_t)(nT - done) : CH;
            if (fread(buf, sizeof *buf, want, tb) != want) die("short table");
            for (size_t r = 0; r < want; r++) {
                const unsigned long long x = done + r;
                if (buf[r].key == ~0ULL) die("reserved key");
                if (x && buf[r].key < K[x - 1]) die("table not sorted");
                if (!x || buf[r].key != K[x - 1]) ndist++;
                K[x] = buf[r].key;
                for (int a = 0; a < 3; a++) E[x].q[a] = buf[r].q[a];
                E[x].sig = buf[r].sig;
            }
            done += want;
        }
        free(buf);
    }
    fclose(tb);
    int full[3];
    for (int a = 0; a < 3; a++) {
        if (ncell[a] < 1) die("bad ncell");
        hw[a] = hw[a] > ~0ULL - 64 ? ~0ULL : hw[a] + 64;  /* the kernel's fixed-point error (< 2 units); |P - Q| <= 2^63 always */
        /* box width 2 hw >= 2^64 - 2^64 / ncell: every cell (the cyclic run lo .. hi may wrap onto its start) */
        const unsigned __int128 W2 = (unsigned __int128)2 * hw[a], lim = ((unsigned __int128)1 << 64) - (((unsigned __int128)1 << 64) / ncell[a]);
        full[a] = W2 >= lim;
    }
    /* cell keys (c0 ncell1 + c1) ncell2 + c2 are taken mod 2^64 as mitm_bottom computes them: cells with equal keys share a
       slot, and the box test separates their entries */
    int hbits = 4; while ((1ULL << hbits) < 2 * ndist + 16) hbits++;
    if (hbits > 40) die("table too large");
    mg_slot *H = malloc(((size_t)1 << hbits) * sizeof *H);
    for (size_t h = 0; h < ((size_t)1 << hbits); h++) { H[h].key = ~0ULL; H[h].first = H[h].count = 0; }
    {   /* each run of equal keys -> one slot (linear probing, as the kernel searches) */
        const unsigned long long mask = (1ULL << hbits) - 1;
        for (unsigned long long r = 0; r < nT;) {
            unsigned long long s = r; while (s < nT && K[s] == K[r]) s++;
            unsigned long long h = (K[r] * 0x9E3779B97F4A7C15ULL) >> (64 - hbits);
            while (H[h].key != ~0ULL) h = (h + 1) & mask;
            H[h].key = K[r]; H[h].first = (unsigned)r; H[h].count = (unsigned)(s - r);
            r = s;
        }
        free(K);
    }

    /* orderings and their digit signatures */
    for (int c = 0; c < MG_NC; c++) rem_[c] = cnt_[c];
    if (j_ > 0) { int tot = 0; for (int c = 0; c < MG_NC; c++) tot += cnt_[c]; if (j_ > tot) die("j exceeds the bag"); }
    gen(0);

    /* Theta_i (the job's decimals) -> round(|Theta_i| 2^64); bits of the largest integer part choose the class */
    mpfr_t th[3], a_, b_, s2; mpz_t z, TH[3]; long tb_bits = 0; int neg[3];
    for (int i = 0; i < 3; i++) {
        mpfr_init2(th[i], prec + 64); if (mpfr_set_str(th[i], tbuf[i], 10, MPFR_RNDN)) die("bad theta value");
        mpz_init(TH[i]); neg[i] = mpfr_sgn(th[i]) < 0;
        mpfr_t x; mpfr_init2(x, prec + 64); mpfr_abs(x, th[i], MPFR_RNDN); mpfr_mul_2ui(x, x, 64, MPFR_RNDN);
        mpfr_get_z(TH[i], x, MPFR_RNDN);                 /* round(|Theta_i| 2^64) */
        long bits = (long)mpz_sizeinbase(TH[i], 2) - 64; if (bits > tb_bits) tb_bits = bits;
        mpfr_clear(x);
    }
    int cls = -1; for (int c = 0; c < NCLS; c++) if (32 * (CLS[c] - 1) >= tb_bits + 80) { cls = c; break; }
    if (cls < 0) {                                       /* beyond the size classes: mitm_top (CPU, MPFR) on the same job */
        static char cpu[4200]; const char *env = getenv("MG_CPU");
        if (env) snprintf(cpu, sizeof cpu, "%s", env);
        else { snprintf(cpu, sizeof cpu, "%s", argv[0]); char *sl = strrchr(cpu, '/'); snprintf(sl ? sl + 1 : cpu, sizeof cpu - (sl ? (size_t)(sl + 1 - cpu) : 0), "mitm_top"); }
        execv(cpu, argv); perror("exec mitm_top"); return 2;
    }
    const int L = CLS[cls], FL = 32 * (L - 1);
    mg_params Pp; memset(&Pp, 0, sizeof Pp);
    Pp.j = j_; Pp.cls = cls; Pp.LW = L; Pp.FL = FL; Pp.Smin = 5; Pp.smax = (int)smax; Pp.hbits = hbits;
    for (int a = 0; a < 3; a++) { Pp.hw[a] = hw[a]; Pp.ncell[a] = ncell[a]; Pp.full[a] = full[a]; }
    for (int c = 0; c < MG_NC; c++) Pp.cntsig |= (unsigned long long)cnt_[c] << (6 * c);
    for (long S = 0; S < 600; S++) {                     /* topgpu: |w| <= 72 10^-S, (n + 1) (S log2 10 - log2 72) >= FL + 1 */
        double bits = S * 3.321928094887362 - 6.17;
        Pp.nterm[S] = bits <= 1 ? 0 : (int)((FL + 1) / bits) + 1;
    }
    int nmax = 1; for (long S = Pp.Smin; S <= smax; S++) if (Pp.nterm[S] > nmax) nmax = Pp.nterm[S];
    const long SS = smax + 1;
    Pp.o_cs = talloc((long long)nord * ((j_ + 3) / 4) + 1);
    Pp.o_sig = talloc(2LL * nord + 2);
    for (int o = 0; o < nord; o++) {
        unsigned long long sg = 0;
        for (int i = 0; i < j_; i++) { TB[Pp.o_cs + (long long)o * ((j_ + 3) / 4) + (i >> 2)] |= (unsigned)ords[o][i] << (8 * (i & 3)); sg += 1ULL << (6 * ords[o][i]); }
        TB[Pp.o_sig + 2LL * o] = (unsigned)sg; TB[Pp.o_sig + 2LL * o + 1] = (unsigned)(sg >> 32);
    }
    Pp.o_lam1 = talloc((long long)MG_NC * SS * L); Pp.o_ia1 = talloc((long long)MG_NC * SS * L); Pp.o_c10 = talloc((long long)MG_NC * SS * L);
    Pp.o_inv = talloc((long long)(nmax + 1) * L);
    mpfr_inits2(FL + 128, a_, b_, s2, (mpfr_ptr)0); mpz_init(z);
    mpfr_set_ui_2exp(s2, 1, FL, MPFR_RNDN);
    for (int c = 0; c < MG_NC; c++) for (long S = 0; S <= smax; S++) {          /* topgpu_host.c's tables, unchanged */
        const long long e_ = ((long long)c * SS + S) * L;
        mpfr_ui_pow_ui(a_, 10, S, MPFR_RNDN); mpfr_ui_div(a_, CVAL[c], a_, MPFR_RNDN);
        mpfr_mul(b_, a_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDD); put_mpz(Pp.o_c10 + e_, z, L);
        mpfr_log1p(b_, a_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_lam1 + e_, z, L);
        mpfr_add_ui(b_, a_, 1, MPFR_RNDN); mpfr_ui_div(b_, 1, b_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_ia1 + e_, z, L);
    }
    for (long n = 1; n <= nmax; n++) { mpz_set_ui(z, 0); mpz_setbit(z, FL); mpz_fdiv_q_ui(z, z, n); put_mpz(Pp.o_inv + n * L, z, L); }
    Pp.o_sl = talloc(4096LL * L); Pp.o_sia = talloc(4096LL * L);
    for (int e_ = 0; e_ < 4096; e_++) {
        mpfr_set_ui(a_, 0, MPFR_RNDN);
        for (int h = 1; h <= 4; h++) {
            const int code = (e_ >> (3 * (h - 1))) & 7;
            if (code) { mpfr_ui_pow_ui(b_, 10, h, MPFR_RNDN); mpfr_ui_div(b_, CVAL[code - 1], b_, MPFR_RNDN); mpfr_add(a_, a_, b_, MPFR_RNDN); }
        }
        mpfr_log1p(b_, a_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_sl + (long long)e_ * L, z, L);
        mpfr_add_ui(b_, a_, 1, MPFR_RNDN); mpfr_ui_div(b_, 1, b_, MPFR_RNDN); mpfr_mul(b_, b_, s2, MPFR_RNDN); mpfr_get_z(z, b_, MPFR_RNDN); put_mpz(Pp.o_sia + (long long)e_ * L, z, L);
    }
    for (int i = 0; i < 3; i++) {                        /* frac_theta returns frac(-Theta' lambda): store Theta' = -Theta_i */
        Pp.o_th[i] = talloc(L + 2);
        TB[Pp.o_th[i]] = neg[i] ? 0u : 1u;               /* -Theta_i < 0 iff Theta_i > 0 */
        mpz_fdiv_r_2exp(z, TH[i], 32 * (L + 1));
        put_mpz(Pp.o_th[i] + 1, z, L + 1);
    }
    Pp.words = TW;
    DUMP = Pp.dump = getenv("MG_DUMP") && atoi(getenv("MG_DUMP")) ? 1 : 0;
    if (DUMP) printf("MGDUMP FL %d L %d theta_bits %ld\n", FL, L, tb_bits);
    if (mg_init(&Pp, TB, E, (long long)nT, H)) die("mg_init");
    const double t_init = now_s();

    /* enumeration: topgpu_host.c's prefix odometer (g0lo = 1, no g0hi), batches within one ordering */
    long batch = getenv("MG_BATCH") ? atol(getenv("MG_BATCH")) : (1L << 18);
    mg_prefix *pf = malloc(sizeof(mg_prefix) * batch); long np = 0;
    long g[MG_MAXJ + 1], ss[MG_MAXJ + 1];
    for (long o = 0; o < nord; o++) {
        long i = 0; ss[0] = 0;
        if (j_ >= 2) g[0] = 0;
        for (;;) {
            if (j_ >= 2) {
                int found = 0;
                while (i >= 0) {
                    g[i]++;
                    long hi = rng[i] - 1 - ss[i];
                    if (g[i] > hi) { i--; continue; }
                    ss[i + 1] = ss[i] + g[i];
                    if (i == j_ - 2) { found = 1; break; }
                    i++; g[i] = 0;
                }
                if (!found) break;
            }
            long acc = 0; for (long q = 0; q < j_ - 1; q++) acc += g[q];
            long ghi = rng[j_ - 1] - 1 - acc, glo = 1;
            if (glo <= ghi) {
                mg_prefix *q = &pf[np++]; memset(q, 0, sizeof *q);
                q->o = (unsigned)o; q->glo = (unsigned short)glo; q->ghi = (unsigned short)ghi;
                for (long s = 0; s < j_ - 1; s++) q->g[s] = (unsigned short)g[s];
                if (np == batch) { run_batch(pf, np); np = 0; }
            }
            if (j_ < 2) break;
        }
        if (np) { run_batch(pf, np); np = 0; }
    }
    printf("DONE %llu %llu %llu\n", NCASE, NPROBE, NCAND);
    if (getenv("MG_TIME"))
        fprintf(stderr, "MGTIME %s class L=%d table %llu entries (%llu cells, 2^%d slots) setup %.3f s, kernel+copies %.3f s, total %.3f s\n",
                mg_device(), L, nT, ndist, hbits, t_init - t_start, T_GPU, now_s() - t_start);
    return 0;
}
