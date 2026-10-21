#define _GNU_SOURCE                                     /* exp10l in gsieve.c: before any header */
/* tsieve_gpu: the tree sieve's main path with the cosets on the GPU (tsieve_gpu.cu); host part.
 * Same region, decisions, counters and dump lines as tsieve6 and `gsieve_k13 main`. The host walks the coset tree
 * (tsieve6.c, node()), computes each coset's constants and exact i0 range, and sends batches to the GPU, which returns
 * counters, the points that passed stage 4 (for gsieve's full_check on the host) and the rare points near 10^XMIN
 * (decided on the host by gsieve's region rule). With dump the GPU returns every s2 point with its decision.
 * usage: tsieve_gpu PARAMFILE G|E M N [n0=N0] [threads=T] [dump] [roots=A:B]    env TG_NT (threads per block, 256),
 *        TG_BATCH (cosets per launch), TG_OUTCAP (points per launch), TG_EMU_THREADS (emulation build)
 * build: see Makefile (tsieve_gpu: CUDA; tsieve_gpu_emu: the same kernel compiled for the CPU)
 */
#include <unistd.h>
#include "tsieve_gpu.h"
#define TSIEVE_NO_MAIN
#define TS_LEAF_HOOK gpu_leaf
static void gpu_leaf(unsigned long long r, unsigned long long ylo, unsigned long long yhi);
#include "tsieve6.c"

static int BATCH = 512; static u64 OUTCAP;
static __thread tg_coset *CB; static __thread int ncb; static __thread tg_ctx *CTX; static __thread unsigned long long gc6[6];

static void run_batch(void) {
    if (!ncb) return;
    tg_point *out; unsigned long long nout;
    if (tg_run(CTX, CB, ncb, DUMP, &out, &nout, gc6)) { fprintf(stderr, "GPU output buffer overflow: raise TG_OUTCAP or lower TG_BATCH\n"); exit(6); }
    for (unsigned long long k = 0; k < nout; k++) {
        u64 a = out[k].a, b = out[k].b, c = out[k].c; unsigned code = out[k].code;
        if (code == TG_CORNER) { survivor((i64)a, (i64)b, (i64)c); continue; }       /* gsieve's region rule on the host */
        if (DUMP) emit("S2 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        if (code == TG_R3) { if (DUMP) emit("R3 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        if (code == TG_R4) { if (DUMP) emit("R4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c); continue; }
        if (code != TG_S4) { fprintf(stderr, "bad decision code %u\n", code); exit(6); }
        if (DUMP) emit("S4 %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, (unsigned long long)c);
        full_check(0, a, b, c, C_EVEN);
    }
    ncb = 0;
}
static void gpu_leaf(unsigned long long r, unsigned long long ylo, unsigned long long yhi) {
    c_leaves++;
    i64 xi[3]; i128 lo, hi;
    if (!bounds(RQ1, n1, 0, xi, (i128)r, &lo, &hi)) return;
    tg_coset *C = &CB[ncb++];
    C->r = r; C->h0 = (ylo << TK) / TM1; C->ybase = yhi;
    C->d0 = mulm1(ylo, TWEPS[0]); C->d1 = mulm1(ylo, TWEPS[1]); C->d2 = mulm1(ylo, TWEPS[2]);
    u64 inv2 = invmod(C->d2 % TMm, TMm);
    C->inv2 = (unsigned)inv2; C->s0 = (unsigned)(inv2 * (yhi % TMm) % TMm);
    C->si0 = (unsigned)(inv2 * (C->d0 % TMm) % TMm); C->si1 = (unsigned)(inv2 * (C->d1 % TMm) % TMm);
    if (lo < -(1LL << 30) || hi > (1LL << 30)) { fprintf(stderr, "i0 range out of bounds\n"); exit(6); }
    C->i0lo = (int)lo; C->i0hi = (int)hi;
    if (ncb == BATCH) run_batch();
}
static void *gworker(void *arg) {
    (void)arg;
    mpz_init(gx); mpz_init(gm); mpz_init(gt);
    OBC = 1 << 16; OB = malloc(OBC); OBN = 0;
    SB = malloc(sizeof(spt) * SBMAX); SB3 = malloc(sizeof(spt) * SBMAX); nsb = 0;
    CB = malloc(sizeof(tg_coset) * BATCH); ncb = 0; CTX = tg_ctx_new(BATCH, OUTCAP);
    for (;;) {
        u64 k = __atomic_fetch_add(&NEXTROOT, 1, __ATOMIC_RELAXED);
        if (k >= RTO) break;
        node(ROOTS[k].j, ROOTS[k].r, ROOTS[k].ylo, ROOTS[k].yhi);
        run_batch(); flush_survivors();
        if (OBN) { pthread_mutex_lock(&OUTMU); fwrite(OB, 1, OBN, stdout); pthread_mutex_unlock(&OUTMU); OBN = 0; }
    }
    cnt[2] += gc6[3]; cnt[3] += gc6[4]; cnt[4] += gc6[5];
    pthread_mutex_lock(&OUTMU);
    G_leaves += c_leaves; G_s1rows += gc6[0]; G_rows += gc6[1]; G_hits += gc6[2];
    for (int k = 0; k < 8; k++) G_cnt[k] += cnt[k];
    pthread_mutex_unlock(&OUTMU);
    free(OB); mpz_clear(gx); mpz_clear(gm); mpz_clear(gt);
    return NULL;
}
static int VINT, VSTOP2; static u64 NVER;
static void *verifier(void *arg) {                     /* every VINT seconds: fail-stop check of all tables (non-ECC) */
    (void)arg;
    for (int t = 0; !__atomic_load_n(&VSTOP2, __ATOMIC_SEQ_CST); t++) {
        sleep(1);
        if ((t + 1) % VINT) continue;
        if (tg_verify()) { fprintf(stderr, "TABLE CHECK FAILED during the run (memory fault?): stopping, discard this run\n"); fflush(stdout); _exit(3); }
        NVER++;
    }
    return NULL;
}
int main(int argc, char **argv) {
    int rc = ts_setup(argc, argv); if (rc) return rc;
    if (!R2OK) { fprintf(stderr, "the projected region does not fit 64-bit arithmetic\n"); return 1; }
    if (TM > 7) { fprintf(stderr, "the GPU build supports M <= 7 (U in shared memory)\n"); return 1; }
    tg_params p; memset(&p, 0, sizeof p);
    p.k = TK; p.m = TM; p.dm = DM; p.nl = NL; p.nt = getenv("TG_NT") ? atoi(getenv("TG_NT")) : 256;
    p.n3 = n3; p.n2 = n2; p.skip = TSKIP; p.nfull = TNFULL; p.nrem = TNREM;
    p.m1 = TM1; p.mm = TMm; p.br1 = (u64)(((u128)1 << 64) / TM1); p.brm = TBRM;
    p.nleaf = 1; for (int i = 0; i < TM; i++) p.nleaf *= DM == 0 ? 4 : 5;
    for (int i = 0; i < TM; i++) { p.p2k[i] = TP2K[i]; p.p5[i] = TL5[i]; }
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) p.v[i][j] = TV[i][j];
    for (int q = 0; q < n3; q++) { for (int m = 0; m < 3; m++) p.a3[q][m] = A3[q][m]; p.b3_0[q] = (i64)RQ3[q].b0; p.b3_r[q] = (i64)RQ3[q].br; }
    for (int t = 0; t < n2; t++) { p.a2[t][0] = R2A[t][0]; p.a2[t][1] = R2A[t][1]; p.b2_0[t] = R2B0[t]; p.b2_r[t] = R2BR[t]; }
    p.cornerlg = XMIN + 2.0; p.l2 = L2; p.l3 = L3; p.l7 = L7; p.tb = TB;
    memcpy(p.ok4, TOK4, 10000); memcpy(p.okr, TOKR, 10000);
    for (int d = 0; d < 10; d++) p.ok1[d] = (unsigned char)okdigits((u64)d, 1);
    if (tg_init(&p, &T2[0][0], &T3[0][0], &T7[0][0], (3 * (u64)NN + 3) * NL, (2 * (u64)NN + 2) * NL, ((u64)NN + 2) * NL,
                (const unsigned long long *)TOP, (1ULL << TB) / 64 + 1)) return 5;
    if (getenv("TG_BATCH")) BATCH = atoi(getenv("TG_BATCH"));
    OUTCAP = getenv("TG_OUTCAP") ? strtoull(getenv("TG_OUTCAP"), 0, 10) : (DUMP ? 1ULL << 22 : 1ULL << 20);
    fprintf(stderr, "tsieve_gpu: device %s, M=%d, block %d threads, batch %d cosets\n", tg_device(), TM, p.nt, BATCH);
    if (tg_verify()) { fprintf(stderr, "TABLE CHECK FAILED before the run\n"); return 3; }
    struct timespec w0, w1; clock_gettime(CLOCK_MONOTONIC, &w0); clock_t c0 = clock();
    pthread_t *tid = malloc(TS_NTH * sizeof *tid), vt;
    VINT = getenv("TG_VERIFY") ? atoi(getenv("TG_VERIFY")) : 60; if (VINT < 1) VINT = 1;
    pthread_create(&vt, NULL, verifier, NULL);
    for (int i = 0; i < TS_NTH; i++) pthread_create(&tid[i], NULL, gworker, NULL);
    for (int i = 0; i < TS_NTH; i++) pthread_join(tid[i], NULL);
    __atomic_store_n(&VSTOP2, 1, __ATOMIC_SEQ_CST); pthread_join(vt, NULL);
    if (tg_verify()) { fprintf(stderr, "TABLE CHECK FAILED at the end of the run (memory fault?): discard this run\n"); return 3; }
    fprintf(stderr, "table checks passed: %llu (device tables and constants against host copies, host copies against their hash)\n",
            (unsigned long long)NVER + 1);
    clock_gettime(CLOCK_MONOTONIC, &w1);
    double wall = (double)(w1.tv_sec - w0.tv_sec) + 1e-9 * (double)(w1.tv_nsec - w0.tv_nsec), cpu = (double)(clock() - c0) / CLOCKS_PER_SEC;
    return ts_report(strcmp(tg_device(), "emulation") ? "treegpu" : "treeemu", cpu, wall);
}
