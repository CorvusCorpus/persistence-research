/* Exact decision of the candidate pairs of the MITM top step: mitm_ref.exact_pairs in C (GMP, MPFR), shared by csrc/mitm_top.c
 * and csrc/mitmgpu_host.c (docs/MITM_STEPS.md section 2).
 * The lattice is rebuilt from the job's exact data (section "exact 1": xGs, xGp, xdps, xU, xW, xb1, xb2, xe7, xphi, xTr):
 *   M rows (b1u/U, b1w/W, l(b1)/eps), (b2u/U, b2w/W, l(b2)/eps), (0, 0, -log10/eps), l(v) = vu log3 + vw log7, eps = 80 10^-Gs;
 *   Mr = Tr M, Minv = Mr^-1, K = ((U - 2)/(2U), 1/2, 1/2), hw_i = sum_r K_r |Minv[r][i]|,
 * at xdps decimal digits (mitm_ref.MitmLattice's precision) plus 64 bits. A pair (top value lam = log(1 + z_T), bottom entry with
 * F and t): v_B = a representative of {(u', w): 3^u' 7^w = F 2^-t mod 5^Gp} (dlog3, then mitm_bottom's rounding; any
 * representative gives the same solutions), tau = ((U + 2 - 2 vBu)/(2U), (W - 2 vBw)/(2W), (lam - t log2 - l(v_B))/eps + 1/2),
 * c = tau Minv; every integer mr with floor(c_i - hw_i - 10^-30) <= mr_i <= ceil(c_i + hw_i + 10^-30): m = mr Tr,
 * u = vBu + m0 b1u + m1 b2u, w = vBw + m0 b1w + m1 b2w, L = m2, kept when 2 <= u <= U, 0 <= w <= W, the value
 * val = t log2 + u log3 + w log7 - L log10 - lam is not outside (0, eps) by more than slack = 10^(digits(U) + 5 - xdps)
 * (mitm_ref's slack; the rounding error is far below it), and 2^t 3^u 7^w = F (mod 10^Gp). */
#ifndef MITM_EXACT_H
#define MITM_EXACT_H
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "dlog5.h"

typedef struct {
    int on, cands, Gs, Gp, dps; mpfr_prec_t prec;
    mpz_t U, W, M10, Tr[3][3];
    long long b1u, b1w, b2u, b2w; unsigned long long e7, phi;
    i128 det2; u64 M5, inv2;
    mpfr_t Minv[3][3], hw[3], eps, l2, l3, l7, l10, margin, slack;
} mx_lat;
static mx_lat MX;

static void mx_die(const char *m) { fprintf(stderr, "mitm_exact: %s\n", m); exit(2); }
static void mx_preinit(void) {
    mpz_inits(MX.U, MX.W, MX.M10, NULL);
    for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) mpz_init(MX.Tr[i][k]);
}
/* one job key of the exact section: 1 consumed, 0 not ours, -1 malformed */
static int mx_key(FILE *f, const char *w) {
    static char b[4096];
    if (!strcmp(w, "exact")) return fscanf(f, "%d", &MX.on) == 1 ? 1 : -1;
    if (!strcmp(w, "xcands")) return fscanf(f, "%d", &MX.cands) == 1 ? 1 : -1;
    if (!strcmp(w, "xGs")) return fscanf(f, "%d", &MX.Gs) == 1 ? 1 : -1;
    if (!strcmp(w, "xGp")) return fscanf(f, "%d", &MX.Gp) == 1 ? 1 : -1;
    if (!strcmp(w, "xdps")) return fscanf(f, "%d", &MX.dps) == 1 ? 1 : -1;
    if (!strcmp(w, "xU")) return fscanf(f, "%4095s", b) == 1 && !mpz_set_str(MX.U, b, 10) ? 1 : -1;
    if (!strcmp(w, "xW")) return fscanf(f, "%4095s", b) == 1 && !mpz_set_str(MX.W, b, 10) ? 1 : -1;
    if (!strcmp(w, "xb1")) return fscanf(f, "%lld %lld", &MX.b1u, &MX.b1w) == 2 ? 1 : -1;
    if (!strcmp(w, "xb2")) return fscanf(f, "%lld %lld", &MX.b2u, &MX.b2w) == 2 ? 1 : -1;
    if (!strcmp(w, "xe7")) return fscanf(f, "%llu", &MX.e7) == 1 ? 1 : -1;
    if (!strcmp(w, "xphi")) return fscanf(f, "%llu", &MX.phi) == 1 ? 1 : -1;
    if (!strcmp(w, "xTr")) {
        for (int i = 0; i < 9; i++) if (fscanf(f, "%4095s", b) != 1 || mpz_set_str(MX.Tr[i / 3][i % 3], b, 10)) return -1;
        return 1;
    }
    return 0;
}
static void mx_init(void) {
    if (MX.Gp < 1 || MX.Gp > 27 || MX.dps < 30 || mpz_sgn(MX.U) <= 0 || mpz_sgn(MX.W) <= 0) mx_die("bad exact section");
    MX.prec = (mpfr_prec_t)(MX.dps * 3.3219280948873623) + 64;
    dlog_init(MX.Gp); dlog_fast_init();
    MX.M5 = P5[MX.Gp]; MX.inv2 = invmod(2, MX.M5);
    MX.det2 = (i128)MX.b1u * MX.b2w - (i128)MX.b1w * MX.b2u;
    if ((MX.det2 < 0 ? -MX.det2 : MX.det2) != (i128)MX.phi || MX.phi != 4 * P5[MX.Gp - 1]) mx_die("basis does not span the coset lattice");
    mpz_ui_pow_ui(MX.M10, 10, (unsigned long)MX.Gp);
    const mpfr_prec_t p = MX.prec;
    mpfr_inits2(p, MX.eps, MX.l2, MX.l3, MX.l7, MX.l10, MX.margin, MX.slack, (mpfr_ptr)0);
    mpfr_t M[3][3], Mr[3][3], x, y, det, K[3];
    for (int i = 0; i < 3; i++) { mpfr_init2(MX.hw[i], p); mpfr_init2(K[i], p); for (int k = 0; k < 3; k++) { mpfr_inits2(p, MX.Minv[i][k], M[i][k], Mr[i][k], (mpfr_ptr)0); } }
    mpfr_inits2(p, x, y, det, (mpfr_ptr)0);
    mpfr_set_ui(x, 2, MPFR_RNDN); mpfr_log(MX.l2, x, MPFR_RNDN); mpfr_set_ui(x, 3, MPFR_RNDN); mpfr_log(MX.l3, x, MPFR_RNDN);
    mpfr_set_ui(x, 7, MPFR_RNDN); mpfr_log(MX.l7, x, MPFR_RNDN); mpfr_set_ui(x, 10, MPFR_RNDN); mpfr_log(MX.l10, x, MPFR_RNDN);
    mpfr_ui_pow_ui(x, 10, (unsigned long)MX.Gs, MPFR_RNDN); mpfr_ui_div(MX.eps, 80, x, MPFR_RNDN);
    mpfr_ui_pow_ui(x, 10, 30, MPFR_RNDN); mpfr_ui_div(MX.margin, 1, x, MPFR_RNDN);
    { char *s = mpz_get_str(NULL, 10, MX.U); long nd = (long)strlen(s); free(s);            /* len(str(U)) */
      long ex = nd + 5 - MX.dps; mpfr_set_ui(x, 10, MPFR_RNDN); mpfr_pow_si(MX.slack, x, ex, MPFR_RNDN); }
    /* M */
    const long long bu[2] = {MX.b1u, MX.b2u}, bw[2] = {MX.b1w, MX.b2w};
    for (int r = 0; r < 2; r++) {
        mpfr_set_si(M[r][0], bu[r], MPFR_RNDN); mpfr_div_z(M[r][0], M[r][0], MX.U, MPFR_RNDN);
        mpfr_set_si(M[r][1], bw[r], MPFR_RNDN); mpfr_div_z(M[r][1], M[r][1], MX.W, MPFR_RNDN);
        mpfr_mul_si(x, MX.l3, bu[r], MPFR_RNDN); mpfr_mul_si(y, MX.l7, bw[r], MPFR_RNDN); mpfr_add(x, x, y, MPFR_RNDN);
        mpfr_div(M[r][2], x, MX.eps, MPFR_RNDN);
    }
    mpfr_set_ui(M[2][0], 0, MPFR_RNDN); mpfr_set_ui(M[2][1], 0, MPFR_RNDN); mpfr_div(M[2][2], MX.l10, MX.eps, MPFR_RNDN); mpfr_neg(M[2][2], M[2][2], MPFR_RNDN);
    /* Mr = Tr M */
    for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) {
        mpfr_set_ui(Mr[i][k], 0, MPFR_RNDN);
        for (int r = 0; r < 3; r++) { mpfr_mul_z(x, M[r][k], MX.Tr[i][r], MPFR_RNDN); mpfr_add(Mr[i][k], Mr[i][k], x, MPFR_RNDN); }
    }
    /* Minv = adj(Mr) / det(Mr): Minv[i][k] = cofactor(k, i) / det */
    #define MR(a, b) Mr[(a) % 3][(b) % 3]
    for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) {          /* cofactor of entry (k, i) */
        mpfr_mul(x, MR(i + 1, k + 1), MR(i + 2, k + 2), MPFR_RNDN); mpfr_mul(y, MR(i + 1, k + 2), MR(i + 2, k + 1), MPFR_RNDN);
        mpfr_sub(MX.Minv[k][i], x, y, MPFR_RNDN);                     /* cyclic indices give the signed cofactor of (i, k) */
    }
    #undef MR
    mpfr_set_ui(det, 0, MPFR_RNDN);
    for (int k = 0; k < 3; k++) { mpfr_mul(x, Mr[0][k], MX.Minv[k][0], MPFR_RNDN); mpfr_add(det, det, x, MPFR_RNDN); }
    if (mpfr_zero_p(det)) mx_die("singular lattice");
    for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) mpfr_div(MX.Minv[i][k], MX.Minv[i][k], det, MPFR_RNDN);
    /* K, hw */
    mpfr_set_z(K[0], MX.U, MPFR_RNDN); mpfr_sub_ui(K[0], K[0], 2, MPFR_RNDN); mpfr_div_z(K[0], K[0], MX.U, MPFR_RNDN); mpfr_div_ui(K[0], K[0], 2, MPFR_RNDN);
    mpfr_set_d(K[1], 0.5, MPFR_RNDN); mpfr_set_d(K[2], 0.5, MPFR_RNDN);
    for (int i = 0; i < 3; i++) {
        mpfr_set_ui(MX.hw[i], 0, MPFR_RNDN);
        for (int r = 0; r < 3; r++) { mpfr_abs(x, MX.Minv[r][i], MPFR_RNDN); mpfr_mul(x, x, K[r], MPFR_RNDN); mpfr_add(MX.hw[i], MX.hw[i], x, MPFR_RNDN); }
    }
    for (int i = 0; i < 3; i++) { mpfr_clear(K[i]); for (int k = 0; k < 3; k++) { mpfr_clear(M[i][k]); mpfr_clear(Mr[i][k]); } }
    mpfr_clears(x, y, det, (mpfr_ptr)0);
}

typedef struct { mpfr_t tau[3], c[3], x, y, lam; mpz_t lo[3], hi[3], mr[3], m[3], u, w, z, F, r1, r2; } mx_scr;
static void mx_scr_init(mx_scr *s) {
    for (int i = 0; i < 3; i++) { mpfr_init2(s->tau[i], MX.prec); mpfr_init2(s->c[i], MX.prec); mpz_inits(s->lo[i], s->hi[i], s->mr[i], s->m[i], NULL); }
    mpfr_inits2(MX.prec, s->x, s->y, s->lam, (mpfr_ptr)0); mpz_inits(s->u, s->w, s->z, s->F, s->r1, s->r2, NULL);
}
/* s->lam = log(1 + sum 9 (d_i - 1) 10^-D_i) at the exact precision */
static void mx_lam(mx_scr *s, int j, const int *digits, const int *D) {
    mpfr_set_ui(s->lam, 0, MPFR_RNDN);
    for (int i = 0; i < j; i++) { mpfr_ui_pow_ui(s->x, 10, (unsigned long)D[i], MPFR_RNDN); mpfr_ui_div(s->x, 9 * (digits[i] - 1), s->x, MPFR_RNDN); mpfr_add(s->lam, s->lam, s->x, MPFR_RNDN); }
    mpfr_log1p(s->lam, s->lam, MPFR_RNDN);
}
/* F = -1 + sum 9 (d_r - 1) 10^h_r from mitm_bottom's codes (hcode: heights as a bit mask, from the top down; dcode: 3-bit
 * indices of the digits 2 3 4 6 7 8 9 in that order) */
static i128 mx_F(u64 hcode, u64 dcode, unsigned b) {
    static const int V[7] = {2, 3, 4, 6, 7, 8, 9};
    i128 F = -1, p; unsigned r = 0;
    for (int h = 63; h >= 0; h--) {
        if (!((hcode >> h) & 1)) continue;
        const unsigned ci = (unsigned)((dcode >> (3 * r)) & 7);
        if (ci > 6 || h > 30) mx_die("window code");
        p = 1; for (int q = 0; q < h; q++) p *= 10;
        F += (i128)(9 * (V[ci] - 1)) * p; r++;
    }
    if (r != b) mx_die("window code: digit count");
    return F;
}
static void mx_set_i128(mpz_t z, i128 v) {
    const int neg = v < 0; u128 a = neg ? (u128)(-v) : (u128)v;
    mpz_set_ui(z, (unsigned long)(u64)(a >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (unsigned long)(u64)a);
    if (neg) mpz_neg(z, z);
}
/* every solution (t, u, w, L) of the pair (s->lam, entry F, t); emit(ctx, t, u, w, L) for each; returns their number */
static int mx_check(mx_scr *s, i128 F, int t, void (*emit)(void *, int, const mpz_t, const mpz_t, const mpz_t), void *ctx) {
    /* v_B */
    i128 Fr = F % (i128)MX.M5; u64 Fm = (u64)(Fr < 0 ? Fr + (i128)MX.M5 : Fr);
    u64 xm = mulmod(Fm, powmod(MX.inv2, (u64)t, MX.M5), MX.M5), dl = dlog3(xm);
    long long c1 = rnd_div((i128)dl * MX.b2w, MX.det2), c2 = rnd_div(-(i128)MX.b1w * dl, MX.det2);
    i128 ru = (i128)dl - (i128)c1 * MX.b1u - (i128)c2 * MX.b2u, rw = -(i128)c1 * MX.b1w - (i128)c2 * MX.b2w;
    if (((ru + (i128)MX.e7 * rw - (i128)dl) % (i128)MX.phi) != 0) mx_die("coset representative");
    const long long vu = (long long)ru, vw = (long long)rw;
    /* tau, c */
    mpz_set_si(s->z, vu); mpz_mul_si(s->z, s->z, -2); mpz_add(s->z, s->z, MX.U); mpz_add_ui(s->z, s->z, 2);
    mpfr_set_z(s->tau[0], s->z, MPFR_RNDN); mpfr_div_z(s->tau[0], s->tau[0], MX.U, MPFR_RNDN); mpfr_div_ui(s->tau[0], s->tau[0], 2, MPFR_RNDN);
    mpz_set_si(s->z, vw); mpz_mul_si(s->z, s->z, -2); mpz_add(s->z, s->z, MX.W);
    mpfr_set_z(s->tau[1], s->z, MPFR_RNDN); mpfr_div_z(s->tau[1], s->tau[1], MX.W, MPFR_RNDN); mpfr_div_ui(s->tau[1], s->tau[1], 2, MPFR_RNDN);
    mpfr_set(s->x, s->lam, MPFR_RNDN);
    mpfr_mul_si(s->y, MX.l2, t, MPFR_RNDN); mpfr_sub(s->x, s->x, s->y, MPFR_RNDN);
    mpfr_mul_si(s->y, MX.l3, vu, MPFR_RNDN); mpfr_sub(s->x, s->x, s->y, MPFR_RNDN);
    mpfr_mul_si(s->y, MX.l7, vw, MPFR_RNDN); mpfr_sub(s->x, s->x, s->y, MPFR_RNDN);
    mpfr_div(s->tau[2], s->x, MX.eps, MPFR_RNDN); mpfr_add_d(s->tau[2], s->tau[2], 0.5, MPFR_RNDN);
    for (int i = 0; i < 3; i++) {
        mpfr_set_ui(s->c[i], 0, MPFR_RNDN);
        for (int r = 0; r < 3; r++) { mpfr_mul(s->x, s->tau[r], MX.Minv[r][i], MPFR_RNDN); mpfr_add(s->c[i], s->c[i], s->x, MPFR_RNDN); }
        mpfr_sub(s->x, s->c[i], MX.hw[i], MPFR_RNDN); mpfr_sub(s->x, s->x, MX.margin, MPFR_RNDN); mpfr_get_z(s->lo[i], s->x, MPFR_RNDD);
        mpfr_add(s->x, s->c[i], MX.hw[i], MPFR_RNDN); mpfr_add(s->x, s->x, MX.margin, MPFR_RNDN); mpfr_get_z(s->hi[i], s->x, MPFR_RNDU);
    }
    mx_set_i128(s->F, F);
    int nsol = 0;
    for (mpz_set(s->mr[0], s->lo[0]); mpz_cmp(s->mr[0], s->hi[0]) <= 0; mpz_add_ui(s->mr[0], s->mr[0], 1))
    for (mpz_set(s->mr[1], s->lo[1]); mpz_cmp(s->mr[1], s->hi[1]) <= 0; mpz_add_ui(s->mr[1], s->mr[1], 1))
    for (mpz_set(s->mr[2], s->lo[2]); mpz_cmp(s->mr[2], s->hi[2]) <= 0; mpz_add_ui(s->mr[2], s->mr[2], 1)) {
        for (int k = 0; k < 3; k++) {
            mpz_set_ui(s->m[k], 0);
            for (int r = 0; r < 3; r++) mpz_addmul(s->m[k], s->mr[r], MX.Tr[r][k]);
        }
        mpz_set_si(s->u, vu); mpz_set_si(s->z, MX.b1u); mpz_addmul(s->u, s->m[0], s->z); mpz_set_si(s->z, MX.b2u); mpz_addmul(s->u, s->m[1], s->z);
        if (mpz_cmp_ui(s->u, 2) < 0 || mpz_cmp(s->u, MX.U) > 0) continue;
        mpz_set_si(s->w, vw); mpz_set_si(s->z, MX.b1w); mpz_addmul(s->w, s->m[0], s->z); mpz_set_si(s->z, MX.b2w); mpz_addmul(s->w, s->m[1], s->z);
        if (mpz_sgn(s->w) < 0 || mpz_cmp(s->w, MX.W) > 0) continue;
        /* val = t log2 + u log3 + w log7 - L log10 - lam, L = m2 */
        mpfr_mul_si(s->x, MX.l2, t, MPFR_RNDN);
        mpfr_mul_z(s->y, MX.l3, s->u, MPFR_RNDN); mpfr_add(s->x, s->x, s->y, MPFR_RNDN);
        mpfr_mul_z(s->y, MX.l7, s->w, MPFR_RNDN); mpfr_add(s->x, s->x, s->y, MPFR_RNDN);
        mpfr_mul_z(s->y, MX.l10, s->m[2], MPFR_RNDN); mpfr_sub(s->x, s->x, s->y, MPFR_RNDN);
        mpfr_sub(s->x, s->x, s->lam, MPFR_RNDN);
        mpfr_neg(s->y, MX.slack, MPFR_RNDN); if (mpfr_lessequal_p(s->x, s->y)) continue;
        mpfr_add(s->y, MX.eps, MX.slack, MPFR_RNDN); if (mpfr_greaterequal_p(s->x, s->y)) continue;
        /* 2^t 3^u 7^w = F (mod 10^Gp) */
        mpz_set_ui(s->z, 2); mpz_powm_ui(s->r1, s->z, (unsigned long)t, MX.M10);
        mpz_set_ui(s->z, 3); mpz_powm(s->r2, s->z, s->u, MX.M10); mpz_mul(s->r1, s->r1, s->r2); mpz_mod(s->r1, s->r1, MX.M10);
        mpz_set_ui(s->z, 7); mpz_powm(s->r2, s->z, s->w, MX.M10); mpz_mul(s->r1, s->r1, s->r2);
        mpz_sub(s->r1, s->r1, s->F); mpz_mod(s->r1, s->r1, MX.M10);
        if (mpz_sgn(s->r1)) continue;
        nsol++; emit(ctx, t, s->u, s->w, s->m[2]);
    }
    return nsol;
}
#endif
