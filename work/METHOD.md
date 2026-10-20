# Complete preimage sets for even-target nodes (multiplicative persistence)

## Setting
f(n) = product of the decimal digits of n. For a final digit (target) δ ∈ {2,4,6,8}, the tree T_δ has
as nodes the digit products whose f-orbit ends at δ; an integer n with f-orbit ending at δ has
persistence depth(f(n)) + 1. Odd targets were settled by Brier–Clavier–Gutsche–Naccache (2021);
target 0 (where persistence 11 lives) is untouched here.

The exhaustive search (this project) lists every zero-free 7-smooth number below 10^5000; from it the
known trees have 33, 9, 84, 51 nodes and heights 5, 4, 7, 5 for δ = 2, 4, 6, 8.

## Lemmas
**Lemma 0 (no factor 5).** Every node of T_δ, and every 7-smooth y with f(y) a node, has the form
2^t 3^u 7^w. *Proof.* If 5 | s and s ends in 0, f(s) = 0. If s ends in 5 and has an even digit,
f(s) ≥ 10 ends in 0, so the next value is 0. If all digits of s are odd, every later value is odd, so
the orbit cannot end at an even digit. A preimage y of a node s with 5 ∤ s has no digit 5 or 0, so it
does not end in 0 or 5. (This closes the step left unproved in Fonga, arXiv 2608.27802, §5.)

**Lemma 1 (shape).** If f(y) = s then y = (1s) + a multiset D of digits from {2,3,4,6,7,8,9} with
product s. With positions a_1 > … > a_k of the non-1 digits, L = number of digits, gaps
g_0 = L − a_1, g_i = a_i − a_{i+1}:  9y = 10^L − 1 + Σ 9(d_i − 1)·10^{a_i}.

**Lemma 2 (bottom).** a_k ≤ 1. *Proof.* If the last two digits were 11, y would be odd, so
y = 3^u 7^w; but 3^u 7^w mod 100 lies in {01,03,07,09,21,…,89} (even tens digit).

**Lemma 3 (2-adic).** t = v_2(y) ≤ t_max(D), computed by the digit-placement recursion
(Fonga 2026, Thm 4.1); re-implemented here, reproduces the four published values and brute force.

**Lemma 4 (linear forms).** For 0 ≤ j < k:
0 < Λ_j < 80·10^{−g_j}, where Λ_j = t log2 + (u+2) log3 + w log7 − L log10 − log(1 + z_j),
z_j = Σ_{i≤j} 9(d_i − 1)·10^{−(g_0+…+g_{i−1})}.

**Lemma 5 (Matveev).** Applying Matveev's theorem (real case, D = 1) level by level bounds every gap,
hence L ≤ X(D). Any gap already proved bounded is substituted, which can shrink X enormously.

## The gap tests (threshold G ≥ 2·log10 X + 12)
*Top test for gap j.* Suppose g_j ≥ G. For each top configuration (d_1..d_j and g_0..g_{j−1} within
their proved ranges), the 4-dimensional lattice spanned by (λ_t,0,0,[C log2]), (0,1,0,[C log3]),
(0,0,λ_w,[C log7]), (0,0,0,−[C log10]) must contain a point within an explicit radius R of
(0,0,0,−[C·log(1+z_j)]). LLL (exact, FLINT) plus the last Gram–Schmidt coordinate shows there is none;
in the rare cases where that bound is not enough, an exact Fincke–Pohst enumeration lists all lattice
points within R and each is checked against the real inequality with interval arithmetic.

*Bottom test for gap j.* Suppose g_j ≥ G. Then the digits below gap j form a configuration with
proved gap ranges and a_k ≤ 1, and y mod 10^{a_{j+1}+G} is known. Exact 2-adic and 5-adic logarithms
turn 2^t 3^{u+2} 7^w ≡ F into a congruence u ≡ αw + β mod Q on a lattice of covolume > (box diameter)²;
at most one lattice point can lie in the box of admissible (u, w), and none does.

*Bootstrapping.* After each proved bound g_j < G, Matveev's chain is recomputed with it; X falls from
~10^160 to ~10^20 and G from ~330 to ~52, and later enumerations use the proved ranges.

*Conclusion per digit multiset.* When every gap has a proved bound and Σ(bound − 1) + 1 ≤ 5000, every
solution has at most 5000 digits and is therefore in the exhaustive list.

## Validation
- Lemma 2 and Lemma 4 checked on all 1,173 linear forms from known numbers.
- Planted solutions at proof-scale parameters: 160/160 detected by the top test, 160/160 by the
  bottom test; fast and slow 10-adic solvers agree on 1,002 configurations.
- Two independent orchestrations (non-bootstrapped and bootstrapped) agree on every node both proved.
- Bugs found and fixed during the work: zero lattice point treated as a solution (fixed by exact
  enumeration with validity checks); non-terminating log series for a digit at position 0 (fixed and
  re-tested); series truncation now uses floor(log_p n) as a safe bound.

## Caveats
Computer-assisted and not peer reviewed. Depends on Matveev's theorem, the exhaustive list below
10^5000, the 2-adic recursion, and the correctness of the code (attached).

---------------------------------------------------------------------------------------------------
# Second session: C inner loops, a precision audit, and the cumulative method

## 1. C inner loops (csrc/top.c, csrc/bottom.c; prove_iter3.py)
*Top.* For each configuration the C loop computes m = round(-C log(1+z)) with MPFR (every operation
correctly rounded; z is a sum of positive terms and log1p preserves relative error, so the computed value
is within C*2.2*(j+12)*2^-prec < 2^-50 of the true one, prec >= bitlen(C)+64), then the exact level-1
test of the Python code (last Gram-Schmidt coordinate, integer arithmetic) and a new exact level-2 test
(the first two levels of the Fincke-Pohst enumeration in exact rational arithmetic). Configurations that
survive both are re-checked in Python exactly (valid_points); the C target is re-derived in Python with
interval arithmetic at both ends of every reported run.
*Bottom.* Exact GMP port of fastbottom.solutions (p-adic logarithms, CRT, Gauss-reduced coset lattice).
*Validation.* Identical targets on 12,000 sampled configurations; identical level-1 failure sets on four
full enumerations; every level-2 kill re-enumerated by Fincke-Pohst has no lattice point within R;
bottom solver agrees with Python on 1,900 queries and recovers 400/400 planted solutions; with the Python
cost model the C prover reproduces four recorded proofs exactly (same bounds, same numbers of top cases
and bottom solves).

## 2. Precision audit (found while porting)
fastbottom.BottomContext precomputes p-adic data modulo 5^Mmax, 2^(Mmax+2) with Mmax = 4G+8, but never
checked that the modulus exponent M of a bottom solve satisfies M <= Mmax. Replaying all 307 recorded
proof chains (they are deterministic) shows the maximum M never exceeded Mmax, so every earlier result
stands; for bags with >= 10 non-1 digits it would have overflowed silently (excess 41 digits at k = 10).
The Python solver now asserts M <= Mmax, and the C solver sizes its precision from the case ranges.

## 3. Cumulative depths (prove_cum.py; csrc/topc, csrc/bottomc)
**Stronger linear-form bound.** With D_j = L - a_{j+1} = g_0 + ... + g_j,
    9y = 10^L (1 + z_j) + R_j,   0 < R_j < 80 * 10^(a_{j+1})   =>   0 < Lam_j < 80 * 10^(-D_j)
(the earlier proofs used the weaker 80*10^(-g_j)). Checked on all 221,467 linear forms of 8,799 known
zero-free 7-smooth numbers: no violation, largest Lam_j * 10^(D_j) = 79.81.
**Facts.** top[j] = T: every solution has D_j < T.   bot[i] = H: every solution has a_i < H (a_0 = L);
bot[k] = 2 by Lemma 2.
**Top test, level j.** Enumerate prefixes (d_1..d_j; g_0..g_{j-1}) allowed by top[0..j-1] (compositions,
not boxes). If D_j >= G then 0 < Lam_j < 80*10^-G; the unchanged lattice test excludes this, so D_j < G.
**Bottom test, level j.** Enumerate digits j+1..k at heights allowed by bot[j+1..k]; with
M = max(G_b, a_{j+1}+1) and a_j >= M, y mod 10^M is fixed by the configuration (all digits between are
1s); the unchanged 10-adic test excludes it, so a_j < M, i.e. a_j < G_b (G_b = max(G, tmax+8), keeping
>= 8 bits of 2-adic information).
**Bootstrapping.** Matveev's bound is applied to D_j directly (D_j < (ln80 + C5 beta Pi (D_{j-1} ln10 + ln9))/ln10);
L <= D_j + a_{j+1} for any split. Stop when the bound on L is <= 5000.
**Effect.** Case counts drop from ~G^j to ~C(G, j): k = 9 bags 10x fewer cases (e.g. 2.95M vs 44.7M top
cases), k = 12 about 75x; final length bounds 101 instead of 358-460.
**Validation.** Prefix enumeration identical to a Python reference (4 cases); bottom enumeration and
verdicts identical to Python on 3 configuration spaces; planted genuine numbers detected 150/150 (top)
and 80/80 (bottom); end-to-end on 8 previously proved bags; every genuine solution satisfies every
proved fact D_j < top[j], a_i < bot[i].
**Survivor resolution.** A lattice survivor is a specific exponent vector; if the last 40 digits of the
corresponding y are not compatible with the digit multiset (a 0, or a digit outside it, or too many of
one digit), it is not a solution and is discarded. Rigorous because the lattice step lists every
exponent vector compatible with the configuration.

## 4. Bases with one free prime (prove_cum_b.py; csrc/topb, csrc/bottomb)
Base b in {3,4,6}: y = B q^c with B a bounded product of the primes dividing b (Fonga's recursion in base b).
Same cumulative structure with eps0 = (b-2) b: 0 < Lam_j < eps0 * b^-D_j; trailing run |Lam_k| < 2 b^-L;
Matveev with n = 3. Top test: one 2-dim lattice, one target per cofactor B, exact level-1/level-2 kills
in C (survivor sets identical to exact enumeration on 3 configuration spaces x all cofactors: 0 mismatches;
planted 20/20 in base 6). Bottom test: the cofactor is forced by the valuations of y mod b^M; the
congruence q^c = rho (mod p^n) is solved with p-adic logarithms and generalized CRT (C port of
basesolver.solve_c; agrees with Python on 7 configuration spaces; planted 160/160 in bases 3, 4, 6).
Bases 3 and 4 re-proved with this method (same trees as the gap-chain proofs).

**Which bases this method reaches.** The bottom test needs the b-adic congruence y = F (mod b^M) to pin down the
free exponents. Modulo b^M the unit group of Z/b^M has one independent logarithm for each prime p | b, so it can
determine at most that many free exponents. Free primes (primes below b not dividing b, after excluding those that
cannot occur) and primes dividing b: base 3 (1 free: 2; p = 3), base 4 (1: 3; p = 2), base 6 (1: 5; p = 2, 3),
base 10 (2: 3, 7, with 5 excluded by Lemma 0; p = 2, 5). In bases 5, 7, 8, 9 there are more free primes than
primes dividing b (base 8: 3, 5, 7 against 2), so only top-side tests remain and the enumeration is far larger.

**Candidate resolution on the bottom side.** When the bottom congruence admits an exponent in the box, that
exponent is unique (base 10: one point per admissible class, base b: c is fixed modulo N > Xc), so the only
possible y is explicit. It is discarded when one of its digits at heights M..M+39 is 0 or not among the
digits left after removing the configuration from the multiset. If N <= Xc (rare, large cofactor
valuations), every exponent of the class up to Xc is checked the same way (at most 200,000, else the case
stays open). Planted tests: the planted exponent is reported in 160/160 cases (bases 3, 4, 6) and genuine
numbers are still detected 80/80 (base 10). Example: for node 80000 in base 6 at an intermediate stage
(Xc ~ 5.7*10^28) one configuration admitted c0 = 2529792095499578018590755839; the digits of 4*5^c0 just
above the configuration contain 0s, so it is not a solution.

## 5. Results of the second session
*Base 6 (all nonzero final digits) is complete:* trees of sizes 1, 10, 3, 14, 1 for final digits 1-5, heights
0, 4, 1, 4, 0, so every positive integer whose base-6 digit-product orbit ends at a nonzero digit has
persistence at most 5 (final digits 2 and 4), 2 (final digit 3) or 1 (final digits 1, 5). The hardest node,
32768 = 2^15, has digit sets with up to 15 non-1 digits; its 15-digit set took 136 proof steps. Bases 3 and 4
were re-proved with the cumulative method (same trees as the gap-chain proofs). See smallbases.md.
*Base 10:* 91 nodes of the even-final-digit trees have complete preimage sets (31 more than before this
session); final-digit trees 2, 4, 6 are determined to depth 3, tree 8 to depth 2. Every open node at the next
depth has 13-19 non-1 digits (estimated 1-2 days each with the present method). See results.md.

Engineering notes. Long proof steps are resumable at the level of (ordering, first gap) on the top side and
(ordering, lowest height) on the bottom side; batched and unbatched runs give identical case counts,
survivor sets and verdicts. The lattice constant's precision is now max(2400, bitlen(C)+200) bits (a safety
assertion caught the fixed 2400 bits being too small for a 13-digit set at the a-priori stage; no result was
affected).

# Third part: the chain sieve (final digit 0 / persistence 12)

See gsieve/CHAIN_SIEVE.md. Result: every 7-smooth number of digit length ≤ 40,000 contains a 0, or a 5
together with an even digit, or is one of 1,117 known numbers (≤ 49 digits, persistence ≤ 10). Hence no
number below 10^40000 has multiplicative persistence ≥ 12 (previous record 2.67·10^30000, Peters 2023).
6.4·10^13 candidates, ≈ 3 h on 2 cores; code, validation suite and logs in gsieve/chain_sieve.zip.
