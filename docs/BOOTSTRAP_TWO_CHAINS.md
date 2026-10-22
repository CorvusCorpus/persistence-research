# Two-sided bootstrapping: a p-adic Baker chain for the bottom digits

Research note (session 5, 22 Oct 2026). Status: **implemented in `work/prove_cum.py` (`twochain_cum`, selected by
`PC_BOOTSTRAP=twochain`) and validated (§6); no committed record uses it yet.** Nothing in `docs/RESULTS.md` depends on it.

## 1. Why the current bootstrapping is expensive
In `prove_cum.matveev_cum` the bound on the length L comes from L < D_j + a_{j+1} for some level j.
- The depth D_j comes from a top fact (D_j < top[j]) or from Matveev's chain.
- The height a_{j+1} comes **only from a bottom fact** (a_{j+1} < bot[j+1]).

Matveev's bound for the next level multiplies the previous level's bound by about 10^18 (C(5)·β·Π). The plans therefore run in
*rounds*: within a round all levels are re-proved at one threshold G until the top frontier (levels with top facts) and the
bottom frontier (levels with bottom facts) gain one level. With `gap` levels between the frontiers,

  XL ≈ 10^(18·gap),  G = ⌈2 log10 XL⌉ + 12 ≈ 36·gap + 15.

Example (bag 2⁹·7⁴, k = 13; `prof_steps.plan_states`): G runs through 486, 458, …, 197, 160, 123, 87, 51. Each round's most
expensive step is the new frontier level j, with about C(G, j)·(orderings) cases. Those steps (top j = 4–6 and bottom j = 7–9 at
G = 87–197) carry almost all the cost.

## 2. A Baker chain from the bottom
**Lemma (bottom chain).** Let 0 ≤ i < k, F = −1 + Σ_{m>i} c_m 10^{a_m} (the digits below digit i; c_m = 9(d_m − 1)), with
a_{i+1} ≤ b_{i+1}. Then

  a_i < (n+1)·min(C₁·Ω₁, C₂·Ω₂)·max{log B, log 5},

where B ≥ max(u′, w, 3) is the coefficient bound of the bootstrapping, and C₁, C₂, Ω₁, Ω₂ are the quantities of Yu's
Theorems 1 and 2 below (n = 3, p = 5). In particular Ω_i ≤ log 3 · log 7 · max(h(F/2^t), small), with
h(F/2^t) ≤ max(log 80 + b_{i+1} log 10, t_max log 2).

*Proof.*
- **F is a 5-adic unit.** Its last digit is never 0 or 5, by Lemma 0. Indeed F ≡ −1 + c_k (mod 10) with c_k ≡ 4(d_k − 1) (mod 5), or F ≡ −1 (mod 10) when a_k = 1.
- **The valuation.** 9y = 2^t 3^{u′} 7^w and 9y − F = 10^L + Σ_{m≤i} c_m 10^{a_m}, every term divisible by 10^{a_i}. So
  a_i ≤ ord₅(9y − F) = ord₅(Ξ − 1) with Ξ = 3^{u′}·7^{w}·(2^t/F)^{1}.
- **Yu's hypotheses.** Ξ ≠ 1, because 9y − F = 10^L + … > 0. The three numbers 3, 7, 2^t/F are 5-adic units, which is (1.16). Hypothesis (1.5) concerns only p = 2. The rank r of ⟨3, 7, 2^t/F⟩ is 2 or 3, and both are covered: take b = {3, 7}, or b = all three.
- **Theorems 1 and 2** (Yu, *p-adic logarithmic forms and a problem of Erdős*, Acta Math. 211 (2013), 315–382, §1.2) give
  ord_p(Ξ − 1) < C_i*(n, d, p, b)·Ω·max{log B, f_p log p}, C_i* = (n+1)·C_i (1.22), (1.26), for K = ℚ (d = 1), p = 5: case (III.2),
  e_p = f_p = 1, q = 2, u = 1 (μ(ℚ) = {±1}, so α₀ = −1). Here B ≥ max{|b₁|, …, |b_n|, 3} (1.23), and h₀ is the absolute
  logarithmic Weil height, so h₀(3) = log 3, h₀(7) = log 7 and h₀(2^t/F) = max(t log 2, log |F|).
  - **δ.** By (1.6), (p^f − 1)/δ(b) = |⟨α₀, ϑ̄₁, …⟩| is the order of a subgroup of 𝔽₅* that contains −1, so δ(b) ∈ {1, 2}.
    δ = 1 gives the larger constant, and the bound below uses it.
  - **ϰ.** (1.2) fixes ϰ by p^{ϰ−1}(p−1) ≤ 2e_p < p^ϰ(p−1); for p = 5, e_p = 1 this gives ϰ = 0.
  - C₁ = c⁽¹⁾ (a⁽¹⁾)ⁿ nⁿ (n+1)^{n+1} / (n! q^u log p) · max{p/(δ (log p)^{n+1}), eⁿ/nⁿ} · max{log(e⁴(n+1)), 1, log p} (1.9), with
    (a⁽¹⁾, ϰ₁) = (7·(p−1)/(p−2), 9·(p−1)/(p−2)) and c⁽¹⁾ = 1790. This gives C₁ = 2.091·10⁹.
  - C₂ = (c⁽²⁾/p^ϰ) (a⁽²⁾ e_p p^ϰ)ⁿ (n+1)^{n+1} / ((n−1)! q^u (log p)³) · max{p/δ, eⁿ/nⁿ·(log p)^{n+1}} · max{log(e⁴(n+1)), 1, log p}
    (1.10), with (a⁽²⁾, ϰ₂) = (7, 25) and c⁽²⁾ = 557·(p−1)/(p−2). This gives C₂ = 1.053·10⁸.
  - Ω = min_b Ω(b), Ω(b) = Π_{α∈b} h₀(α) · Π_{α∈a∖b} h⁽ⁿ⁾(α) (1.21), where b runs over the multiplicatively independent subsets of
    a = {3, 7, 2^t/F} of size r = rank⟨a⟩, and h⁽ⁿ⁾(α) = max{h₀(α), max{n, log p}/(ϰ₁(n+5))} (1.20, Theorem 1) or
    max{h₀(α), max{n/log p, 1}/(ϰ₂ p^ϰ)} (1.25, Theorem 2). If r = 3 then b = a; if r = 2 then b = {3, 7} is admissible. In both
    cases Ω ≤ log 3 · log 7 · h⁽ⁿ⁾(2^t/F). C_i(n, d, p, b) depends on b only through δ(b), and δ = 1 gives its largest value. ∎

The values above were checked against the typeset pages of the paper ((1.2), (1.6), (1.9), (1.10), (1.20)–(1.27) and the table
of §1.3 for case (III.2)); the PDF's text layer drops ϰ, which first gave a C₂ five times too small.

Per level the bound multiplies h(F) ≈ b_{i+1} log 10 by about (n+1)·C₂·log 3·log 7·log B·log 10 ≈ 10^11.5 (log B ≈ 170),
against about 10^17.9 for the top chain.

## 3. A 4-logarithm top chain
Matveev's chain uses the linear form in five logarithms t log 2 + u′ log 3 + w log 7 − L log 10 − log(1 + z_j).
- With N_j = 10^{S_j}(1 + z_j) ∈ ℤ (S_j = D_{j−1}, the depth of digit j) it is u′ log 3 + w log 7 − (L − S_j) log 10 − log(N_j / 2^t):
  four logarithms, because t ≤ t_max is small.
- The height of N_j/2^t is at most max(log N_j, t_max log 2), and log N_j ≤ D_{j−1} log 10 + log 9.
- Matveev with n = 4: C(4) = 1.4·30⁷·4^4.5 replaces C(5) = 1.4·30⁸·5^4.5. That is about 80× smaller, partly offset because the factor
  log 2 < 1 drops out of the product of heights. Per level: about 10^16.5 instead of 10^17.9.

## 4. Combined bound
L < min_j (D_j bound from the top + a_{j+1} bound from the bottom). Each bound comes from a fact where one exists, else from
its chain. With per-level factors F_top ≈ 10^16.5 and F_bot ≈ 10^11.5 the gap is shared:
log10 XL ≈ gap · 16.5·11.5/(16.5 + 11.5) ≈ 6.8·gap instead of 18·gap, so G ≈ 13.6·gap + 15.

## 5. Effect (dry-run plans in prove_cum's cost model)
Per bag, today → with both chains (Yu 2013, Theorem 2; 4-logarithm top chain):

| bag | k | today | two chains | gain | first G |
|---|---|---|---|---|---|
| 2⁷·3³·7 | 11 | 3.05·10⁹ | 2.51·10⁸ | 12.2× | 411 → 173 |
| 2⁹·7⁴ | 13 | 3.51·10¹⁰ | 2.57·10⁹ | 13.6× | 486 → 201 |
| 2¹¹·3³·4 | 15 | 4.39·10¹² | 1.15·10¹¹ | 38.1× | 561 → 226 |
| 2¹³·3³ | 16 | 5.17·10¹² | 1.10·10¹¹ | 47.2× | 599 → 240 |

Whole open nodes (`PC_BOOTSTRAP=twochain python3 node_costs.py --json logs/node_costs_twochain.json 221184 1161216 11239424
33718272`, log `work/logs/node_costs_twochain.out`; "today" from `work/logs/node_costs.out`). Times use the pool rate measured
on the idle machine with 10 workers (9.5·10⁷ units/s, `work/logs/pool_bench_221184.out`):

| node | bags | k ≤ | today (units) | two chains (units) | gain | time |
|---|---|---|---|---|---|---|
| 221184 | 110 | 16 | 2.23·10¹³ | 8.03·10¹¹ | 28× | ≈ 2.3 h |
| 1161216 | 118 | 16 | 1.39·10¹⁴ | 6.27·10¹² | 22× | ≈ 18 h |
| 11239424 | 27 | 18 | 1.32·10¹⁴ + 2 bags without estimate | 7.05·10¹² | — | ≈ 21 h |
| 33718272 | 51 | 19 | 1.58·10¹⁵ + 6 bags without estimate | 1.64·10¹⁴ | — | ≈ 20 days |

Today the bags with k ≥ 18 have no estimate, because the a-priori bound overflows the floats. With two chains the a-priori bound
is about the square root of today's, so every bag can be planned.

## 6. Implementation and validation
- **Code.** `prove_cum.twochain_cum(k, tm, top, bot)` computes the bound. `matveev_cum` dispatches to it when
  `PC_BOOTSTRAP=twochain`; the default (`matveev`) is unchanged. Only the bound on L changes. The facts, the tests, their
  thresholds G = ⌈2 log10 XL⌉ + 12 and the final L ≤ 5000 criterion stay as they were.
- **Records.** A bag proved with two chains gets `"bootstrap": "twochain"` in its `bagresC` record, and a fifth element
  `"twochain"` in its entry of the node's `proofC` file. Records without the field use Matveev's chain.
- **Facts do not depend on the bootstrap.** Every fact is proved for all (t, u′, w, L) in the box of a valid bound, so
  checkpoints and records of either kind can be mixed.
- **`work/val_twochain.py`** (ALL OK, 22 Oct 2026):
  - C: the constants re-evaluated with mpmath at 60 digits from the formulas. prove_cum's floats are at most 1.3·10⁻¹⁶ below
    them in relative terms (C4: 1.4 is not a binary fraction); the final bound carries a 1% margin.
  - R: the default bound equals a frozen copy of the original code on 3,000 random fact states.
  - D: the two-chain bound never exceeds Matveev's, on those states and on every state of three bags' plans.
  - G: 1,003 genuine zero-free 7-smooth numbers (no digit 5, k ≤ 20), in exact arithmetic:
    - the 5-adic units, the valuations a_i ≤ ord₅(9y − F_i) and the Yu bound at their exact heights (10,040 inequalities);
    - 0 < Λ_j < 80·10^(−D_j), N_j < 9·10^(D_(j−1)), and the 4-logarithm Matveev bound;
    - twochain_cum ≥ L for every partial set of their exact facts.
- **`work/val_twochain_e2e.sh NODE`** re-proves a proved node from scratch with two chains on the GPU, in an empty scratch
  directory. It checks the preimages against the committed proofC file, every bound ≤ 5000 and every bag's tag, and every
  fact against the genuine solutions. Results (`validation_pc/twochain/`, topgpu v2 `5ac4892e…`, 8 workers):

  | node | bags | k ≤ | preimages | two chains (units) | committed records (units) | gain | wall |
  |---|---|---|---|---|---|---|---|
  | 24192 | 39 | 11 | the 6 committed ones | 1.18·10⁹ | 8.2·10⁹ | 7.0× | 323 s |
  | 1229312 | 12 | 13 | none, as committed | 5.8·10⁹ | 6.7·10¹⁰ | 11.6× | 1,859 s |

  In both runs every bound is ≤ 5000 and every bag carries the tag. The 6 genuine solutions of 24192's bags satisfy every fact.
  Most of 1229312's wall time went to the final top steps at G = 38, which topgpu v2 handed to topc3 (|θ| ≥ 2⁻⁶⁰). In the
  two-chain plans these steps are about 30% of the work (221184: 2.4·10¹¹ of 8.0·10¹¹ units). topgpu v3 covers them
  (`docs/GPU_EVEN_TARGETS.md` §3).
