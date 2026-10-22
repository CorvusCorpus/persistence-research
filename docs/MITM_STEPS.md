# Meet-in-the-middle steps for the even-target prover (design and cost model)

Research note (session 5, 22 Oct 2026). Status: **cost model only; nothing implemented, nothing in `docs/RESULTS.md` depends
on it.** It builds on the two-chain bootstrapping (`docs/BOOTSTRAP_TWO_CHAINS.md`), which is implemented.

## 1. Where the cost is after two chains
`prove_cum` proves facts D_j < G (top) and a_i < G (bottom) at G = ⌈2 log10 XL⌉ + 12, where XL is the current bound on
the length L. A top step at level j tests every configuration of the top j digits, about C(G, j)·(orderings) of them. The test
is one lattice computation per configuration. The margin 12 makes a lattice point in the box (|u'| ≤ U ≈ 2.1·XL, |w| ≤ W ≈ 1.2·XL)
improbable (about 10⁻¹² per configuration).

With two chains the costly steps sit in the middle and last rounds: thresholds 40–110, 6–10 digits on one side. Example,
2¹³·3³ (the largest bag of 221184, k = 16): the last four rounds (top levels 6–9 at G = 80, 61, 47, 38) are 7·10¹⁰ of the
bag's 1.0·10¹¹ units.

## 2. The MITM step: top configurations against a table of bottom windows
**Statement.** Fix a top level j, a threshold Gs and a window Gp. Let T run over the configurations of the top j digits
(ranges from the top facts), and let B run over the bottom windows: all ways the digits of height < Gp can look (the lowest
b digits, any b, the next digit at height ≥ Gp). For a solution with D_j ≥ Gs and L ≥ (depth of digit j) + Gp, both of these
hold, with T and B the solution's own configurations:
- **top:** 0 < Λ_T = t log 2 + u' log 3 + w log 7 − L log 10 − log(1 + z_T) < 80·10^(−Gs);
- **bottom:** 2^t 3^u' 7^w ≡ F_B (mod 10^Gp), F_B = −1 + Σ_{m∈B} c_m 10^(a_m). So t = v₂(F_B) when v₂(F_B) < Gp, and
  3^u' 7^w ≡ F_B 2^(−t) (mod 5^Gp).

If no pair (T, B) admits a common (t, u', w, L) in the box, then D_j < Gs, or L is below (top window) + Gp. The second
alternative is a few dozen digits and enters the bootstrapping as max(·, XL).

**Why it is cheaper.** Each relation alone needs about 2 log10 XL + 12 digits of precision. Together they need only
Gs + 0.699·Gp ≈ 2 log10 XL + (margin), so the top threshold drops by about 0.7·Gp. Pairs are not enumerated:
- Λ_B = {(u', w) : 3^u' 7^w ≡ 1 (mod 5^Gp)} is a lattice of index 0.8·5^Gp (3 is a primitive root mod 5^n).
- Each B fixes a coset v_B + Λ_B. Write (u', w) = v_B + λ, with v_B reduced into a fundamental domain, so that λ lies in a
  slightly enlarged box.
- The top condition is linear in λ: ℓ(λ) ∈ [x_T − ℓ(v_B), x_T − ℓ(v_B) + ε) mod 1, where ℓ(u', w) = (u' log 3 + w log 7)/log 10,
  x_T = (log(1 + z_T) − t log 2)/log 10 and ε = 80·10^(−Gs)/log 10.
- That is a closest-vector problem in a 3-dimensional lattice. After LLL it reads: P_T − Q_B ∈ E (mod ℤ³), where
  P_T = κ + x_T·θ and Q_B = (ℓ(v_B) + t log2/log10)·θ, with θ, κ ∈ ℝ³ fixed per step. E is a small parallelepiped of volume
  ≈ 108·10^(2 log10 XL − Gs − 0.699 Gp) (about 1/det).
- A hash join on a grid of cells the size of E finds every close pair in N_T + N_B operations. Close pairs are rare false
  pairs (expected N_T·N_B·vol E) and are decided exactly: a 3-dimensional Fincke–Pohst enumeration, then the candidate's digits.

The bottom MITM step (bottom configurations of level i against a table of top windows, depth < Gt) is symmetric. Its own
side carries 5-adic precision 5^Gs. Only the 5-adic congruence is used: bags with many 2s allow t ≤ 55 (2¹³·3³), which leaves
the 2-adic part almost nothing at Gp < 55.

**Soundness.** The relations are prove_cum's (F3') and the bottom congruence. The pair set found by the hash join must be a
superset of the true pairs: cells are searched with E enlarged by the rigorous error of P_T and Q_B, computed in fixed point
with the same error accounting as topgpu. Each candidate pair is then decided by exact arithmetic. Validation would follow
the existing pattern:
- planted genuine configurations from Z.json must be found;
- on small parameters, the MITM step's facts must agree with brute force over all (t, u', w) of a short length range;
- independent re-proofs on the CPU path.

## 3. Deep-digit split inside one step
In a standard top step at threshold G, the digits deeper than h = ⌈(G + 5)/2⌉ enter linearly:
log(1 + z_A + z_B) = log(1 + z_A) + z_B/(1 + z_A) + e, where z_A comes from the digits shallower than h. Since
|z_B| < 80·10^(−h), |e| ≤ 3.2·10^(3−2h) ≤ 3.2·10^(−G−2). That is 1/2500 of the test's window 80·10^(−G), and it is added to
the test radius.
- For each configuration A of the shallow digits, the deep digits' contributions to the torus point are vectors
  v(c, d) = θ·c·10^(−d)/(1 + z_A), one per position and value.
- A configuration B = B₁ ∪ B₂ of the deep digits passes only if θ·log(1 + z_A) + Σ_{B₁} v + Σ_{B₂} v ∈ E. Hash the B₂ sums,
  look up each B₁ sum, and check the position order (last of B₁ < first of B₂) on the matches.
- Cost per A: (G − h)·(values) + |B₁| + |B₂| instead of |B₁|·|B₂|.
- The 5-adic analogue works for bottom steps: high digits (height ≥ ⌈Gs/2⌉) enter log F linearly mod 5^Gs.
- The torus test must use at least three coordinates: the two-level filter of topgpu lets 0.02–0.05% of cases through, too
  many false matches for a join.

## 4. Cost model (`work/plan_mitm.py`; dry runs, round by round)
Each round re-proves one side's levels up to its frontier at one threshold. It chooses the side, the threshold and the kind
of step with the smallest cost per digit of log10 XL gained. The costs are in `node_costs.py`'s units:

| Item | Cost (units) |
|---|---|
| top configuration | 1 |
| bottom configuration or bottom-window entry | 4 |
| top-window entry | 1 |
| false pair | 1,000, and false pairs may cost at most as much as the step itself |

The largest bag of each open node (two chains throughout; "standard" is the same planner with prove_cum's steps only):

| bag (node) | k | today | two chains, standard | + MITM | + deep split | + either, per round | vs today |
|---|---|---|---|---|---|---|---|
| 2¹³·3³ (221184) | 16 | 5.16·10¹² | 1.02·10¹¹ | 2.30·10¹⁰ | 2.26·10¹⁰ | 1.31·10¹⁰ (7.8×) | 394× |
| 2¹¹·3⁴·7 (1161216) | 16 | 3.49·10¹³ | 8.84·10¹¹ | 1.55·10¹¹ | 2.10·10¹¹ | 9.83·10¹⁰ (9.0×) | 355× |
| 2¹¹·4²·7³ (11239424) | 16 | 6.17·10¹³ | 1.62·10¹² | 2.83·10¹¹ | 3.85·10¹¹ | 1.78·10¹¹ (9.1×) | 347× |
| 2¹³·3·4·7³ (33718272) | 18 | no estimate | 3.54·10¹³ | 1.01·10¹³ | 6.72·10¹² | 4.51·10¹² (7.9×) | — |

The two kinds of step help in different rounds:
- The MITM step pays off where the other side's windows are cheap, i.e. in the last rounds. For 2¹³·3³ the last round drops
  from 2.1·10¹⁰ to 1.6·10⁸.
- The deep split pays off where many digits sit deep, i.e. in the middle rounds, on both sides.
- With both, the middle rounds still dominate. For 2¹³·3·4·7³ the planner uses deep-split rounds 5–15 and MITM rounds 16–18;
  rounds 12–16 (thresholds 61–112) are 4.1·10¹² of the 4.5·10¹² units.

## 5. Other things tried
- **Closing step** (every digit in one of two windows, L < Gt + Gb in one step): never competitive. Each window must hold
  configurations with any number of digits.
- **Bottom threshold max(G, tm + 8):** dropping the tm + 8 floor (per-configuration t instead) gains 1.7× on 2¹³·3³.
- **Margin 12:** with candidate handling at about 10⁴ units each, the optimum is about 9 for tm = 55, worth 1.3–1.5×.
  It is worth more for bags with small tm.
- **Scaling.** None of this changes the exponential growth in k. The cost of a frontier step is still about C(G, j), with G
  falling about 14 per round. Bags 3^(k/2)·7^(k/2) (`plan_mitm.py`, odd digits, so tm = 0):

  | k | standard (two chains) | MITM or deep |
  |---|---|---|
  | 16 | 1.1·10¹¹ | 2.3·10¹⁰ |
  | 20 | 3.2·10¹³ | 3.6·10¹² |
  | 24 | 9.6·10¹⁵ | 5.3·10¹⁴ |

  So four more digits cost about 150× more, even with both steps. Nodes deeper in the trees (up to 44 non-1 digits) need a
  different idea.

## 6. Order of implementation (if pursued)
1. The bottom-window table (C, 64-bit arithmetic for Gp ≤ 27: 5-adic discrete log, coset reduction, Q_B, cell index).
2. A top kernel variant that outputs P_T (three coordinates) and probes the table on the GPU, with the existing prefix
   machinery and error accounting.
3. Exact decision of candidate pairs (Python or C Fincke–Pohst in 3 dimensions), then the existing could_be_solution check.
4. A planner in prove_cum that chooses (Gs, Gp) per round; facts and records as now, plus the step type.
5. Validation as in §2, then re-proofs of proved nodes and comparison of facts.
