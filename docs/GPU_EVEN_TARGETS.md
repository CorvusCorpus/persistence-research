# Even-target tree nodes on the GPU

Design notes for moving the case loops of the base-10 even-target prover (`work/prove_cum.py`, C loops `csrc/topc3`
and `csrc/bottomc3`) to the RTX 3060 Ti. Status: **both implemented and validated on the GPU** (bottom solver `csrc/bottomgpu`,
top filter `csrc/topgpu`; NOTEBOOK E.3–E.5): byte-identical output to bottomc3/topc3 on whole configuration spaces, planted
solutions found, kernel tests at production scale. Nothing here changes a claim in `docs/RESULTS.md`.

## 1. Where the time goes

### 1.1 Cost of the open nodes (`work/node_costs.py`, log `work/logs/node_costs.out`)
`prove_cum.prove_bag(bag, dry=True)` replays the bootstrapping plan of a bag (every test assumed to pass at its
threshold) and adds up its cost model: one unit per top case, `BOTTOM_COST = 4` units per bottom solve. The split by
side is recovered from the facts before and after each step and must add up to `res['stats']['cost']`.

| node | target | bags | max k | cost (units) | top | bottom | already proved |
|---|---|---|---|---|---|---|---|
| 131712 | 8 | 15 | 11 | 5.66·10⁹ | 1.46·10⁹ | 4.21·10⁹ | 14 bags (2.61·10⁹) |
| 122472 | 6 | 26 | 11 | 7.88·10⁹ | 2.06·10⁹ | 5.82·10⁹ | 23 bags (2.94·10⁹) |
| 24192 | 6 | 39 | 11 | 8.22·10⁹ | 2.13·10⁹ | 6.08·10⁹ | 0 |
| 1229312 | 2 | 12 | 13 | 6.70·10¹⁰ | 2.31·10¹⁰ | 4.39·10¹⁰ | 0 |
| 314928 | 2 | 50 | 13 | 1.33·10¹¹ | 4.73·10¹⁰ | 8.55·10¹⁰ | 0 |
| 41472 | 6 | 84 | 13 | 1.40·10¹¹ | 4.99·10¹⁰ | 8.97·10¹⁰ | 0 |
| 93312 | 2 | 85 | 13 | 1.72·10¹¹ | 6.11·10¹⁰ | 1.11·10¹¹ | 0 |
| 221184 | 6 | 110 | 16 | 2.23·10¹³ | 7.18·10¹² | 1.51·10¹³ | 0 |
| 1161216 | 4 | 118 | 16 | 1.39·10¹⁴ | 4.10·10¹³ | 9.84·10¹³ | 0 |
| 11239424 | 2 | 27 | 18 | 1.32·10¹⁴ + 2 bags | 3.84·10¹³ | 9.35·10¹³ | 0 |
| 33718272 | 8 | 51 | 19 | 1.58·10¹⁵ + 6 bags | 4.27·10¹⁴ | 1.15·10¹⁵ | 0 |

"+ n bags": `prove_cum.matveev_cum` cannot even form the a-priori Matveev bound for bags with k ≥ 18 (its fixed-point
iteration starts at 10³⁰⁰ and the bound exceeds the float range; `assert nb <= B`). Those bags need the bound computed in
log space before any estimate is possible.

### 1.2 Profile (`work/prof_steps.py`; cycle counters in instrumented copies of top3.c and bottom3.c)
Representative bag: 2·2·2·2·2·2·2·3·3·3·7 (k = 11, the largest bag of 24192, 3.05·10⁹ units), the eight costliest steps of
its plan. TSC cycles per case / per solve (the machine was loaded by another job at the time, so absolute times are
pessimistic; a re-measurement on the idle machine is in NOTEBOOK E):

| step | prefix | per-case MPFR log | level 1 | level 2 |
|---|---|---|---|---|
| top j=3, G=197 (C 956 bits) | 105 | 5,179 | 2,215 | 487 |
| top j=4, G=123 (588 bits) | 193 | 3,470 | 1,242 | 128 |
| top j=4, G=87 (408 bits) | 193 | 2,689 | 761 | 92 |
| top j=5, G=51 (228 bits) | 2,150 | 2,303 | 458 | 36 |

| step | 5-adic series | 2-adic | mt + β (CRT) | lattice | low(Flow), per solve |
|---|---|---|---|---|---|
| bottom j=7, G=160 | 2,479 | 807 | 1,090 | 2,206 | 5,390 |
| bottom j=7, G=123 | 2,491 | 893 | 1,157 | 2,378 | 5,918 |
| bottom j=6, G=87 | 1,990 | 827 | 1,076 | 2,147 | 5,766 |
| bottom j=6, G=51 | 1,985 | 955 | 1,057 | 2,320 | 6,386 |

A bottom solve costs about 12,000 cycles, a top case 3,700–8,000, so in time the bottom side is 52–63% of every node
(more than the cost model's units suggest, fewer than 4×). Inside a solve, `low(Flow)` — two p-adic logarithms and two
inverses per lower configuration, computed at precision Mmax = 4G + 8 although the solves need only 5^G — is about half.
Inside a top case the MPFR logarithm is about two thirds and the exact level-1 test a quarter.

## 2. The bottom solver on the GPU (`work/csrc/bottomgpu*.{cu,c,h}`, `bigfix.h`) — implemented

Drop-in replacement for `bottomc3` in `cum` mode: same job file, same output (CAND lines in bottomc3's order and format,
PROG per ordering, DONE nsolve ndirect5 ndirect2). The host part includes `bottom3.c` unchanged.

**Method.** One thread per lower configuration (Flow = −1 + Σ c_q 10^{a_q}), as in bottom3.c's `cum_rec`: low(Flow) once,
then every height p of the top digit. Every quantity is the residue bottom3.c computes, mod 5^M (M = G for every solve of
a covered job) and mod 2^n (n = M − t), obtained without exact divisions:
- log1p(p^v y) = Σ_n (−1)^{n+1} p^{nv−e} yⁿ/n′ (n = pᵉn′) by Horner in y with tabulated coefficients, over n with
  nv − e < M; bottom3.c sums xⁿ/pᵉ mod p^{M+12} instead. Both are the exact value of the convergent series mod p^M.
- plog(x) = log(x^{p^s})/p^s = Σ (−1)^{n+1} p^{n·vz−e−s} z′ⁿ/n′ with x^{p^s} = 1 + p^{vz} z′ (powering mod 5^{G+s}).
- Flow⁻¹ mod 5^G and mod 2^K by Newton; T/5 and z/5^{vz} as products with 5⁻¹, 5^{−vz} mod 2^{32L} (exact quotients).
- CRT by Garner: β = β₂ + 2^{n−2}((β₅ − β₂)·2^{−(n−2)} mod 5^{G−1}).
- 5-adic products are Montgomery products (CIOS) with one R = 2^{32 L} for the three moduli 5^{G+s}, 5^G, 5^{G−1};
  tables hold one operand in Montgomery form so that products come out as plain residues. Every loop bound is a
  compile-time constant (five size classes: L = 6, 8, 11, 13, 16 limbs for G ≤ 60, 90, 125, 165, 200), so limbs stay in
  registers; runtime shifts use a barrel of compile-time limb offsets.
- Lattice step: the coset point nearest the box centre has coordinates ⌊K − off0·R⌋ / 2^F in the Gauss-reduced basis
  (R = b/det, K tabulated; F = bits(Q) + 64, error < 2^−64). When the fraction is within 2^−60 of an integer the solve is
  handed back; the candidate coordinates are then exact two's-complement products.

**Handed back to the host** (bottom3.c's own `low()` and `solutions()`): ambiguous lattice estimates, v2(Flow) beyond the
2-adic width, and any t for which bottom3.c's `mt()` would stop ("lattice too fine"). **Handed to bottomc3 itself**
(exec): other modes, jobs with heff[0] > G (M would vary), G > 200 (only cheap steps in practice).

**Exactness argument.** Each step computes the same residue of the same exact expression as bottom3.c (the series and
plog identities above, with all truncated terms ≡ 0 mod p^M), so β is identical; the class tests and the lattice step are
the same integer computations (the floor via a certified estimate, or the host). `work/bottom_ref.py` holds a Python
transcription of bottom3.c next to the restructured form; its self-check compares them.

**Validation** (`work/val_bottomgpu.py`; first on the CPU emulation build `bottomgpu_emu`, then on the GPU build with the same
results and 3,900 kernel queries, NOTEBOOK E.5):
- K, kernel tests at production scale: G = 14 … 200 in all five size classes, boxes U, W ≈ 2·XL up to 10⁹⁴ with
  G = ⌈2 log10 XL⌉ + 12; planted F = 2^t 3^{u′} 7^w mod 10^G split as c·10^p + Flow (all p, including p = 0 and small p,
  which take the direct 5-adic and 2-adic paths) and arbitrary Flows. 2,600 queries: every planted point found
  (1,949/1,949), the kernel's candidate lists equal to bottom_ref's exact lists, no deferrals.
- C, configuration spaces: 57 whole steps (and slices) of the bootstrapping plans of five bags, G = 51 … 383: output
  byte-identical to bottomc3 (47 through the kernel, 10 handed to bottomc3 for G > 200).
- P, genuine numbers: 60 configurations built from the digits of genuine zero-free numbers (work/Z.json): byte-identical
  to bottomc3, the genuine solution among the CAND lines every time.

**Integration.** `work/cbottom5.py` runs a whole step in one process (one per (ordering, height group) would be dominated
by CUDA start-up), resumable at ordering boundaries with the same checkpoint format as cbottom4. `work/prove_cum_gpu.py`
runs prove_cum's prove_bag unchanged with cbottom5 as its bottom module and proves the bags of a node in parallel
processes; `work/val_prove_gpu.sh` proves a node from scratch with both bottom paths and compares every bag's facts,
bound and counters and the proofC file.

## 3. The top side on the GPU (`work/csrc/topgpu*.{cu,c,h}`) — implemented

Implemented as designed below, except that the prefix logarithms come from host tables at level 1 and from the series at
levels ≥ 2, with prefixes whose heights fall below 5 handed to the host. Validation (`work/val_topgpu.py`, NOTEBOOK E.5): ψ
errors at most 2 units of 2^−64 against arb (the margins assume 4), 79 jobs byte-identical to topc3, 150/150 planted genuine
configurations kept. Speed: 100× one core at G = 123, 30× at G = 51; steps with Θ > 2^592 (G ≥ ~180) go to topc3 because
nvcc 13.3 crashes on the L = 24 instantiation.

**Exact test (top3.c).** For a case (ordering c_1..c_j, gaps), m = round(v) with v = −C·log(1+z) (MPFR); level 1 kills
when dist(m·θ, ℤ) > ρ, θ = qa/qb = (B⁻¹)ₙₙ, ρ = √(R²/g_{n−1}); level 2 is the next Fincke–Pohst level with
θ₂ = a2/D and μ = P/Q.

**Key observation.** θ and θ₂ are tiny (B is LLL-reduced with entries ~C^{1/3}, so B⁻¹ entries ~C^{−1/3}: e.g. θ ≈ 2^−197
at G = 123). With m = v + δ, |δ| ≤ 1/2 + 2^−50, m·θ = −Θ·log(1+z) + δθ with Θ = C·θ (≈ 2^390 at G = 123) and |δθ| < 2^−196.
So both levels depend only on frac(Θ·log(1+z)) and frac(Θ₂·log(1+z)), up to a negligible, bounded term — m itself is not
needed to kill a case.

**Design: GPU filter + exact host path, identical decisions by construction.**
- The GPU computes ψ₁ ≈ frac(−Θ·λ), ψ₂ ≈ frac(−Θ₂·λ), λ = log(1+z), with a rigorous bound ε (fixed point, F = bits(Θ)
  + 64 + guard fraction bits; every truncation and series tail counted). It kills a case only when the exact test is
  certain to kill: |s̃| > ρ + ε + |θ|, or (d̃₂ − ε₂)²g₂ + (|s̃| − ε)²g₃ > R² for level 2 (s̃ the signed distance of ψ₁, d̃₂
  that of −ψ₂ − μ s̃, ε₂ = ε(1 + |μ|)).
- Everything not killed (0.02–0.05% of cases survive level 2 at the profiled steps, 1.2–2.4% level 1) goes through
  top3.c's own per-case code on the host (MPFR m, exact level 1 and 2, the same run merging and SURV lines). The output
  and the DONE counters are then those of topc3: a killed case is one topc3 kills at level 1 or 2, and all others take
  the identical code path.
- λ: per prefix (all digits but the last) λ_A = log A and 1/A, computed level by level down the prefix tree
  (λ_i = λ_{i−1} + log1p(c_i 10^{−S_i}/A_{i−1}), 1/A_i by the matching geometric series); per case
  λ = λ_A + log1p(c_l·10^{−e}/A), whose series needs only N(e) ≈ (bits(Θ) + 64)/(3.32 e) terms (1–3 for the typical
  e ≈ G). Then two "middle products" with Θ and Θ₂ (only the product bits around 2⁰ are needed).
- Estimated GPU cost: a few thousand 32-bit multiply-adds per case including the amortised prefix work, against about
  5,000 cycles per case for topc3 on one core.

**Validation plan.** Kernel tests: ψ₁, ψ₂ against MPFR at 2× precision on random and boundary cases, the claimed ε checked
on every one; SURV output and DONE counters byte-identical to topc3 on whole steps; the planted top tests
(`planted_top_cum3.py`: genuine numbers must survive) through the new path.

## 4. Run plan
1. GPU validation of bottomgpu (val_bottomgpu.py K, C, P with the CUDA build), benchmark, `val_prove_gpu.sh` on proved
   nodes (1344, 6272): identical facts and proofC.
2. The three cheap open nodes (131712: one bag left; 122472: three; 24192: 39) with prove_cum_gpu (GPU bottom, CPU top,
   bags in parallel): about an hour. Claims only after the proofs complete and the outputs are checked.
3. The k = 13 nodes need the top filter (their largest bags have ~10¹⁰ top cases each, about 5 core-hours per bag on the
   CPU even with the bottom on the GPU).
