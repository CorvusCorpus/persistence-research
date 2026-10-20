# Session summary (19–20 Oct 2026): multiplicative persistence

A long working session, run in a 2-core cloud container (Xeon 2.1 GHz, AVX-512, 7 GB RAM), split over
several context windows. This file is the narrative baseline; `RESULTS.md` is the precise ledger and
`ROADMAP.md` the plan.

## 0. Background

- p(n) = product of decimal digits; persistence(n) = number of steps of p to reach a single digit.
- Records (OEIS A003001): a(11) = 277777788888899 (15 digits). Record lengths for a(1)…a(11):
  2, 2, 2, 2, 3, 4, 5, 7, 8, 10, 15 digits. Conjecture: 11 is the maximum.
- Every p(n) is 7-smooth; for persistence ≥ 3 the intermediate numbers are zero-free 7-smooth numbers
  not divisible by 10, i.e. of the forms 2^a 3^b 7^c or 3^a 5^b 7^c (b ≥ 1).
- The orbit of 277777788888899: → 4996238671872 → 438939648 → 4478976 → 338688 → 27648 → 2688 → 768 →
  336 → 54 → 20 → 0. Final digit 0, which is where all long orbits end.
- The orbits split by final digit. Odd final digits were settled by Brier, Clavier, Gutsche and Naccache
  (2021). Even nonzero digits reduce to finitely many S-unit equations, made effective via 2-adic bounds
  (Fonga 2026, arXiv 2608.27802). Final digit 0 is the hard case: it is equivalent to understanding
  zero-free 7-smooth numbers, which is in the same family as Erdős-type digit problems.

## 1. How the session went (questions asked, in order)

1. **"Disprove this tweet with a counterexample"** (a tweet about persistence 11). No counterexample
   exists in any searched range; the tweet was accurate.
2. **"Find a persistence-12 number or prove 11 is the maximum."** Neither is achievable with known
   mathematics. The final-digit-0 case contains an Erdős-type open problem. Work pivoted to rigorous
   partial results.
3. **"Be clever / prove E is finite."** E is the set of zero-free 7-smooth numbers. Its finiteness is open
   and out of reach; in base 3 the analogue is the open question whether 2^15 = 1122221122₃ is the last
   power of 2 with no ternary 0. Explored: density/counting bounds, digit-lifting automata, 2-adic trees.
   None gives finiteness.
4. **"Chip away with smaller proofs."** Built a rigorous prover for the even-final-digit trees in base 10
   (Section 3).
5. **C port + small bases.** Ported the inner loops to C, then proved complete theorems for bases 3, 4
   and 6 for all nonzero final digits (Section 4).
6. **Literature check.** Results appear new. The cumulative depth inequality is standard (Section 7).
7. **Base-10 compute here.** More nodes were proved. The streams are paused (Section 3.4).
8. **"Novel ways to dent final digit 0."** The chain sieve gave a new record, a(12) > 10^40000
   (Section 5).
9. **Chance, scaling, significance, GPU questions** (Sections 6, 8).

## 2. Data produced

- `work/Z.json`: all zero-free 7-smooth numbers below 10^5000.
  - Source: `zf7.c` (exhaustive OpenMP sieve on the low 36/72/144 digits, survivors in
    `logs/run5000.txt`) plus `verify.py` (exact).
  - 12,606 numbers ≥ 10. The largest is 2^25·3^227·7^28 (140 digits). None between 10^140 and 10^5000;
    Chaffin (OEIS) reports none up to 10^20000.
  - Persistence distribution: 1:13, 2:11994, 3:387, 4:142, 5:41, 6:12, 7:8, 8:5, 9:2, 10:2.
  - Zero-free depth D (the number of consecutive zero-free numbers in the orbit) is at most 9. The only
    numbers with D = 9 are 4996238671872 = 2^19·3^4·7^6 and 937638166841712 = 2^4·3^20·7^5.
- `work/gsieve/G_small.json`: the 1,117 chain numbers.
  - A chain number is in E and does not contain both a 5 and an even digit.
  - Largest: 49 digits. Persistence ≤ 10.
  - Every intermediate c₀…c₈ of a persistence-12 orbit must be a chain number.

## 3. Base 10, even final digits (2, 4, 6, 8)

### 3.1 Problem
A node z is a 7-smooth number in the tree of final digit d. Its children are all 7-smooth y with
p(y) = z. Writing y = (10^L − 1)/9 + Σ (d_i − 1)·10^(a_i) gives an S-unit equation
9y = 10^L − 1 + Σ 9(d_i − 1)·10^(a_i), y = 2^α 3^β 7^γ, for each multiset ("bag") of non-1 digits
with product z.

### 3.2 Method (`docs/METHOD_EVEN_TARGETS.md` has the details)
- **Upper bound:** Matveev's bound on L, via a cumulative chain of linear forms.
- **Cumulative inequality:** 0 < Λ_j < 80·10^(−D_j), with D_j = L − a_(j+1). This is standard in the
  Baker-method literature.
- **Top gaps:** reduced with lattice reduction on real logarithms (MPFR targets with rigorous error).
- **Bottom gaps:** reduced with p-adic logarithms (2-adic/5-adic) and Gauss-reduced coset lattices.
- **Bookkeeping:** "facts" top[j], bot[i] over compositions of gap patterns. Once the length bound is at
  most 5000 digits, `Z.json` finishes the job.
- **Survivor handling:** lattice survivors and congruence candidates are explicit numbers, discarded by
  their trailing or leading digits. WEAK classes are enumerated.
- **Code:**
  - Python drivers: `prove_cum.py` (base 10), `prove_cum_b.py` (bases b).
  - C inner loops in `csrc/` (topc, topc3, bottomc2, bottomc3, topb, bottomb).
  - Checkpoint/resume at bag, proof and partial level.
- **Validation:** planted-solution tests (`planted_*.py`), C-vs-Python agreement (`val_*.py`), a
  precision audit (`audit_M.py`), and end-to-end re-derivations (`val_cum_e2e.py`).

### 3.3 Results
- **92 nodes proved** (`work/certificate.json`; table in `docs/results_base10.md`).
- **Tree levels determined exactly:** targets 2, 4, 6 to depth 3; target 8 to depth 2.
- **Open blocking nodes:**
  - Target 2: 93312, 314928, 1229312, 11239424.
  - Target 4: 1161216. Estimated cost ≈ 1.4·10^14 units, about 40,000 core-hours.
  - Target 6: 24192, 41472, 122472, 221184 (≥ 1.4·10^13).
  - Target 8: 131712, 33718272. Deeper nodes there have up to 44 non-1 digits.
- **Node 8192** (target 6, depth 3) was completed on the morning of 20 Oct: no preimages.

### 3.4 Paused runs
- `runner4.py 8192 122472`: 8192 done, 122472 in progress.
- `runner4.py 131712 24192`: 131712 in progress.
- Checkpoints are in `ckptC/` and `bagresC/`. Resume by rerunning the same commands from `work/`.

## 4. Small bases 3, 4, 6 (nonzero final digits)

- **Theorems** (`docs/results_smallbases.md`):
  - Base 3: maximum persistence 3 (final digit 1) and 2 (digit 2).
  - Base 4: 3 (digit 2) and 1 (digits 1, 3).
  - Base 6: 5 (digits 2, 4), 2 (digit 3) and 1 (digits 1, 5).
  - All maxima are attained, and all agree with the previously observed numerical maxima, which had no proof.
- **Method:** "one free prime." In these bases y = B·q^c with the cofactor forced by valuations, so the
  same cumulative machinery applies. Base 10 also works; other bases fail on p-adic rank.
- **Validation:** a lattice survivor and a congruence hit came up in base 6, and both were resolved by
  explicit digits. WEAK enumerations kept the prover sound.

## 5. Final digit 0: the chain sieve and a new record

### 5.1 Key observation (chain lemma)
If n has persistence ≥ 12, then c₀ = p(n) is a chain number (otherwise p(c₀) would be ≡ 0 mod 10).
Also m(c₀) ≤ len(n), where m(x) = a/3 + b/2 + c + d for x = 2^a 3^b 5^d 7^c. So to show that there is
no persistence-12 number with at most N digits, it suffices to show that every chain number with
m ≤ N is known.

### 5.2 Algorithmic ideas
- **Suffix collapse.** If 2^K divides x = 2^a 3^b 7^c, then the last K digits of x equal those of 2^e,
  with e ≡ a + b·log₂3 + c·log₂7 (mod 4·5^(K−1)). So one bitmap of 4·5^12 bits tests the last 13 digits
  of every candidate, scanned 64 exponents per word along each row (b, c).
- **Exact branching.** With digits restricted to exclude 0 and 5 (even numbers), every lift keeps exactly
  4 of 5 options, so there are exactly 4^k admissible k-digit endings (4^13 = 67,108,864 confirmed).
- **Affine lift.** Digits 14–26 are h₀ + 2^13·((A_hi + A_lo·w·j) mod 5^13), with j constant along a row.
  This is computed with AVX-512 (Montgomery multiplication mod 5^13, IFMA divisions, gathers from small
  tables).
- **Leading digits.** A 64-bit fixed-point frac(log10 x) and a 2^18-bin conservative table cover the top
  5 digits.
- **Later stages.** A base-10^9 limb engine checks digits 27–126; GMP then checks the last 600 digits
  and the full number.
- **Odd branch (3^a 5^b 7^c).** All digits must be odd, so exactly one lift survives at each step and
  there is a unique 2-adic thread. For b ≥ 40 this gives O(1) work per (a, c); smaller b is stepped
  directly.

### 5.3 Validation
- Every stage decision was compared with exact big-integer arithmetic on full regions (N = 150, G and
  zero-free modes, K = 12 and 13).
- In zero-free mode, the run end to end (N = 160) rediscovers exactly the 314 known large zero-free numbers.
- The direct and odd paths agree with brute force.
- The n0 split partitions the region exactly.
- Below 10^131, exhaustive enumeration gives exactly the 1,117 known chain numbers.
- Production-scale slices (c = 35,500 and c = 34,200 at N = 36,000) were re-checked against exact
  arithmetic.
- Stage survival rates over the full runs match theory to 4 significant figures (0.06872 and 0.05497).
- The production binary rebuilds byte-identically from `gsieve.c`.

### 5.4 Production
- N = 36,000: 4.671·10^13 candidates, about 3.9 core-hours.
- Shell 36,000–40,000: 1.736·10^13 candidates, about 1.6 core-hours (short rows make it slower per point).
- No new chain numbers. **Result: no number below 10^40000 has persistence ≥ 12.** This beats
  2.67·10^30000 (Peters 2023). Earlier bounds: Chaffin 2016, zero-free products up to 10^20000;
  Radcliffe 2019, 10^20585.
- Throughput: 0.286 ns per candidate per core (3.5·10^9 candidates/s/core) on the main path.

### 5.5 Structure above the bound (bridge lemma)
Any persistence-12 orbit contains either nine chain numbers with m > 40,000, or a bridge. A bridge is a
chain number with more than 27,950 digits whose digit product is one of the 1,117 known chain numbers
(< 2^163), so over 99.4% of its digits are 1s. Bridges lead to S-unit equations with up to 163 terms,
far beyond Baker-method reach.

## 6. Heuristics (random-digit model)

- **Calibration** against the complete list, predicted / actual:

  | Digits | Zero-free numbers | Chain numbers |
  |---|---|---|
  | 20–39 | 5,622 / 5,581 | 189.3 / 186 |
  | 40–59 | 2,059 / 2,111 | 7.3 / 7 |
  | 120–139 | 3.3 / 4 | – |

- **Persistence 12 exists** if and only if some 7-smooth number has digit product 4996238671872 or
  937638166841712. The expected number of such numbers under the model is 0.0186, about a 2% chance
  a priori.
  - That mass is concentrated at 18–22 digits (peak at 20; 99.99% within 30). So a(12) "should" have had
    about 21–25 digits, and that range is searched and empty.
  - Beyond 30 digits the expected count is about 10^−6; beyond 100 digits, about 10^−58.
- **Expected new chain numbers** beyond 49 digits: about 1.2. Beyond 25,000 digits: 10^−2428.
- **Conclusion:** larger bounds certify an exponentially thin tail; they do not change the odds.

## 7. Literature (`docs/LITERATURE.md`)

- **Small bases:** the maxima were numerical observations (Perez–Styer 2015; de Faria–Tresser 2014;
  Lamont-Smith 2021; OEIS A064868). No proofs were found.
- **Base 10:** BCGN 2021 settles odd targets. Fonga 2026 gives a 2-adic bound making even targets finite,
  but the algorithm is not executed there. Mondot's OEIS sequences A350180–A350187 are "postulated
  complete."
- **a(12) bounds:** Peters 2023 > 2.67·10^30000 (enumerated about 3.6·10^13 candidates, checking the last
  306 digits for zeros); Chaffin 2016; Radcliffe 2019.
- **Not new:** the cumulative depth inequality (Bennett–Bugeaud–Mignotte style).
- **New as far as found:** the chain-number reduction, suffix collapse, exact 4^k branching, odd thread,
  bridge lemma and calibrated heuristic.

## 8. Scaling and next-step estimates

- **Chain sieve cost** grows like N³: about 3.9 core-hours × (N/36,000)³.
  - N = 10^5 takes about 84 core-hours (1.3 h on 64 cores). N = 10^6 takes about 84,000 core-hours.
  - The floor is the bitmap pass (≈0.03 ns per candidate). Sub-cubic would need a fundamentally new idea.
- **RTX 3060 Ti** (8 GB, 448 GB/s) estimate: 5·10^10–2·10^11 candidates/s once ported to CUDA. That would
  put 10^5 digits at about 2–10 hours. This is unbenchmarked.
- **Memory** is constant: about 1.1 GB per process for the K = 13 tables.

## 9. Dead ends and limits (do not repeat blindly)

- **Finiteness of E:** no approach found. Digit-lifting trees branch about 4.5× per digit (supercritical);
  only rigid branches (the odd thread) are decidable.
- **Bridges via Baker:** too many terms (k ≤ 163; even the D = 9 bridges have 15–29 non-1 digits).
  Congruence sieves on mostly-ones numbers do not bite.
- **Parity forcing:** the bound on v₂ of bridges is only exponential in k when made rigorous.
- **Sub-cubic sieve:** lattice directions, φ-binning and multi-level bitmaps give constant factors only.
- **Zero-free-mode sieve (extending Chaffin's 10^20000):** about 4× costlier per candidate than chain
  mode. Not run.
- **Node 1161216 and k ≥ 13 mixed nodes:** cost explodes combinatorially with gap patterns.

## 10. Engineering pitfalls hit (and fixed)

- Cloud background jobs die at 2 hours, so every long job needs checkpoint/resume. Do not grep appended
  logs for stale completion lines.
- `pkill -f <pattern>` killed the calling shell. Use precise PIDs.
- C division truncates toward 0. Guard negative numerators in region bounds, otherwise extra
  out-of-region points are counted.
- Fixed-width digit windows count leading zeros as digits. Only apply them to numbers longer than the
  window (the chain sieve requires x ≥ 10^130 for the 126-digit engine).
- `exp10l` needs `_GNU_SOURCE`. An implicit declaration silently returns int.
- Lattice precision must scale with the coefficient size (an assertion fired at k = 13).
  Python bottom-test precision was audited against Mmax.
- Write JSON atomically (temp file plus rename).
- Avoid permutations of large bags: use a multiset generator. Group small bags to avoid per-call overhead.
- Large random-access tables suffer TLB/L2 misses. Use huge pages (madvise) and keep hot tables small.
