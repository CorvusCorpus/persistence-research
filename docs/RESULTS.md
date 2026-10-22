# Results ledger

Status legend: **THEOREM** (complete proof, computer-assisted where stated) · **COMPUTATIONAL** (exhaustive
computation with validation; no hand proof) · **HEURISTIC** (model-based, not a proof) · **OPEN**.
All computer-assisted claims are unreviewed; see "Trust base" at the end.

Notation: p(n) = product of the decimal digits of n; persistence = number of applications of p needed to reach
one digit. 7-smooth = all prime factors ≤ 7. E = zero-free 7-smooth numbers ≥ 10. For x = 2^a 3^b 5^d 7^c,
digit length m(x) = a/3 + b/2 + c + d. Chain number = 7-smooth x ≥ 10 with no digit 0 that does not contain
both a 5 and an even digit (set G).

## 1. Final digit 0 / persistence 12 (this session's headline)

| # | Claim | Status | Evidence |
|---|---|---|---|
| 1.1 | Every chain number x with m(x) ≤ 100,000 is one of 1,117 known numbers (largest: 4613824623662469877463718941499166849983276122112, 49 digits). | COMPUTATIONAL | `work/gsieve/`: 1.0004·10^15 candidates. **m ≤ 40,000:** 6.4067·10^13 candidates (cloud, ECC memory), logs in `work/gsieve/logs/`; re-executed in full twice (session 3) with the independent GPU implementation of stages 1–2 on another machine, every counter and every per-c count equal (`validation_pc/task3_cuda/`). **40,000 < m ≤ 100,000:** six shells, 9.3635·10^14 candidates (20 Oct 2026, session 3, CUDA variant with fail-stop table checks, binaries validated by `val_cuda.sh` and `val_all.sh`), logs in `work/gsieve/logs/prod_N*_n0*.tar.xz`. `scripts/check_coverage.py` shows full coverage and no members for every run, and the survival rates match 0.8^12 and 0.8^13. **Second executions:** the whole main path for m ≤ 100,000 (1.0·10^15 candidates) by the GPU tree sieve, an independent stage-1–2 algorithm (coset enumeration) with its own stages 3–4, every counter equal to the production totals (22 Oct 2026, NOTEBOOK T.7, `validation_pc/tsieve/v2/`); all paths of the shells 40,000–60,000 by the CPU variant on cloud hardware, every counter and per-c count equal (NOTEBOOK S.4); the direct and odd paths of the shells 60,000–100,000 by the CPU binary on the desktop, every counter and per-c count equal (22 Oct 2026, NOTEBOOK E.8, `validation_pc/task4_check/check_do_N70000-100000.tar.xz`). Every path for m ≤ 100,000 has now been executed twice |
| 1.2 | No number below 10^100000 has multiplicative persistence ≥ 12 (any final digit), so a(12) > 2.67·10^100000. Previous record: a(12) > 2.67·10^30000 (T. Peters, OEIS A003001, 2023). | THEOREM given 1.1 | proof in `docs/CHAIN_SIEVE.md` (corollary) |
| 1.3 | The 1,117 chain numbers have persistence ≤ 10; the only ones with persistence 10 are 4996238671872 and 937638166841712. | COMPUTATIONAL | `work/gsieve/G_small.json` + check in summary |
| 1.4 | Exact branching: for even 7-smooth x divisible by 2^k, the number of admissible k-digit endings (no 0, no 5) is exactly 4^k (out of 4·5^(k−1) residues). | THEOREM (short proof) + checked k ≤ 8 and k = 13 | `docs/CHAIN_SIEVE.md` L2 |
| 1.5 | Suffix collapse: for a ≥ K, 2^a 3^b 7^c ≡ 2^e (mod 10^K) with e ≡ a + b·log₂3 + c·log₂7 (mod 4·5^(K−1)). | THEOREM | L1 |
| 1.6 | Affine lift formula for digits K+1..2K. | THEOREM + randomized checks | L3, `gconst.py` |
| 1.7 | Odd branch: a 7-smooth multiple of 5^b (b ≥ 3) with only odd digits has its last b digits equal to those of one fixed 10-adic number …9313975971735999337397515773773193359375 (unique 2-adic thread, sign −1). | THEOREM (exhaustive lifting) | L4, `gconst.py` |
| 1.8 | Bridge structure: any persistence-12 orbit contains either nine chain numbers with m > 100,000, or a "bridge": a chain number with > 69,897 digits, > 99.77% of them 1s, whose digit product is a known chain number. | THEOREM given 1.1 | `docs/CHAIN_SIEVE.md` |
| 1.9 | Random-digit model: matches the complete zero-free list well (e.g. 20–39 digits: 5,622 predicted / 5,581 actual). Expected number of 7-smooth numbers with a persistence-10 digit product: 0.019 (≈2% a priori), peaking at 20 digits; beyond 25,000 digits ≤ 10^−2428. | HEURISTIC | `docs/SESSION_SUMMARY.md` §6 |
| 1.10 | Finiteness of E (equivalently: does the final-digit-0 tree terminate?) | OPEN (Erdős-type; same family as "2^86 is the last zero-free power of 2") | — |

## 2. Base 10, nonzero even final digits (2, 4, 6, 8)

| # | Claim | Status | Evidence |
|---|---|---|---|
| 2.1 | For 99 nodes of the even-final-digit trees, the complete set of 7-smooth preimages (all sizes) is determined. | THEOREM (computer-assisted: Matveev + lattice/p-adic reduction + exhaustive list) | `work/certificate.json`, `work/proofC_*.json`, `work/proofJ_*.json`, table `docs/results_base10.md`. Seven nodes (131712, 122472, 24192, 1229312, 314928, 41472, 93312; 22 Oct 2026) were proved with the GPU case loops (`work/prove_cum_gpu.py`, validated against the CPU programs, `docs/GPU_EVEN_TARGETS.md`); one bag of each of six of them was re-proved from scratch on the CPU path with identical facts, bounds and counters (NOTEBOOK E.6–E.8) |
| 2.2 | Tree levels determined exactly: target 2 to depth 3; target 4 to depth 3; target 6 to depth 3; target 8 to depth 2. | THEOREM (from 2.1) | `docs/results_base10.md` |
| 2.3 | Open (blocking) nodes, one per target: target 2: 11239424 (1.3·10^14 work units, plus two bags with k ≥ 17 whose a-priori bound overflows the current code); target 4: 1161216 (1.4·10^14); target 6: 221184 (2.2·10^13); target 8: 33718272 (1.6·10^15, plus six bags with k ≥ 17). | OPEN | costs `work/logs/node_costs.out`; checkpoints `work/ckptC/`, `work/bagresC/` |
| 2.4 | Odd final digits in base 10 are settled in the literature (Brier–Clavier–Gutsche–Naccache 2021). | KNOWN | `docs/LITERATURE.md` |

## 3. Small bases (nonzero final digits)

| # | Claim | Status | Evidence |
|---|---|---|---|
| 3.1 | Base 3: orbits ending at 1 have persistence ≤ 3; ending at 2: ≤ 2 (attained). | THEOREM (computer-assisted) | `work/base3_cum.json`, `docs/results_smallbases.md` |
| 3.2 | Base 4: ending at 2: ≤ 3; ending at 1 or 3: ≤ 1 (attained). | THEOREM | `work/base4_cum.json` |
| 3.3 | Base 6: ending at 2 or 4: ≤ 5; ending at 3: ≤ 2; ending at 1 or 5: ≤ 1 (attained). | THEOREM | `work/base6_cum.json`, `work/trees_b6_D800.json` |
| 3.4 | These match the numerically observed maxima (Perez–Styer 2015; Lamont-Smith 2021) for all orbits that end at a nonzero digit; no prior proof found. | literature check | `docs/LITERATURE.md` |
| 3.5 | The method works only in bases 3, 4, 6, 10 (p-adic rank vs number of free primes). | observation | `docs/METHOD_EVEN_TARGETS.md` §4 |

## 4. Data

| Item | Content | Status |
|---|---|---|
| `work/Z.json` | all zero-free 7-smooth numbers below 10^5000 (12,615 incl. 1–9; largest 2^25·3^227·7^28, 140 digits) | COMPUTATIONAL (zf7.c sieve + verify.py) |
| `work/trees_b6_D800.json` | base 6 analogue below 6^800 | COMPUTATIONAL |
| `work/gsieve/G_small.json` | the 1,117 chain numbers (all < 10^131 exhaustively) | COMPUTATIONAL |

## Trust base
Matveev's lower bound for linear forms in logarithms (base-10 and small-base trees); correctness of the
exhaustive lists; correctness of the C/Python code (validated by planted tests, brute-force comparisons and
independent re-derivations as documented); GMP/MPFR. Nothing here has been peer reviewed. Independent
re-running (ideally re-implementation) is the main missing step before publication.
