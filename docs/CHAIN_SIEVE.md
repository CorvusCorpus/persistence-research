# The chain sieve: a new lower bound for multiplicative persistence 12

## Statement

Let p(n) be the product of the decimal digits of n. Call a 7-smooth number x ≥ 10 a **chain number**
if x has no digit 0 and does not contain both a digit 5 and an even digit. Write G for the set of
chain numbers. For x = 2^a 3^b 5^d 7^c put m(x) = a/3 + b/2 + c + d ("digit length").

**Theorem (computational).** Every chain number x with m(x) ≤ 40,000 is one of 1,117 known numbers,
the largest being 4613824623662469877463718941499166849983276122112 (49 digits).

**Corollary.** No number below 10^40000 has multiplicative persistence ≥ 12 (any final digit).

Proof of the corollary. Let n have persistence ≥ 12 and c_i = p^(i+1)(n). Since p(x) < x for x ≥ 10,
c_0 > c_1 > c_2 ≥ p^11(n) ≥ 10. c_1 and c_0 contain no 0 (their digit products are nonzero), so c_1 is a
zero-free number ≥ 10, hence not divisible by 10; therefore c_0 cannot contain both a 5 and an even
digit, i.e. c_0 ∈ G. Every digit d of n contributes at most 1 to m(p(n)) (2: 1/3, 3: 1/2, 4: 2/3, 5: 1,
6: 5/6, 7, 8, 9: 1), so m(c_0) ≤ len(n). If len(n) ≤ 40,000, c_0 is one of the 1,117 known chain numbers;
all of them have persistence ≤ 10, so n has persistence ≤ 11. ∎

Previous record: Tim Peters (OEIS A003001, Sep 2023): if a(12) exists it is > 2.67·10^30000.

## Why chain numbers

Every intermediate number c_0, …, c_8 of a persistence-12 orbit is a chain number. G is much thinner
than the zero-free 7-smooth numbers E: E has 12,606 known members up to 140 digits, G has 1,117 up to
49 digits. Per digit a random number survives the chain condition with probability 0.8 (vs 0.9).

Forms: x ∈ G is never divisible by 10, so x = 2^a 3^b 7^c (form i) or x = 3^a 5^b 7^c with b ≥ 1
(form iii). In form iii the last digit is 5, so every digit must be odd.

## Lemmas used by the sieve

**L1 (suffix collapse).** For a ≥ K, 2^a 3^b 7^c ≡ 2^e (mod 10^K) where e ≥ K and
e ≡ a + b·μ3 + c·μ7 (mod 4·5^(K−1)), μp = log₂ p in (Z/5^K)^×.
(x ≡ 0 mod 2^K and x ≡ 2^e mod 5^K since 2 is a primitive root mod 5^K.)
So one bitmap indexed by e mod 4·5^12 (976,562,500 bits) answers "are the last 13 digits allowed" for
every 7-smooth number with a ≥ 13, and along a row (b, c fixed) the index runs through consecutive
residues, 64 exponents per machine word.

**L2 (exact branching).** For even numbers with the chain condition (digits ∉ {0,5}) the set S_K of
admissible residues has exactly 4^K elements. Lifting from K to K+1 digits, the 5 lifts realise the 5
digits of one fixed parity; exactly one of them is 0 or 5. (Build check: 4^13 = 67,108,864 ✓.)

**L3 (affine lift).** Let e = e13 + P13·j (P13 = 4·5^12, j < 5^13). With y0 = 2^(e13−26) mod 5^26 =
A_lo + 5^13·A_hi, w = (2^P13 mod 5^26 − 1)/5^13:
digits 14–26 of x equal h0 + 2^13·((A_hi + A_lo·w·j) mod 5^13), h0 = ⌊2^13·A_lo/5^13⌋.
(Because (1 + 5^13 w)^j ≡ 1 + 5^13 w j mod 5^26.) j is constant along a row, so the next 13 digits cost
one modular multiply; implemented with AVX-512 (Montgomery multiply, IFMA divisions).

**L4 (odd-branch rigidity).** For x = 3^a 5^b 7^c, x mod 2^K = (−1)^(a+c)·5^(b + a·f3 + c·f7).
If b ≥ k ≥ 3, the next digit has two candidates differing by 5 — one odd, one even — so exactly one
lift survives the all-odd condition. Exhaustive lifting from k = 3 shows there is exactly one 2-adic
thread (sign −1). Hence for b ≥ 40: all last 40 digits odd ⟺ a + c odd and
b ≡ t − a·f3 − c·f7 (mod 2^38). Every such number ends in
…9313975971735999337397515773773193359375. This branch costs O(1) per (a, c).

**L5 (leading digits).** t = a·L2 + b·L3 + c·L7 mod 2^64 with L_p = round(2^64 log10 p) approximates
2^64·frac(log10 x) within (a+b+c)/2 units. A 2^18-bin table marks bins whose every point (with margin
2^−40) has a 0 or 5 among the first 5 significant digits; only those are rejected.

## Pipeline (form i, a ≥ 26, x ≥ 10^130)

1. bitmap (last 13 digits) → 6.87% survive
2. affine lift, digits 14–26 → 0.38% survive
3. leading-digit table (5 digits) → ~40% of those survive
4. exact base-10^9 limb product, digits 27–126
5. GMP: last 600 digits, then the full number

Form i with a < 26: three 18-digit limbs stepped along b. Form iii: b ≤ 39 stepped along a with an
18-digit filter; b ≥ 40 via L4. Numbers below 10^131: exhaustive (4,483,467 numbers).

## Validation

- Stage decisions (bitmap+lift, leading-digit rejections, limb engine) compared one by one with exact
  big-integer arithmetic on full regions (N = 150; G and zero-free modes; K = 12 and 13 builds).
- Zero-free mode end to end (main + direct + odd paths, N = 160, numbers ≥ 10^67): reproduces exactly
  the 314 zero-free 7-smooth numbers of that size from the independent list (both builds).
- Direct and odd paths: survivor sets equal to brute force (N = 160).
- Region split by digit length (n0) partitions points and members exactly.
- Small range (< 10^131): exhaustive enumeration gives exactly the 1,117 known chain numbers.
- Production scale: complete slices c = 35,500 and c = 34,200 of the N = 36,000 run re-checked against
  exact arithmetic (val_prodrow.py): stage-2 survivor sets identical, all limb-engine decisions and a
  sample of leading-digit rejections verified. Odd-thread candidates verified to end in the 40-digit thread.
- Survival rates over the full runs match theory to 4 significant figures (0.06872, 0.05497).

## Production runs

Base run N ≤ 36,000, then the shell 36,000 < m ≤ 40,000 (flag n0=36000). Both: every c covered for
every path (aggregate.py), no chain numbers outside the known list.

| run | main | direct | odd | total |
|---|---|---|---|---|
| m ≤ 36,000 | 46,626,195,781,376 | 33,689,253,108 | 50,490,005,607 | 46,710,375,040,091 |
| 36,000 < m ≤ 40,000 | 17,337,008,704,000 | 7,903,308,000 | 11,850,072,000 | 17,356,762,084,000 |
| **m ≤ 40,000** | | | | **64,067,137,124,091** (+ 4,483,467 below 10^131) |

### Base run detail (N = 36,000)

| path | region | candidates | exceptions |
|---|---|---|---|
| main (2^a3^b7^c, a ≥ 26, ≥ 130 digits) | 2a+3b+6c ≤ 216000 | 46,626,195,781,376 | 0 |
| direct (2^a3^b7^c, a < 26, ≥ 130 digits) | same | 33,689,253,108 | 0 |
| odd (3^a5^b7^c, ≥ 130 digits) | a+2b+2c ≤ 72000 | 50,490,005,607 | 0 |
| tiny (all forms, < 10^131) | all | 4,483,467 | 1,117 known |
| **total** | | **46,710,379,523,558** | |

Stage survival on the main path: 0.06872 after the bitmap (theory 0.8^12 = 0.068719),
0.05497 after the affine lift (theory 0.8^13 = 0.054976). Every c in [0, 36000] is covered by
checkpointed logs for each path (aggregate.py). Wall time ≈ 2 h on 2 cores (Xeon, 2.1 GHz, AVX-512).

## Structure of any counterexample above the bound

Let n have persistence ≥ 12 and c_i = p^(i+1)(n). Then c_0, …, c_8 are chain numbers and
m(c_0) > 40,000. Let t be the largest index ≤ 8 with m(c_t) > 40,000. Either t = 8 (nine chain numbers
beyond the bound), or c_t is a **bridge**: a chain number with more than 27,950 digits (m ≤ 1.431·log10 x)
whose digit product is one of the 1,117 known chain numbers (< 2^163), so fewer than 163 of its digits
differ from 1 — over 99.4% ones.

## Heuristic (random-digit model)

Model: digits of 7-smooth numbers behave like independent uniform digits. Calibration against the
complete list below 10^5000: zero-free numbers with 20–39 digits predicted 5,622 / actual 5,581;
40–59: 2,059 / 2,111; 120–139: 3.3 / 4; chain numbers with 40–59 digits: 7.3 / 7.

A persistence-12 number exists iff some 7-smooth number has digit product 4996238671872 or
937638166841712. Expected number of such 7-smooth numbers: 0.019, peaking at 20 digits
(86% within 22 digits, 99.99% within 30). Expected number of new chain numbers beyond 25,000
digits: 10^-2428.
