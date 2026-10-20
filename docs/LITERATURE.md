# Literature check (20 Oct 2026): are the small-base results known?

Searched: persistence papers (Sloane 1973 as cited; Perez-Styer 2015; de Faria-Tresser 2014; Lamont-Smith 2021;
Brier-Clavier-Gutsche-Naccache 2021; Fonga 2026), OEIS A003001, A064867-A064872, A245760, A350180-A350187, web/MathOverflow.

Findings
- Perez & Styer, Involve 8 (2015): numerical evidence only; maximal persistence 3, 3, 5 in bases 3, 4, 6.
- de Faria & Tresser, arXiv 1307.1188: base 3 bound nu_3(n) <= 2(1 + log_3 log_3 n); bases 3 and 4 reduced to
  conjectures that 2^k (base 3), and 3^k and 2*3^k (base 4), eventually always contain a 0 digit; no base-6 result;
  no treatment of orbits ending at a nonzero digit.
- Lamont-Smith, JIS 24 (2021): proves only two-digit results; all maximal-persistence values per base are
  computational/conjectural (Table 5: base 3: 3, base 4: 3, base 6: 5).
- OEIS A064868: "a(3) and a(4) do not seem to exist" (no persistence-4 number in bases 3, 4): observation, no proof.
- BCGN 2021: base 10 only (odd final digits).
- Fonga 2026 (arXiv 2608.27802): 2-adic valuation bound for fixed digit multiplicities, valid in every base >= 3
  and every prime p | b; says it yields a finite algorithm in base 10 when combined with BCGN; does not run it.
- OEIS A350180-A350187 (Mondot): 7-smooth digit products of each persistence, "postulated finite and complete".

Conclusion
- No source found proves the nonzero-final-digit trees (or persistence bounds) in bases 3, 4 or 6. The values we
  prove coincide with the numerically observed maxima, so those maxima are now proved for every orbit that ends at
  a nonzero digit; orbits ending at 0 remain open (in bases 3, 4 they would follow from de Faria-Tresser's digit conjectures).
- The base-10 even-digit trees are an execution of the algorithm Fonga's bound makes possible (with our cumulative
  refinement and C implementation); our complete node sets settle the nonzero-final-digit parts of Mondot's
  postulated-complete sequences up to the depths reached.
- The cumulative depth inequality 0 < Lam_j < 80*10^(-(L - a_(j+1))) is the natural form used in the Baker-method
  literature on integers with few nonzero digits (e.g. Bennett-Bugeaud-Mignotte); it improved on our own earlier
  weaker form, and should not be claimed as new.
- Caveat: not exhaustive (theses, unindexed preprints and Kurz's original web pages were not checked).

## Re-check of the persistence-12 bound (20 Oct 2026, session 2)
Source: OEIS A003001 as of 20 Oct 2026, 14:57 EDT. The page text was supplied by the user, because oeis.org is
blocked from the cloud container.

- **Latest bound.** Tim Peters, 19 Sep 2023: if a(12) exists, it is > 2.67·10^30000. No later bound is recorded.
  - His candidates are the digit strings matching eight regular expressions: the smallest integer of each
    digit-product class.
  - Every candidate with at most 30000 digits (a bit over 36·10^12 of them) has a digit product with a zero
    among its last 306 digits.
  - That method differs from the chain sieve: it enumerates digit strings, while the sieve enumerates digit
    products.
- **Chaffin (Sep 2016) agrees with our data.**
  - The largest digit product with persistence 10 is 2^4·3^20·7^5; the only other known one is 2^19·3^4·7^6.
  - The largest with persistence 2 is 2^25·3^227·7^28, with 140 digits.
  - Every p(n) in (10^140, 10^20000) contains a 0.
  - Radcliffe (Mar 2019): no terms up to 10^20585.
- **OEIS-style statement of our result.** Every n with at most 40,000 digits has persistence ≤ 11. If a(12)
  exists it is a candidate, so it is at least the smallest 40,001-digit candidate: 26 followed by 39,999 sevens.
  Hence a(12) > 2.67·10^40000.
- **Other sources.** A web search finds the same 2.67·10^30000 (30,000 digits) bound cited by Wikipedia and by
  J. D. Cook's blog (Jun 2025). The search is not exhaustive: unindexed or private work cannot be excluded.

## Re-check before the write-up (20 Oct 2026, session 3, desktop)
- **OEIS A003001**, fetched directly from oeis.org (the page was last modified Oct 20 2026, 17:14 EDT).
  - The newest bound is still Peters's: a(12) > 2.67·10^30000 (Sep 19 2023).
  - No comment dated 2024–2026 is about bounds, and no comment mentions chain numbers.
- **Web and arXiv searches** for bounds on a(12), zero-free 7-smooth searches, and digit lifting for powers of 2
  found nothing that supersedes the above.
  - The google-deepmind/formal-conjectures issue #7184 (Oct 2026) also lists persistence 12 as open.
- **Prior work that the write-up must credit:**
  - K. McElwee, "An algorithm for multiplicative persistence research", blog post, 13 Jul 2019.
    - It enumerates zero-free products separately in the families {2,3,7} and {3,5,7}, up to 10^525, and finds
      2^25·3^227·7^28 as the largest.
    - It remarks that a number containing both 2s and 5s maps to a multiple of 10.
    - So the form split, and the elementary 5-with-even observation, are not new.
  - B. C. Dumitru, arXiv:2503.23177 (2025), on powers of two with all even digits. This and the surrounding
    discussions use residue lifting of 2^n mod 10^k over exponents mod 4·5^(k−1).
    - Lemmas L1 and L2 are in that standard style. Only the exact count 4^k for the excluded set {0,5}, and its use
      as a bitmap sieve, are claimed here, with "not found stated elsewhere" as the wording.
  - P. N. Fonga, arXiv:2608.27802 (28 Aug 2026), as before: even nonzero targets, with no bound on a(12).
- **Novelty, as worded in `paper/chain_sieve.tex`:**
  - the bound itself;
  - the chain-number sieve weighted by digit length m;
  - the affine lift (L3) and the odd thread (L4) as sieve tools;
  - the bridge lemma.

  Everything else is credited.
