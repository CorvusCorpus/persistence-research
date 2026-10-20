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
