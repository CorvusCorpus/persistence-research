# Draft comment for OEIS A003001 (not submitted)

Submitting is the author's decision, and it is done by the author on oeis.org (an OEIS account is needed). The current
newest bound in the entry is: "New lower bound: if a(12) exists, it must be > 2.67*10^30000." (Tim Peters, Sep 19 2023).
Literature and OEIS were re-checked on 20 Oct 2026 (oeis.org fetched directly: last modified Oct 20 2026, no bound newer
than Peters's); see `docs/LITERATURE.md`.

Only the bound proved and validated in `docs/RESULTS.md` may go in. Today that is 10^40000; replace the numbers if the
extension to 10^100000 is completed and recorded as COMPUTATIONAL.

## COMMENTS entry

```
If a(12) exists, it must be > 2.67*10^40000. If n has persistence >= 12, then its digit product c = p(n) has no digit 0
and does not contain both a 5 and an even digit (otherwise p(c) would end in 0), and c = 2^i*3^j*5^k*7^l satisfies
i/3 + j/2 + k + l <= (number of digits of n), since each digit of n contributes at most 1 to that sum. An exhaustive
search of all such c with i/3 + j/2 + k + l <= 40000 (6.4*10^13 candidates, using a bitmap of the last 13 digits of
powers of 2 and an exact lift of the next 13 digits) found only the 1117 such numbers below 10^49, none of which has
persistence 11. So every number with at most 40000 digits has persistence <= 11, and the smallest candidate with 40001
digits is 2677...7 = 2.677...*10^40000. - [Name], Oct 20 2026
```

## LINKS entry (optional)

```
[Name], <a href="https://github.com/CorvusCorpus/persistence-research">Chain sieve for multiplicative persistence</a>,
code, logs and validation (2026).
```

## Notes for the submitter
- The persistence-11 test of the 1117 numbers and the candidate argument (`docs/LITERATURE.md`, "OEIS-style statement")
  are short and can be checked by hand or by a few lines of code.
- The OEIS asks for comments that are verifiable. The repository has the logs (`work/gsieve/logs/`), the coverage check
  (`make records-check`) and the validation records (`docs/NOTEBOOK.md`, `validation_pc/`).
- Do not mention unvalidated extensions.
