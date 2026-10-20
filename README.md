# Multiplicative persistence: computer-assisted results

Baseline snapshot of a research session (19–20 Oct 2026). Start with `docs/SESSION_SUMMARY.md`. Then read
`docs/RESULTS.md`, which lists every claim with its status and evidence, and `docs/ROADMAP.md`.

## Headline results
1. **No number below 10^40000 has multiplicative persistence ≥ 12** (any final digit). The previous record
   was 2.67·10^30000 (T. Peters, OEIS A003001, 2023). This rests on the "chain sieve": every 7-smooth
   number of digit length ≤ 40,000 contains a 0, or contains a 5 together with an even digit, or is one of
   1,117 known numbers. 6.4·10^13 candidates were checked. See `docs/CHAIN_SIEVE.md`.
2. **Bases 3, 4, 6:** complete theorems for every nonzero final digit. Maximum persistence for orbits
   ending at each digit:
   - base 3: 3 (digit 1), 2 (digit 2);
   - base 4: 3 (digit 2), 1 (digits 1 and 3);
   - base 6: 5 (digits 2 and 4), 2 (digit 3), 1 (digits 1 and 5).

   These were previously only numerical observations.
3. **Base 10, even final digits:** complete 7-smooth preimage sets for 92 tree nodes, which determines the
   trees of final digits 2, 4 and 6 to depth 3 and of final digit 8 to depth 2.
4. **Heuristic:** a random-digit model, checked against the complete data, gives a persistence-12 number
   about a 2% chance a priori. That chance sits at about 20 digits, a range already excluded. Beyond
   25,000 digits the chance is below 10^−2400.

## Layout
```
README.md            this file
CLAUDE.md            working rules for Claude Code sessions (read before changing anything)
HANDOFF_PROMPT.md    prompt used to start the follow-up session
docs/                summary, results ledger, roadmap, method notes, literature, result tables
work/                original working directory (flat layout; run scripts from here)
  *.py               provers, drivers, validation scripts (base 10 and small bases)
  csrc/*.c           C inner loops (MPFR/GMP); `make` builds them
  zf7.c, verify.py   exhaustive zero-free 7-smooth sieve -> Z.json
  Z.json             all zero-free 7-smooth numbers below 10^5000
  proofC_*.json, proofJ_*.json, proofI_*.json, certificate.json   base-10 node proofs
  base{3,4,6}_cum.json, trees_b6_D800.json                        small-base proofs and data
  ckpt*/, bagres*/   checkpoints of paused base-10 runs (122472, 131712, 24192)
  logs/              run logs, sieve outputs
  gsieve/            chain sieve: gsieve.c (AVX-512), constants, drivers, validation, logs/
```

## Quick start (Linux; gcc, GMP, MPFR, Python 3 with python-flint and mpmath)
```
cd work && make                    # builds csrc helpers, zf7 and gsieve
make -C gsieve test                # chain-sieve validation suite (needs AVX-512 CPU)
python3 gsieve/aggregate.py 36000 -1 gsieve/logs         # base run coverage
python3 gsieve/aggregate.py 40000 36000 gsieve/logs      # extension shell coverage
```
Re-deriving a base-10 node from scratch: run in a copy whose `bagres*/` and `ckpt*/` are empty, and move
that node's `proofC_*.json` aside, e.g. `python3 prove_cum.py 98 147`.

## Caveats
Computer-assisted and not peer reviewed. The trees rely on Matveev's theorem; everything relies on the
exhaustive lists and on the correctness of the code. Validation is extensive (see each method document),
but independent re-running and re-implementation are still to do.
