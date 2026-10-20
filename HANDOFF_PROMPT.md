You are taking over a research codebase on multiplicative persistence: the digit-product map n → p(n),
and how many steps it takes to reach a single digit. The previous session produced validated results and a
working toolchain. Your job is to turn this directory into a proper repository and continue the research
rigorously.

Read first, in this order: README.md, CLAUDE.md, docs/SESSION_SUMMARY.md, docs/RESULTS.md, docs/ROADMAP.md,
docs/CHAIN_SIEVE.md. The rules in CLAUDE.md are binding.

Already established (re-run only to validate, not to re-derive):
- No number below 10^40000 has persistence ≥ 12. This comes from the "chain sieve": 6.4·10^13 candidates,
  logs in work/gsieve/logs/. The previous record was 2.67·10^30000 (Peters 2023).
- Bases 3, 4 and 6: complete theorems for every nonzero final digit.
- Base 10, even final digits: 92 tree nodes proved. The open nodes are listed in docs/RESULTS.md, and the
  paused runs have checkpoints in work/ckptC and work/bagresC.
- Heuristic: a persistence-12 number had about a 2% chance a priori, concentrated around 20 digits, a range
  that is already excluded.

Where you are running matters:
- In a cloud container with no GPU (e.g. a Claude Code web session): do task 1, plus the CPU-only parts of
  task 2 if the CPU has AVX-512. Then stop and summarise, because tasks 3+ need the desktop.
- On the desktop (RTX 3060 Ti): continue with tasks 2–5.

Tasks, in order. Report after each one before moving on.

1. Set up and reproduce.
   - The repository already contains the baseline commit (the bundle as produced). Keep .gitignore; the
     production logs are committed as plain files (all small). Add a top-level Makefile or test script.
   - Check the environment: gcc, GMP, MPFR, Python with python-flint and mpmath. CPU flags:
     `grep -o 'avx512[a-z_]*\|bmi2' /proc/cpuinfo | sort -u`. GPU: `nvidia-smi`, `nvcc --version`.
   - Build with `cd work && make`. If AVX-512 is available, run `make -C work/gsieve test`. Run the
     base-10 smoke test: in a fresh copy with empty bagres*/ and ckpt*/ directories, move proofC_98.json
     and proofC_147.json aside and run `python3 prove_cum.py 98 147`; the output must match
     docs/results_base10.md.
   - Commit the result as the baseline, and start docs/NOTEBOOK.md (dated entries: command, binary/flags,
     commit, timing, outcome).
2. Make the chain sieve portable and faster on CPU.
   - If the CPU lacks AVX-512, add an AVX2 path for the main-path stage 2 (affine lift and digit checks)
     behind a compile flag. Prove it makes identical decisions using work/gsieve/val_dump.py and
     val_prodrow.py, extended as needed.
   - Add pthreads with shared tables. Benchmark candidates per second.
3. CUDA port for an RTX 3060 Ti (8 GB).
   - Port stage 1 (bitmap scan of the last 13 digits) and stage 2 (affine lift, digits 14–26): one row
     (b, c) segment per block, tables (≈1 GB) in device memory. Stream stage-2 survivor positions back to
     the host and keep stages 3+ (leading-digit table, limb engine, GMP) on the CPU, unchanged.
   - Validate exact equality of the stage-2 survivor sets against the CPU build on several c-slices at
     N = 36000 and on E-mode at N = 160 (the E-mode end-to-end test must still find exactly the 314 known
     numbers).
   - Then benchmark. The pre-port estimate is 5·10^10–2·10^11 candidates/s; report the actual figure.
4. Extend the record in shells with the n0 option: 40000 → 50000 → … → 100000.
   - Use checkpointed long runs, run `aggregate.py N N0` after each shell, and archive the logs.
   - Claim a new bound only when coverage is complete, there are no members, the stage survival rates
     match theory (0.8^12 and 0.8^13 on the main path), and the validation suite passes on the exact
     binaries used.
5. Write-up.
   - A short LaTeX paper covering: the chain lemma, suffix collapse, exact 4^k branching, the affine lift,
     the odd 2-adic thread, the bridge lemma, the computation and validation, and the calibrated heuristic.
   - A draft OEIS A003001 comment.
   - A slow, simple independent reference implementation of stages 1–2 for audits.
6. Optional, if CPU time is spare: resume the base-10 even-target nodes (`python3 runner4.py 122472`,
   `python3 runner4.py 131712 24192` from work/).

Constraints:
- Be rigorous and conservative. A filter may only discard a candidate when it is certain.
- Never hand-edit proof or log artefacts.
- Keep heuristics separate from proofs, and don't overclaim novelty: re-check the literature (OEIS
  A003001 comments, arXiv) before calling anything new.
- Ask me before deleting anything or starting runs longer than about 12 hours.
