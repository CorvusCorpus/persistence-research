# CLAUDE.md: working rules for this repository

## What this project is
Computer-assisted research on multiplicative persistence. Current state is summarised in
`docs/SESSION_SUMMARY.md`; every claim and its status is listed in `docs/RESULTS.md`; next steps are in
`docs/ROADMAP.md`. Read those three before starting work.

## Non-negotiable rules
1. **No claim without validation.** A computational result counts only when:
   - every stage is validated against exact arithmetic on test regions;
   - a positive control passes, e.g. E-mode rediscovering the known zero-free numbers;
   - the production logs show full coverage (`work/gsieve/aggregate.py`).

   State each result as THEOREM, COMPUTATIONAL, HEURISTIC or OPEN, as in `docs/RESULTS.md`.
2. **Never weaken soundness for speed.** Filters may only reject a candidate when the rejection is
   certain. Floating point is used only with explicit error margins; ambiguous cases must pass through to
   exact checks.
3. **Do not edit proof artefacts by hand.**
   - `work/proofC_*.json`, `work/proofJ_*.json`, `work/certificate.json`, `work/base*_cum.json` and the
     production logs are records.
   - Regenerate tables with `work/results_table.py`, `work/certificate.py` and `work/smallbases_table.py`.
4. **Keep the heuristics separate from the proofs.** Do not overclaim novelty; check `docs/LITERATURE.md`
   and search again before calling something new.
5. **Reproducibility.**
   - Record the exact command, binary build flags and git commit for every production run.
   - Keep the per-c checkpoint logs.
   - Binaries must rebuild from source (the production gsieve binary was byte-identical to its source build).

## Environment
- Linux, gcc ≥ 11, GMP, MPFR (`libgmp-dev`, `libmpfr-dev`), Python ≥ 3.10 with `python-flint` and `mpmath`.
- The chain sieve (`work/gsieve/gsieve.c`) needs AVX-512 (F, DQ, BW, VL, IFMA, VBMI2) and BMI2. On other
  CPUs an AVX2 port is needed (ROADMAP P1); a CUDA port is ROADMAP P2.
- Build with `cd work && make`. Test with `make -C work/gsieve test`. Run scripts from `work/` (base-10
  and small-base code) or `work/gsieve/` (chain sieve), because they use relative paths.

## Key commands
- **Chain sieve production.** `cd work/gsieve && python3 prod.py N A|B BUDGET_SECONDS [N0]`. Worker A and
  worker B split the work; `N0` restricts to the shell N0 < m ≤ N.
  - Use `JOBS=mode:c0:c1,...` to override the job list (keep the original c-ranges when resuming).
  - Coverage check: `python3 aggregate.py N N0 [logdir]`.
- **Single run.** `./gsieve_k13 main|direct|odd G N c0 c1 [prog] [n0=M] [dump]`. G = chain mode,
  E = zero-free mode (validation), P = permissive.
- **Base-10 node.** `python3 prove_cum.py NODE ...`. Resumable via `ckptC/` and `bagresC/`.
  `runner4.py` runs node lists.
- **Small bases.** `python3 prove_cum_b.py BASE DEPTH [digits...]`.
- **Zero-free list.** `zf7 D K` (OpenMP) plus `verify.py`/`struct.py` produce `Z.json`.

## Pitfalls already hit (see SESSION_SUMMARY §10)
- **Signed division in C.** Integer division truncates toward 0; guard negative numerators in region bounds.
- **Leading zeros.** Fixed-width digit windows treat leading zeros as zero digits. Only apply a window to
  numbers longer than it (the sieve requires x ≥ 10^130).
- **`exp10l`** needs `#define _GNU_SOURCE`.
- **`pkill -f`** can kill your own shell. Use PIDs.
- **Long jobs.** Use checkpointed drivers. When waiting on appended logs, never match old completion lines.
- **Memory locality.** Large random-access tables cost TLB/L2 misses. Use huge pages and keep hot tables
  small; measure before optimising (the stage breakdown benchmarks used `-DBENCH*` flags).
