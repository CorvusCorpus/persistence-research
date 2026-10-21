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
- The chain sieve (`work/gsieve/gsieve.c`) builds on any x86-64. The build picks the best stage-2 variant
  the target supports: `ifma` (AVX-512 IFMA, the original production build), `avx512` (AVX-512 F/DQ/BW/VL),
  `avx2`, then `scalar`. All four make identical decisions (`docs/CHAIN_SIEVE.md`).
- The CUDA variant (stages 1–2 of the main path on the GPU) builds with `make -C work/gsieve cuda`
  (nvcc, CUDA 13). It makes the same decisions as the CPU variants and is validated by `val_cuda.sh`.
- Desktop PC (Windows 11, RTX 3060 Ti): use the WSL2 distribution `persistence` (`scripts/setup_wsl_cuda.sh`).
  - Run commands as `wsl -d persistence -u root --cd <repo> --exec bash -lc '<command>'`.
  - The login shell puts the venv and CUDA on `PATH`. `--exec` stops an outer shell from expanding `$VARIABLES`.
- Build with `make` (repository root). Test with `make test`, or `make -C work/gsieve test [V=variant]`.
  Run scripts from `work/` (base-10 and small-base code) or `work/gsieve/` (chain sieve), because they use
  relative paths. Log every run in `docs/NOTEBOOK.md`.

## Key commands
- **Chain sieve production.** `cd work/gsieve && python3 prod.py N A|B|ALL BUDGET_SECONDS [N0]`.
  - Worker A and worker B split the work between two processes; `ALL` runs every path in one process with
    `THREADS` threads (default: all CPUs).
  - `N0` restricts to the shell N0 < m ≤ N.
  - `EXE` selects the binary.
  - Use `JOBS=mode:c0:c1,...` to override the job list (keep the original c-ranges when resuming).
  - Coverage check: `python3 aggregate.py N N0 [logdir]`.
- **Single run.** `./gsieve_k13 main|direct|odd G N c0 c1 [prog] [n0=M] [dump] [threads=T]`. G = chain mode,
  E = zero-free mode (validation), P = permissive. Output does not depend on T.
- **Production-scale regression.** `python3 val_proddone.py EXE THREADS [--only mode:c0:c1] LOG...` re-runs
  complete logged invocations and compares every counter with the production DONE line.
  `compare_cdone.py N N0 DIR_A DIR_B` compares per-c counts of two log sets.
  `rerun_production.sh DIR EXE THREADS shell|base|both` re-runs the whole logged 10^40000 computation.
- **CUDA variant.** `./val_cuda.sh LOGDIR` runs its whole validation: kernel tests, suite, production slices,
  cross-variant dumps and a logged invocation.
  - `./gsieve_k13_cuda kerneltest G|E NJ SEED` checks the GPU's stage-1/2 decisions on the whole table for NJ lift
    indices, and on 2^20 random points against 2^E mod 10^26.
  - `EXE=./gsieve_k13_cuda prod.py ...` runs production with it.
- **Shells with checks.** `./run_shells.sh RUNDIR "N0:N ..." gpu|cpu|both` runs production shells with coverage,
  survival-rate checks and log archives (gpu), and the independent second pass (cpu).
  - A MEMBER tripwire stops everything at once.
  - Coverage of archived runs: `scripts/check_coverage.py` (RUNS list).
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
- **Windows checkouts.** Git for Windows checks out CRLF. `.gitattributes` forces LF: the record hashes and the
  shell scripts need it.
- **WSL `/tmp` is wiped** when the distribution stops between commands. Keep scratch files on disk (e.g.
  `/root/scratch`), and never compare files without checking that they exist.
- **Non-ECC memory.** GeForce memory and desktop RAM have no ECC.
  - The CUDA variant re-checks its tables every minute and writes results only after a passed check (fail-stop).
  - Double-check new production with an independent run (another variant or a second pass), and compare the
    counters (`compare_cdone.py`, DONE lines).
- **CUDA copies.** A pageable `cudaMemcpy` may return before its DMA lands. Synchronize before reading the data back
  on another stream.
