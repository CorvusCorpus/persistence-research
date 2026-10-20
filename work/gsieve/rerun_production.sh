#!/bin/bash
# Re-run the logged chain-sieve production (base run N=36000 and/or shell 36000 < m <= 40000) with another
# binary and compare it with the production logs:
#   - every complete production invocation is re-run over the same c-range (val_proddone.py) and all of its
#     counters (points, s1, s2, s3, s4, gmp, full, members) must equal the production DONE line;
#   - the ranges of killed production invocations (no DONE line) are re-run too, so the re-run covers every c;
#   - aggregate.py must show full coverage and no members, and compare_cdone.py must find identical per-c counts.
# Re-run logs go to DIR (aggregate.py format); complete logs already there are reused, so the script resumes.
# usage: rerun_production.sh DIR EXE THREADS shell|base|both      (run from work/gsieve)
set -u
R=$1; EXE=$2; T=$3; WHAT=$4
mkdir -p "$R"; st=0
if [ "$WHAT" = shell ] || [ "$WHAT" = both ]; then
  python3 val_proddone.py "$EXE" "$T" --logdir "$R" logs/prod_N40000_n036000_*.log || st=1
  python3 aggregate.py 40000 36000 "$R" || st=1
  python3 compare_cdone.py 40000 36000 logs "$R" || st=1
fi
if [ "$WHAT" = base ] || [ "$WHAT" = both ]; then
  python3 val_proddone.py "$EXE" "$T" --logdir "$R" logs/prod_N36000_n0-1_*.log || st=1
  for job in "main 0 7540" "direct 0 11381"; do         # killed production invocations: re-run their ranges
    set -- $job; f="$R/prod_N36000_n0-1_$1_$2_$3.log"
    if ! grep -q '^DONE' "$f" 2>/dev/null; then
      echo "# invocation start c=$2 at $(date) exe=$EXE threads=$T (range of a killed production invocation)" > "$f"
      "$EXE" "$1" G 36000 "$2" "$3" prog "threads=$T" >> "$f" 2>&1 || st=1
    fi
    grep '^DONE' "$f"
  done
  python3 aggregate.py 36000 -1 "$R" || st=1
  python3 compare_cdone.py 36000 -1 logs "$R" || st=1
fi
exit $st
