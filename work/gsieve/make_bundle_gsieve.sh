#!/bin/sh
# Package the chain-sieve code, validation scripts, write-up and production logs.
cd "$(dirname "$0")"
rm -f chain_sieve.zip
tar czf prod_logs.tar.gz logs/prod_N*.log logs/drv*.out 2>/dev/null
zip -q chain_sieve.zip CHAIN_SIEVE.md gsieve.c gconst.py gconst.h prod.py aggregate.py tiny.py G_small.json \
  val_dump.py val_e2e.py val_direct_odd.py val_prodrow.py thread_sample.txt prod_logs.tar.gz
ls -la chain_sieve.zip
