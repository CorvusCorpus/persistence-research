#!/bin/sh
# Print SHA-256 sums of the record files (CLAUDE.md rule 3): production and run logs, proof records,
# exhaustive lists. `make manifest` writes them to MANIFEST.sha256; `make records-check` verifies them.
# Checkpoints (ckpt*/, bagres*/) are resumable state, not records, and are not listed.
cd "$(dirname "$0")/.." || exit 1
LC_ALL=C
export LC_ALL
{
  ls work/gsieve/logs/* work/logs/*
  ls work/proof*.json work/certificate.json work/base*_cum.json work/base*_proof.json
  ls work/Z.json work/trees_b6_D800.json work/gsieve/G_small.json
} | sort | xargs sha256sum
