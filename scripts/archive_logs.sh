#!/bin/sh
# Reproducible .tar.xz archive of a directory of chain-sieve logs (same tar flags as logs_tarball.sh), written next to
# the directory, with the SHA-256 of every archived file listed in <dir>.sha256 (check after unpacking with
# `sha256sum -c`). The directory itself is left in place.
# usage: scripts/archive_logs.sh DIR      (e.g. validation_pc/task3_cuda/rerun_cuda_final)
set -e
d=${1%/}; [ -d "$d" ] || { echo "no directory $d"; exit 1; }
parent=$(dirname "$d"); base=$(basename "$d")
( cd "$parent" && LC_ALL=C find "$base" -type f | LC_ALL=C sort | xargs sha256sum > "$base.sha256" )
tar -C "$parent" --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-01-01 00:00Z' -cJf "$d.tar.xz" "$base"
echo "$d.tar.xz: $(wc -c < "$d.tar.xz") bytes, $(wc -l < "$d.sha256") files; SHA-256 $(sha256sum "$d.tar.xz" | cut -c1-16)"
