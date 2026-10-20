#!/bin/sh
# Release tarball of all run logs plus the record manifest, for attaching to a GitHub release.
# Output: dist/persistence-logs-<commit>.tar.xz and its .sha256 (dist/ is git-ignored).
cd "$(dirname "$0")/.." || exit 1
set -e
sha256sum --quiet -c MANIFEST.sha256
rev=$(git rev-parse --short=12 HEAD)
mkdir -p dist
out="dist/persistence-logs-$rev.tar.xz"
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-01-01 00:00Z' \
    -cJf "$out" MANIFEST.sha256 work/gsieve/logs work/logs
sha256sum "$out" > "$out.sha256"
ls -l "$out"; cat "$out.sha256"
