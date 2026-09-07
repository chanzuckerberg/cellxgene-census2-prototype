#!/usr/bin/env bash
# Build the workflow bundle that AWS HealthOmics runs.
#
# HealthOmics has no way to pass `-profile` to Nextflow, so the healthomics
# profile's settings have to be the defaults inside the bundle. This script
# copies the pipeline and appends that include, rather than the repository
# carrying a HealthOmics-only default that would surprise a local run.
#
# Usage: tools/build-healthomics-bundle.sh [outfile.zip]
set -euo pipefail

cd "$(dirname "$0")/.."
OUT="${1:-$PWD/healthomics-bundle.zip}"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# Only what the workflow needs at run time. No tests, no venv, no docker context.
for item in main.nf nextflow.config conf modules; do
    cp -R "$item" "$STAGE/"
done

cat >> "$STAGE/nextflow.config" <<'CFG'

// Appended by tools/build-healthomics-bundle.sh. HealthOmics cannot select a
// Nextflow profile, so the platform settings are applied unconditionally here.
includeConfig 'conf/healthomics.config'
CFG

rm -f "$OUT"
( cd "$STAGE" && zip -qr "$OUT" . )
echo "bundle: $OUT ($(du -h "$OUT" | cut -f1))"
unzip -l "$OUT" | tail -n +4 | head -12
