#!/usr/bin/env bash
# Build the workflow bundle that AWS HealthOmics runs.
#
# The healthomics profile's settings are the defaults inside the bundle: this
# script copies the pipeline and appends the include. The repository therefore
# carries no HealthOmics-only default that would change a local run. StartRun
# can also select a profile (`--engine-settings`), which would remove this step;
# the pipeline does not use that yet.
#
# Usage: tools/build-healthomics-bundle.sh [outfile.zip]
set -euo pipefail

cd "$(dirname "$0")/.."
# Absolute: zip runs from the staging directory, so a relative path would land
# there and be removed with it.
ARG="${1:-healthomics-bundle.zip}"
OUT="$(cd "$(dirname "$ARG")" && pwd)/$(basename "$ARG")"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# Only what the workflow needs at run time. No tests, no venv, no docker context.
for item in main.nf nextflow.config conf modules; do
    cp -R "$item" "$STAGE/"
done

# Never ship bytecode, even if the working tree holds some.
find "$STAGE" -name __pycache__ -type d -prune -exec rm -rf {} +

cat >> "$STAGE/nextflow.config" <<'CFG'

// Appended by tools/build-healthomics-bundle.sh, so that the platform settings
// apply without a profile.
includeConfig 'conf/healthomics.config'
CFG

rm -f "$OUT"
( cd "$STAGE" && zip -qr "$OUT" . )
echo "bundle: $OUT ($(du -h "$OUT" | cut -f1))"
unzip -l "$OUT" | tail -n +4 | head -12
