#!/usr/bin/env bash
# Stage 4: human variants at the anchor and its shell (needs stage1 results)
#   01  ClinVar + gnomAD v4 (gnomAD GraphQL API; rate-limited, ~10 min)
#   02  somatic cancer mutations from all public cBioPortal studies
#   03  depletion / enrichment tests at the anchor and its shell
#   04  AlphaMissense at the anchor vs other buried aromatics; aromatic swaps
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
mkdir -p data results
$PY scripts/01_variants.py
$PY scripts/02_cbioportal.py
$PY scripts/03_constraint.py
$PY scripts/04_alphamissense.py
