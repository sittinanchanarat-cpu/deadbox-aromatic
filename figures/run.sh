#!/usr/bin/env bash
# Main-figure data panels (SC Lab style), written to figures/out/ as PDF, TIFF and PNG.
# Then the structure renderings (PyMOL, ray-traced; ~10 min on a laptop CPU).
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
for f in fig1_definition fig2_contact_pivot fig3_aromatic_anchor fig4_evolution fig5_variants; do
  $PY "$f.py"
done
$PY structures.py
