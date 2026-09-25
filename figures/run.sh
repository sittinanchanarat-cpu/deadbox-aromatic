#!/usr/bin/env bash
# Main-figure data panels (SC Lab style), written to figures/out/ as PDF, TIFF and PNG.
# Structure renderings (anchor, stacker, Gln, adenine) are made separately.
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
for f in fig1_definition fig2_contact_pivot fig3_trp_census fig4_evolution fig5_variants; do
  $PY "$f.py"
done
