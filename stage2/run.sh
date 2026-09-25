#!/usr/bin/env bash
# Stage 2: AlphaFold check of the anchor (needs stage1 results)
#   01  pick orthologs (Trp subfamilies, Rok1, Phe controls) and fetch AFDB models
#   02  transplant adenine from 2DB3 (Vasa-AMPPNP) and measure every aromatic at Q-20..Q-30
#   03  packing-shell covariation, W vs F (between and within subfamilies)
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
mkdir -p data/af results
$PY scripts/01_select_fetch.py
$PY scripts/02_af_geometry.py
$PY scripts/03_shell_covariation.py
