#!/usr/bin/env bash
# Stage 3: phylogeny of the anchor (needs stage1 results; ~2 h on 8 threads)
#   01  anchor identity by eukaryotic supergroup
#   02  per-subfamily gene trees (MAFFT, IQ-TREE 3 LG+G4 -fast, -asr, --rate)
#   03  ancestral state at the anchor, switches, site rate
#   04  family-level tree (UFBoot 1000): are the Trp subfamilies related?
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
mkdir -p work results
$PY scripts/01_supergroups.py
$PY scripts/02_build_trees.py
$PY scripts/03_asr.py
$PY scripts/04_family_tree.py
