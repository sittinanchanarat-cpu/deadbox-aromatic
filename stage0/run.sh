#!/usr/bin/env bash
# Stage 0: PDB structure census of the upstream anchor (UA)
#   01-03  nucleotide-bound SF2 structures -> stacking geometry -> UA -> family summary
#          (family labels are hand-curated in data/family_curation.tsv; keep that file)
#   04-05  partner-complex census (DEAD-box chains in complexes; no symmetry mates)
#   06     UA -- adenine N6 geometry
#   07     per-residue partner contacts and loop conformation vs free structures
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
PY="$HOME/.local/bin/micromamba run -n helicase python"
$PY scripts/01_fetch_entries.py "$@"
$PY scripts/02_stacking_geometry.py
$PY scripts/02b_upstream_anchor.py
$PY scripts/03_families_and_summary.py
$PY scripts/04_partner_fetch.py
$PY scripts/05_partner_contacts.py
$PY scripts/06_n6_geometry.py
$PY scripts/07_loop_detail.py
