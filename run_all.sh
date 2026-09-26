#!/usr/bin/env bash
# Rebuild every result of the DEAD-box upstream-anchor project from public data.
#   ./run_all.sh              all stages (several hours; stage 3 trees take ~2 h on 8 threads)
#   ./run_all.sh 2 3          only the listed stages (each needs the ones before it to have run)
# Stages
#   0  PDB structure census: anchor geometry, N6 contact, partner complexes, loop conformation
#   1  sequence census: UniProt reference proteomes, anchor identity per subfamily and genus
#   2  AlphaFold models: is the called aromatic inserted at adenine N6?
#   3  phylogeny: ancestral anchor state, switches, independent origins of Trp
#   4  human variants: ClinVar, gnomAD v4, cBioPortal
#   figures  main-figure data panels (figures/out/)
# Environment: micromamba env "helicase" (env/helicase.yml). Data versions: VERSIONS.md.
set -euo pipefail
cd "$(dirname "$0")"
stages=("$@")
[ ${#stages[@]} -eq 0 ] && stages=(0 1 2 3 4)
for s in "${stages[@]}"; do
  echo "=== stage $s: $(date '+%F %T')"
  "stage$s/run.sh"
done
echo "=== figures: $(date '+%F %T')"
figures/run.sh
echo "=== supplementary tables: $(date '+%F %T')"
MAMBA_ROOT_PREFIX=~/micromamba ~/.local/bin/micromamba run -n helicase python supplement/build_tables.py
./record_versions.sh
echo "=== done: $(date '+%F %T')"
