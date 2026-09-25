#!/usr/bin/env bash
# Write VERSIONS.md (data releases and access dates) and env/helicase.yml (software).
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
MM="$HOME/.local/bin/micromamba"
mkdir -p env
$MM env export -n helicase > env/helicase.yml
uniprot=$(curl -sI "https://rest.uniprot.org/uniprotkb/search?query=accession:P21372&size=1" | grep -i '^x-uniprot-release:' | cut -d' ' -f2 | tr -d '\r')
cbio=$(curl -s https://www.cbioportal.org/api/info | python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('portalVersion','?'), d.get('dbVersion','?'))" 2>/dev/null || echo "?")
mtime() { [ -e "$1" ] && date -r "$1" '+%F' || echo "not downloaded"; }
cat > VERSIONS.md <<MD
# Data and software versions

Recorded $(date '+%F %T') by record_versions.sh.

| Source | Version / access date | Used in |
|---|---|---|
| UniProtKB (current release) | ${uniprot:-?} | stage 1, 2, 4 |
| UniProtKB download (IPR014014, reference proteomes) | $(mtime stage1/data/dead_refprot.tsv.gz) | stage 1 |
| RCSB PDB entries (Stage 0 set) | $(mtime stage0/data/entries.tsv) | stage 0 |
| RCSB PDB partner complexes | $(mtime stage0/data/partner_entities.tsv) | stage 0 |
| AlphaFold DB models (v6) | $(mtime stage2/data/selection.tsv) | stage 2 |
| gnomAD v4 (gnomad_r4) + ClinVar via gnomAD API | $(mtime stage4/data/gnomad) | stage 4 |
| cBioPortal public portal (portal, db) | ${cbio} / $(mtime stage4/data/cbio_mutations.tsv.gz) | stage 4 |

Software: see env/helicase.yml (Python 3.12, gemmi, Biopython, pandas, SciPy, MAFFT,
MMseqs2, IQ-TREE 3, HMMER).
MD
echo "wrote VERSIONS.md and env/helicase.yml"
