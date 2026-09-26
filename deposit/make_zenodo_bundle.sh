#!/usr/bin/env bash
# Build the Zenodo deposit: a code snapshot plus the data too large or too
# time-dependent for git, each group as one .tar.gz, with a SHA-256 manifest.
#
#   deposit/make_zenodo_bundle.sh [--with-cbio] [version]
#
#   --with-cbio  also archive the raw cBioPortal dump (every mutation returned for
#                the 37 genes, with patient IDs; each study has its own reuse terms).
#                By default only the filtered mutations (Table S6B, in the code
#                snapshot) are deposited, and stage4/scripts/02_cbioportal.py repeats
#                the query.
#   version      defaults to the version in CITATION.cff
#
# Output: deposit/out/deadbox-aromatic-<version>/  (upload every file in it)
# Left out on purpose: PDB coordinates (re-fetched by ID from the tracked entry lists),
# MMseqs2 working files and IQ-TREE checkpoints/distance matrices (rebuilt from the
# archived inputs), manuscript drafts and notes (manuscript/, docs/).
set -euo pipefail
cd "$(dirname "$0")/.."

WITH_CBIO=0
if [[ ${1:-} == --with-cbio ]]; then WITH_CBIO=1; shift; fi
VERSION=${1:-$(awk -F'"' '/^version:/{print $2}' CITATION.cff)}
DOI=$(awk -F'"' '/^doi:/{print $2}' CITATION.cff)
NAME="deadbox-aromatic-$VERSION"
OUT="deposit/out/$NAME"

if [[ -n $(git status --porcelain --untracked-files=no) ]]; then
  echo "Uncommitted changes to tracked files: commit first, so the code snapshot matches HEAD." >&2
  exit 1
fi
for f in stage1/data/dead_refprot.tsv.gz stage1/results/census_all.tsv stage3/work/prp5/t.state \
         stage2/data/af stage4/data/gnomad supplement/Supplementary_Tables_S1-S6.xlsx; do
  [[ -e $f ]] || { echo "Missing $f: run run_all.sh (or the stage) first." >&2; exit 1; }
done

rm -rf "$OUT"
mkdir -p "$OUT"
COMMIT=$(git rev-parse --short HEAD)
# deterministic archives: fixed order, owner and timestamp
TAR=(tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="$(git log -1 --format=%cI)")
pack() {  # pack <archive name> <paths...>
  local name=$1; shift
  "${TAR[@]}" -cf - "$@" | gzip -n -9 > "$OUT/$name.tar.gz"
  printf '  %-28s %s\n' "$name.tar.gz" "$(du -h "$OUT/$name.tar.gz" | cut -f1)"
}

echo "Building $OUT from commit $COMMIT"

# 1. code, curated inputs and small result tables at HEAD (manuscript/ and docs/ are not tracked)
git archive --format=tar --prefix="$NAME/" HEAD \
  | gzip -n -9 > "$OUT/code.tar.gz"
printf '  %-28s %s\n' code.tar.gz "$(du -h "$OUT/code.tar.gz" | cut -f1)"

# 2. source snapshots that cannot be downloaded again in the same form
pack uniprot_2026_03 stage1/data
pack variant_databases stage4/data/gnomad stage4/data/alphamissense \
     $([[ $WITH_CBIO == 1 ]] && echo stage4/data/cbio_mutations.tsv.gz)
pack alphafold_models stage2/data/selection.tsv stage2/data/af stage4/data/af

# 3. per-sequence and per-tree results
pack sequence_census stage1/results/census_all.tsv stage1/results/census_balanced_assigned.tsv
mapfile -t TREES < <(find stage3/work -mindepth 2 -maxdepth 2 -type f \
  \( -name 'aln.fa' -o -name 'aln.mask.fa' -o -name 'in.fa' -o -name 'meta.tsv' -o -name 'anchor_col.txt' \
     -o -name 't.treefile' -o -name 't.contree' -o -name 't.iqtree' -o -name 't.state' -o -name 't.rate' \) | sort)
pack gene_trees "${TREES[@]}" stage3/results/family_tree_trp.txt

# 4. supplementary tables, licence, citation, description
cp supplement/Supplementary_Tables_S1-S6.xlsx LICENSE CITATION.cff "$OUT/"
cat > "$OUT/README.md" <<EOF
# DEAD-box aromatic anchor: code and data ($VERSION)

Deposit for "An aromatic anchor, not the Q motif, marks the DEAD-box ATP site"
(Chanarat). DOI: https://doi.org/$DOI. Code: https://github.com/sittinanchanarat-cpu/deadbox-aromatic
(commit $COMMIT). Built on $(date +%F).

| File | Contents |
|---|---|
| code.tar.gz | Scripts for every stage, environment (env/helicase.yml), curated family labels, PDB entry lists and all small result tables. \`run_all.sh\` rebuilds everything; see README.md inside. |
| uniprot_2026_03.tar.gz | UniProtKB release 2026_03 download: DEAD-box Q-motif entries of reference proteomes, proteome list, Hub1/UBL5 entries (stage1/data). |
| variant_databases.tar.gz | gnomAD v4 / ClinVar API responses per gene and AlphaMissense tables per protein, as retrieved (stage4/data).$([[ $WITH_CBIO == 1 ]] && echo " Includes the raw cBioPortal mutation dump.") |
| alphafold_models.tar.gz | AlphaFold DB (version 6) models analysed: 261 census orthologues (stage2/data/af, with selection.tsv) and 35 human proteins (stage4/data/af). |
| sequence_census.tar.gz | Anchor call for all 156,510 sequences (census_all.tsv) and the genus-balanced set with subfamily assignments (census_balanced_assigned.tsv). |
| gene_trees.tar.gz | Per gene tree (stage3/work/<subfamily>/): input sequences, alignment, masked alignment, metadata, anchor column, IQ-TREE tree and report, marginal ancestral states (t.state) and site rates (t.rate); family-level tree (stage3/work/family/). |
| Supplementary_Tables_S1-S6.xlsx | Supplementary Tables S1-S6 of the article. |
| MANIFEST.sha256 | SHA-256 checksums of every file above. |

To reproduce: unpack code.tar.gz, then unpack the data archives inside the resulting
directory (paths are relative to the repository root), create the environment from
env/helicase.yml and run the stage scripts. PDB coordinates are fetched by ID.

Licences: code under MIT (LICENSE). Tables and trees produced here: CC BY 4.0.
Third-party data (UniProt, AlphaFold DB, AlphaMissense, gnomAD, ClinVar, cBioPortal)
remain under the terms of their providers.
EOF

( cd "$OUT" && sha256sum -- * | grep -v ' MANIFEST.sha256$' > MANIFEST.sha256 )
echo "Total: $(du -sh "$OUT" | cut -f1) in $OUT"
