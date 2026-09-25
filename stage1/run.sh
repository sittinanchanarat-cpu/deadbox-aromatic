#!/usr/bin/env bash
# Stage 1: sequence census of the anchor across DEAD-box proteins
#   downloads  UniProtKB reference proteomes: IPR014014 (Q motif) sequences, Hub1/UBL5,
#              eukaryotic reference proteome list with BUSCO
#   01         anchor call, taxonomy, one proteome per genus
#   mmseqs     unlabelled eukaryotic sequences vs labelled ones (nearest homolog)
#   02-03      Trp vs Hub1 table; subfamily assignment of unlabelled sequences
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=~/micromamba
MM="$HOME/.local/bin/micromamba run -n helicase"
UNIPROT=https://rest.uniprot.org
mkdir -p data results work
if [ ! -s data/dead_refprot.tsv.gz ]; then
  curl -sf -D data/uniprot_release.headers -o data/dead_refprot.tsv.gz \
    "$UNIPROT/uniprotkb/stream?compressed=true&format=tsv&fields=accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,lineage_ids,lineage,xref_proteomes,length,ft_motif,ft_domain,protein_families,sequence&query=xref:interpro-IPR014014+AND+keyword:KW-1185"
fi
if [ ! -s data/hub1_refprot.tsv.gz ]; then
  curl -sf -o data/hub1_refprot.tsv.gz \
    "$UNIPROT/uniprotkb/stream?compressed=true&format=tsv&fields=accession,reviewed,protein_name,gene_names,organism_name,organism_id,xref_proteomes,length,sequence&query=(xref:interpro-IPR039732+OR+xref:pfam-PF00240+AND+(protein_name:hub1+OR+protein_name:ubl5+OR+gene:ubl5+OR+gene:hub1))+AND+keyword:KW-1185+AND+taxonomy_id:2759"
fi
if [ ! -s data/euk_proteomes.tsv ]; then
  curl -sf -o data/euk_proteomes.tsv \
    "$UNIPROT/proteomes/stream?format=tsv&fields=upid,organism,organism_id,protein_count,busco,cpd&query=(taxonomy_id:2759)+AND+(reference:true)"
fi
$MM python scripts/01_census.py
$MM python scripts/02_trp_hub1.py
# nearest labelled homolog for unlabelled eukaryotic sequences
$MM python - <<'PY'
import pandas as pd
s = pd.read_csv("data/dead_refprot.tsv.gz", sep="\t", usecols=["Entry", "Sequence"]).set_index("Entry").Sequence
c = pd.read_csv("results/census_balanced.tsv", sep="\t")
c = c[c.domain == "Eukaryota"]
for name, sub in (("query", c[c.subfamily == "unassigned"]), ("ref", c[c.subfamily != "unassigned"])):
    with open(f"work/{name}.fa", "w") as fh:
        for a in sub.acc:
            fh.write(f">{a}\n{s[a]}\n")
PY
$MM mmseqs easy-search work/query.fa work/ref.fa work/hits.m8 work/tmp --max-seqs 50 -e 1e-10 \
  --threads 8 --format-output query,target,pident,bits,qcov,tcov > work/mmseqs.log 2>&1
$MM python scripts/03_assign_unlabelled.py
