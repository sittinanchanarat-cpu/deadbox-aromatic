#!/usr/bin/env python
"""Somatic cancer mutations at the anchor and its shell, from the public cBioPortal.

Genes and positions come from data/genes.tsv (01_variants.py). All public
MUTATION_EXTENDED profiles are queried (batches of 40) for the genes' Entrez IDs.
A mutation is kept when its protein position is one of the listed sites and its
reference residue (first letter of proteinChange) matches UniProt. The same
patient + protein change seen in several studies (overlapping cohorts, e.g. TCGA
legacy and PanCancer Atlas) is counted once. Cell-line collections (CCLE, NCI-60,
other cell-line studies) re-profile the same lines across releases under different
IDs, so they are listed but left out of the counts (tumour_only = False).

Output
  data/cbio_mutations.tsv.gz        every mutation returned for these genes (cached)
  results/cancer_positions.tsv      unique patient-level mutations at the listed sites
  results/cancer_missense_counts.tsv  per gene: missense patients at UA vs per-residue average
"""
import gzip
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
API = "https://www.cbioportal.org/api"


def get(path):
    with urllib.request.urlopen(f"{API}{path}", timeout=300) as r:
        return json.load(r)


def post(path, body):
    req = urllib.request.Request(f"{API}{path}", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=600) as r:
        return json.load(r)


def entrez(acc):
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.tsv?fields=xref_geneid,sequence", timeout=60) as r:
        _, row = r.read().decode().strip().split("\n")
    gid, seq = row.split("\t")
    return int(gid.strip(";").split(";")[0]), seq


def main():
    g = pd.read_csv(ROOT / "data" / "genes.tsv", sep="\t")
    info = {r.acc: entrez(r.acc) for r in g.itertuples()}
    g["entrez"] = g.acc.map(lambda a: info[a][0])
    seqs = {r.gene: info[r.acc][1] for r in g.itertuples()}
    cache = ROOT / "data" / "cbio_mutations.tsv.gz"
    if cache.exists():
        m = pd.read_csv(cache, sep="\t")
    else:
        profiles = [p["molecularProfileId"] for p in get("/molecular-profiles?projection=SUMMARY&pageSize=100000")
                    if p.get("molecularAlterationType") == "MUTATION_EXTENDED"]
        rows = []
        for i in range(0, len(profiles), 40):
            batch = profiles[i:i + 40]
            try:
                rows += post("/mutations/fetch?projection=SUMMARY",
                             {"entrezGeneIds": [int(x) for x in g.entrez], "molecularProfileIds": batch})
            except Exception as ex:
                print(f"batch {i}: {ex}", file=sys.stderr)
            print(f"{i + len(batch)}/{len(profiles)} profiles, {len(rows)} mutations", file=sys.stderr, flush=True)
        m = pd.DataFrame(rows)[["molecularProfileId", "studyId", "patientId", "sampleId", "entrezGeneId",
                                "proteinChange", "mutationType", "proteinPosStart"]]
        m.to_csv(cache, sep="\t", index=False, compression="gzip")
    gene_of = dict(zip(g.entrez, g.gene))
    m["gene"] = m.entrezGeneId.map(gene_of)
    m = m.dropna(subset=["proteinPosStart"])
    m["pos"] = m.proteinPosStart.astype(int)
    m["ref"] = m.proteinChange.str[0]
    m["uniprot_match"] = [p <= len(seqs[gn]) and seqs[gn][p - 1] == ref for gn, p, ref in zip(m.gene, m.pos, m.ref)]
    m = m.drop_duplicates(["gene", "patientId", "proteinChange"])
    m["tumour_only"] = ~m.studyId.str.contains(r"ccle|cellline|nci60|_cell_line|celline", case=False)
    hits, counts = [], []
    for r in g.itertuples():
        pos = json.loads(r.positions)
        where = {v: k for k, v in pos.items()}
        x = m[(m.gene == r.gene) & m.pos.isin(where)]
        for y in x.itertuples():
            hits.append({"gene": r.gene, "anchor": r.ua, "site": where[y.pos], "pos": y.pos,
                         "change": y.proteinChange, "type": y.mutationType, "study": y.studyId,
                         "patient": y.patientId, "uniprot_match": y.uniprot_match,
                         "tumour_only": y.tumour_only})
        mis = m[(m.gene == r.gene) & (m.mutationType == "Missense_Mutation") & m.uniprot_match & m.tumour_only]
        L = len(seqs[r.gene])
        ua = pos["UA"]
        counts.append({"gene": r.gene, "anchor": r.ua, "missense_patients": len(mis),
                       "per_residue_mean": round(len(mis) / L, 3), "at_UA": int((mis.pos == ua).sum()),
                       "in_loop": int(mis.pos.between(ua - 8, ua + 3).sum())})
    pd.DataFrame(hits).to_csv(ROOT / "results" / "cancer_positions.tsv", sep="\t", index=False)
    pd.DataFrame(counts).to_csv(ROOT / "results" / "cancer_missense_counts.tsv", sep="\t", index=False)
    print(f"{len(hits)} cancer mutations at listed sites", file=sys.stderr)


if __name__ == "__main__":
    main()
