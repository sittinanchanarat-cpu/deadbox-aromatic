#!/usr/bin/env python
"""AlphaMissense at the anchor and its shell, human DEAD-box proteins.

For each gene in data/genes.tsv the AlphaMissense substitution table (AlphaFold DB,
per UniProt accession) and the AlphaFold model are fetched.
  * per position: mean AlphaMissense pathogenicity over the 19 substitutions
  * site classes UA, UA+1, Q-stacker, Q-4, Q: that mean and its percentile among
    (a) all residues and (b) buried aromatics of the same protein (F/W/Y with
    relative SASA < 0.1 in the AlphaFold model, the anchor excluded) -- the fair
    control for "is the anchor more critical than any buried ring?"
  * aromatic swaps at the anchor: F->W, F->Y, W->F, W->Y, and F/W->L, ->A as references

Output
  data/alphamissense/<acc>.csv, data/af/<acc>.cif.gz   (cached downloads)
  results/alphamissense_sites.tsv   one row per gene x site class
  results/alphamissense_swaps.tsv   anchor substitutions per gene
"""
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parents[1]
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
           "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
           "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}  # Tien et al. 2013
SITES = ["UA", "UA+1", "Q-stacker", "Q-4", "Q"]


def fetch(acc):
    am, af = ROOT / "data" / "alphamissense" / f"{acc}.csv", ROOT / "data" / "af" / f"{acc}.cif.gz"
    if not (am.exists() and af.exists()):
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60) as r:
            meta = json.load(r)[0]
        am.parent.mkdir(parents=True, exist_ok=True)
        af.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(meta["amAnnotationsUrl"], am)
        with urllib.request.urlopen(meta["cifUrl"], timeout=120) as r:
            af.write_bytes(gzip.compress(r.read()))
    return pd.read_csv(am), af


def rsasa(af):
    with gzip.open(af, "rt") as fh:
        s = MMCIFParser(QUIET=True).get_structure("m", io.StringIO(fh.read()))
    ShrakeRupley().compute(s[0], level="R")
    return {r.id[1]: (r.get_resname(), r.sasa / MAX_ASA.get(r.get_resname(), 200)) for r in s[0].get_residues()}


def main():
    g = pd.read_csv(ROOT / "data" / "genes.tsv", sep="\t")
    sites, swaps = [], []
    for r in g.itertuples():
        try:
            am, af = fetch(r.acc)
        except Exception as ex:
            print(f"{r.gene}: {ex}", file=sys.stderr)
            continue
        am["pos"] = am.protein_variant.str[1:-1].astype(int)
        am["ref"], am["alt"] = am.protein_variant.str[0], am.protein_variant.str[-1]
        per = am.groupby("pos").am_pathogenicity.mean()
        sa = rsasa(af)
        pos = json.loads(r.positions)
        ua = pos["UA"]
        buried_arom = [p for p, (rn, x) in sa.items() if rn in ("PHE", "TRP", "TYR") and x < 0.1 and p != ua and p in per]
        for s in SITES:
            p = pos[s]
            if p not in per:
                continue
            v = per[p]
            sites.append({"gene": r.gene, "anchor": r.ua, "site": s, "pos": p, "residue": sa.get(p, ("?",))[0],
                          "rsasa": round(sa.get(p, ("", float("nan")))[1], 3), "am_mean": round(v, 3),
                          "pct_all": round(100 * (per < v).mean(), 1),
                          "pct_buried_aromatics": round(100 * (per[buried_arom] < v).mean(), 1) if buried_arom else None,
                          "n_buried_aromatics": len(buried_arom)})
        x = am[am.pos == ua].set_index("alt").am_pathogenicity
        ref = r.ua[0]
        swaps.append({"gene": r.gene, "anchor": r.ua, **{f"{ref}->{a}": round(x.get(a, float("nan")), 3) for a in "FWYLA" if a != ref}})
        print(f"{r.gene} done", file=sys.stderr, flush=True)
    s = pd.DataFrame(sites)
    s.to_csv(ROOT / "results" / "alphamissense_sites.tsv", sep="\t", index=False)
    w = pd.DataFrame(swaps)
    w.to_csv(ROOT / "results" / "alphamissense_swaps.tsv", sep="\t", index=False)
    print(s.groupby("site")[["am_mean", "pct_all", "pct_buried_aromatics"]].median().round(2).to_string())


if __name__ == "__main__":
    main()
