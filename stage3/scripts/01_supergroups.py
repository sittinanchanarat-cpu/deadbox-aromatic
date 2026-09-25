#!/usr/bin/env python
"""Anchor residue by eukaryotic supergroup, per subfamily.

A residue found across both sides of the eukaryotic root (Amorphea vs
Diaphoretickes, plus Discoba/Metamonada wherever the root is placed) is most
parsimoniously a LECA state. Supergroups come from the UniProt lineage names:
  Amorphea        Opisthokonta, Amoebozoa, Apusozoa
  Archaeplastida  Viridiplantae, Rhodophyta
  SAR             Stramenopiles, Alveolata, Rhizaria (Sar)
  Haptista, Cryptophyceae, Discoba, Metamonada

Output
  results/supergroups.tsv   subfamily x supergroup: n genera, % W, % F
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT.parent / "stage1"
GROUPS = [  # first match wins
    ("Opisthokonta", "Amorphea"), ("Amoebozoa", "Amorphea"), ("Apusozoa", "Amorphea"),
    ("Viridiplantae", "Archaeplastida"), ("Rhodophyta", "Archaeplastida"),
    ("Sar (clade)", "SAR"), ("Haptista", "Haptista"), ("Cryptophyceae", "Cryptophyceae"),
    ("Discoba", "Discoba"), ("Metamonada", "Metamonada"),
]
ORDER = ["Amorphea", "Archaeplastida", "SAR", "Haptista", "Cryptophyceae", "Discoba", "Metamonada"]


def supergroup(lineage):
    return next((g for key, g in GROUPS if key in str(lineage)), "other")


def main():
    lin = pd.read_csv(S1 / "data" / "dead_refprot.tsv.gz", sep="\t",
                      usecols=["Entry", "Taxonomic lineage"]).set_index("Entry")["Taxonomic lineage"]
    c = pd.read_csv(S1 / "results" / "census_balanced_assigned.tsv", sep="\t")
    c = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool)
          & ~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    c = c.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    c["supergroup"] = c.acc.map(lin).map(supergroup)
    c["anchor"] = c.ua_aa.where(c.ua_found.fillna(False).astype(bool), "-")
    rows = []
    for (sf, sg), g in c.groupby(["subfamily_final", "supergroup"]):
        rows.append({"subfamily": sf, "supergroup": sg, "n": len(g),
                     "pct_W": round(100 * (g.anchor == "W").mean()), "pct_F": round(100 * (g.anchor == "F").mean())})
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "results" / "supergroups.tsv", sep="\t", index=False)
    cell = d.assign(cell=d.pct_W.astype(str) + "%W/" + d.n.astype(str)).pivot(index="subfamily", columns="supergroup", values="cell")
    with pd.option_context("display.width", 250):
        print(cell.reindex(columns=[g for g in ORDER if g in cell.columns]).fillna("").to_string())


if __name__ == "__main__":
    main()
