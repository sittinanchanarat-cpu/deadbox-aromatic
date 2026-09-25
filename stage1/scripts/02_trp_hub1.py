#!/usr/bin/env python
"""Trp anchors vs Hub1: which lineages break the subfamily rule, and does Hub1 follow?

Per genus (the proteome chosen in 01_census.py):
  hub1        a Hub1/UBL5 sequence (IPR039732, or a Ubl named Hub1/UBL5) in that proteome
  busco_c     BUSCO completeness of the proteome; Hub1 absence is only called
              when busco_c >= 90 (otherwise 'unknown')
  <subfamily> anchor residue of that genus's representative (F/W/Y, '-' = no aromatic,
              blank = subfamily not found in the proteome)

Output
  results/genus_table.tsv        one row per eukaryotic genus
  results/trp_exceptions.tsv     genera whose Prp5/Prp28/Mak5/Spb4 anchor is not W
  printed: anchor identity x Hub1 status for each Trp subfamily, by major clade
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRP_SF = ["DDX46/PRP5", "DDX23/PRP28", "DDX24/MAK5", "DDX55/SPB4"]
BUSCO_MIN = 90.0


def main():
    b = pd.read_csv(ROOT / "results" / "census_balanced.tsv", sep="\t")
    b = b[(b.domain == "Eukaryota") & b.q_ok]
    b["anchor"] = b.ua_aa.where(b.ua_found.astype(bool), "-")

    prot = pd.read_csv(ROOT / "data" / "euk_proteomes.tsv", sep="\t")
    prot["busco_c"] = prot.BUSCO.str.extract(r"C:([\d.]+)%")[0].astype(float)
    hub = pd.read_csv(ROOT / "data" / "hub1_refprot.tsv.gz", sep="\t")
    hub_prot = set(hub.Proteomes.dropna().str.split(":").str[0])

    g = (b.groupby("genus")
          .agg(proteome=("proteome", "first"), organism=("organism", "first"), kingdom=("kingdom", "first"),
               phylum=("phylum", "first"), cls=("class", "first"), order=("order", "first")))
    g["group"] = g.kingdom.fillna(g.phylum).fillna("unclassified")
    for sf in TRP_SF + ["DDX3/DED1", "eIF4A"]:
        a = b[b.subfamily == sf].drop_duplicates("genus").set_index("genus").anchor
        g[sf] = a.reindex(g.index).fillna("")
    g = g.join(prot.set_index("Proteome Id")[["busco_c"]], on="proteome")
    g["hub1"] = [("yes" if p in hub_prot else "no" if c >= BUSCO_MIN else "unknown")
                 for p, c in zip(g.proteome, g.busco_c.fillna(0))]
    g.to_csv(ROOT / "results" / "genus_table.tsv", sep="\t")

    print(f"{len(g)} eukaryotic genera; Hub1 yes/no/unknown:", g.hub1.value_counts().to_dict())
    for sf in TRP_SF:
        x = g[g[sf] != ""]
        print(f"\n== {sf}: {len(x)} genera; anchor x Hub1")
        print(pd.crosstab(x[sf], x.hub1, margins=True).to_string())
    ex = g[[any(g.at[i, sf] not in ("", "W") for sf in TRP_SF) for i in g.index]]
    ex.to_csv(ROOT / "results" / "trp_exceptions.tsv", sep="\t")
    for sf in TRP_SF:
        x = g[(g[sf] != "") & (g[sf] != "W")]
        print(f"\n-- non-W {sf} by group/phylum/class ({len(x)}):")
        print(x.groupby(["group", "phylum", "cls", sf, "hub1"], dropna=False).size().to_string())


if __name__ == "__main__":
    main()
