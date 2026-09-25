#!/usr/bin/env python
"""Stage 2, step 1: pick orthologs for the AlphaFold anchor check and fetch their models.

From stage1/results/census_balanced_assigned.tsv (eukaryotes, Q-motif Gln found):
  Trp subfamilies (Prp5, Prp28, Mak5, Spb4)  up to 30 each: one per taxonomic class
      first (reviewed entries preferred), then filled from other orders; every
      non-W variant is kept up to 10 per subfamily (Prp5-F, Spb4-F in plants)
  Rok1 (anchor ambiguous: Y at Q-22 vs F at Q-29)   up to 25
  Phe controls (eIF4A, DDX3, DDX6, DDX5, Rrp3, DDX54)  12 each, plus up to 8 of the
      lineage-specific W variants of Rrp3 (insects) and DDX54 (nematodes)
  named yeast/human references (REFS), whether or not they were sampled.

Models: AlphaFold DB, current version, fetched by accession via the API.

Output
  data/selection.tsv   acc, subfamily, organism, class, census anchor call
  data/af/AF-<acc>.cif.gz
"""
import gzip
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT.parent / "stage1" / "results" / "census_balanced_assigned.tsv"
TRP = ["DDX46/PRP5", "DDX23/PRP28", "DDX24/MAK5", "DDX55/SPB4"]
CONTROLS = ["eIF4A", "DDX3/DED1", "DDX6/DHH1", "DDX5/DBP2", "DDX47/RRP3", "DDX54/DBP10"]
REFS = {  # accession: (subfamily, label)
    "P21372": ("DDX46/PRP5", "yeast Prp5"), "Q7L014": ("DDX46/PRP5", "human DDX46"),
    "P23394": ("DDX23/PRP28", "yeast Prp28"), "Q9BUQ8": ("DDX23/PRP28", "human DDX23"),
    "P38112": ("DDX24/MAK5", "yeast Mak5"), "Q9GZR7": ("DDX24/MAK5", "human DDX24"),
    "P25808": ("DDX55/SPB4", "yeast Spb4"), "Q8NHQ9": ("DDX55/SPB4", "human DDX55"),
    "P45818": ("DDX52/ROK1", "yeast Rok1"), "Q9Y2R4": ("DDX52/ROK1", "human DDX52"),
    "P10081": ("eIF4A", "yeast eIF4A"), "P60842": ("eIF4A", "human eIF4A1"),
    "O00571": ("DDX3/DED1", "human DDX3X"), "P26196": ("DDX6/DHH1", "human DDX6"),
}


def pick(g, n, seed=0):
    g = g.sample(frac=1, random_state=seed)
    g = g.sort_values("Reviewed", key=lambda s: s.ne("reviewed"), kind="stable")
    first = g.drop_duplicates("class")
    rest = g.drop(first.index).drop_duplicates("order")
    return pd.concat([first, rest]).head(n)


def fetch(acc):
    out = ROOT / "data" / "af" / f"AF-{acc}.cif.gz"
    if out.exists():
        return "cached"
    try:
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60) as r:
            meta = json.load(r)
        url = next((m["cifUrl"] for m in meta if m.get("uniprotAccession", acc) == acc), meta[0]["cifUrl"])
        with urllib.request.urlopen(url, timeout=120) as r:
            out.write_bytes(gzip.compress(r.read()))
        return "ok"
    except Exception as ex:  # no model (e.g. too long, or not in AFDB)
        return f"fail:{type(ex).__name__}"


def main():
    c = pd.read_csv(CENSUS, sep="\t")
    c = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool) & ~c.short_n.fillna(True).astype(bool)]
    c["anchor"] = c.ua_aa.where(c.ua_found.fillna(False).astype(bool), "-")
    sel = []
    for sf in TRP:
        g = c[c.subfamily_final == sf]
        sel += [pick(g[g.anchor != "W"], 10), pick(g[g.anchor == "W"], 30)]
    sel.append(pick(c[c.subfamily_final == "DDX52/ROK1"], 25))
    for sf in CONTROLS:
        g = c[c.subfamily_final == sf]
        sel.append(pick(g[g.anchor == "F"], 12))
        if sf in ("DDX47/RRP3", "DDX54/DBP10"):
            sel.append(pick(g[g.anchor == "W"], 8))
    s = pd.concat(sel).drop_duplicates("acc")
    s = s[["acc", "subfamily_final", "organism", "kingdom", "phylum", "class", "order", "anchor",
           "ua_pos", "ua_spacing", "q_pos", "ua_context"]].rename(columns={"subfamily_final": "subfamily"})
    s["label"] = ""
    extra = pd.DataFrame([{"acc": a, "subfamily": sf, "label": lab} for a, (sf, lab) in REFS.items()
                          if a not in set(s.acc)])
    s = pd.concat([s, extra], ignore_index=True)
    s.loc[s.acc.isin(REFS), "label"] = s.acc.map(lambda a: REFS.get(a, ("", ""))[1])
    status = []
    for i, a in enumerate(s.acc):
        status.append(fetch(a))
        if i % 25 == 0:
            print(i, a, status[-1], file=sys.stderr)
    s["af"] = status
    s.to_csv(ROOT / "data" / "selection.tsv", sep="\t", index=False)
    print(s.groupby(["subfamily", "anchor"], dropna=False).size().to_string())
    print(s.af.str.split(":").str[0].value_counts().to_dict())


if __name__ == "__main__":
    main()
