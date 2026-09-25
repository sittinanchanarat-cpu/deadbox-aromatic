#!/usr/bin/env python
"""Stage 1 sequence census: upstream-anchor (UA) identity across DEAD-box proteins.

Input  data/dead_refprot.tsv.gz  UniProtKB, InterPro IPR014014 (PROSITE Q_MOTIF
       profile, PRU00552) AND keyword KW-1185 (reference proteome).

UA call. Anchored on the Q-motif Gln (span_end - 1 in the ProRule Q_MOTIF
span, else span_end / span_end - 2). In 11 structurally verified proteins the
anchor sits at Q-24..Q-27, and the profile column (span start + 2) is Q-25 in 98%
of sequences. So
  ua_pos    = F/W at Q-25; otherwise the F/W nearest Q-25 in Q-22..Q-30 (ties ->
              the upstream one); only then a Tyr, by the same order. No aromatic ->
              ua_aa = the residue at Q-25 and ua_found = False.
  ua_spacing = Q - ua_pos (the spacer length; Spb4 W is at Q-27, confirmed in 7NAC)
  short_n   = Q-30 falls before residue 1 (an N-terminal truncation could hide it)

Balancing. One organism per genus (the proteome with most DEAD proteins), then
one sequence per (genus, subfamily): prefer reviewed, q_ok, not short_n, and
the median length of that subfamily.

Output
  results/census_all.tsv       every sequence with its calls
  results/census_balanced.tsv  the genus-balanced set
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AROM = set("FWY")
RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]


def lineage(s):
    out = {}
    for part in str(s).split(", "):
        m = re.match(r"(.+) \((\w+)\)$", part)
        if m and m.group(2) in RANKS and m.group(2) not in out:
            out[m.group(2)] = m.group(1)
    return out


def calls(row):
    m = re.search(r'MOTIF (\d+)\.\.(\d+); /note="Q motif"', str(row.Motif))
    seq = row.Sequence
    if not m:
        return pd.Series({"q_start": np.nan})
    s, e = int(m.group(1)), int(m.group(2))
    q = next((p for p in (e - 1, e, e - 2) if 1 <= p <= len(seq) and seq[p - 1] == "Q"), None)
    out = {"q_start": s, "q_end": e, "q_pos": q, "q_ok": q is not None}
    if q is None:
        return pd.Series(out)
    # F/W anywhere in the window beats Y: AlphaFold models put Rok1's anchor at
    # the Phe at Q-29, not the Tyr at Q-22 (stage2/results/af_anchor.tsv)
    order = (25, 26, 24, 27, 23, 28, 22, 29, 30)
    ua = next((q - d for aa in ("FW", "Y") for d in order if q - d >= 1 and seq[q - d - 1] in aa), None)
    found = ua is not None
    if not found:
        ua = q - 25
    aa = seq[ua - 1] if ua >= 1 else ""
    out.update({"ua_pos": ua, "ua_aa": aa, "ua_found": found, "ua_spacing": q - ua,
                "ua_context": seq[max(0, ua - 5):max(0, ua - 1)] + "[" + aa + "]" + seq[ua:ua + 4] if ua >= 1 else "",
                "short_n": q - 30 < 1})
    return pd.Series(out)


def main():
    d = pd.read_csv(ROOT / "data" / "dead_refprot.tsv.gz", sep="\t")
    d = d.rename(columns={"Entry": "acc", "Organism": "organism", "Organism (ID)": "taxid",
                          "Proteomes": "proteome", "Protein families": "families", "Length": "length"})
    d["subfamily"] = (d.families.fillna("").str.extract(r"DEAD box helicase family, ([^;,]+ subfamily)")[0]
                      .str.replace(" subfamily", "").fillna("unassigned"))
    d = pd.concat([d, d.apply(calls, axis=1)], axis=1)
    tax = d["Taxonomic lineage"].map(lineage).apply(pd.Series)
    d = pd.concat([d, tax.reindex(columns=RANKS)], axis=1)
    d["genus"] = d.genus.fillna(d.organism.str.split().str[0])
    d["proteome"] = d.proteome.fillna("").str.split(":").str[0]
    keep = ["acc", "Reviewed", "Gene Names", "Protein names", "organism", "taxid", "proteome", "length",
            "subfamily", "q_start", "q_end", "q_pos", "q_ok", "ua_pos", "ua_aa", "ua_found", "ua_spacing",
            "ua_context", "short_n"] + RANKS
    d = d[keep]
    d.to_csv(ROOT / "results" / "census_all.tsv", sep="\t", index=False)

    # one proteome per genus: the one with most DEAD proteins
    per_prot = d.groupby(["genus", "proteome"]).size().rename("n").reset_index()
    best = per_prot.sort_values("n", ascending=False).drop_duplicates("genus")
    b = d.merge(best[["genus", "proteome"]], on=["genus", "proteome"])
    med = b.groupby("subfamily").length.transform("median")
    b = b.assign(_rev=b.Reviewed.eq("reviewed"), _dev=(b.length - med).abs())
    b = b.sort_values(["_rev", "q_ok", "short_n", "_dev"], ascending=[False, False, True, True])
    # 'unassigned' sequences are kept individually (they are not one subfamily)
    assigned = b[b.subfamily != "unassigned"].drop_duplicates(["genus", "subfamily"])
    b = pd.concat([assigned, b[b.subfamily == "unassigned"]]).drop(columns=["_rev", "_dev"])
    b.to_csv(ROOT / "results" / "census_balanced.tsv", sep="\t", index=False)
    print(f"all {len(d)} seqs, {d.proteome.nunique()} proteomes; balanced {len(b)} seqs, {b.genus.nunique()} genera")


if __name__ == "__main__":
    main()
