#!/usr/bin/env python
"""Give unlabelled eukaryotic DEAD-box sequences a subfamily by nearest labelled homolog.

Input  work/hits.m8   mmseqs easy-search of the unlabelled balanced sequences
                      (query) against the labelled ones (target); columns
                      query,target,pident,bits,qcov,tcov
Rule   the subfamily of the best-scoring target, if its bit score beats the best
       hit from any other subfamily by >= 10% and the hit covers >= 50% of the
       query; otherwise 'ambiguous'. No hit -> 'no_hit'.

Output results/census_balanced_assigned.tsv   census_balanced + subfamily_final,
                                              assign_how (label / homolog / ambiguous / no_hit)
       printed: the per-subfamily anchor table redone with the assigned sequences,
                deduplicated again to one sequence per (genus, subfamily)
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MARGIN, MIN_QCOV = 1.10, 0.5


def main():
    c = pd.read_csv(ROOT / "results" / "census_balanced.tsv", sep="\t")
    lab = c.set_index("acc").subfamily
    h = pd.read_csv(ROOT / "work" / "hits.m8", sep="\t", names=["q", "t", "pident", "bits", "qcov", "tcov"])
    h = h[h.qcov >= MIN_QCOV]
    h["sf"] = h.t.map(lab)
    best = h.sort_values("bits", ascending=False).drop_duplicates(["q", "sf"])
    rows = {}
    for q, g in best.groupby("q", sort=False):
        g = g.sort_values("bits", ascending=False)
        top = g.iloc[0]
        second = g.bits.iloc[1] if len(g) > 1 else 0
        rows[q] = (top.sf, "homolog") if top.bits >= MARGIN * second else ("ambiguous", "ambiguous")
    c["subfamily_final"] = c.subfamily
    c["assign_how"] = "label"
    un = c.subfamily == "unassigned"
    c.loc[un, "subfamily_final"] = [rows.get(a, ("no_hit",))[0] for a in c.loc[un, "acc"]]
    c.loc[un, "assign_how"] = [rows.get(a, (None, "no_hit"))[1] for a in c.loc[un, "acc"]]
    c.to_csv(ROOT / "results" / "census_balanced_assigned.tsv", sep="\t", index=False)
    e = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool)]
    print("unlabelled eukaryotic:", e[e.subfamily == "unassigned"].assign_how.value_counts().to_dict())

    e = e[~e.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    e = e.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    e["anchor"] = e.ua_aa.where(e.ua_found.fillna(False).astype(bool), "-").where(lambda s: s.isin(list("FWY-")), "?")
    t = pd.crosstab(e.subfamily_final, e.anchor)
    t["n"] = t.sum(axis=1)
    t["%W"] = (100 * t.get("W", 0) / t.n).round(1)
    t["genera_added"] = e[e.assign_how == "homolog"].groupby("subfamily_final").size().reindex(t.index).fillna(0).astype(int)
    print(t[t.n >= 30].sort_values("%W", ascending=False).to_string())
    left = c[(c.domain == "Eukaryota") & c.subfamily_final.isin(["ambiguous", "no_hit"])]
    print("\nstill unassigned:", len(left), "anchor:", left.ua_aa.value_counts().head(5).to_dict())


if __name__ == "__main__":
    main()
