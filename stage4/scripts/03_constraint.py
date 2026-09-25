#!/usr/bin/env python
"""Is the anchor site depleted of missense variation in people, or hit in cancer?

For each site class (UA, UA+1, Q-stacker, Q-4, Q) the observed number of distinct
missense variants (gnomAD v4) or tumour patients (cBioPortal, cell lines excluded)
is compared with the expectation for an average residue of the same gene
(total missense / protein length), summed over genes; one-sided Poisson tests.
Observed and expected both count only variants whose reference residue matches
UniProt at that position, so a transcript-numbering offset affects both equally.

Output
  results/constraint.tsv
"""
from pathlib import Path

import pandas as pd
from scipy.stats import poisson

ROOT = Path(__file__).resolve().parents[1]
SITES = ["UA", "UA+1", "Q-stacker", "Q-4", "Q"]


def test(obs, exp):
    return {"observed": int(obs), "expected": round(exp, 1), "ratio": round(obs / exp, 2),
            "p_depleted": float(f"{poisson.cdf(obs, exp):.2g}"),
            "p_enriched": float(f"{poisson.sf(obs - 1, exp):.2g}")}


def main():
    v = pd.read_csv(ROOT / "results" / "variants_positions.tsv", sep="\t")
    gn = pd.read_csv(ROOT / "results" / "gnomad_missense_counts.tsv", sep="\t")
    g = v[(v.source == "gnomAD") & v.uniprot_match & (v.consequence == "missense_variant")]
    exp_g = (gn.missense_variants_total / gn.length).sum()

    cm = pd.read_csv(ROOT / "results" / "cancer_missense_counts.tsv", sep="\t")
    cp = pd.read_csv(ROOT / "results" / "cancer_positions.tsv", sep="\t")
    cp = cp[cp.uniprot_match & cp.tumour_only & (cp.type == "Missense_Mutation")]
    exp_c = cm.per_residue_mean.sum()

    rows = []
    for s in SITES:
        rows.append({"dataset": "gnomAD v4 missense variants", "site": s, "genes": len(gn),
                     **test((g.site == s).sum(), exp_g)})
    for s in SITES:
        rows.append({"dataset": "cancer missense patients (tumours)", "site": s, "genes": len(cm),
                     **test((cp.site == s).sum(), exp_c)})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "constraint.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
