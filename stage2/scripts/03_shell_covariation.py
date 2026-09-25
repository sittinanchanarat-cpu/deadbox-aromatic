#!/usr/bin/env python
"""Does the anchor's packing shell change when the anchor is Trp?

Shell positions come from the AlphaFold models (residues within 4.5 A of the
anchor side chain in >= 45% of 237 models): Q-7 (Q-motif stacker), Q-6, Q-5,
Q-4 (Q-motif Pro), Q-12, Q+4, Q+29, and the anchor's own loop UA+3, UA+5
(written relative to the anchor, because the anchor-Q spacing varies).

Two comparisons on the genus-balanced census (stage1, one seq per genus and
subfamily, anchor found):
  between   Trp subfamilies (Prp5, Prp28, Mak5, Spb4) vs all Phe subfamilies
  within    each subfamily that switched: W-anchored vs F-anchored members
            (the switch lineages: Rrp3 insects, DDX54 nematodes, DDX56 budding
            yeasts, DDX6 Sordariomycetes, Spb4 plants, Prp5 F genera)
For each position: residue frequencies in each group, and Fisher's exact test on
the most common residue of the F group (present vs absent) between groups.

Output
  results/shell_covariation.tsv
"""
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT.parent / "stage1"
Q_POS = {"Q-7": -7, "Q-6": -6, "Q-5": -5, "Q-4": -4, "Q-12": -12, "Q+4": 4, "Q+29": 29}
UA_POS = {"UA+1": 1, "UA+3": 3, "UA+5": 5}
TRP_SF = ["DDX46/PRP5", "DDX23/PRP28", "DDX24/MAK5", "DDX55/SPB4"]
SWITCH_SF = ["DDX47/RRP3", "DDX54/DBP10", "DDX56/DBP9", "DDX6/DHH1", "DDX55/SPB4", "DDX46/PRP5"]


def residues(c, seqs):
    out = {}
    for name, d in Q_POS.items():
        out[name] = [s[int(q) - 1 + d] if 0 <= int(q) - 1 + d < len(s) else "-" for s, q in zip(seqs, c.q_pos)]
    for name, d in UA_POS.items():
        out[name] = [s[int(u) - 1 + d] if 0 <= int(u) - 1 + d < len(s) else "-" for s, u in zip(seqs, c.ua_pos)]
    return pd.DataFrame(out, index=c.index)


def compare(a, b, label_a, label_b, context):
    rows = []
    for pos in a.columns:
        fa, fb = a[pos].value_counts(normalize=True), b[pos].value_counts(normalize=True)
        top = fb.index[0]  # dominant residue in the reference (F) group
        na, nb = (a[pos] == top).sum(), (b[pos] == top).sum()
        _, p = fisher_exact([[na, len(a) - na], [nb, len(b) - nb]])
        fmt = lambda f: " ".join(f"{k}{100 * v:.0f}" for k, v in f.head(4).items())
        rows.append({"context": context, "pos": pos, "group_a": label_a, "n_a": len(a), "top_a": fmt(fa),
                     "group_b": label_b, "n_b": len(b), "top_b": fmt(fb), "ref_residue": top,
                     "frac_ref_a": round(na / len(a), 3), "frac_ref_b": round(nb / len(b), 3), "p": p})
    return rows


def main():
    c = pd.read_csv(S1 / "results" / "census_balanced_assigned.tsv", sep="\t")
    c = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool) & c.ua_found.fillna(False).astype(bool)]
    c = c[~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    c = c.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    seq = pd.read_csv(S1 / "data" / "dead_refprot.tsv.gz", sep="\t", usecols=["Entry", "Sequence"]).set_index("Entry").Sequence
    r = residues(c, c.acc.map(seq))
    rows = []
    trp = c.subfamily_final.isin(TRP_SF) & c.ua_aa.eq("W")
    phe = ~c.subfamily_final.isin(TRP_SF) & c.ua_aa.eq("F")
    rows += compare(r[trp], r[phe], "Trp subfamilies (W)", "Phe subfamilies (F)", "between")
    for sf in SWITCH_SF:
        m = c.subfamily_final == sf
        w, f = m & c.ua_aa.eq("W"), m & c.ua_aa.eq("F")
        if w.sum() >= 5 and f.sum() >= 5:
            rows += compare(r[w], r[f], f"{sf} W", f"{sf} F", "within")
    out = pd.DataFrame(rows)
    out["p"] = out.p.map(lambda x: float(f"{x:.2g}"))
    out.to_csv(ROOT / "results" / "shell_covariation.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 40):
        print(out.drop(columns=["context", "ref_residue"]).to_string())


if __name__ == "__main__":
    main()
