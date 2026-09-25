#!/usr/bin/env python
"""Read the ancestral states at the anchor column and count switches.

For each subfamily tree (work/<tag>/t.*):
  * root the tree on the outgroup; the ingroup root is the MRCA of the ingroup
  * state of every internal node at the anchor column = IQ-TREE marginal MAP
    state (with its posterior); tips = observed residue
  * a switch = an edge whose parent and child states differ; each switch is
    reported with the clade below it (size, supergroups, lowest shared taxon)
  * site rate: the anchor column's empirical Bayes rate (--rate) and its
    percentile among all columns (0 = slowest)

Output
  results/asr_summary.tsv    one row per subfamily
  results/asr_switches.tsv   one row per switch
"""
import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
TAGS = ["prp5", "prp28", "mak5", "spb4", "rrp3", "dbp10", "dbp9", "dhh1", "dbp5", "drs1", "eif4a", "ded1"]
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus"]


def read_states(path, col):
    st = pd.read_csv(path, sep="\t", comment="#")
    st = st[st.Site == col]
    probs = st.filter(like="p_")
    return {n: (s, float(p.max())) for n, s, (_, p) in zip(st.Node, st.State, probs.iterrows())}


def aln_column(path, col):
    out = {}
    for block in open(path).read().split(">")[1:]:
        name, *s = block.split("\n")
        out[name.strip()] = "".join(s)[col - 1]
    return out


def lowest_shared(meta, tips):
    m = meta.loc[[t for t in tips if t in meta.index]]
    for rank in reversed(RANKS):
        v = m[rank].dropna().unique()
        if len(v) == 1 and m[rank].notna().all():
            return f"{rank}:{v[0]}"
    for rank in RANKS:
        v = m[rank].dropna().unique()
        if len(v) > 1:
            return "mixed " + rank + ":" + ",".join(sorted(v)[:4]) + ("..." if len(v) > 4 else "")
    return "?"


def analyse(tag):
    d = ROOT / "work" / tag
    col = int((d / "anchor_col.txt").read_text())
    meta = pd.read_csv(d / "meta.tsv", sep="\t").set_index("acc")
    tree = Phylo.read(d / "t.treefile", "newick")
    out = [t for t in tree.get_terminals() if meta.at[t.name, "role"] == "out"]
    tree.root_with_outgroup(out)
    ing = [t for t in tree.get_terminals() if meta.at[t.name, "role"] == "in"]
    iroot = tree.common_ancestor(ing)
    anc = read_states(d / "t.state", col)
    tips = aln_column(d / "aln.mask.fa", col)

    def state(clade):
        if clade.is_terminal():
            return tips[clade.name], 1.0
        return anc.get(clade.name, ("?", 0.0))

    sw = []
    for parent in iroot.find_clades(order="preorder"):
        for child in parent.clades:
            a, b = state(parent)[0], state(child)[0]
            if a != b:
                leaves = [t.name for t in child.get_terminals()]
                sub = meta.loc[leaves]
                sw.append({"subfamily": meta.subfamily_final.iloc[0], "from": a, "to": b,
                           "p_parent": round(state(parent)[1], 2), "p_child": round(state(child)[1], 2),
                           "n_tips": len(leaves), "supergroups": ",".join(sorted(sub.supergroup.unique())),
                           "clade": lowest_shared(meta, leaves),
                           "tip_states": "".join(f"{k}{v}" for k, v in sub.anchor.value_counts().items())})
    rates = pd.read_csv(d / "t.rate", sep="\t", comment="#")
    r = rates.loc[rates.Site == col, "Rate"].iloc[0]
    ing_meta = meta[meta.role == "in"]
    summary = {"subfamily": ing_meta.subfamily_final.iloc[0], "n_in": len(ing_meta),
               "tips": "".join(f"{k}{v}" for k, v in ing_meta.anchor.value_counts().items()),
               "root_state": state(iroot)[0], "root_p": round(state(iroot)[1], 3),
               "n_switches": len(sw),
               "to_W": sum(s["to"] == "W" for s in sw), "from_W": sum(s["from"] == "W" for s in sw),
               "anchor_rate": round(r, 3),
               "rate_percentile": round(100 * (rates.Rate < r).mean(), 1),
               "frac_sites_slower_or_equal": round((rates.Rate <= r).mean(), 3)}
    return summary, sw


def main():
    rows, sws = [], []
    for tag in (sys.argv[1:] or TAGS):
        if not (ROOT / "work" / tag / "t.state").exists():
            continue
        s, sw = analyse(tag)
        rows.append(s)
        sws += sw
    pd.DataFrame(rows).to_csv(ROOT / "results" / "asr_summary.tsv", sep="\t", index=False)
    pd.DataFrame(sws).to_csv(ROOT / "results" / "asr_switches.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 70):
        print(pd.DataFrame(rows).to_string(index=False))
        if sws:
            print(pd.DataFrame(sws).to_string(index=False))


if __name__ == "__main__":
    main()
