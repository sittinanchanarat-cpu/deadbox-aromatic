#!/usr/bin/env python
"""Family-level tree: did the Trp anchor arise once or several times before LECA?

Six sequences per eukaryotic subfamily (>= 100 genera; one per kingdom/phylum,
Swiss-Prot first) plus 30 bacterial DEAD-box sequences as the outgroup (one per
bacterial subfamily/phylum). Region anchor-10 .. Q+360. MAFFT --auto, columns
> 50% gaps dropped, IQ-TREE 3 LG+G4 with 1000 ultrafast bootstraps.

Each subfamily is then collapsed to its majority anchor (W/F) and the Trp
subfamilies are checked for monophyly (the smallest clade holding all their
sequences, its bootstrap, and which other subfamilies it contains).

Output
  work/family/*   alignment and IQ-TREE outputs
  results/family_tree_trp.txt
"""
import subprocess
from pathlib import Path

import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT.parent / "stage1"
MAMBA = [str(Path.home() / ".local/bin/micromamba"), "run", "-n", "helicase"]
TRP = ["DDX46/PRP5", "DDX23/PRP28", "DDX24/MAK5", "DDX55/SPB4"]
PER_SF, N_BAC = 6, 30


def main():
    d = ROOT / "work" / "family"
    d.mkdir(parents=True, exist_ok=True)
    seq = pd.read_csv(S1 / "data" / "dead_refprot.tsv.gz", sep="\t", usecols=["Entry", "Sequence"]).set_index("Entry").Sequence
    c = pd.read_csv(S1 / "results" / "census_balanced_assigned.tsv", sep="\t")
    c = c[c.q_ok.fillna(False).astype(bool) & c.ua_found.fillna(False).astype(bool)]
    euk = c[(c.domain == "Eukaryota") & ~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    big = euk.subfamily_final.value_counts()
    big = big[big >= 100].index
    picks = []
    for sf in big:
        g = euk[euk.subfamily_final == sf].sample(frac=1, random_state=3)
        g = g.sort_values("Reviewed", key=lambda s: s.ne("reviewed"), kind="stable")
        g = g.assign(k=g.kingdom.fillna(g.phylum))
        picks.append(pd.concat([g.drop_duplicates("k"), g]).drop_duplicates("acc").head(PER_SF))
    bac = c[c.domain == "Bacteria"].sample(frac=1, random_state=3)
    bac = pd.concat([bac.drop_duplicates("subfamily"), bac.drop_duplicates("phylum"), bac]).drop_duplicates("acc").head(N_BAC)
    s = pd.concat(picks + [bac.assign(subfamily_final="bacterial:" + bac.subfamily)])
    s["label"] = s.subfamily_final.str.replace("/", "_").str.replace(":", "_") + "|" + s.ua_aa + "|" + s.acc
    with open(d / "in.fa", "w") as fh:
        for r in s.itertuples():
            q, u = int(r.q_pos), int(r.ua_pos)
            fh.write(f">{r.label}\n{seq[r.acc][max(0, u - 11):q + 360]}\n")
    with open(d / "aln.fa", "w") as fh:
        subprocess.run(MAMBA + ["mafft", "--auto", "--thread", "8", "--quiet", str(d / "in.fa")], stdout=fh, check=True)
    aln = {}
    for block in open(d / "aln.fa").read().split(">")[1:]:
        name, *x = block.split("\n")
        aln[name.strip()] = "".join(x)
    L = len(next(iter(aln.values())))
    keep = [j for j in range(L) if sum(v[j] == "-" for v in aln.values()) / len(aln) <= 0.5]
    with open(d / "aln.mask.fa", "w") as fh:
        for k, v in aln.items():
            fh.write(f">{k}\n{''.join(v[j] for j in keep)}\n")
    outg = ",".join(l for l in s.label if l.startswith("bacterial"))
    subprocess.run(MAMBA + ["iqtree3", "-s", str(d / "aln.mask.fa"), "-m", "LG+G4", "-B", "1000", "-o", outg,
                            "-T", "8", "-pre", str(d / "t"), "-redo", "-quiet"], check=True)

    tree = Phylo.read(d / "t.treefile", "newick")
    tree.root_with_outgroup([t for t in tree.get_terminals() if t.name.startswith("bacterial")])
    lines = []
    sf_of = lambda n: n.split("|")[0]
    for sf in TRP:
        tips = [t for t in tree.get_terminals() if sf_of(t.name) == sf.replace("/", "_")]
        ca = tree.common_ancestor(tips)
        inside = sorted({sf_of(t.name) for t in ca.get_terminals()})
        lines.append(f"{sf}: {len(tips)} seqs; monophyletic={inside == [sf.replace('/', '_')]}; "
                     f"support={ca.confidence}")
    all_trp = [t for t in tree.get_terminals() if sf_of(t.name) in {x.replace("/", "_") for x in TRP}]
    ca = tree.common_ancestor(all_trp)
    inside = pd.Series([sf_of(t.name) for t in ca.get_terminals()]).value_counts()
    lines.append(f"\nSmallest clade with all Trp-subfamily sequences: {ca.count_terminals()} tips, support={ca.confidence}")
    lines.append(inside.to_string())
    # which subfamily is the sister of each Trp subfamily?
    for sf in TRP:
        tips = [t for t in tree.get_terminals() if sf_of(t.name) == sf.replace("/", "_")]
        node = tree.common_ancestor(tips)
        path = tree.get_path(node)
        parent = path[-2] if len(path) > 1 else tree.root
        sis = sorted({sf_of(t.name) for t in parent.get_terminals()} - {sf.replace("/", "_")})
        lines.append(f"sister of {sf}: {', '.join(sis)} (parent support={parent.confidence})")
    (ROOT / "results" / "family_tree_trp.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
