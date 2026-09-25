#!/usr/bin/env python
"""Gene trees with ancestral states for the anchor, one per subfamily.

Sampling (genus-balanced census, anchor found, Q-motif Gln found):
  * every genus outside Amorphea and Archaeplastida (rare supergroups), capped at 120
  * one genus per order within Amorphea and Archaeplastida
  * one genus per family for the minority anchor state of that subfamily (the
    switch lineages), capped at 60
  total capped at 400; plus 8 outgroup sequences (eIF4A; DDX6 for the eIF4A tree)
Region: anchor-40 .. Q+360 (upstream loop, RecA1, RecA2), so N/C-terminal
extensions do not dominate the alignment.
Alignment: MAFFT --auto; columns with > 50% gaps are dropped except the anchor
column (the modal alignment column of the anchor residues).
Tree: IQ-TREE 3, LG+G4, -fast, rooted on the outgroup, with -asr (marginal
ancestral states) and --rate (per-site empirical Bayes rates).

Output (per subfamily, work/<tag>/)
  in.fa, aln.fa, aln.mask.fa, meta.tsv (acc, supergroup, anchor, lineage),
  anchor_col.txt (1-based column in aln.mask.fa), iqtree outputs (t.*)
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S1 = ROOT.parent / "stage1"
sys.path.insert(0, str(ROOT / "scripts"))
from importlib import import_module  # noqa: E402

supergroup = import_module("01_supergroups").supergroup
MAMBA = [str(Path.home() / ".local/bin/micromamba"), "run", "-n", "helicase"]
SUBFAMILIES = ["DDX46/PRP5", "DDX23/PRP28", "DDX24/MAK5", "DDX55/SPB4",
               "DDX47/RRP3", "DDX54/DBP10", "DDX56/DBP9", "DDX6/DHH1", "DDX19/DBP5", "DDX27/DRS1",
               "eIF4A", "DDX3/DED1"]
UP, DOWN, CAP, N_OUT = 40, 360, 400, 8


def tag(sf):
    return sf.split("/")[-1].lower()


def load():
    raw = pd.read_csv(S1 / "data" / "dead_refprot.tsv.gz", sep="\t",
                      usecols=["Entry", "Sequence", "Taxonomic lineage"]).set_index("Entry")
    c = pd.read_csv(S1 / "results" / "census_balanced_assigned.tsv", sep="\t")
    c = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool) & c.ua_found.fillna(False).astype(bool)
          & ~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    c = c.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    c["seq"] = c.acc.map(raw.Sequence)
    c["supergroup"] = c.acc.map(raw["Taxonomic lineage"]).map(supergroup)
    c["anchor"] = c.ua_aa
    return c


def sample(g):
    g = g.sample(frac=1, random_state=1)
    rare = g[~g.supergroup.isin(["Amorphea", "Archaeplastida"])].head(120)
    common = g[g.supergroup.isin(["Amorphea", "Archaeplastida"])].drop_duplicates("order")
    minority_state = g.anchor.value_counts().index[1] if g.anchor.nunique() > 1 else None
    minority = g[g.anchor == minority_state].drop_duplicates("family").head(60) if minority_state else g.head(0)
    s = pd.concat([minority, rare, common]).drop_duplicates("acc").head(CAP)
    return s


def segment(r):
    start = max(0, int(r.ua_pos) - 1 - UP)
    return r.seq[start:int(r.q_pos) + DOWN], int(r.ua_pos) - 1 - start   # (segment, anchor index in segment)


def run(cmd, **kw):
    subprocess.run(MAMBA + cmd, check=True, **kw)


def build(sf, c):
    d = ROOT / "work" / tag(sf)
    d.mkdir(parents=True, exist_ok=True)
    ing = sample(c[c.subfamily_final == sf])
    out_sf = "DDX6/DHH1" if sf == "eIF4A" else "eIF4A"
    og = c[c.subfamily_final == out_sf].sample(frac=1, random_state=2).drop_duplicates("kingdom").head(N_OUT)
    meta = pd.concat([ing.assign(role="in"), og.assign(role="out")])
    seg = meta.apply(segment, axis=1)
    meta["segment"], meta["ua_idx"] = seg.str[0], seg.str[1]
    with open(d / "in.fa", "w") as fh:
        for r in meta.itertuples():
            fh.write(f">{r.acc}\n{r.segment}\n")
    with open(d / "aln.fa", "w") as fh:
        run(["mafft", "--auto", "--thread", "8", "--quiet", str(d / "in.fa")], stdout=fh)
    aln = {}
    for block in open(d / "aln.fa").read().split(">")[1:]:
        name, *s = block.split("\n")
        aln[name.strip()] = "".join(s).upper()
    # anchor column: modal alignment column of the ingroup anchors
    cols = []
    for r in meta[meta.role == "in"].itertuples():
        s, k = aln[r.acc], -1
        for j, ch in enumerate(s):
            if ch != "-":
                k += 1
                if k == r.ua_idx:
                    cols.append(j)
                    break
    anchor = pd.Series(cols).mode()[0]
    n = len(aln)
    keep = [j for j in range(len(next(iter(aln.values()))))
            if j == anchor or sum(s[j] == "-" for s in aln.values()) / n <= 0.5]
    with open(d / "aln.mask.fa", "w") as fh:
        for name, s in aln.items():
            fh.write(f">{name}\n{''.join(s[j] for j in keep)}\n")
    (d / "anchor_col.txt").write_text(f"{keep.index(anchor) + 1}\n")
    meta[["acc", "role", "subfamily_final", "supergroup", "kingdom", "phylum", "class", "order", "family",
          "genus", "organism", "anchor", "ua_spacing"]].to_csv(d / "meta.tsv", sep="\t", index=False)
    outg = ",".join(meta[meta.role == "out"].acc)
    run(["iqtree3", "-s", str(d / "aln.mask.fa"), "-m", "LG+G4", "-fast", "-asr", "--rate", "-o", outg,
         "-T", "8", "-pre", str(d / "t"), "-redo", "-quiet"])
    print(f"{sf}: {len(ing)} ingroup + {len(og)} outgroup, {len(keep)} columns, anchor column "
          f"{keep.index(anchor) + 1} ({100 * len(cols) / len(ing):.0f}% of anchors located)", flush=True)


def main():
    c = load()
    for sf in (sys.argv[1:] or SUBFAMILIES):
        build(sf, c)


if __name__ == "__main__":
    main()
