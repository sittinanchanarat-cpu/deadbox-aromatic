#!/usr/bin/env python
"""Figure 5: human variants at the anchor and its shell.

A  ClinVar missense and in-frame variants at the anchor (UA), the residue after it
   (UA+1), the Q-motif stacker, the Q-motif Pro (Q-4) and the Q-motif Gln, for human
   DEAD-box genes with at least one such variant; the most severe class per cell is
   shown, with the change named for pathogenic/likely pathogenic calls
   (stage4/results/variants_positions.tsv; UniProt-matched positions only)
B  AlphaMissense: mean predicted pathogenicity over the 19 substitutions at each site,
   one dot per gene; grey = the gene's median buried aromatic (F/W/Y, relative SASA < 0.1
   in the AlphaFold model, anchor excluded) (stage4/results/alphamissense_sites.tsv)
C  AlphaMissense at the anchor by substitution: aromatic swaps vs aliphatic/Ala, for Phe
   anchors (31 genes) and Trp anchors (4 genes); dashed lines = AlphaMissense class
   thresholds (0.34 benign, 0.564 pathogenic) (stage4/results/alphamissense_swaps.tsv)
Supplementary Figure S1 (figS1_constraint): observed/expected missense counts per site
class, summed over 37 genes: gnomAD v4 variants and tumour patients (cBioPortal, cell
lines excluded); exact Poisson 95% CIs (stage4/results/constraint.tsv)
"""
import pandas as pd
from scipy.stats import chi2

from sc_style import *  # noqa: F403

S4 = ROOT / "stage4" / "results"
SITES = ["UA", "UA+1", "Q-stacker", "Q-4", "Q"]
SITE_LABEL = {"UA": "Anchor", "UA+1": "UA+1", "Q-stacker": "Q-motif\nstacker", "Q-4": "Q-motif\nPro", "Q": "Q-motif\nGln"}
AA1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Ter": "*"}


def cls(sig):
    s = str(sig).lower()
    if "pathogenic" in s and "conflicting" not in s and "benign" not in s:
        return "pathogenic"
    if "benign" in s and "conflicting" not in s:
        return "benign"
    return "uncertain"


def short(h):
    for k, v in AA1.items():
        h = h.replace(k, v)
    return h.replace("p.", "")


def panel_a(ax):
    v = pd.read_csv(S4 / "variants_positions.tsv", sep="\t")
    v = v[(v.source == "ClinVar") & v.uniprot_match & v.site.isin(SITES)
          & v.consequence.isin(["missense_variant", "inframe_deletion", "inframe_insertion"])]
    v["cls"] = v.clinvar.map(cls)
    rank = {"pathogenic": 0, "uncertain": 1, "benign": 2}
    genes = sorted(v.gene.unique(), key=lambda g: (v[v.gene == g].cls.map(rank).min(), g))
    for i, g in enumerate(genes):
        for j, s in enumerate(SITES):
            x = v[(v.gene == g) & (v.site == s)]
            if not len(x):
                continue
            best = min(x.cls, key=rank.get)
            ax.scatter(j, i, s=26 + 8 * (len(x) - 1), color=CLINVAR_COLOR[best], lw=0, alpha=0.9, zorder=3)
            if best == "pathogenic":
                names = ", ".join(sorted({short(h) for h in x[x.cls == "pathogenic"].hgvsp}))
                ax.text(j, i + 0.3, names.replace(", ", "\n"), fontsize=5, ha="center", va="top", color=CRIMSON)
    anchors = v.drop_duplicates("gene").set_index("gene").anchor
    ax.set_yticks(range(len(genes)), [f"{g} ({anchors[g]})" for g in genes], fontsize=6)
    ax.set_xticks(range(len(SITES)), [SITE_LABEL[s] for s in SITES], fontsize=6)
    ax.set_xlim(-0.5, len(SITES) - 0.2)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    ax.grid(axis="both", color=SC_GRID, lw=0.5, zorder=0)
    ax.tick_params(length=0)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", ms=4, color=CLINVAR_COLOR[k], label=l)
                       for k, l in (("pathogenic", "Pathogenic / likely pathogenic"),
                                    ("uncertain", "Uncertain / conflicting"), ("benign", "Benign / likely benign"))],
              loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=6)


def poisson_ci(k, alpha=0.05):
    lo = 0.0 if k == 0 else chi2.ppf(alpha / 2, 2 * k) / 2
    hi = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    return lo, hi


def panel_b(ax):
    s = pd.read_csv(S4 / "alphamissense_sites.tsv", sep="\t")
    ref = s.drop_duplicates("gene").buried_aromatics_median.dropna()
    ax.axhspan(ref.quantile(0.25), ref.quantile(0.75), color=SC_GRID, lw=0, zorder=0)
    ax.axhline(ref.median(), color=SC_ANNOT, lw=0.8, ls="--", zorder=1)
    ax.text(len(SITES) - 0.45, ref.median(), "buried\naromatics", fontsize=5.5, color=SC_ANNOT, va="center")
    for i, site in enumerate(SITES):
        y = s[s.site == site].am_mean
        col = FOREST if site == "UA" else SAGE if site == "UA+1" else MOSS
        ax.scatter(i + np.random.default_rng(i).uniform(-0.15, 0.15, len(y)), y, s=7, color=col, lw=0, alpha=ALPHA, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [y.median()] * 2, color=SC_AXIS, lw=1.2, zorder=4)
    ax.set_xticks(range(len(SITES)), [SITE_LABEL[x] for x in SITES], fontsize=6)
    ax.set_xlim(-0.5, len(SITES) - 0.1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("AlphaMissense (mean of 19 substitutions)")


def panel_c(ax):
    w = pd.read_csv(S4 / "alphamissense_swaps.tsv", sep="\t")
    cols = [("F", "F->W", "F\u2192W", CRIMSON), ("F", "F->Y", "F\u2192Y", SAND), ("F", "F->L", "F\u2192L", SC_ANNOT),
            ("F", "F->A", "F\u2192A", SC_ANNOT), ("W", "W->F", "W\u2192F", FOREST), ("W", "W->Y", "W\u2192Y", SAND),
            ("W", "W->L", "W\u2192L", SC_ANNOT), ("W", "W->A", "W\u2192A", SC_ANNOT)]
    xs = [0, 1, 2, 3, 4.8, 5.8, 6.8, 7.8]
    for x, (aa, key, lab, col) in zip(xs, cols):
        y = w[w.anchor.str[0] == aa][key].dropna()
        ax.scatter(x + np.random.default_rng(int(x * 10)).uniform(-0.15, 0.15, len(y)), y, s=7, color=col,
                   lw=0, alpha=ALPHA, zorder=3)
        ax.plot([x - 0.28, x + 0.28], [y.median()] * 2, color=SC_AXIS, lw=1.2, zorder=4)
    for t in (0.34, 0.564):
        ax.axhline(t, color=SC_ANNOT, lw=0.5, ls="--", zorder=1)
    ax.set_xticks(xs, [c[2] for c in cols], fontsize=6)
    ax.text(1.5, -0.2, f"Phe anchors (n = {(w.anchor.str[0] == 'F').sum()})", ha="center", fontsize=6,
            transform=ax.get_xaxis_transform())
    ax.text(6.3, -0.2, f"Trp anchors (n = {(w.anchor.str[0] == 'W').sum()})", ha="center", fontsize=6,
            transform=ax.get_xaxis_transform())
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("AlphaMissense pathogenicity")


def constraint_panel(ax):
    c = pd.read_csv(S4 / "constraint.tsv", sep="\t")
    for k, (ds, col, dy, lab) in enumerate((("gnomAD", FOREST, -0.14, "gnomAD v4 variants"),
                                            ("cancer", TERRACOTTA, 0.14, "Tumour patients"))):
        x = c[c.dataset.str.startswith(ds)].set_index("site").reindex(SITES)
        for i, s in enumerate(SITES):
            r = x.loc[s]
            lo, hi = poisson_ci(int(r.observed))
            ax.errorbar(r.ratio, i + dy, xerr=[[r.ratio - lo / r.expected], [hi / r.expected - r.ratio]],
                        fmt="o", ms=3.5, color=col, elinewidth=1.0, capsize=0, label=lab if i == 0 else None)
            p = min(r.p_depleted, r.p_enriched)
            if p < 0.05:
                ax.text(hi / r.expected + 0.05, i + dy, pval_to_stars(p), va="center", fontsize=6, color=col)
    ax.axvline(1, color=SC_ANNOT, lw=0.6, ls="--")
    ax.set_yticks(range(len(SITES)), [SITE_LABEL[s].replace("\n", " ") for s in SITES], fontsize=6)
    ax.set_ylim(len(SITES) - 0.5, -0.5)
    ax.set_xlim(0, 2.2)
    ax.set_xlabel("Observed / expected missense")
    ax.legend(loc="lower right", fontsize=6)


def main():
    fig = plt.figure(figsize=(JOURNAL_2COL, 6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], height_ratios=[1.15, 1])
    a, b, c = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    panel_a(a)
    panel_b(b)
    panel_c(c)
    add_panel_label(a, "A", x=-0.3)
    add_panel_label(b, "B", x=-0.25)
    add_panel_label(c, "C", x=-0.25)
    save_figure(fig, "fig5_variants")
    fig, ax = plt.subplots(figsize=(JOURNAL_1COL, 2.8), constrained_layout=True)
    constraint_panel(ax)
    save_figure(fig, "figS1_constraint")
    print("Fig 5 and Fig S1 written")


if __name__ == "__main__":
    main()
