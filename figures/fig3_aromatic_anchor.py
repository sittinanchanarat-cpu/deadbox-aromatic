#!/usr/bin/env python
"""Figure 3: the anchor is aromatic in every subfamily; which aromatic varies.

A  anchor residue per eukaryotic DEAD-box subfamily (one sequence per genus; subfamilies
   with >= 100 genera), plus all bacterial DEAD-box proteins; right column = % aromatic
   (stage1/results/census_balanced_assigned.tsv)
B  AlphaFold models with adenine transplanted from 2DB3: closest anchor ring atom -> N6
   vs ring centroid -> Q-motif stacker centroid; Phe and Trp anchors make the same contact
   (stage2/results/af_anchor.tsv)
"""
import pandas as pd

from sc_style import *  # noqa: F403

TRP_SF = {"DDX23/PRP28", "DDX46/PRP5", "DDX24/MAK5", "DDX55/SPB4"}


def panel_a(ax):
    c = pd.read_csv(ROOT / "stage1" / "results" / "census_balanced_assigned.tsv", sep="\t", low_memory=False)
    c = c[c.q_ok.fillna(False).astype(bool)]
    c["anchor"] = c.ua_aa.where(c.ua_found.fillna(False).astype(bool), "-").where(lambda s: s.isin(list("FWY-")), "-")
    e = c[(c.domain == "Eukaryota") & ~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    e = e.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    t = pd.crosstab(e.subfamily_final, e.anchor)
    t = t[t.sum(axis=1) >= 100]
    bac = c[c.domain == "Bacteria"].anchor.value_counts()
    t.loc["Bacteria (all)"] = bac.reindex(t.columns).fillna(0)
    frac = t.div(t.sum(axis=1), axis=0) * 100
    frac["arom"] = frac[["F", "W", "Y"]].sum(axis=1)
    order = ["Bacteria (all)"] + list(frac.drop("Bacteria (all)").sort_values(["arom", "W"], ascending=[True, True]).index)
    y = np.arange(len(order))
    left = np.zeros(len(order))
    for aa, lab in (("W", "Trp"), ("F", "Phe"), ("Y", "Tyr"), ("-", "none")):
        v = frac.reindex(order)[aa].fillna(0).values if aa in frac else np.zeros(len(order))
        ax.barh(y, v, 0.78, left=left, color=AA_COLOR[aa], alpha=ALPHA, lw=0, label=lab)
        left += v
    names = [s.replace("/", " / ") for s in order]
    ax.set_yticks(y, names, fontsize=6)
    for i, s in enumerate(order):
        ax.text(101.5, i, f"{frac.loc[s, 'arom']:.0f}", va="center", fontsize=5.5, color=SC_ANNOT)
        ax.text(111, i, f"{int(t.loc[s].sum()):,}", va="center", fontsize=5.5, color=SC_ANNOT)
    ax.text(101.5, len(order) - 0.2, "arom.\n(%)", fontsize=5.5, color=SC_ANNOT, va="bottom")
    ax.text(111, len(order) - 0.2, "genera", fontsize=5.5, color=SC_ANNOT, va="bottom")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("Anchor residue (% of genera)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4)
    return t


def panel_b(ax):
    p = pd.read_csv(ROOT / "stage2" / "results" / "af_anchor.tsv", sep="\t")
    p = p[(p.status == "ok") & p.struct_aa.isin(["W", "F"])]
    for aa, col, lab in (("F", FOREST, "Phe anchor"), ("W", CRIMSON, "Trp anchor")):
        x = p[p.struct_aa == aa]
        ax.scatter(x.to_qarom, x.n6_closest, s=8, color=col, lw=0, alpha=ALPHA, label=f"{lab} (n = {len(x)})")
    ax.axhspan(3.3, 4.1, color=SC_GRID, lw=0, zorder=0)
    ax.text(5.6, 4.13, "PDB range", fontsize=5.5, color=SC_ANNOT, va="bottom")
    ax.set_xlabel("Anchor ring – Q-motif stacker ring (Å)")
    ax.set_ylabel("Closest anchor ring atom – N6 (Å)")
    ax.set_xlim(5.5, 9.0)
    ax.set_ylim(2.9, 4.8)
    ax.legend(loc="upper left", fontsize=6, bbox_to_anchor=(0, 1.02))
    for aa, col in (("F", FOREST), ("W", CRIMSON)):
        m = p[p.struct_aa == aa].n6_closest.median()
        ax.axhline(m, color=col, lw=0.6, ls=":")
    return len(p)


def main():
    fig, (a, b) = plt.subplots(1, 2, figsize=(JOURNAL_2COL, 4.2), gridspec_kw={"width_ratios": [1.15, 1]},
                               constrained_layout=True)
    t = panel_a(a)
    n = panel_b(b)
    add_panel_label(a, "A", x=-0.32)
    add_panel_label(b, "B", x=-0.16)
    save_figure(fig, "fig3_aromatic_anchor")
    print(f"Fig 3: {len(t)} rows; {n} AlphaFold models with an inserted anchor")


if __name__ == "__main__":
    main()
