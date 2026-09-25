#!/usr/bin/env python
"""Figure 3: Trp is a fixed trait of four subfamilies and makes the same N6 contact as Phe.

A  anchor residue per eukaryotic DEAD-box subfamily (one sequence per genus; subfamilies
   with >= 100 genera), plus all bacterial DEAD-box proteins
   (stage1/results/census_balanced_assigned.tsv)
B  AlphaFold models with adenine transplanted from 2DB3: closest anchor ring atom -> N6
   vs ring centroid -> Q-motif stacker centroid, Trp vs Phe anchors
   (stage2/results/af_anchor.tsv)
"""
import pandas as pd

from sc_style import *  # noqa: F403

TRP_ROLE = {"DDX23/PRP28": "splicing", "DDX46/PRP5": "splicing", "DDX24/MAK5": "60S assembly",
            "DDX55/SPB4": "60S assembly"}


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
    order = ["Bacteria (all)"] + list(frac.drop("Bacteria (all)").sort_values(["W", "F"], ascending=[True, False]).index)
    y = np.arange(len(order))
    left = np.zeros(len(order))
    for aa, lab in (("W", "Trp"), ("F", "Phe"), ("Y", "Tyr"), ("-", "none")):
        v = frac.reindex(order)[aa].fillna(0).values if aa in frac else np.zeros(len(order))
        ax.barh(y, v, 0.78, left=left, color=AA_COLOR[aa], alpha=ALPHA, lw=0, label=lab)
        left += v
    names = [s.replace("/", " / ") for s in order]
    ax.set_yticks(y, names, fontsize=6)
    for tick, s in zip(ax.get_yticklabels(), order):
        if s in TRP_ROLE:
            tick.set_color(CRIMSON)
            tick.set_fontweight("bold")
    for i, s in enumerate(order):
        n = int(t.loc[s].sum())
        ax.text(101, i, f"{n:,}", va="center", fontsize=5.5, color=SC_ANNOT)
    ax.text(101, len(order) - 0.2, "genera", fontsize=5.5, color=SC_ANNOT, va="bottom")
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
    trp_sf = p.subfamily.isin(TRP_ROLE)
    groups = [("F", ~trp_sf, FOREST, "Phe, Phe subfamilies"), ("F", trp_sf, SAND, "Phe, Trp subfamilies"),
              ("W", trp_sf, CRIMSON, "Trp, Trp subfamilies"), ("W", ~trp_sf, TERRACOTTA, "Trp, other subfamilies")]
    for aa, m, col, lab in groups:
        x = p[m & (p.struct_aa == aa)]
        ax.scatter(x.to_qarom, x.n6_closest, s=8, color=col, lw=0, alpha=ALPHA, label=f"{lab} (n = {len(x)})")
    ax.axhspan(3.3, 4.1, color=SC_GRID, lw=0, zorder=0)
    ax.text(5.6, 4.13, "PDB range", fontsize=5.5, color=SC_ANNOT, va="bottom")
    ax.set_xlabel("Anchor ring – Q-motif stacker ring (Å)")
    ax.set_ylabel("Closest anchor ring atom – N6 (Å)")
    ax.set_xlim(5.5, 9.0)
    ax.set_ylim(2.9, 4.8)
    ax.legend(loc="upper left", fontsize=6, bbox_to_anchor=(0, 1.02))
    return len(p)


def main():
    fig, (a, b) = plt.subplots(1, 2, figsize=(JOURNAL_2COL, 4.2), gridspec_kw={"width_ratios": [1.15, 1]},
                               constrained_layout=True)
    t = panel_a(a)
    n = panel_b(b)
    add_panel_label(a, "A", x=-0.32)
    add_panel_label(b, "B", x=-0.16)
    save_figure(fig, "fig3_trp_census")
    print(f"Fig 3: {len(t)} rows; {n} AlphaFold models with an inserted anchor")


if __name__ == "__main__":
    main()
