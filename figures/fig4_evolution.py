#!/usr/bin/env python
"""Figure 4: Trp arose independently several times before LECA, then was kept.

A  % Trp anchor per subfamily across eukaryotic supergroups (one sequence per genus;
   cells with < 3 genera left blank) (stage3/results/supergroups.tsv)
B  anchor switches reconstructed on per-subfamily gene trees (IQ-TREE marginal ASR):
   Trp losses (W->F) and gains (F->W). Switches on the edges next to the outgroup
   root (clades holding > 50% of the ingroup) are rooting artefacts and are left out
   (stage3/results/asr_switches.tsv, asr_summary.tsv)
C  family-level tree (LG+G4, 1000 ultrafast bootstraps) reduced to one sequence per
   subfamily, midpoint-rooted for display (bacterial sequences fall among eukaryotic
   subfamilies, so they are not a clean outgroup and are not shown); tips coloured by
   the subfamily's anchor, UFBoot >= 90 marked (stage3/work/family/t.treefile)
"""
import pandas as pd
from Bio import Phylo

from sc_style import *  # noqa: F403

S3 = ROOT / "stage3"
GROUPS = ["Amorphea", "Archaeplastida", "SAR", "Haptista", "Discoba", "Metamonada"]  # Cryptophyceae: < 3 genera
SHOW = ["DDX23/PRP28", "DDX24/MAK5", "DDX46/PRP5", "DDX55/SPB4", "DDX47/RRP3", "DDX54/DBP10",
        "DDX56/DBP9", "DDX6/DHH1", "DDX19/DBP5", "DDX27/DRS1", "eIF4A", "DDX3/DED1"]
TRP_SF = {"DDX23_PRP28", "DDX24_MAK5", "DDX46_PRP5", "DDX55_SPB4"}


def panel_a(ax):
    d = pd.read_csv(S3 / "results" / "supergroups.tsv", sep="\t")
    m = np.full((len(SHOW), len(GROUPS)), np.nan)
    for i, sf in enumerate(SHOW):
        for j, g in enumerate(GROUPS):
            r = d[(d.subfamily == sf) & (d.supergroup == g)]
            if len(r) and r.n.iloc[0] >= 3:
                m[i, j] = r.pct_W.iloc[0]
                ax.text(j, i, f"{r.pct_W.iloc[0]:.0f}", ha="center", va="center", fontsize=5,
                        color="white" if r.pct_W.iloc[0] > 55 else SC_AXIS)
    im = ax.imshow(np.ma.masked_invalid(m), cmap=SC_WARM_CMAP, vmin=0, vmax=100, aspect="auto")
    ax.set_facecolor("#F7F7F7")
    ax.set_xticks(range(len(GROUPS)), GROUPS, rotation=40, ha="right", fontsize=6)
    ax.set_yticks(range(len(SHOW)), [s.replace("/", " / ") for s in SHOW], fontsize=6)
    for t, s in zip(ax.get_yticklabels(), SHOW):
        if s.replace("/", "_") in TRP_SF:
            t.set_color(CRIMSON)
            t.set_fontweight("bold")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Trp anchor (% of genera)", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_visible(False)


def panel_b(ax):
    sw = pd.read_csv(S3 / "results" / "asr_switches.tsv", sep="\t")
    sm = pd.read_csv(S3 / "results" / "asr_summary.tsv", sep="\t").set_index("subfamily")
    sw = sw[sw["from"].isin(["W", "F"]) & sw["to"].isin(["W", "F"])]
    sw = sw[sw.n_tips <= 0.5 * sw.subfamily.map(sm.n_in)]
    y = np.arange(len(SHOW))
    loss = [((sw.subfamily == s) & (sw["from"] == "W")).sum() for s in SHOW]
    gain = [((sw.subfamily == s) & (sw["to"] == "W")).sum() for s in SHOW]
    ax.barh(y, [-v for v in loss], 0.7, color=FOREST, alpha=ALPHA, lw=0, label="Trp → Phe (loss)")
    ax.barh(y, gain, 0.7, color=CRIMSON, alpha=ALPHA, lw=0, label="Phe → Trp (gain)")
    for i, s in enumerate(SHOW):
        ax.text(21, i, f"{int(sm.loc[s, 'n_in'])}", va="center", fontsize=5.5, color=SC_ANNOT)
    ax.text(21, -1.1, "seqs", fontsize=5.5, color=SC_ANNOT, va="center")
    ax.axvline(0, color=SC_AXIS, lw=0.6)
    ax.set_yticks(y, [x.split("/")[-1] for x in SHOW], fontsize=5.5)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(SHOW) - 0.5, -1.5)
    ax.set_xlim(-22, 20)
    ax.set_xticks([-20, -10, 0, 10, 20], ["20", "10", "0", "10", "20"])
    ax.set_xlabel("Anchor switches on the gene tree")
    ax.legend(loc="lower right", fontsize=6, bbox_to_anchor=(1.0, 1.0), ncol=1)


def panel_c(ax):
    tree = Phylo.read(S3 / "work" / "family" / "t.treefile", "newick")
    keep, seen = [], set()
    for t in tree.get_terminals():
        sf = t.name.split("|")[0]
        if not sf.startswith("bacterial") and sf not in seen:
            seen.add(sf)
            keep.append(t)
            t.name = sf.replace("_", " / ")
    for t in [t for t in tree.get_terminals() if t not in keep]:
        tree.prune(t)
    tree.root_at_midpoint()
    tree.ladderize()
    anchor = {}
    for t in tree.get_terminals():
        key = t.name.replace(" / ", "_")
        anchor[t.name] = "W" if key in TRP_SF else "F"
    for c in tree.find_clades():
        c.branch_length = c.branch_length or 0.0
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
               label_func=lambda c: c.name if c.is_terminal() else "",
               label_colors=lambda n: CRIMSON if anchor.get(n) == "W" else SC_AXIS)
    for t in ax.texts:
        t.set_fontsize(5.5)
        if anchor.get(t.get_text().strip()) == "W":
            t.set_fontweight("bold")
    depths = tree.depths()
    ys = {}
    for i, t in enumerate(tree.get_terminals()):
        ys[t] = i + 1
    for c in tree.get_nonterminals(order="postorder"):
        ys[c] = np.mean([ys[k] for k in c.clades])
        if c.confidence is not None and c.confidence >= 90 and c is not tree.root:
            ax.plot(depths[c], ys[c], "o", ms=2.2, color=SC_ANNOT)
    for t in tree.get_terminals():
        if anchor.get(t.name) == "W":
            ax.plot(depths[t], ys[t], "s", ms=3, color=CRIMSON)
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Substitutions per site")


def main():
    fig = plt.figure(figsize=(JOURNAL_2COL, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 0.8], height_ratios=[1, 1.25])
    a, b, c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    panel_a(a)
    panel_b(b)
    panel_c(c)
    add_panel_label(a, "A", x=-0.3)
    add_panel_label(b, "B", x=-0.08)
    add_panel_label(c, "C", x=-0.02, y=1.02)
    save_figure(fig, "fig4_evolution")
    print("Fig 4 written")


if __name__ == "__main__":
    main()
