#!/usr/bin/env python
"""Figure 1 (data panels): the anchor defines DEAD-box helicases; the Q motif does not.

A  per SF2 family (PDB, nucleotide-bound): % proteins with a Q-motif Gln on adenine
   vs % with an inserted upstream anchor (stage0/results/per_protein.tsv)
B  anchor -> Q-motif Gln spacing across eukaryotic DEAD-box proteins (log scale), one sequence
   per genus and subfamily (stage1/results/census_balanced_assigned.tsv)
C-E  the adenine pocket in the same orientation (superposed on Vasa RecA1): Vasa-AMPPNP
   (2DB3, Phe anchor), Prp5-ADP (4LJY, Trp anchor), Mtr4-ADP (2XGJ, Ski2-like, no anchor)
   (renderings from structures.py)
"""
import pandas as pd
from scipy.stats import fisher_exact

from sc_style import *  # noqa: F403

FAMILIES = [("DEAD", ["DEAD"]), ("RLR/Dicer", ["RLR/Dicer"]), ("RecQ", ["RecQ"]), ("Ski2-like", ["Ski2-like"]),
            ("SNF2", ["SWI2/SNF2"]), ("DEAH", ["DEAH"]), ("other SF2", ["other", "Rad3/XPD", "viral"])]


def panel_a(ax):
    p = pd.read_csv(ROOT / "stage0" / "results" / "per_protein.tsv", sep="\t")
    rows = []
    for name, fams in FAMILIES:
        g = p[p.family.isin(fams)]
        ev = g[g.ua_evaluable.astype(bool)]
        rows.append((name, len(g), (g.frac_q >= 0.5).mean() * 100, len(ev),
                     (ev.ua_inserted > 0).mean() * 100 if len(ev) else np.nan, (ev.ua_inserted > 0).sum()))
    x = np.arange(len(rows))
    w = 0.38
    for i, (name, n, q, nev, ua, nua) in enumerate(rows):
        col = TERRACOTTA if name == "DEAD" else SAGE
        ax.bar(i - w / 2, q, w, color="white", edgecolor=col, lw=0.9)
        if nev:
            ax.bar(i + w / 2, ua, w, color=col, alpha=ALPHA, lw=0)
            ax.text(i + w / 2, ua + 2, f"{nua}/{nev}", ha="center", va="bottom", fontsize=5.5, color=SC_ANNOT)
        else:
            ax.text(i + w / 2, 2, "n.e.", ha="center", va="bottom", fontsize=5.5, color=SC_ANNOT)
    ax.set_xticks(x, [r[0] for r in rows], rotation=35, ha="right")
    ax.set_ylabel("Proteins (%)")
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    handles = [mpl.patches.Patch(facecolor="white", edgecolor=SC_AXIS, lw=0.9, label="Q-motif Gln on adenine"),
               mpl.patches.Patch(facecolor=SC_AXIS, alpha=ALPHA, label="Inserted anchor")]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.08))
    dead = p[p.family == "DEAD"]
    rest = p[(p.family != "DEAD") & (p.family != "?") & p.ua_evaluable.astype(bool)]
    dead = dead[dead.ua_evaluable.astype(bool)]
    _, pv = fisher_exact([[(dead.ua_inserted > 0).sum(), (dead.ua_inserted == 0).sum()],
                          [(rest.ua_inserted > 0).sum(), (rest.ua_inserted == 0).sum()]])
    return pv


def panel_b(ax):
    c = pd.read_csv(ROOT / "stage1" / "results" / "census_balanced_assigned.tsv", sep="\t", low_memory=False)
    c = c[(c.domain == "Eukaryota") & c.q_ok.fillna(False).astype(bool) & c.ua_found.fillna(False).astype(bool)
          & ~c.subfamily_final.isin(["ambiguous", "no_hit", "unassigned"])]
    c = c.sort_values("assign_how", key=lambda s: s.ne("label")).drop_duplicates(["genus", "subfamily_final"])
    spacings = range(22, 31)
    bottom = np.zeros(len(spacings))
    for aa in ["F", "W", "Y"]:
        n = np.array([((c.ua_spacing == s) & (c.ua_aa == aa)).sum() for s in spacings])
        ax.bar(list(spacings), n, 0.75, bottom=bottom, color=AA_COLOR[aa], alpha=ALPHA, lw=0,
               label={"F": "Phe", "W": "Trp", "Y": "Tyr"}[aa])
        bottom += n
    ax.set_yscale("log")
    ax.set_ylim(1, bottom.max() * 3)
    ax.set_xticks(list(spacings))
    ax.set_xlabel("Anchor position (residues before the Q-motif Gln)")
    ax.set_ylabel("Sequences")
    ax.legend(loc="upper right")
    total = bottom.sum()
    ax.text(25, bottom[3] * 1.4, f"{100 * bottom[3] / total:.0f}%", ha="center", fontsize=6, color=SC_ANNOT)
    return int(total)


def main():
    fig = plt.figure(figsize=(JOURNAL_2COL, 5.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.05])
    a, b = fig.add_subplot(gs[0, :3]), fig.add_subplot(gs[0, 3:])
    c, d, e = fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:])
    pv = panel_a(a)
    n = panel_b(b)
    show_structure(c, "struct_vasa_phe", "Vasa (DEAD), Phe anchor")
    show_structure(d, "struct_prp5_trp", "Prp5 (DEAD), Trp anchor")
    show_structure(e, "struct_mtr4_none", "Mtr4 (Ski2-like), no anchor")
    add_panel_label(a, "A", x=-0.16)
    add_panel_label(b, "B", x=-0.16)
    for ax, lab in zip((c, d, e), "CDE"):
        add_panel_label(ax, lab, x=0.0, y=1.02)
    save_figure(fig, "fig1_definition")
    print(f"Fig 1: DEAD vs other SF2 anchor, Fisher p = {pv:.1e}; spacing panel n = {n}")


if __name__ == "__main__":
    main()
