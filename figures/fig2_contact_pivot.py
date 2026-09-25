#!/usr/bin/env python
"""Figure 2: the anchor reads adenine N6 edge-on, is never touched by partners, and stays put.

A  N6 -> ring centroid distance vs approach angle (angle between ring normal and the
   centroid->N6 vector; 90 deg = N6 in the ring plane) for 119 inserted-anchor sites
   (stage0/results/n6_geometry.tsv)
B  closest anchor ring atom -> N6 and Q-motif Gln -> N6 distances by nucleotide state
C  partner contacts along the anchor loop: fraction of partner-bound chains in which
   each position is within 4.5 A of a partner (deposited coordinates only; symmetry
   mates excluded), and mean buried area (stage0/results/loop_contacts.tsv)
D  anchor-loop C-alpha RMSD to the nearest free structure of the same protein (another
   entry), free chains vs partner-bound chains (stage0/results/loop_summary.tsv)
"""
import pandas as pd

from sc_style import *  # noqa: F403

S0 = ROOT / "stage0" / "results"
STATES = ["ATP-like", "ADP+TS-mimic", "ADP", "AMP"]


def jitter(n, w=0.12, seed=0):
    return np.random.default_rng(seed).uniform(-w, w, n)


def panel_a(ax, n):
    for aa, col, lab in (("PHE", FOREST, "Phe anchor"), ("TRP", CRIMSON, "Trp anchor")):
        x = n[n.ua.str.startswith(aa)]
        if len(x):
            ax.scatter(x.approach_angle, x.n6_to_ring, s=9, color=col, lw=0, alpha=ALPHA, label=f"{lab} (n = {len(x)})")
    ax.axvspan(0, 30, color=SC_GRID, lw=0)
    ax.text(15, 6.35, "face-on\n(NH-π)", ha="center", va="top", fontsize=5.5, color=SC_ANNOT)
    ax.set_xlim(0, 95)
    ax.set_ylim(4.0, 6.5)
    ax.set_xlabel("Approach angle to ring normal (°)")
    ax.set_ylabel("N6 – ring centroid (Å)")
    ax.legend(loc="lower left")


def panel_b(ax, n):
    for i, st in enumerate(STATES):
        x = n[n.state == st]
        for dx, col, key in ((-0.18, FOREST, "min_ring_atom_dist"), (0.18, SAND, "q_hbond_n6")):
            y = x[key].dropna()
            ax.scatter(i + dx + jitter(len(y), 0.08, i), y, s=5, color=col, lw=0, alpha=ALPHA)
            ax.plot([i + dx - 0.13, i + dx + 0.13], [y.median()] * 2, color=SC_AXIS, lw=1.0)
    ax.set_xticks(range(len(STATES)), [f"{s}\n(n = {(n.state == s).sum()})" for s in STATES], fontsize=6)
    ax.set_ylabel("Distance to adenine N6 (Å)")
    ax.set_ylim(2.4, 4.6)
    ax.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=FOREST, ms=3, label="Anchor ring edge"),
                       mpl.lines.Line2D([], [], marker="o", ls="", color=SAND, ms=3, label="Q-motif Gln")],
              loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.1))


def panel_c(ax):
    c = pd.read_csv(S0 / "loop_contacts.tsv", sep="\t")
    chains = c[["pdb", "chain"]].drop_duplicates()
    touched = c.assign(hit=c.partner_residues.notna()).groupby("pos").hit.mean() * 100
    buried = c.groupby("pos").d_sasa.mean()
    pos = touched.index
    cols = [CRIMSON if p == 0 else TERRACOTTA if p == 1 else SAGE for p in pos]
    ax.bar(pos, touched.values, 0.8, color=cols, alpha=ALPHA, lw=0)
    ax.set_xlabel("Position relative to the anchor")
    ax.set_ylabel(f"Partner-bound chains in contact (%)\n(n = {len(chains)} chains)")
    ax.set_xticks([-12, -8, -4, 0, 4])
    ax.set_xticklabels(["−12", "−8", "−4", "UA", "+4"])
    ax.text(1, touched.loc[1] + 1, "+1", ha="center", va="bottom", fontsize=6, color=TERRACOTTA)
    ax.annotate("anchor: never\ncontacted", xy=(0, 1), xytext=(-5.5, touched.max() * 0.8), fontsize=6,
                color=SC_ANNOT, arrowprops={"arrowstyle": "-", "color": SC_ANNOT, "lw": 0.6})
    ax2 = ax.twinx()
    ax2.plot(buried.index, buried.values, color=SC_ANNOT, lw=0.9, marker="o", ms=2)
    ax2.set_ylabel("Mean buried area (Å²)", color=SC_ANNOT)
    ax2.spines["right"].set_visible(True)
    ax2.spines["top"].set_visible(False)


def panel_d(ax):
    s = pd.read_csv(S0 / "loop_summary.tsv", sep="\t").dropna(subset=["loop_rmsd_nearest_free"])
    prots = [p for p in s.protein.unique() if (s[s.protein == p].label != "free").any()
             and (s[s.protein == p].label == "free").any()]
    for i, p in enumerate(prots):
        x = s[s.protein == p]
        f, b = x[x.label == "free"], x[x.label != "free"]
        ax.scatter(i - 0.15 + jitter(len(f), 0.06, i), f.loop_rmsd_nearest_free, s=7, color=SAGE, lw=0)
        ax.scatter(i + 0.15 + jitter(len(b), 0.06, i + 50), b.loop_rmsd_nearest_free, s=7, color=TERRACOTTA, lw=0)
        last = -1.0
        for lab, y in sorted(b.groupby("label").loop_rmsd_nearest_free.max().items(), key=lambda t: t[1]):
            y = max(y, last + 0.09)  # keep stacked labels from overlapping
            ax.text(i + 0.27, y, lab, fontsize=5, va="center", color=SC_ANNOT)
            last = y
    ax.set_xticks(range(len(prots)), [p.replace(" (", "\n(") for p in prots], fontsize=6)
    ax.set_ylabel("Anchor-loop Cα RMSD to\nnearest free structure (Å)")
    ax.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=SAGE, ms=3, label="free"),
                       mpl.lines.Line2D([], [], marker="o", ls="", color=TERRACOTTA, ms=3, label="with partner")],
              loc="upper left")


def main():
    n = pd.read_csv(S0 / "n6_geometry.tsv", sep="\t")
    fig = plt.figure(figsize=(JOURNAL_2COL, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15])
    a, b, c, d = (fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2))
    panel_a(a, n)
    panel_b(b, n)
    panel_c(c)
    panel_d(d)
    for ax, lab in zip((a, b, c, d), "ABCD"):
        add_panel_label(ax, lab, x=-0.2)
    save_figure(fig, "fig2_contact_pivot")
    print(f"Fig 2: {len(n)} sites; edge-on median angle {n.approach_angle.median():.0f} deg")


if __name__ == "__main__":
    main()
