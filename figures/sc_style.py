"""SC Lab figure style (Desert & Moss) plus this project's colour key.

Import in every figure script:  from sc_style import *
Colour key used across all figures of the paper:
  anchor residue  Trp = Crimson, Phe = Deep Forest, Tyr = Sand, none = Sage
  families        DEAD = Terracotta, other SF2 = Sage
  ClinVar         pathogenic/likely = Crimson, uncertain = Sand, benign/likely = Moss
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt  # noqa: F401  (re-exported)
import numpy as np  # noqa: F401

SC_COLORS = ["#9B2335", "#C96A3A", "#DFB166", "#B8C4A2", "#5E7A52", "#2F4A3E"]
CRIMSON, TERRACOTTA, SAND, SAGE, MOSS, FOREST = SC_COLORS
SC_ANNOT, SC_GRID, SC_AXIS = "#4A4A4A", "#E8E8E8", "#2C2C2C"
SC_WARM_CMAP = mpl.colors.LinearSegmentedColormap.from_list("sc_warm", ["#F5ECD7", SAND, TERRACOTTA, CRIMSON])
SC_DIV_CMAP = mpl.colors.LinearSegmentedColormap.from_list("sc_div", [FOREST, "#F5F5F5", CRIMSON])

AA_COLOR = {"W": CRIMSON, "F": FOREST, "Y": SAND, "-": SAGE}
FAMILY_COLOR = {"DEAD": TERRACOTTA, "other": SAGE}
CLINVAR_COLOR = {"pathogenic": CRIMSON, "uncertain": SAND, "benign": MOSS}
ALPHA = 0.8  # filled areas and bars

mpl.rcParams.update({
    # Arial where installed; Liberation Sans is metric-identical to Arial
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75,
    "axes.edgecolor": SC_AXIS, "axes.facecolor": "white", "axes.grid": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.major.width": 0.75, "ytick.major.width": 0.75,
    "xtick.color": SC_AXIS, "ytick.color": SC_AXIS,
    "lines.linewidth": 1.5, "lines.markersize": 5,
    "legend.frameon": False, "legend.handlelength": 1.2,
    "figure.facecolor": "white", "figure.dpi": 300,
    "savefig.facecolor": "white", "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

JOURNAL_1COL, JOURNAL_15COL, JOURNAL_2COL = 3.5, 5.0, 7.0
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "out"


def pval_to_stars(p):
    return "ns" if p > 0.05 else "*" if p > 0.01 else "**" if p > 0.001 else "***" if p > 0.0001 else "****"


def add_panel_label(ax, label, x=-0.18, y=1.05):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top", ha="left", color=SC_AXIS)


def save_figure(fig, stem):
    """PDF (600 dpi) + TIFF (300 dpi) for the journal, PNG (150 dpi) for quick viewing."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", dpi=600)
    fig.savefig(OUT / f"{stem}.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / f"{stem}.png", dpi=150)
