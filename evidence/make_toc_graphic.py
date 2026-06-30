#!/usr/bin/env python3
"""Regenerate the ACS Omega graphical abstract / Table-of-Contents graphic
(toc_graphic.png and toc_graphic.pdf) as a clean funnel sized to the ACS TOC box
(3.25 in wide x 1.75 in tall). GLP-1R headline survival counts and metrics.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"

NAVY = "#1f3b73"
TITLE = "#15315f"
GREY = "#444444"
# (top half-width, bottom half-width, y0, y1, fill, num_color, count, label, stage, sub)
TIERS = [
    (0.30, 0.235, 0.60, 0.80, "#9ecae1", TITLE, "1,000,000", "ZINC compounds",
     "Stage 0–2", "topological gates", "#2c7fb8"),
    (0.225, 0.16, 0.375, 0.565, "#3f8fc4", "white", "47,812", "shortlist",
     "Stage 3", "3D conformer rerank", "#2c7fb8"),
    (0.15, 0.095, 0.15, 0.34, "#08519c", "white", "1,000", "ranked candidates",
     "Native rerank", "peptide-contact fit", "#2a7d8c"),
]
XC = 0.32
XLAB = 0.66


def trapezoid(ax, hw_top, hw_bot, y0, y1, color):
    pts = [(XC - hw_top, y1), (XC + hw_top, y1), (XC + hw_bot, y0), (XC - hw_bot, y0)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor="white", linewidth=1.5))


def main() -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.955, "Topology-to-Native Pharmacophore Cascade",
            ha="center", va="center", fontsize=15, fontweight="bold", color=TITLE)
    ax.text(0.5, 0.875, "million-scale screening at the GLP-1R interface",
            ha="center", va="center", fontsize=10.5, style="italic", color=GREY)

    for hw_top, hw_bot, y0, y1, fill, numc, count, label, stage, sub, stagec in TIERS:
        trapezoid(ax, hw_top, hw_bot, y0, y1, fill)
        yc = (y0 + y1) / 2
        ax.text(XC, yc + 0.028, count, ha="center", va="center", fontsize=19,
                fontweight="bold", color=numc)
        ax.text(XC, yc - 0.045, label, ha="center", va="center", fontsize=10, color=numc)
        ax.text(XLAB, yc + 0.03, stage, ha="left", va="center", fontsize=13,
                fontweight="bold", color=stagec)
        ax.text(XLAB, yc - 0.04, sub, ha="left", va="center", fontsize=10, color=GREY)

    box = FancyBboxPatch((0.03, 0.005), 0.94, 0.085,
                         boxstyle="round,pad=0.004,rounding_size=0.02",
                         linewidth=1.3, edgecolor=NAVY, facecolor="#eef4fb")
    ax.add_patch(box)
    ax.text(0.5, 0.048, "ROC-AUC 0.800    ·    EF1% 30    ·    5/10 actives recovered",
            ha="center", va="center", fontsize=12, fontweight="bold", color=NAVY)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(OUT / "toc_graphic.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / "toc_graphic.pdf", bbox_inches="tight", pad_inches=0.02)
    print("wrote toc_graphic.png and toc_graphic.pdf")


if __name__ == "__main__":
    main()
