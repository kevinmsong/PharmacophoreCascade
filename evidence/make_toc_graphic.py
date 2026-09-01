#!/usr/bin/env python3
"""Regenerate the ACS Omega graphical abstract / Table-of-Contents graphic
(toc_graphic.png and toc_graphic.pdf), sized to the ACS TOC box (3.25 in wide
x 1.75 in tall).

The graphic states the paper's actual conclusion: the terminal native
peptide-contact scoring supplies the enrichment, and the staged cascade
upstream supplies tractability. An earlier version headlined the cascade with
the full-cascade ROC-AUC, which advertised a claim the component analysis does
not support.

Colors come from evidence/figstyle.py so the TOC graphic matches the figures.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_resubmission"

SCREEN = fs.OKABE_ITO["blue"]
NATIVE = fs.OKABE_ITO["green"]
GREY = "#4D4D4D"
INK = "#1A1A1A"

# (top half-width, bottom half-width, y0, y1, fill, text color, count, label)
TIERS = [
    (0.235, 0.185, 0.615, 0.790, "#BBD6E8", INK, "1,000,000", "ZINC compounds"),
    (0.178, 0.126, 0.415, 0.590, "#7FB2D4", "white", "47,812", "shortlist"),
    (0.119, 0.070, 0.215, 0.390, SCREEN, "white", "5,000", "native scored"),
]
XC = 0.235


def trapezoid(ax, hw_top, hw_bot, y0, y1, color):
    pts = [(XC - hw_top, y1), (XC + hw_top, y1), (XC + hw_bot, y0), (XC - hw_bot, y0)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color,
                         edgecolor="white", linewidth=1.4))


def main() -> None:
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.955, "Which stage supplies the enrichment?",
            ha="center", va="center", fontsize=16, fontweight="bold", color=INK)
    ax.text(0.5, 0.878, "peptide-contact scoring at four receptor interfaces",
            ha="center", va="center", fontsize=10.5, style="italic", color=GREY)

    # Left: the funnel, labelled as what it actually buys.
    for hw_top, hw_bot, y0, y1, fill, numc, count, label in TIERS:
        trapezoid(ax, hw_top, hw_bot, y0, y1, fill)
        yc = (y0 + y1) / 2
        ax.text(XC, yc + 0.024, count, ha="center", va="center",
                fontsize=15, fontweight="bold", color=numc)
        ax.text(XC, yc - 0.042, label, ha="center", va="center",
                fontsize=9, color=numc)

    ax.text(XC, 0.155, "staged cascade", ha="center", va="center",
            fontsize=12, fontweight="bold", color=SCREEN)
    ax.text(XC, 0.092, "buys tractability: 13.3 h, not 31 CPU-days",
            ha="center", va="center", fontsize=9, color=GREY)

    # Right: the terminal stage, named as the source of the retrieval.
    box = FancyBboxPatch((0.50, 0.36), 0.46, 0.40,
                         boxstyle="round,pad=0.012,rounding_size=0.03",
                         linewidth=1.6, edgecolor=NATIVE, facecolor="#E6F4EF")
    ax.add_patch(box)
    ax.text(0.73, 0.700, "native peptide-contact", ha="center", va="center",
            fontsize=12, fontweight="bold", color=NATIVE)
    ax.text(0.73, 0.635, "pharmacophore scoring", ha="center", va="center",
            fontsize=12, fontweight="bold", color=NATIVE)
    ax.text(0.73, 0.545, "supplies the enrichment", ha="center", va="center",
            fontsize=10, color=INK)
    ax.text(0.73, 0.462, "beats a single-pass 3D pharmacophore", ha="center",
            va="center", fontsize=8.6, color=GREY)
    ax.text(0.73, 0.410, "on all four interfaces tested", ha="center",
            va="center", fontsize=8.6, color=GREY)

    ax.add_patch(FancyArrowPatch((0.335, 0.30), (0.545, 0.345),
                                 arrowstyle="-|>", mutation_scale=13,
                                 linewidth=1.4, color=NATIVE,
                                 connectionstyle="arc3,rad=-0.2"))

    strip = FancyBboxPatch((0.50, 0.135), 0.46, 0.165,
                           boxstyle="round,pad=0.008,rounding_size=0.025",
                           linewidth=1.2, edgecolor=GREY, facecolor="#F2F2F2")
    ax.add_patch(strip)
    ax.text(0.73, 0.255, "native-only matches or beats the cascade",
            ha="center", va="center", fontsize=8.8, color=INK)
    ax.text(0.73, 0.190, "on 3 of 4 systems",
            ha="center", va="center", fontsize=8.8, fontweight="bold", color=INK)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "toc_graphic.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / "toc_graphic.pdf", bbox_inches="tight", pad_inches=0.02)
    print(f"wrote toc_graphic.png and toc_graphic.pdf to {OUT.name}")


if __name__ == "__main__":
    main()
