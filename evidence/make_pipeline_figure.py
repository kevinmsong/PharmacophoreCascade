#!/usr/bin/env python3
"""Regenerate the cascade architecture schematic (Scheme 1, figure1_pipeline.png)
as a clean, de-jargoned funnel: the five screening stages plus the terminal
native branch, with the real million-compound survival counts. Replaces the
earlier flowchart that carried implementation jargon ("bitmask", "heap").
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"

# (title, subtitle, count, color) top -> bottom
SCREEN = "#2c7fb8"
NATIVE = "#2a9d8f"
INPUT = "#5a6473"
FINAL = "#1d3f5e"

STAGES = [
    ("Input library", "One million ZINC compounds (tranches H17-H20)", "1,000,000", INPUT, 0.92),
    ("Stage 0 · Standardization and property/chemistry gate",
     "Drug-likeness and structural-alert filters", "956,240", SCREEN, 0.85),
    ("Stage 1 · Receptor hotspot compatibility",
     "Weighted match to GLP-1R contact features", "948,971", SCREEN, 0.78),
    ("Stage 2 · Pairwise feature-distance agreement",
     "Blended with Stage 1 into the cascade score", "blended", SCREEN, 0.71),
    ("Cascade-score shortlist (top 5%)",
     "Only the shortlist proceeds to 3D scoring", "47,812", SCREEN, 0.64),
    ("Stage 3 · 3D conformer-level geometric rerank",
     "ETKDG conformers aligned to the query", "47,689", SCREEN, 0.57),
]
NATIVE_STAGES = [
    ("Native branch · diversified candidate pool",
     "Scaffold-aware caps preserve chemotype diversity", "20,000 → 5,000", NATIVE, 0.46),
    ("Native peptide-contact rescoring",
     "Scored against the active GLP-1 receptor complex", "4,997 scored", NATIVE, 0.39),
    ("Final ranked output",
     "Peptide-mimicking small molecules", "1,000", FINAL, 0.32),
]


def box(ax, y, title, subtitle, count, color, width):
    x0 = 0.5 - width / 2
    ax.add_patch(FancyBboxPatch((x0, y - 0.024), width, 0.048,
                                boxstyle="round,pad=0.006,rounding_size=0.012",
                                linewidth=1.1, edgecolor=color, facecolor=color + "22"))
    ax.text(0.5, y + 0.009, title, ha="center", va="center", fontsize=11.5,
            fontweight="bold", color="#1a1a1a")
    ax.text(0.5, y - 0.012, subtitle, ha="center", va="center", fontsize=9.3, color="#444")
    if count:
        ax.annotate(count, xy=(x0 + width + 0.018, y), xycoords="data", ha="left", va="center",
                    fontsize=10, color=color, fontweight="bold")


def main() -> None:
    fig, ax = plt.subplots(figsize=(8.4, 10.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0.26, 1.0)
    ax.axis("off")

    rows = STAGES + NATIVE_STAGES
    # funnel: width narrows down the screening stages
    widths = [0.66, 0.78, 0.78, 0.78, 0.70, 0.70, 0.72, 0.66, 0.50]
    ys = [r[4] for r in rows]
    for (title, sub, cnt, col, y), w in zip(rows, widths):
        box(ax, y, title, sub, cnt, col, w)
    # arrows between consecutive boxes
    for i in range(len(ys) - 1):
        y_top = ys[i] - 0.024
        y_bot = ys[i + 1] + 0.024
        col = "#2a9d8f" if i >= len(STAGES) - 1 else "#2c7fb8"
        ax.add_patch(FancyArrowPatch((0.5, y_top), (0.5, y_bot),
                                     arrowstyle="-|>", mutation_scale=14,
                                     linewidth=1.3, color=col))
    ax.text(0.5, 0.975, "Topological-to-3D Native Pharmacophore Cascade",
            ha="center", va="center", fontsize=13.5, fontweight="bold", color="#1a1a1a")
    fig.tight_layout()
    fig.savefig(OUT / "figure1_pipeline.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "figure1_pipeline.pdf", bbox_inches="tight")
    print("wrote figure1_pipeline.png")


if __name__ == "__main__":
    main()
