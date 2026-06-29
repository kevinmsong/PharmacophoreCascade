#!/usr/bin/env python3
"""Regenerate the GLP-1R benchmark summary figure (Figure 2, benchmark_plots.pdf)
with professional, publication-ready labels, matching the style of the
cross-system generalization figure. Reads the benchmark summary CSV and plots
ROC-AUC, PR-AUC, EF1%, and BEDROC for the cascade and the two conventional
baselines, with stratified-bootstrap 95% CI error bars.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"
CSV = ROOT / "evidence/outputs/benchmark_glp1r_full/benchmark_summary.csv"

METHODS = ["full_cascade", "stage3_only", "standard_3d_pharmacophore"]
METHOD_LABEL = {"full_cascade": "Full cascade", "stage3_only": "Stage-3 only",
                "standard_3d_pharmacophore": "Single-pass 3D"}
COLOR = {"full_cascade": "#2c7fb8", "stage3_only": "#7fb8d8",
         "standard_3d_pharmacophore": "#b0b0b0"}
PANELS = [("roc_auc", "ROC-AUC"), ("pr_auc", "PR-AUC"),
          ("ef_1pct", "EF1%"), ("bedroc", r"BEDROC ($\alpha = 20$)")]


def main() -> None:
    df = pd.read_csv(CSV).set_index("method")
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
                         "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    labels = [METHOD_LABEL[m] for m in METHODS]
    colors = [COLOR[m] for m in METHODS]
    for ax, (col, title) in zip(axes.flat, PANELS):
        vals = [df.loc[m, col] for m in METHODS]
        lo = [df.loc[m, col] - df.loc[m, f"{col}_ci_low"] for m in METHODS]
        hi = [df.loc[m, f"{col}_ci_high"] - df.loc[m, col] for m in METHODS]
        x = range(len(METHODS))
        ax.bar(x, vals, color=colors, width=0.65,
               yerr=[lo, hi], capsize=4, error_kw={"elinewidth": 1.2, "ecolor": "#444"})
        ax.set_title(title)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels)
        ax.spines[["top", "right"]].set_visible(False)
        if col == "ef_1pct":
            ax.set_ylabel("Fold enrichment")
        elif col == "bedroc":
            ax.set_ylabel("Score")
        else:
            ax.set_ylabel("Area under curve")
    fig.suptitle("Retrospective GLP-1R benchmark: cascade versus conventional baselines",
                 fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "benchmark_plots.pdf")
    fig.savefig(OUT / "benchmark_plots.png", dpi=150)
    print("wrote benchmark_plots.pdf")


if __name__ == "__main__":
    main()
