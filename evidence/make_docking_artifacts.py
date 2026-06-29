#!/usr/bin/env python3
"""Build the redocking table (Table tab:top10docking) and Figure S3 from the
current top-10 final-ranked GLP-1R redocking summary (reviewer points 9-11).

Inputs:
  evidence/outputs/revision/docking_top10_summary.csv
    (ligand_id, final_rank, best_active, best_inactive, active_pref)
Outputs:
  ACS_Omega_submission/top10_docking_table.tex
  ACS_Omega_submission/figure9_orthogonal_docking.pdf
and prints summary statistics for the main text / response letter.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"
SUMMARY = ROOT / "evidence/outputs/revision/docking_top10_summary.csv"


def _fmt(v: float) -> str:
    return f"{v:.1f}" if pd.notna(v) else "--"


def make_table(df: pd.DataFrame) -> None:
    lines = [
        r"\begin{table}[tbp]", r"\centering",
        r"\caption{\textbf{Redocking of the current top-10 final-ranked GLP-1R ligands "
        r"against active- and inactive-state receptors.} Each ligand was docked with "
        r"AutoDock Vina (1.2.7) against three active-state structures (PDB 6X18, 7KI0, "
        r"7LCJ) and two inactive-state structures (PDB 5VEW, 6LN2) using three random "
        r"seeds per structure (exhaustiveness 16, 9 modes). Entries are the best "
        r"(most negative) affinity over the structures of each state. The active-state "
        r"preference ($\Delta = $ best active $-$ best inactive; negative favors the "
        r"active state) is reported only for ligands that produced a valid pose against "
        r"both states.}",
        r"\label{tab:top10docking}", r"\footnotesize",
        r"\begin{tabular}{@{}clrrr@{}}", r"\toprule",
        r"Final rank & Ligand ID & Best active & Best inactive & $\Delta$ (act.\ pref.) \\",
        r" & & (kcal/mol) & (kcal/mol) & (kcal/mol) \\", r"\midrule",
    ]
    for _, r in df.iterrows():
        dpref = f"{r['active_pref']:+.1f}" if pd.notna(r["active_pref"]) else "--"
        lines.append(
            f"{int(r['final_rank'])} & {r['ligand_id']} & {_fmt(r['best_active'])} & "
            f"{_fmt(r['best_inactive'])} & {dpref} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (OUT / "top10_docking_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote top10_docking_table.tex")


def make_figure(df: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 15,
                         "legend.fontsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 12})
    d = df.sort_values("final_rank")
    x = np.arange(len(d))
    w = 0.4
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(x - w / 2, -d["best_active"].values, w, label="Active state (best of 6X18/7KI0/7LCJ)",
           color="#2c7fb8")
    inact = -d["best_inactive"].values
    ax.bar(x + w / 2, np.nan_to_num(inact, nan=0.0), w,
           label="Inactive state (best of 5VEW/6LN2)", color="#f0a202")
    for xi, vi in zip(x, inact):
        if np.isnan(vi):
            ax.text(xi + w / 2, 0.1, "n/a", ha="center", va="bottom", fontsize=9, color="#888")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{int(r)}\n{l[:12]}" for r, l in zip(d["final_rank"], d["ligand_id"])],
                       rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Best Vina affinity\n($-$kcal/mol, higher = stronger)")
    ax.set_xlabel("Final-ranked ligand (rank, ID)")
    ax.set_title("Redocking of top-10 final-ranked GLP-1R ligands by receptor state")
    ax.legend(loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "figure9_orthogonal_docking.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figure9_orthogonal_docking.png", dpi=150, bbox_inches="tight")
    print("wrote figure9_orthogonal_docking.pdf")


def main() -> None:
    df = pd.read_csv(SUMMARY)
    df = df.sort_values("final_rank").reset_index(drop=True)
    make_table(df)
    make_figure(df)
    both = df.dropna(subset=["active_pref"])
    print("\n=== docking summary stats (for text / letter) ===")
    print(f"ligands docked in both states: {len(both)}/{len(df)}")
    print(f"median best active affinity (all): {df['best_active'].median():.2f} kcal/mol")
    print(f"best (most negative) active affinity: {df['best_active'].min():.2f} kcal/mol")
    if len(both):
        print(f"median active-state preference (both-state): {both['active_pref'].median():+.2f} kcal/mol")
        print(f"prefer active state: {(both['active_pref'] < 0).sum()}/{len(both)}")


if __name__ == "__main__":
    main()
