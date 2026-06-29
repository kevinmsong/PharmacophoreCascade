#!/usr/bin/env python3
"""
Generate revision figures and LaTeX tables from benchmark output CSVs.

Reads the per-system ``benchmark_summary.csv`` files plus the docking and
reproducibility artifacts under ``evidence/outputs/revision`` and writes:

* ``ACS_Omega_submission/figure8_generalization.pdf`` - cross-system enrichment of the
  full cascade vs the conventional single-pass 3D pharmacophore.
* ``ACS_Omega_submission/figure9_orthogonal_docking.pdf`` - AutoDock Vina orthogonal
  validation (active vs inactive state).
* LaTeX tables under ``ACS_Omega_submission/`` for the generalization summary,
  the no-manual-curation ablation, and (SI) the docking table.

It is tolerant of which systems are present; pass the systems to include.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Larger fonts for readability at journal-column scale (reviewer point 17).
plt.rcParams.update({
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 15,
})

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"
REV = ROOT / "evidence" / "outputs" / "revision"

# system_key -> (display label, benchmark_summary.csv path)
SUMMARY_PATHS = {
    "glp1r": ROOT / "evidence" / "outputs" / "benchmark_summary.csv",
    "glp1r_automated": ROOT / "evidence" / "outputs" / "benchmark_glp1r_automated" / "benchmark_summary.csv",
    "ghsr": ROOT / "evidence" / "outputs" / "benchmark_ghsr" / "benchmark_summary.csv",
    "ntsr1": ROOT / "evidence" / "outputs" / "benchmark_ntsr1" / "benchmark_summary.csv",
    "mdm2": ROOT / "evidence" / "outputs" / "benchmark_mdm2" / "benchmark_summary.csv",
}
LABELS = {
    "glp1r": "GLP-1R",
    "glp1r_automated": "GLP-1R (automated)",
    "ghsr": "GHSR\\\\(ghrelin)",
    "ntsr1": "NTSR1\\\\(neurotensin)",
    "mdm2": "MDM2--p53",
}
PLAIN = {"glp1r": "GLP-1R", "ghsr": "GHSR (ghrelin)", "ntsr1": "NTSR1 (neurotensin)", "mdm2": "MDM2-p53"}


def load_summary(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path).set_index("method")


def fmt(x, nd=3):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return "--"


def generalization_figure(systems: list[str]) -> None:
    metrics = [("roc_auc", "ROC-AUC"), ("ef_1pct", "EF1%"), ("bedroc", "BEDROC ($\\alpha=20$)")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    x = np.arange(len(systems))
    width = 0.38
    for ax, (m, title) in zip(axes, metrics):
        casc, base = [], []
        for s in systems:
            df = load_summary(SUMMARY_PATHS[s])
            casc.append(float(df.loc["full_cascade", m]) if df is not None else np.nan)
            base.append(float(df.loc["standard_3d_pharmacophore", m]) if df is not None else np.nan)
        ax.bar(x - width / 2, casc, width, label="Full cascade", color="#2b6cb0")
        ax.bar(x + width / 2, base, width, label="Single-pass 3D", color="#a0aec0")
        ax.set_title(title, fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels([PLAIN[s].split(" ")[0] for s in systems], rotation=20, fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Cross-system enrichment: staged cascade vs conventional single-pass pharmacophore", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "figure8_generalization.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "figure8_generalization.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote figure8_generalization")


def docking_figure() -> None:
    p = REV / "docking_orthogonal_summary.csv"
    if not p.exists():
        print("skip docking figure (no data)")
        return
    agg = pd.read_csv(p).sort_values("best_active")
    fig, ax = plt.subplots(figsize=(7, 3.6))
    y = np.arange(len(agg))
    ax.barh(y, agg["best_active"], color="#2b6cb0", label="best active-state")
    mask = agg["best_inactive"].notna()
    ax.scatter(agg.loc[mask, "best_inactive"], y[mask.values], color="#dd6b20", zorder=3, label="best inactive-state")
    ax.set_yticks(y)
    ax.set_yticklabels(agg["ligand_id"], fontsize=7)
    ax.set_xlabel("AutoDock Vina affinity (kcal/mol)")
    ax.set_title("Top native-ranked ligands: active vs inactive GLP-1R", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "figure9_orthogonal_docking.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "figure9_orthogonal_docking.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote figure9_orthogonal_docking")


def generalization_table(systems: list[str]) -> None:
    rows = []
    for s in systems:
        df = load_summary(SUMMARY_PATHS[s])
        if df is None:
            continue
        fc = df.loc["full_cascade"]
        b = df.loc["standard_3d_pharmacophore"]
        n_act = int(round(float(fc.get("top10_recovery", 0)) * 10)) if "top10_recovery" in fc else None
        rows.append(
            {
                "system": PLAIN[s],
                "fc_roc": fmt(fc["roc_auc"]),
                "fc_pr": fmt(fc["pr_auc"]),
                "fc_ef1": fmt(fc["ef_1pct"], 0),
                "fc_bedroc": fmt(fc["bedroc"]),
                "b_roc": fmt(b["roc_auc"]),
                "b_ef1": fmt(b["ef_1pct"], 0),
            }
        )
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{Generalization of the cascade across peptide-receptor interfaces.} "
        "Retrospective active-vs-decoy enrichment of the full cascade and the conventional "
        "single-pass 3D pharmacophore baseline on each system's labeled benchmark, using a "
        "fully automated interface pharmacophore (no manual curation) for the three additional "
        "systems. Larger values indicate better enrichment.}",
        "\\label{tab:generalization}",
        "\\footnotesize",
        "\\begin{tabular}{@{}lrrrrrr@{}}",
        "\\toprule",
        " & \\multicolumn{4}{c}{Full cascade} & \\multicolumn{2}{c}{Single-pass 3D} \\\\",
        "\\cmidrule(lr){2-5}\\cmidrule(lr){6-7}",
        "System & ROC-AUC & PR-AUC & EF1\\% & BEDROC & ROC-AUC & EF1\\% \\\\",
        "\\midrule",
    ]
    for r in rows:
        lines.append(
            f"{r['system']} & {r['fc_roc']} & {r['fc_pr']} & {r['fc_ef1']} & {r['fc_bedroc']} & {r['b_roc']} & {r['b_ef1']} \\\\"
        )
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    (OUT / "generalization_summary_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote generalization_summary_table.tex")


def ablation_table() -> None:
    cur = load_summary(SUMMARY_PATHS["glp1r"])
    auto = load_summary(SUMMARY_PATHS["glp1r_automated"])
    if cur is None or auto is None:
        print("skip ablation table")
        return
    mets = [("roc_auc", "ROC-AUC"), ("pr_auc", "PR-AUC"), ("ef_1pct", "EF1\\%"),
            ("bedroc", "BEDROC ($\\alpha=20$)"), ("top10_recovery", "Actives in top 10")]
    lines = [
        "\\begin{table}[tbp]",
        "\\centering",
        "\\caption{\\textbf{Cascade performance with and without manual curation (GLP-1R).} "
        "The fully automated interface pharmacophore removes all hand-picked contact groups, "
        "weight bonuses, native-supported expansion, and curated rescue; priority and contact "
        "groups are assigned automatically. The staged architecture retains strong global "
        "discrimination without any manual curation (and still far exceeds the single-pass "
        "baseline), while manual curation specifically improves early recovery.}",
        "\\label{tab:ablation_curation}",
        "\\footnotesize",
        "\\begin{tabular}{@{}lrr@{}}",
        "\\toprule",
        "Metric & Curated & Automated (no curation) \\\\",
        "\\midrule",
    ]
    for m, lab in mets:
        cv = float(cur.loc["full_cascade", m])
        av = float(auto.loc["full_cascade", m])
        if m == "top10_recovery":
            lines.append(f"{lab} & {int(round(cv*10))} of 10 & {int(round(av*10))} of 10 \\\\")
        elif m == "ef_1pct":
            lines.append(f"{lab} & {cv:.0f} & {av:.0f} \\\\")
        else:
            lines.append(f"{lab} & {cv:.3f} & {av:.3f} \\\\")
    # single-pass baseline ROC for context
    lines.append("\\midrule")
    lines.append(
        f"Single-pass 3D ROC-AUC & {float(cur.loc['standard_3d_pharmacophore','roc_auc']):.3f} & "
        f"{float(auto.loc['standard_3d_pharmacophore','roc_auc']):.3f} \\\\"
    )
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    (OUT / "ablation_curation_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote ablation_curation_table.tex")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="+", default=["glp1r", "ghsr", "ntsr1", "mdm2"])
    args = ap.parse_args()
    present = [s for s in args.systems if SUMMARY_PATHS[s].exists()]
    print("systems present:", present)
    generalization_figure(present)
    generalization_table(present)
    ablation_table()
    docking_figure()


if __name__ == "__main__":
    main()
