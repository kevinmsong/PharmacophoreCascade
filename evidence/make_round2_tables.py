#!/usr/bin/env python3
"""Generate the second-round revision LaTeX tables from benchmark output CSVs.

Tolerant of missing inputs: writes whichever tables have data available.
Outputs \\input-able .tex files into ACS_Omega_submission/.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"
OUTPUTS = ROOT / "evidence" / "outputs"
SYS_LABEL = {"glp1r": "GLP-1R", "ghsr": "GHSR", "ntsr1": "NTSR1", "mdm2": "MDM2--p53"}
SYSTEMS = ["glp1r", "ghsr", "ntsr1", "mdm2"]
# In-domain actives per system: GLP-1R is data-limited to 10; the others expanded to 50.
N_ACTIVES = {"glp1r": 10, "ghsr": 50, "ntsr1": 50, "mdm2": 50}


def _summary(system: str):
    p = OUTPUTS / f"benchmark_{system}_full" / "benchmark_summary.csv"
    return pd.read_csv(p).set_index("method") if p.exists() else None


def native_only_table() -> None:
    rows = []
    for s in SYSTEMS:
        d = _summary(s)
        if d is None or "native_only" not in d.index:
            continue
        fc, no = d.loc["full_cascade"], d.loc["native_only"]
        rows.append((SYS_LABEL[s], N_ACTIVES[s], fc["roc_auc"], no["roc_auc"], fc["bedroc"], no["bedroc"]))
    if not rows:
        return
    lines = [r"\begin{table}[tbp]", r"\centering",
             r"\caption{\textbf{Native-only baseline versus the full cascade.} "
             r"Each molecule was scored directly with the terminal native pharmacophore (native-only), "
             r"using the same preparation and tie-breaking as the cascade's native branch but bypassing "
             r"Stages 0--3. The terminal native scoring supplies most of the retrospective enrichment, "
             r"matching or exceeding the full cascade on GLP-1R, GHSR, and MDM2--p53; on NTSR1 the staged "
             r"topological gates instead add discriminative power, so the full cascade leads there. The "
             r"upstream cascade therefore contributes scalability across systems and, in some, additional "
             r"enrichment. $n$ is the number of in-domain actives (GLP-1R is data-limited to 10; the others "
             r"were expanded to 50), each with 30 matched decoys.}",
             r"\label{tab:nativeonly}", r"\footnotesize",
             r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
             r" & & \multicolumn{2}{c}{ROC-AUC} & \multicolumn{2}{c}{BEDROC ($\alpha=20$)} \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
             r"System & $n$ & Full cascade & Native-only & Full cascade & Native-only \\", r"\midrule"]
    for name, n, fr, nr, fb, nb in rows:
        lines.append(f"{name} & {n} & {fr:.3f} & {nr:.3f} & {fb:.3f} & {nb:.3f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (OUT / "native_only_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote native_only_table.tex")


def production_table() -> None:
    rows_by_sys = {}
    for s in SYSTEMS:
        p = OUTPUTS / f"production_{s}" / "stage_survival.csv"
        if p.exists():
            rows_by_sys[s] = pd.read_csv(p)
    if not rows_by_sys:
        return
    stages = ["input", "stage0_pass", "stage1_pass", "stage3_shortlist", "native_pool", "final_ranked"]
    stage_disp = {"input": "Input", "stage0_pass": "Stage 0", "stage1_pass": "Stage 1",
                  "stage3_shortlist": "Stage-3 shortlist (5\\%)", "native_pool": "Native pool",
                  "final_ranked": "Final ranking"}
    lines = [r"\begin{table}[tbp]", r"\centering",
             r"\caption{\textbf{Active survival under production constraints (stage by stage).} "
             r"Each system's actives (10 for GLP-1R, which is data-limited; 50 for the other systems) and their "
             r"matched decoys (30 per active) were embedded in a 30,000-molecule ZINC background and screened "
             r"with the real production pipeline (the Stage 0--2 gate, the 5\% Stage-3 shortlist, and the "
             r"native-pool caps of the million-compound run). Entries give surviving actives; the Input row "
             r"gives each system's active count.}",
             r"\label{tab:production}", r"\footnotesize", r"\begin{tabular}{@{}l" + "r" * len(rows_by_sys) + "@{}}",
             r"\toprule", "Stage & " + " & ".join(SYS_LABEL[s] for s in rows_by_sys) + r" \\", r"\midrule"]
    for st in stages:
        cells = []
        for s in rows_by_sys:
            df = rows_by_sys[s].set_index("stage")
            cells.append(str(int(df.loc[st, "actives"])) if st in df.index else "--")
        lines.append(f"{stage_disp[st]} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (OUT / "production_survival_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote production_survival_table.tex")


def robustness_table() -> None:
    """Decoy-selection robustness: n=10 five-replicate mean +/- SD, plus the single
    expanded-benchmark (n=50, 1 set) value where available (GLP-1R is data-limited)."""
    p = OUTPUTS / "revision" / "decoy_replicate_robustness.csv"
    if not p.exists():
        return
    agg = pd.read_csv(p)
    fc = agg[agg["method"] == "full_cascade"].set_index("system")
    if fc.empty:
        return

    def ms(r, m, prec=3):
        if (m + "_mean") not in r or pd.isna(r[m + "_mean"]):
            return "--"
        sd = r.get(m + "_std")
        return (f"{r[m+'_mean']:.{prec}f} $\\pm$ {sd:.{prec}f}" if pd.notna(sd)
                else f"{r[m+'_mean']:.{prec}f}")

    lines = [r"\begin{table}[tbp]", r"\centering",
             r"\caption{\textbf{Robustness to decoy selection (full cascade).} "
             r"For each system, the original 10-active benchmark is resampled across five independently "
             r"sampled, equally property-matched decoy sets (mean $\pm$ SD); the expanded 50-active benchmark "
             r"(GHSR, NTSR1, MDM2--p53) is reported for a single matched decoy set. The small five-set spread "
             r"shows the enrichment does not depend on one favorable decoy draw, and it persists at the "
             r"expanded active sets. GLP-1R is data-limited to 10 in-domain actives.}",
             r"\label{tab:robustness}", r"\footnotesize",
             r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
             r" & \multicolumn{2}{c}{$n=10$ (5 decoy sets, mean $\pm$ SD)} & \multicolumn{2}{c}{$n=50$ (1 decoy set)} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r"System & ROC-AUC & BEDROC & ROC-AUC & BEDROC \\", r"\midrule"]
    for s in SYSTEMS:
        if s not in fc.index:
            continue
        r10 = fc.loc[s]
        roc10, bed10 = ms(r10, "roc_auc"), ms(r10, "bedroc")
        d50 = _summary(s)
        if s != "glp1r" and d50 is not None and "full_cascade" in d50.index:
            roc50 = f"{d50.loc['full_cascade','roc_auc']:.3f}"
            bed50 = f"{d50.loc['full_cascade','bedroc']:.3f}"
        else:
            roc50 = bed50 = "--"
        lines.append(f"{SYS_LABEL[s]} & {roc10} & {bed10} & {roc50} & {bed50} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (OUT / "decoy_robustness_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote decoy_robustness_table.tex")


if __name__ == "__main__":
    native_only_table()
    production_table()
    robustness_table()
