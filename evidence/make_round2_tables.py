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
        rows.append((SYS_LABEL[s], fc["roc_auc"], no["roc_auc"], int(round(fc["top10_recovery"] * 10)),
                     int(round(no["top10_recovery"] * 10))))
    if not rows:
        return
    lines = [r"\begin{table}[tbp]", r"\centering",
             r"\caption{\textbf{A native-only baseline matches or exceeds the full cascade.} "
             r"Each molecule was scored directly with the terminal native pharmacophore (native-only), "
             r"using the same preparation and tie-breaking as the cascade's native branch but bypassing "
             r"Stages 0--3. The retrospective enrichment is supplied by the terminal native stage; the "
             r"upstream cascade contributes scalability, not enrichment. For the three automated-pharmacophore "
             r"systems the strict gate filtered some actives before native scoring, so the native-only "
             r"baseline (which scores every molecule) recovered more of them.}",
             r"\label{tab:nativeonly}", r"\footnotesize",
             r"\begin{tabular}{@{}lrrcc@{}}", r"\toprule",
             r" & \multicolumn{2}{c}{ROC-AUC} & \multicolumn{2}{c}{Actives in top 10} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r"System & Full cascade & Native-only & Full cascade & Native-only \\", r"\midrule"]
    for name, fr, nr, ft, nt in rows:
        lines.append(f"{name} & {fr:.3f} & {nr:.3f} & {ft} & {nt} \\\\")
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
             r"Each system's 10 actives and 300 matched decoys were embedded in a large ZINC background and "
             r"screened with the real production pipeline (the Stage 0--2 gate, the 5\% Stage-3 shortlist, and "
             r"the native-pool caps of the million-compound run). Entries give surviving actives (of 10).}",
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
    p = OUTPUTS / "revision" / "decoy_replicate_robustness.csv"
    if not p.exists():
        return
    agg = pd.read_csv(p)
    fc = agg[agg["method"] == "full_cascade"]
    if fc.empty:
        return
    lines = [r"\begin{table}[tbp]", r"\centering",
             r"\caption{\textbf{Robustness across five independent matched-decoy sets (full cascade).} "
             r"Mean $\pm$ SD across five independently sampled, equally property-matched decoy sets per system. "
             r"The modest spread indicates the reported enrichment does not depend on one favorable decoy selection.}",
             r"\label{tab:robustness}", r"\footnotesize", r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r"System & ROC-AUC & PR-AUC & BEDROC & EF1\% \\", r"\midrule"]
    for _, r in fc.iterrows():
        def ms(m, prec=3):
            if (m + "_mean") not in r or pd.isna(r[m + "_mean"]):
                return "--"
            sd = r.get(m + "_std")
            return (f"{r[m+'_mean']:.{prec}f} $\\pm$ {sd:.{prec}f}" if pd.notna(sd)
                    else f"{r[m+'_mean']:.{prec}f}")
        lines.append(f"{SYS_LABEL.get(r['system'], r['system'])} & {ms('roc_auc')} & {ms('pr_auc')} & {ms('bedroc')} & {ms('ef_1pct', 1)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (OUT / "decoy_robustness_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote decoy_robustness_table.tex")


if __name__ == "__main__":
    native_only_table()
    production_table()
    robustness_table()
