#!/usr/bin/env python3
"""Emit SI LaTeX tables for the additional peptide-receptor case studies."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "ACS_Omega_submission"
DATA = ROOT / "evidence" / "data"

SYSTEMS = {
    "ghsr": ("Ghrelin receptor (GHSR)", "CHEMBL4616", "7F9Y", "EC$_{50}$"),
    "ntsr1": ("Neurotensin receptor 1 (NTSR1)", "CHEMBL4123", "4GRV", "EC$_{50}$"),
    "mdm2": ("MDM2--p53 interface", "CHEMBL5023", "1YCR", "IC$_{50}$/$K_i$"),
}


def actives_table(systems: list[str]) -> None:
    lines = [
        "{\\footnotesize",
        "\\setlength{\\tabcolsep}{4pt}",
        "\\begin{longtable}{@{}l l l r r r@{}}",
        "\\caption{\\textbf{In-domain actives used in the additional retrospective benchmarks.} "
        "Actives were curated from ChEMBL with the same Stage-0 design-space, canonical-SMILES, and "
        "Murcko-scaffold deduplication filters used for GLP-1R; peptide-GPCR systems use functional "
        "agonist potencies and the MDM2--p53 interface uses p53-pocket binding/inhibition potencies. "
        "Each active has a distinct Murcko scaffold.}",
        "\\label{tab:si_new_actives}\\\\",
        "\\toprule",
        "System & Ligand (ChEMBL) & Potency type & Potency (nM) & MW & cLogP \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        "System & Ligand (ChEMBL) & Potency type & Potency (nM) & MW & cLogP \\\\",
        "\\midrule",
        "\\endhead",
    ]
    for key in systems:
        name, _, _, _ = SYSTEMS[key]
        df = pd.read_csv(DATA / f"{key}_external_benchmark_library.csv")
        act = df[df["label"] == "active"].copy()
        act = act.sort_values("potency_value")
        first = True
        for _, r in act.iterrows():
            sys_label = name if first else ""
            first = False
            pv = r["potency_value"]
            pv_s = f"{pv:.3g}" if pd.notna(pv) else "--"
            lines.append(
                f"{sys_label} & {r['ligand_id']} & {r.get('potency_type','')} & {pv_s} & "
                f"{r['mw']:.0f} & {r['logp']:.2f} \\\\"
            )
        lines.append("\\midrule")
    if lines[-1] == "\\midrule":
        lines[-1] = "\\bottomrule"
    else:
        lines.append("\\bottomrule")
    lines += ["\\end{longtable}", "}", ""]
    (OUT / "si_new_systems_actives_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote si_new_systems_actives_table.tex")


def counts_table(systems: list[str]) -> None:
    summ = json.loads((DATA / "additional_benchmarks_summary.json").read_text())
    lines = [
        "\\begin{table}[H]",
        "\\centering",
        "\\footnotesize",
        "\\caption{\\textbf{Composition of the additional retrospective benchmarks.} "
        "Each labeled universe contains 10 in-domain actives (distinct Murcko scaffolds) and 300 "
        "property-matched, scaffold-distinct, ECFP4-dissimilar decoys drawn from the local ZINC pool.}",
        "\\label{tab:si_new_counts}",
        "\\begin{tabular}{@{}l l l c c@{}}",
        "\\toprule",
        "System & ChEMBL target & Reference complex & Actives & Decoys \\\\",
        "\\midrule",
    ]
    for key in systems:
        name, chembl, pdb, _ = SYSTEMS[key]
        info = summ.get(key, {})
        lines.append(
            f"{name} & {chembl} & PDB {pdb} & {info.get('actives','--')} & {info.get('decoys','--')} \\\\"
        )
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    (OUT / "si_new_systems_counts_table.tex").write_text("\n".join(lines), encoding="utf-8")
    print("wrote si_new_systems_counts_table.tex")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="+", default=["ghsr", "ntsr1", "mdm2"])
    args = ap.parse_args()
    present = [s for s in args.systems if (DATA / f"{s}_external_benchmark_library.csv").exists()]
    actives_table(present)
    counts_table(present)
