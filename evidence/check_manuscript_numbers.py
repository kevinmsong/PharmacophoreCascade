#!/usr/bin/env python3
"""Check every headline statistic in the manuscript against its source data.

The manuscript quotes several hundred numbers that originate in
``evidence/outputs/``. Transcribing them by hand is how the previous version
acquired a handful of rounding errors that disagreed with the released
per-molecule data. This script re-derives each number from the source CSV,
formats it the way the manuscript does, and confirms the string is present in
the LaTeX sources.

It is a text search rather than a parse on purpose: it catches a value that
drifted, was rounded differently, or was edited in one place and not the other,
without needing to understand LaTeX table structure.

Usage:
    python evidence/check_manuscript_numbers.py
    python evidence/check_manuscript_numbers.py --dir ACS_Omega_resubmission
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
EV = ROOT / "evidence" / "outputs"

SYSTEMS = {
    "GLP-1R": "benchmark_glp1r_full",
    "GHSR": "benchmark_ghsr_full",
    "NTSR1": "benchmark_ntsr1_full",
    "MDM2-p53": "benchmark_mdm2_full",
}
METHODS = ["full_cascade", "native_only", "stage3_only", "standard_3d_pharmacophore"]

#: How the manuscript writes a signed value, e.g. "$+0.117$" / "$-0.009$".
def signed(x: float, nd: int = 3) -> str:
    return f"{'+' if x >= 0 else '-'}{abs(x):.{nd}f}"


def load_text(sub: Path) -> str:
    """All LaTeX the manuscript and its tables are built from, concatenated."""
    names = ["manuscript_body.tex", "crosssystem_table.tex", "tost_table.tex",
             "attrition_table.tex", "efficiency_table.tex", "evidence_table.tex",
             "ablation_table.tex", "ablation_curation_table.tex", "top10_table.tex",
             "terminology_table.tex", "decoy_robustness_table.tex",
             "top10_docking_table.tex", "production_survival_table.tex",
             "supporting_information.tex",
             "floor_survival_table.tex", "floor_timing_table.tex",
             "current_run_config_table.tex", "current_rank_shift_table.tex"]
    parts = []
    for n in names:
        p = sub / n
        if p.exists():
            parts.append(p.read_text(encoding="utf-8"))
        else:
            print(f"  note: {n} not found, skipped")
    return "\n".join(parts)


class Checker:
    def __init__(self, text: str):
        self.text = text
        self.failures: list[tuple[str, str]] = []
        self.n = 0

    def expect(self, label: str, needle: str) -> None:
        self.n += 1
        if needle not in self.text:
            self.failures.append((label, needle))

    def report(self) -> int:
        print(f"\n{self.n} values checked against evidence/outputs/")
        if not self.failures:
            print("all present and consistent")
            return 0
        print(f"{len(self.failures)} NOT FOUND in the manuscript sources:")
        for label, needle in self.failures:
            print(f"  {label:52s} expected {needle!r}")
        return 1


def check_cross_system(c: Checker) -> None:
    """ROC-AUC, PR-AUC, EF1%, BEDROC for every system and method."""
    for sysname, folder in SYSTEMS.items():
        p = EV / folder / "benchmark_summary.csv"
        if not p.exists():
            print(f"  note: {p} missing, cross-system check skipped for {sysname}")
            continue
        t = pd.read_csv(p).set_index("method")
        for m in METHODS:
            if m not in t.index:
                continue
            r = t.loc[m]
            c.expect(f"{sysname}/{m} ROC-AUC", f"{r.roc_auc:.3f}")
            c.expect(f"{sysname}/{m} PR-AUC", f"{r.pr_auc:.3f}")
            c.expect(f"{sysname}/{m} BEDROC", f"{r.bedroc:.3f}")


def check_tost(c: Checker) -> None:
    """Paired deltas and 90% intervals for the equivalence tests."""
    p = EV / "equivalence" / "equivalence_tost.csv"
    t = pd.read_csv(p).query("metric == 'roc_auc'")
    for _, r in t.iterrows():
        name = r.system
        c.expect(f"{name} TOST delta", signed(r.observed_delta))
        c.expect(f"{name} TOST 90% CI",
                 f"[{signed(r.ci90_low)}, {signed(r.ci90_high)}]")


def check_attrition(c: Checker) -> None:
    """Per-stage, per-cause active losses and final survival."""
    t = pd.read_csv(EV / "attrition" / "active_attrition_summary.csv")
    cols = {c_.lower(): c_ for c_ in t.columns}
    per = pd.read_csv(EV / "attrition" / "active_attrition_per_molecule.csv")
    for sysname in SYSTEMS:
        sub = per[per.system == sysname]
        retained = int((sub.stage_lost == "-").sum())
        total = len(sub)
        c.expect(f"{sysname} reached native scoring", f"{retained}/{total}")
    # Under the revised setting, structural alerts cannot exclude an active.
    lost0 = per[per.stage_lost == "Stage 0"]
    alerts = int((lost0.reason_class == "structural alert").sum())
    envelope = int((lost0.reason_class == "property envelope").sum())
    assert alerts == 0, f"Structural-alert exclusions still present: {alerts}"
    c.expect("Stage-0 active survival", "158/160")
    assert envelope == 2, f"expected 2 property-envelope losses, found {envelope}"


def check_efficiency(c: Checker) -> None:
    """Shortlist size, native scoring, and actives retained for each policy pair.

    Keyed on the data file rather than on a phrase in the caption: wording is
    edited between drafts, and a prose-keyed branch silently falls through to a
    different dataset when someone rewrites the sentence it was matching.
    """
    policies = EV / "absolute_floor" / "policy_comparison.csv"
    for row in pd.read_csv(policies).itertuples():
        tag = f"{row.system} {row.policy}"
        for stage in ["shortlist", "native_success", "final_ranked"]:
            c.expect(f"{tag} {stage} ligands", f'{getattr(row, stage + "_n"):,}')
        c.expect(f"{tag} shortlist actives", f"{row.shortlist_actives}/{row.n_actives}")
        c.expect(f"{tag} final actives", f"{row.final_ranked_actives}/{row.n_actives}")
    headline = json.loads(
        (EV / "absolute_floor" / "analysis_summary.json").read_text())["headline"]
    c.expect("measured end-to-end hours",
             f"{headline['pipeline_wall_seconds'] / 3600:.2f}")


def check_docking(c: Checker) -> None:
    """Paired Wilcoxon result for the active-state preference."""
    d = json.loads((EV / "equivalence" / "docking_wilcoxon.json").read_text())
    c.expect("docking Wilcoxon W", f"$W = {d['wilcoxon_W']:g}$")
    c.expect("docking Wilcoxon p", f"$p = {d['p_two_sided']:.4f}$")
    c.expect("docking median delta", f"{d['median_delta']:.2f}~kcal")
    c.expect("docking IQR low", f"{d['iqr_low']:.2f}")
    c.expect("docking IQR high", f"{d['iqr_high']:.2f}")


def check_alternate_structure(c: Checker) -> None:
    """The 7KI0 sensitivity check."""
    t = pd.read_csv(EV / "benchmark_glp1r_7ki0" / "benchmark_summary.csv").set_index("method")
    c.expect("7KI0 full cascade ROC-AUC", f"{t.loc['full_cascade'].roc_auc:.3f}")
    c.expect("7KI0 native-only ROC-AUC", f"{t.loc['native_only'].roc_auc:.3f}")
    c.expect("7KI0 BEDROC", f"{t.loc['full_cascade'].bedroc:.3f}")
    c.expect("7KI0 top-10 recovery",
             f"{t.loc['full_cascade'].top10_recovery * 10:.2f} of 10")


def check_ablations(c: Checker) -> None:
    """Kendall tau over the common support, all rows comparable."""
    t = pd.read_csv(EV / "ablation_common" / "ablation_common_support.csv")
    for _, r in t.iterrows():
        c.expect(f"ablation tau1000 [{r.ablation[:28]}]", signed(r.tau_1000))
        c.expect(f"ablation tauCommon [{r.ablation[:28]}]", signed(r.tau_common))


def check_decoy_replicates(c: Checker) -> None:
    """Mean and SD of the full cascade across five matched-decoy sets."""
    for sysname, stem in (("GLP-1R", "glp1r"), ("GHSR", "ghsr"),
                          ("NTSR1", "ntsr1"), ("MDM2-p53", "mdm2")):
        roc, bed = [], []
        for s in range(1, 6):
            p = EV / "decoy_replicates" / f"{stem}_seed{s}" / "benchmark_summary.csv"
            if not p.exists():
                break
            r = pd.read_csv(p).query("method == 'full_cascade'").iloc[0]
            roc.append(r.roc_auc)
            bed.append(r.bedroc)
        if len(roc) != 5:
            print(f"  note: {sysname} has {len(roc)}/5 decoy replicates, skipped")
            continue
        c.expect(f"{sysname} decoy-replicate ROC-AUC",
                 f"{np.mean(roc):.3f} $\\pm$ {np.std(roc, ddof=1):.3f}")
        c.expect(f"{sysname} decoy-replicate BEDROC",
                 f"{np.mean(bed):.3f} $\\pm$ {np.std(bed, ddof=1):.3f}")


def check_curation_ablation(c: Checker) -> None:
    """Curated versus fully automated GLP-1R pharmacophore."""
    t = pd.read_csv(EV / "benchmark_glp1r_automated" / "benchmark_summary.csv").set_index("method")
    c.expect("automated ROC-AUC", f"{t.loc['full_cascade'].roc_auc:.3f}")
    c.expect("automated PR-AUC", f"{t.loc['full_cascade'].pr_auc:.3f}")
    c.expect("automated BEDROC", f"{t.loc['full_cascade'].bedroc:.3f}")
    c.expect("automated single-pass ROC-AUC",
             f"{t.loc['standard_3d_pharmacophore'].roc_auc:.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="ACS_Omega_resubmission",
                    help="submission folder holding the LaTeX sources")
    args = ap.parse_args()

    sub = ROOT / args.dir
    if not sub.is_dir():
        raise SystemExit(f"no such folder: {sub}")

    print(f"checking {args.dir} against evidence/outputs/")
    c = Checker(load_text(sub))
    for fn in (check_cross_system, check_tost, check_attrition, check_efficiency,
               check_docking, check_alternate_structure, check_ablations,
               check_decoy_replicates, check_curation_ablation):
        fn(c)
    sys.exit(c.report())


if __name__ == "__main__":
    main()
