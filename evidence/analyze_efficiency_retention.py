#!/usr/bin/env python3
"""Efficiency-retention trade-off of the staged cascade.

The cascade's stated purpose is to make terminal native peptide-contact
scoring affordable at library scale. Whether that trade is favorable depends
on two quantities that the manuscript previously reported only at a single
operating point: how much native scoring the cascade avoids, and how many true
actives it discards on the way.

This script sweeps the Stage-3 shortlist fraction from 0.5% to 100% for every
system and reports, at each setting, the actives retained and the number of
molecules that would reach native scoring. The result is a retention-versus-
compute curve on which the production 5% setting can be located, so a reader
can judge the trade rather than take it on assertion.

It also quantifies the retention that the strict chemistry gate costs, by
recomputing survival with structural-alert rejections readmitted. Those
rejections are a configuration choice, not a property of the topological
front end, and they turn out to account for every Stage-0 active loss.

Reads ``evidence/outputs/production_<system>/stage012_evaluation.csv`` and the
per-molecule attrition table written by ``analyze_active_attrition.py``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "evidence" / "outputs" / "efficiency"

SYSTEMS = {"glp1r": "GLP-1R", "ghsr": "GHSR", "ntsr1": "NTSR1", "mdm2": "MDM2-p53"}
PRODUCTION_FRACTION = 0.05


def sweep_system(key: str, fractions: np.ndarray) -> pd.DataFrame:
    df = pd.read_csv(ROOT / "evidence" / "outputs" / f"production_{key}" / "stage012_evaluation.csv")
    n_input = len(df)
    n_actives = int((df.label == "active").sum())
    cand = df[df.gate == "candidate"].sort_values("cascade_score_pct", ascending=False)

    rows = []
    for frac in fractions:
        k = min(len(cand), int(np.ceil(frac * len(cand))))
        top = cand.head(k)
        n_act = int((top.label == "active").sum())
        rows.append(
            {
                "system": SYSTEMS[key],
                "shortlist_fraction": frac,
                "n_native_scored": k,
                "n_actives_retained": n_act,
                "active_retention": n_act / n_actives,
                # How much native scoring the cascade avoids relative to
                # scoring every input molecule.
                "compute_reduction": n_input / max(k, 1),
                "n_input": n_input,
                "n_actives_input": n_actives,
                "n_candidates": len(cand),
            }
        )
    return pd.DataFrame(rows)


def chemistry_gate_cost() -> pd.DataFrame:
    """Actives recoverable by relaxing the strict structural-alert gate."""
    path = ROOT / "evidence" / "outputs" / "attrition" / "active_attrition_per_molecule.csv"
    if not path.exists():
        raise SystemExit("Run evidence/analyze_active_attrition.py first.")
    att = pd.read_csv(path)
    rows = []
    for sysname in SYSTEMS.values():
        sub = att[att.system == sysname]
        total = len(sub)
        retained = int((sub.stage_lost == "-").sum())
        alerts = int((sub.reason_class == "structural alert").sum())
        rows.append(
            {
                "system": sysname,
                "n_actives": total,
                "retained_strict": retained,
                "lost_to_structural_alerts": alerts,
                "retention_strict": retained / total,
                # Upper bound: alert-flagged actives would still have to clear
                # Stages 1-3, so this brackets rather than predicts the gain.
                "max_retained_if_alerts_readmitted": retained + alerts,
            }
        )
    return pd.DataFrame(rows)


def latex_operating_points(sweep: pd.DataFrame) -> str:
    """Table contrasting the production setting with a relaxed shortlist."""
    lines = [
        r"\begin{table}[!t]",
        r"\caption{Efficiency-retention trade-off of the Stage-3 shortlist. For each system the",
        r"table contrasts the production setting (top 5\% of Stage-1 survivors) with a fully",
        r"permissive shortlist. ``Native scored'' is the number of molecules reaching the",
        r"expensive terminal stage out of the 30\,000-molecule screen; ``retained'' counts",
        r"in-domain actives. The trade is favorable where the two systems differ little in",
        r"retention and greatly in cost, and unfavorable for GHSR, where the 5\% cut discards",
        r"two thirds of the actives that the gate had already admitted.}",
        r"\label{tab:efficiency}",
        r"\centering",
        r"\begin{tabular}{lrrrr}",
        r"\hline",
        r" & \multicolumn{2}{c}{Production (5\%)} & \multicolumn{2}{c}{Permissive (100\%)} \\",
        r"System & Native scored & Retained & Native scored & Retained \\",
        r"\hline",
    ]
    for sysname in SYSTEMS.values():
        sub = sweep[sweep.system == sysname]
        prod = sub.iloc[(sub.shortlist_fraction - PRODUCTION_FRACTION).abs().argmin()]
        full = sub.iloc[sub.shortlist_fraction.argmax()]
        lines.append(
            f"{sysname} & {int(prod.n_native_scored):,} & "
            f"{int(prod.n_actives_retained)}/{int(prod.n_actives_input)} & "
            f"{int(full.n_native_scored):,} & "
            f"{int(full.n_actives_retained)}/{int(full.n_actives_input)} \\\\"
        )
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--points", type=int, default=60)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    fractions = np.unique(
        np.concatenate([np.geomspace(0.005, 1.0, args.points), [PRODUCTION_FRACTION]])
    )
    sweep = pd.concat([sweep_system(k, fractions) for k in SYSTEMS], ignore_index=True)
    sweep.to_csv(OUT / "shortlist_sweep.csv", index=False)

    cost = chemistry_gate_cost()
    cost.to_csv(OUT / "chemistry_gate_cost.csv", index=False)
    (OUT / "efficiency_table.tex").write_text(latex_operating_points(sweep), encoding="utf-8")

    print("=== Shortlist sweep (actives retained / native-scoring calls) ===")
    for sysname in SYSTEMS.values():
        sub = sweep[sweep.system == sysname]
        print(f"\n{sysname} (candidates after Stages 0-2: {int(sub.n_candidates.iloc[0]):,})")
        for frac in (0.01, 0.05, 0.10, 0.25, 1.00):
            r = sub.iloc[(sub.shortlist_fraction - frac).abs().argmin()]
            mark = "  <- production" if abs(frac - PRODUCTION_FRACTION) < 1e-9 else ""
            print(
                f"   {frac * 100:6.1f}%  native-scored {int(r.n_native_scored):6,}  "
                f"actives {int(r.n_actives_retained):3d}/{int(r.n_actives_input)}"
                f"  ({r.active_retention * 100:5.1f}%){mark}"
            )

    print("\n=== Cost of the strict structural-alert gate ===")
    print(cost.to_string(index=False))
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    # Replaced inferred native admission with executed paired-policy results.
    from analyze_absolute_floor import main as current_main
    current_main()
