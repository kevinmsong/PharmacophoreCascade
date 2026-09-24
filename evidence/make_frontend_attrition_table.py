#!/usr/bin/env python3
"""Build the Supporting Information table of front-end attrition in the
retrospective benchmarks.

The four retrospective methods form two pairs that share a scoring function and
differ only in which molecules reach it: the full cascade against native-only
over the terminal native score, and Stage-3-only against single-pass 3D over the
three-dimensional pharmacophore score. This table reports, per system, how many
labeled molecules each front end removes and how closely the two members of each
pair then agree. It is the evidence behind the statement that the matching
MDM2--p53 rows of the cross-system table are structural rather than coincidental.

Counts come from the deposited per-molecule records; the rank agreements are the
published pairwise statistics for the same runs.

Writes ``frontend_attrition_table.tex`` next to the other table sources.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCORED = ROOT / "evidence" / "data" / "machine_readable"
PAIRWISE = ROOT / "evidence" / "outputs"

SYSTEMS = [
    ("glp1r", "GLP-1R"),
    ("ghsr", "GHSR"),
    ("ntsr1", "NTSR1"),
    ("mdm2", "MDM2--p53"),
]

# The two comparators in each pair, ordered as they appear in the deposited
# pairwise files.
PAIRS = [
    ("full_cascade", "native_only"),
    ("stage3_only", "standard_3d_pharmacophore"),
]


def attrition(system: str) -> dict[str, int]:
    """Count Stage-0 and Stage-1 removals from the full-cascade status column."""
    with open(SCORED / f"{system}_benchmark_scored.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    status = [r["full_cascade_status"] for r in rows]
    return {
        "n": len(rows),
        "stage0": status.count("property_filtered"),
        "stage1": status.count("hotspot_filtered"),
        "reaching": status.count("scored"),
    }


def rank_agreement(system: str) -> dict[tuple[str, str], float]:
    """Read the published Kendall tau for each comparator pair."""
    path = PAIRWISE / f"benchmark_{system}_full" / "benchmark_pairwise.csv"
    with open(path, newline="") as handle:
        table = {(r["method_a"], r["method_b"]): float(r["kendall_tau"])
                 for r in csv.DictReader(handle)}
    out = {}
    for a, b in PAIRS:
        if (a, b) in table:
            out[(a, b)] = table[(a, b)]
        else:
            out[(a, b)] = table[(b, a)]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "ACS_Omega_resubmission_R2"
                                                / "frontend_attrition_table.tex"))
    args = parser.parse_args()

    body = []
    for key, label in SYSTEMS:
        counts = attrition(key)
        taus = rank_agreement(key)
        body.append(
            "{label} & {n:,} & {s0} & {s1:,} & {reach:,} & {t1:.3f} & {t2:.3f} \\\\".format(
                label=label,
                n=counts["n"],
                s0=counts["stage0"],
                s1=counts["stage1"],
                reach=counts["reaching"],
                t1=taus[PAIRS[0]],
                t2=taus[PAIRS[1]],
            )
        )

    tex = r"""\begin{table}[H]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{\textbf{Front-end attrition in the retrospective benchmarks and its effect on the
paired comparators.} Stage 0 removes molecules outside the property limits after
standardization and Stage 1 removes molecules failing the hotspot gate; Stage 2 removes no
labeled molecule in any of these runs. Removed molecules are retained at the bottom of the
ranking under a recorded status rather than dropped. The two rightmost columns give the
Kendall rank correlation between the members of each comparator pair, which share a scoring
function and differ only in which molecules reach it. At MDM2--p53 the front end removes
nothing, both correlations are exactly 1.000, and the paired metrics in Table~2 of the
manuscript therefore coincide exactly.}
\label{tab:s16_frontend}
\begin{tabular}{@{}lrrrrrr@{}}
\toprule
 & Labeled & \multicolumn{2}{c}{Removed by} & Reaching & \multicolumn{2}{c}{Kendall $\tau$ within pair} \\
\cmidrule(lr){3-4}\cmidrule(lr){6-7}
System & molecules & Stage 0 & Stage 1 & scoring & Cascade/native & Stage 3/single-pass \\
\midrule
__BODY__
\bottomrule
\end{tabular}
\end{table}
""".replace("__BODY__", "\n".join(body))

    out_path = Path(args.output)
    out_path.write_text(tex, encoding="utf-8")
    print(f"wrote {out_path}")
    for line in body:
        print("  " + line)


if __name__ == "__main__":
    main()
