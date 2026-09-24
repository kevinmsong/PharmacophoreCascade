#!/usr/bin/env python3
"""Rebuild the decoy-replicate robustness table from the per-replicate runs.

Reads ``evidence/outputs/decoy_replicates/<system>_seed<n>/benchmark_summary.csv``
for the five seeds of each system and writes the mean and sample standard
deviation of the full cascade's ROC-AUC, PR-AUC, and BEDROC.

The replicate design holds a system's active set fixed and resamples only the
decoys, which is what isolates sensitivity to decoy selection. GLP-1R carries
10 actives, the ceiling for that system; the other three carry the same 50
actives as their primary benchmarks.

Writes ``decoy_robustness_table.tex`` next to the other table sources.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REP = ROOT / "evidence" / "outputs" / "decoy_replicates"

SYSTEMS = [("glp1r", "GLP-1R"), ("ghsr", "GHSR"),
           ("ntsr1", "NTSR1"), ("mdm2", "MDM2--p53")]
METRICS = ["roc_auc", "pr_auc", "bedroc"]
SEEDS = [1, 2, 3, 4, 5]

CAPTION_TEMPLATE = (
    r"\textbf{Robustness of the full cascade across five independent matched-decoy sets.} "
    r"Values are mean $\pm$ sample standard deviation over five runs per system. Within a "
    r"system the active set is held fixed and only the decoys are resampled, so the spread "
    r"reflects decoy selection alone. Each replicate pairs a system's actives with up to 30 "
    r"matched decoys per active: 10 actives and 300 decoys for GLP-1R, which is data-limited, "
    r"and 50 actives with __MDM2NOTE__ for the other three, the same active sets as their "
    r"primary benchmarks. Across those three systems no metric varies by more than "
    r"__SPREAD__ between draws, so the standard deviations shown round to 0.000. These "
    r"replicates draw decoys from the property-matched neighborhood rather than taking the "
    r"strictly nearest neighbors used for the primary libraries, so their decoys are less "
    r"tightly matched and their values sit slightly above the corresponding entries of "
    r"Table~2; the comparison here is between draws, not against the primary benchmarks. "
    r"Structural alerts are recorded rather than excluding."
)


def load(system: str) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        path = REP / f"{system}_seed{seed}" / "benchmark_summary.csv"
        if not path.exists():
            raise SystemExit(f"missing replicate summary: {path}")
        frame = pd.read_csv(path).query("method == 'full_cascade'")
        if len(frame) != 1:
            raise SystemExit(f"expected one full_cascade row in {path}")
        record = {"seed": seed, "n_ranked": int(frame.iloc[0]["n_ranked"])}
        record.update({m: float(frame.iloc[0][m]) for m in METRICS})
        rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "ACS_Omega_resubmission_R2"
                                                / "decoy_robustness_table.tex"))
    args = parser.parse_args()

    body = []
    spread = 0.0
    decoy_counts: list[int] = []
    for key, label in SYSTEMS:
        frame = load(key)
        sizes = sorted(set(frame["n_ranked"]))
        cells = [f"{frame[m].mean():.3f} $\\pm$ {frame[m].std(ddof=1):.3f}" for m in METRICS]
        body.append(f"{label} & " + " & ".join(cells) + r" \\")
        if key != "glp1r":
            # The 50-active systems are the ones the caption characterizes.
            spread = max(spread, max(frame[m].max() - frame[m].min() for m in METRICS))
            decoy_counts.extend(int(n) - 50 for n in sizes)
        print(f"{label:10s} n={sizes}  " +
              "  ".join(f"{m}={frame[m].mean():.3f}+/-{frame[m].std(ddof=1):.3f}"
                        for m in METRICS))

    # MDM2--p53 exhausts the matched pool for a few actives, so its replicates
    # carry slightly fewer than 1,500 decoys. Say so rather than rounding.
    low, high = min(decoy_counts), max(decoy_counts)
    mdm2_note = (f"{low:,} decoys" if low == high
                 else f"{low:,} to {high:,} decoys")
    # Round the bound up, never down: "no more than 0.015" would be false for
    # an observed spread of 0.0151.
    bound = math.ceil(spread * 1000) / 1000
    caption = (CAPTION_TEMPLATE
               .replace("__SPREAD__", f"{bound:.3f}")
               .replace("__MDM2NOTE__", mdm2_note))
    print(f"\nlargest spread across the 50-active systems: {spread:.4f}")
    print(f"decoys per 50-active replicate: {low:,} to {high:,}")

    tex = r"""\begin{table}[tbp]
\centering
\caption{__CAPTION__}
\label{tab:decoy_robustness}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lrrr@{}}
\toprule
System & ROC-AUC & PR-AUC & BEDROC \\
\midrule
__BODY__
\bottomrule
\end{tabular}
\end{table}
""".replace("__CAPTION__", caption).replace("__BODY__", "\n".join(body))

    out = Path(args.output)
    out.write_text(tex, encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
