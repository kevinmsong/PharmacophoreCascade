#!/usr/bin/env python3
"""Reproduce the retrospective enrichment metrics reported in the manuscript
directly from the released per-molecule benchmark scores.

This script is fully self-contained: it depends only on numpy, pandas, and
scikit-learn (all pinned in requirements.txt) and reads only the public
machine-readable benchmark tables under ``evidence/data/machine_readable/``.
It does NOT require the production screening engine. Each table lists, for
every benchmarked molecule, its active/decoy label and the rank and score
assigned by each method (full cascade, native-only, Stage-3-only, and the
single-pass 3D pharmacophore baseline), so every ROC-AUC, PR-AUC, EF, BEDROC,
and top-k recovery value in the paper can be recomputed from these columns.

Usage:
    python reproduce/reproduce_benchmarks.py
    python reproduce/reproduce_benchmarks.py --system glp1r

Metric definitions match the Methods section:
  * ROC-AUC / PR-AUC: global ranking quality (scikit-learn).
  * EF k%: fold enrichment of actives in the top ceil(k% * N) ranked molecules
    relative to the overall active prevalence.
  * BEDROC (alpha = 20): early-recognition-weighted score (Truchon & Bayly,
    J. Chem. Inf. Model. 2007, 47, 488-508).
  * top-k recovery: number of actives ranked within the top k positions.
Molecules that failed an intermediate stage are ranked below scored molecules
(as in the manuscript), so ordering by the released rank reproduces that.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "evidence" / "data" / "machine_readable"

SYSTEMS = ["glp1r", "ghsr", "ntsr1", "mdm2"]
METHODS = ["full_cascade", "native_only", "stage3_only", "standard_3d_pharmacophore"]


def bedroc(rank_positions: np.ndarray, n_total: int, n_actives: int, alpha: float = 20.0) -> float:
    """Truchon-Bayly BEDROC from 1-indexed active rank positions."""
    if n_actives == 0 or n_actives == n_total:
        return float("nan")
    ra = n_actives / n_total
    s = np.sum(np.exp(-alpha * rank_positions / n_total))
    rie = s * (1.0 / n_actives) * n_total * (math.exp(alpha / n_total) - 1.0) / (1.0 - math.exp(-alpha))
    factor = ra * math.sinh(alpha / 2.0) / (math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra))
    return rie * factor + 1.0 / (1.0 - math.exp(alpha * (1.0 - ra)))


def enrichment_factor(ranks: np.ndarray, labels: np.ndarray, pct: float) -> float:
    """EF k% = (fraction of actives recovered in the top ceil(k% * N)) / (k/100),
    matching the convention used in the manuscript and evidence summaries."""
    n = len(labels)
    n_act = int(labels.sum())
    m = max(1, math.ceil(pct / 100.0 * n))
    recovered = labels[np.argsort(ranks)][:m].sum()
    return (recovered / n_act) / (pct / 100.0)


def metrics_for_method(df: pd.DataFrame, method: str) -> dict:
    labels = (df["label"] == "active").to_numpy().astype(int)
    ranks = df[f"{method}_rank"].to_numpy().astype(float)
    score = -ranks  # rank 1 (top) -> highest score; ties impossible (unique ranks)
    order = np.argsort(ranks)
    sorted_labels = labels[order]
    active_positions = np.where(sorted_labels == 1)[0] + 1  # 1-indexed
    n, n_act = len(labels), int(labels.sum())
    return {
        "method": method,
        "ROC-AUC": roc_auc_score(labels, score),
        "PR-AUC": average_precision_score(labels, score),
        "EF1%": enrichment_factor(ranks, labels, 1.0),
        "EF5%": enrichment_factor(ranks, labels, 5.0),
        "BEDROC": bedroc(active_positions, n, n_act, 20.0),
        "top10_rec": int((ranks[labels == 1] <= 10).sum()),
        "top25_rec": int((ranks[labels == 1] <= 25).sum()),
    }


def run_system(system: str) -> pd.DataFrame:
    path = DATA / f"{system}_benchmark_scored.csv"
    df = pd.read_csv(path)
    rows = [metrics_for_method(df, m) for m in METHODS if f"{m}_rank" in df.columns]
    out = pd.DataFrame(rows).set_index("method")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--system", choices=SYSTEMS, help="single system (default: all)")
    args = ap.parse_args()
    systems = [args.system] if args.system else SYSTEMS
    for system in systems:
        print(f"\n=== {system.upper()} (recomputed from {system}_benchmark_scored.csv) ===")
        out = run_system(system)
        with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
            print(out.to_string())


if __name__ == "__main__":
    main()
