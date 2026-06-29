#!/usr/bin/env python3
"""Aggregate the per-run decoy-replicate benchmark summaries into mean +/- SD."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REP = ROOT / "evidence" / "outputs" / "decoy_replicates"
OUT = ROOT / "evidence" / "outputs" / "revision"
METRICS = ["roc_auc", "pr_auc", "bedroc", "ef_1pct", "top10_recovery"]
SYSTEMS = ["glp1r", "ghsr", "ntsr1", "mdm2"]
SEEDS = [1, 2, 3, 4, 5]


def main() -> None:
    rows = []
    for s in SYSTEMS:
        for seed in SEEDS:
            p = REP / f"{s}_seed{seed}" / "benchmark_summary.csv"
            if not p.exists():
                continue
            d = pd.read_csv(p).set_index("method")
            for method in d.index:
                rec = {"system": s, "seed": seed, "method": method}
                for m in METRICS:
                    rec[m] = float(d.loc[method, m]) if m in d.columns else float("nan")
                rows.append(rec)
    raw = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    raw.to_csv(OUT / "decoy_replicate_raw.csv", index=False)
    agg = raw.groupby(["system", "method"])[METRICS].agg(["mean", "std"]).round(4)
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()
    agg.to_csv(OUT / "decoy_replicate_robustness.csv", index=False)
    n = raw.groupby(["system"])["seed"].nunique().to_dict()
    print("replicate counts per system:", n)
    fc = agg[agg["method"] == "full_cascade"]
    print(fc[["system", "roc_auc_mean", "roc_auc_std", "ef_1pct_mean", "ef_1pct_std",
              "top10_recovery_mean", "top10_recovery_std"]].to_string(index=False))


if __name__ == "__main__":
    main()
