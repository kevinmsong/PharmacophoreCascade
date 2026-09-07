#!/usr/bin/env python3
"""Equivalence testing and paired significance tests for the revision.

Two analyses that the manuscript previously reported informally:

1. **Cascade versus native-only equivalence (TOST).** The manuscript
   previously read a large two-sided p-value as evidence that the two methods
   perform equally. A non-significant difference is not a demonstration of
   equivalence, so this script runs a proper two one-sided tests (TOST)
   procedure on the grouped bootstrap distribution of the paired ROC-AUC
   difference, against a pre-specified margin. Equivalence is declared only if
   the 90% confidence interval of the difference lies entirely inside
   +/- margin, which is the interval-inclusion form of TOST.

2. **Docking active-state preference (Wilcoxon signed-rank).** The manuscript
   previously reported only a median active-minus-inactive affinity
   difference. This script runs the paired test across the ten ligands.

The bootstrap resamples each active together with its matched decoys, matching
the procedure used for the manuscript's confidence intervals, and reuses the
metric definitions from ``reproduce/reproduce_benchmarks.py``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from itertools import product

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "reproduce"))
from reproduce_benchmarks import bedroc, enrichment_factor  # noqa: E402

from sklearn.metrics import roc_auc_score  # noqa: E402
from tie_aware_metrics import evaluate as evaluate_ties, clustered_counts

DATA = ROOT / "evidence" / "data" / "machine_readable"
OUT = ROOT / "evidence" / "outputs" / "equivalence"

SYSTEMS = {"glp1r": "GLP-1R", "ghsr": "GHSR", "ntsr1": "NTSR1", "mdm2": "MDM2-p53"}

#: Pre-specified equivalence margin on ROC-AUC. 0.05 is a conventional choice
#: for "practically indistinguishable ranking quality" and is fixed here before
#: the tests are run.
MARGIN = 0.05


def _metric(sub: pd.DataFrame, method: str, name: str) -> float:
    frame=sub[['ligand_id','label',f'{method}_rank',f'{method}_score',f'{method}_status']].rename(
        columns={f'{method}_rank':'rank',f'{method}_score':'score',f'{method}_status':'status'})
    return evaluate_ties(frame)[name]


def grouped_bootstrap_delta(
    df: pd.DataFrame,
    method_a: str,
    method_b: str,
    metric: str,
    n_boot: int = 5000,
    seed: int = 42,
) -> np.ndarray:
    """Bootstrap distribution of ``metric(a) - metric(b)``.

    Each active is resampled together with the decoys matched to it, which
    respects the benchmark's matching structure and keeps the resulting
    intervals honest about the small number of actives.
    """
    results=[]
    for method in [method_a,method_b]:
        frame=df[['ligand_id','label',f'{method}_rank',f'{method}_score',f'{method}_status']].rename(
            columns={f'{method}_rank':'rank',f'{method}_score':'score',f'{method}_status':'status'})
        _,boot=clustered_counts(frame,df,draws=n_boot,seed=seed)
        results.append(boot[metric])
    return results[0]-results[1]


def tost(deltas: np.ndarray, margin: float) -> dict:
    """Interval-inclusion TOST on a bootstrap delta distribution.

    Equivalence at level alpha = 0.05 holds when the 90% interval lies wholly
    inside (-margin, +margin). The two one-sided bootstrap p-values are the
    tail masses beyond each equivalence bound.
    """
    lo90, hi90 = np.percentile(deltas, [5, 95])
    lo95, hi95 = np.percentile(deltas, [2.5, 97.5])
    p_lower = float(np.mean(deltas <= -margin))  # H0: delta <= -margin
    p_upper = float(np.mean(deltas >= margin))  # H0: delta >= +margin
    return {
        "delta_mean": float(np.mean(deltas)),
        "ci90_low": float(lo90),
        "ci90_high": float(hi90),
        "ci95_low": float(lo95),
        "ci95_high": float(hi95),
        "margin": margin,
        "p_tost": float(max(p_lower, p_upper)),
        "equivalent": bool(lo90 > -margin and hi90 < margin),
    }


def run_equivalence(n_boot: int, seed: int) -> pd.DataFrame:
    rows = []
    for key, label in SYSTEMS.items():
        df = pd.read_csv(DATA / f"{key}_benchmark_scored.csv")
        if "native_only_rank" not in df.columns:
            continue
        for metric in ("roc_auc", "bedroc"):
            deltas = grouped_bootstrap_delta(
                df, "full_cascade", "native_only", metric, n_boot=n_boot, seed=seed
            )
            res = tost(deltas, MARGIN)
            obs_a = _metric(df, "full_cascade", metric)
            obs_b = _metric(df, "native_only", metric)
            rows.append(
                {
                    "system": label,
                    "n_actives": int((df.label == "active").sum()),
                    "metric": metric,
                    "full_cascade": obs_a,
                    "native_only": obs_b,
                    "observed_delta": obs_a - obs_b,
                    **res,
                }
            )
    return pd.DataFrame(rows)


def run_docking_test() -> dict:
    """Paired Wilcoxon signed-rank on the active-minus-inactive Vina deltas."""
    path = ROOT / "evidence" / "outputs" / "revision" / "docking_top10_summary.csv"
    all_rows = pd.read_csv(path).sort_values("final_rank")
    df=all_rows.replace([np.inf,-np.inf],np.nan).dropna(subset=['best_active','best_inactive','active_pref'])
    d = df["active_pref"].to_numpy(dtype=float)
    if len(d)==0:raise RuntimeError('No ligands have successful docking against both states')
    # Exhaustive conditional sign flips give the exact signed-rank null even
    # when rounded Vina scores produce ties. Drop zero differences explicitly.
    nonzero=d[np.abs(d)>1e-12]
    ranks=rankdata(np.round(np.abs(nonzero),10),method='average')
    positive=float(ranks[nonzero>0].sum());total=float(ranks.sum())
    stat=min(positive,total-positive)
    null_positive=np.array([np.dot(bits,ranks) for bits in product([0,1],repeat=len(ranks))])
    p_two=float((np.minimum(null_positive,total-null_positive)<=stat+1e-12).mean())
    p_less=float((null_positive<=positive+1e-12).mean())
    return {
        "n": int(len(d)),
        "n_selected":int(len(all_rows)),"n_zero_differences":int(len(d)-len(nonzero)),
        "test_method":"exact signed-rank permutation; exhaustive sign flips, average ranks for ties, zero differences excluded",
        "median_delta": float(np.median(d)),
        "iqr_low": float(np.percentile(d, 25)),
        "iqr_high": float(np.percentile(d, 75)),
        "n_favoring_active": int((d < 0).sum()),
        "wilcoxon_W": float(stat),
        "p_two_sided": float(p_two),
        "p_one_sided_less": float(p_less),
        "source": str(path.relative_to(ROOT)),
        "ligands": df["ligand_id"].tolist(),
    }


def latex_tost(eq: pd.DataFrame) -> str:
    sub = eq[eq.metric == "roc_auc"]
    lines = [
        r"\begin{table}[tbp]",
        r"\centering",
        r"\caption{\textbf{Equivalence testing of the full cascade against the native-only",
        r"baseline (ROC-AUC).} Two one-sided tests were run on the grouped-bootstrap",
        r"distribution of the paired difference against a pre-specified margin of $\pm$"
        + f"{MARGIN:.2f}" + r".",
        r"Equivalence requires the 90\% interval to lie wholly inside that margin. A",
        r"non-significant difference alone must not be read as equivalence.}",
        r"\label{tab:tost}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{6pt}",
        r"\begin{tabular}{@{}lrrrl@{}}",
        r"\toprule",
        r"System & $n$ & $\Delta$ROC-AUC & 90\% CI & Equivalent? \\",
        r"\midrule",
    ]
    for _, r in sub.iterrows():
        eqv = "yes" if r.equivalent else "no"
        # The manuscript sets the protein-protein pair with an en dash throughout.
        system = r.system.replace("MDM2-p53", "MDM2--p53")
        lines.append(
            f"{system} & {int(r.n_actives)} & {r.observed_delta:+.3f} & "
            f"[{r.ci90_low:+.3f}, {r.ci90_high:+.3f}] & {eqv} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bootstrap", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    eq = run_equivalence(args.bootstrap, args.seed)
    eq.to_csv(OUT / "equivalence_tost.csv", index=False)
    (OUT / "tost_table.tex").write_text(latex_tost(eq), encoding="utf-8")

    print("=== TOST: full cascade vs native-only ===")
    print(f"pre-specified margin: +/- {MARGIN}")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(
            eq[
                ["system", "n_actives", "metric", "observed_delta", "ci90_low", "ci90_high",
                 "p_tost", "equivalent"]
            ].to_string(index=False, float_format=lambda v: f"{v:.4f}")
        )

    dock = run_docking_test()
    (OUT / "docking_wilcoxon.json").write_text(json.dumps(dock, indent=2), encoding="utf-8")
    print("\n=== Docking: paired Wilcoxon signed-rank (active - inactive) ===")
    print(f"n = {dock['n']}, favoring active = {dock['n_favoring_active']}/{dock['n']}")
    print(f"median delta = {dock['median_delta']:.2f} kcal/mol "
          f"(IQR {dock['iqr_low']:.2f} to {dock['iqr_high']:.2f})")
    print(f"W = {dock['wilcoxon_W']:.1f}, exact two-sided p = {dock['p_two_sided']:.5f}")
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    main()
