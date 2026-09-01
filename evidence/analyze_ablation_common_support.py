#!/usr/bin/env python3
"""Recompute the component ablations on a common support.

The published ablation table reported Kendall tau over whatever ligands each
ablated ranking happened to share with the canonical output. Because that
shared set ranged from 52 to 1,000 ligands depending on the ablation, the tau
values were not comparable with one another: an ablation restricted to the
native pool was scored on a different, much smaller population than one that
reranks the whole Stage-3 table.

This script recomputes every comparison twice:

* ``tau_shared``   -- the original convention, over the intersection only.
* ``tau_common``   -- over a fixed common support (the canonical top-1,000
  plus every ligand any ablation ranks), with ligands an ablation does not
  rank assigned a tied rank below all ranked ligands.

``tau_common`` uses the same population for every row, so the rows can be
compared. Both are reported so the change in convention is visible.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "evidence" / "outputs" / "cache"
OUT = ROOT / "evidence" / "outputs" / "ablation_common"

LABELS = {
    "ablate_remove_stage2": "Remove Stage 2 (hotspot only)",
    "ablate_remove_stage1": "Remove Stage 1 (pair-hash only)",
    "ablate_remove_stage3": "Cascade score only within native pool",
    "ablate_native_top_stage3_only": "Native pool: top Stage-3 only",
    "ablate_alternative_tiebreak": "Alternative final tie-break",
}
ORDER = [
    "ablate_remove_stage2", "ablate_remove_stage1", "ablate_remove_stage3",
    "ablate_native_top_stage3_only", "ablate_alternative_tiebreak",
]


def _ranked(df: pd.DataFrame, id_col: str, score_col: str, ascending: bool) -> pd.Series:
    d = df[[id_col, score_col]].dropna()
    d = d.sort_values(score_col, ascending=ascending, kind="stable")
    return pd.Series(np.arange(1, len(d) + 1), index=d[id_col].to_numpy())


RESULTS = ROOT / "results"
BUNDLE = RESULTS / "screening_full_1M_topological_hashed_native_terminal_bundle" / "analysis"


def build_rankings() -> tuple[pd.Series, dict[str, pd.Series]]:
    """Rebuild the canonical and ablated rankings from the released tables.

    The cached ``native_scored.parquet`` holds only the final top-1,000
    ligands, so the full per-ligand native summary is read instead: all 4,997
    successfully native-scored ligands, which is the population the ablations
    are defined over.
    """
    stage3 = pd.read_parquet(CACHE / "stage3.parquet")
    native = pd.read_csv(BUNDLE / "ligand_best_native_mapping_summary.csv")
    final = pd.read_csv(RESULTS / "top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv")

    canonical = pd.Series(
        np.arange(1, len(final) + 1),
        index=final.sort_values("final_rank")["zinc_id"].to_numpy(),
    )

    stage3 = stage3.assign(ligand_id=stage3["zinc_id"])
    ns = native.assign(ligand_id=native["zinc_id"])
    ns = ns.sort_values("native_weighted_coverage_pct", ascending=False)
    ns = ns.drop_duplicates("ligand_id")
    rank_col = "source_input_rank" if "source_input_rank" in ns.columns else "shortlist_rank"

    ablations = {
        "ablate_remove_stage1": _ranked(stage3, "ligand_id", "pair_hash_overlap_pct", False),
        "ablate_remove_stage2": _ranked(stage3, "ligand_id", "hotspot_weighted_pct", False),
        "ablate_remove_stage3": _ranked(
            ns.sort_values("cascade_score_pct", ascending=False).drop_duplicates("ligand_id"),
            "ligand_id", "cascade_score_pct", False),
        "ablate_native_top_stage3_only": _ranked(
            ns.nsmallest(5000, rank_col), "ligand_id", "native_weighted_coverage_pct", False),
        "ablate_alternative_tiebreak": _ranked(ns, "ligand_id", "fit_rmsd_angstrom", True),
    }
    return canonical, ablations


def compare(canonical: pd.Series, ablations: dict[str, pd.Series]) -> pd.DataFrame:
    # Common support: every ligand the canonical output or any ablation ranks.
    support = set(canonical.index)
    for s in ablations.values():
        support |= set(s.index)
    support = sorted(support)
    bottom = len(support) + 1

    ref_common = canonical.reindex(support).fillna(bottom).to_numpy()

    rows = []
    for name in ORDER:
        cand = ablations[name]
        shared = canonical.index.intersection(cand.index)
        tau_shared = kendalltau(canonical.reindex(shared), cand.reindex(shared)).correlation
        rho_shared = spearmanr(canonical.reindex(shared), cand.reindex(shared)).correlation

        cand_common = cand.reindex(support).fillna(bottom).to_numpy()
        tau_common = kendalltau(ref_common, cand_common).correlation
        rho_common = spearmanr(ref_common, cand_common).correlation

        ref_ids = canonical.sort_values().index.tolist()
        cand_ids = cand.sort_values().index.tolist()

        def jac(k: int) -> float:
            a, b = set(ref_ids[:k]), set(cand_ids[:k])
            return len(a & b) / len(a | b) if a | b else 0.0

        rows.append({
            "ablation": LABELS[name],
            "n_shared": len(shared),
            "tau_shared": tau_shared,
            "rho_shared": rho_shared,
            "n_common": len(support),
            "tau_common": tau_common,
            "rho_common": rho_common,
            "J10": jac(10), "J50": jac(50), "J100": jac(100), "J500": jac(500),
        })
    return pd.DataFrame(rows)


def latex(df: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[!t]",
        r"\caption{Component ablations against the canonical final 1\,000-ligand ranking.",
        r"Each ablation replaces one pipeline element with a simplified alternative. All",
        r"comparisons are computed over the same population, so $\tau$ is comparable across",
        r"rows: $\tau_{1000}$ correlates the two orderings over the canonical top 1\,000, and",
        r"$\tau_\mathrm{common}$ repeats the comparison over every ligand either ranking",
        r"reaches, with unranked ligands tied below all ranked ones. $J_k$ is the Jaccard",
        r"overlap of the top-$k$ sets. An earlier version of this table reported $\tau$ over",
        r"shared sets of only 52--80 ligands for the three native-branch rows; those small",
        r"sets came from an incomplete export of the native-scored population and are",
        r"superseded here.}",
        r"\label{tab:ablation}",
        r"\centering",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular}{lrrrrr}",
        r"\hline",
        r"Ablation & $n$ & $\tau_{1000}$ & $\tau_\mathrm{common}$"
        r" & $J_{10}$ & $J_{100}$ \\",
        r"\hline",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"{r.ablation} & {int(r.n_shared):,} & {r.tau_shared:+.3f} & {r.tau_common:+.3f}"
            f" & {r.J10:.3f} & {r.J100:.3f} \\\\"
        )
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    canonical, ablations = build_rankings()
    df = compare(canonical, ablations)
    df.to_csv(OUT / "ablation_common_support.csv", index=False)
    (OUT / "ablation_table.tex").write_text(latex(df), encoding="utf-8")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nCommon support: {int(df.n_common.iloc[0]):,} ligands")
    print(f"Wrote outputs to {OUT}")


if __name__ == "__main__":
    main()
