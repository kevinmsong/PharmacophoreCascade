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
    "ablate_remove_stage2": "Hotspot-score ordering",
    "ablate_remove_stage1": "Pair-overlap ordering",
    "ablate_remove_stage3": "Cascade score only within native pool",
    "ablate_native_top_stage3_only": "Observed Stage-3 top-5,000 subset",
    "ablate_alternative_tiebreak": "Native-fit RMSD ordering",
}
ORDER = [
    "ablate_remove_stage2", "ablate_remove_stage1", "ablate_remove_stage3",
    "ablate_native_top_stage3_only", "ablate_alternative_tiebreak",
]


def _ranked(df: pd.DataFrame, id_col: str, score_col: str, ascending: bool) -> pd.Series:
    d = df[[id_col, score_col]].dropna()
    d = d.sort_values(score_col, ascending=ascending, kind="stable")
    return pd.Series(np.arange(1, len(d) + 1), index=d[id_col].to_numpy())


RESULTS = ROOT / "results" / "absolute_floor_1000"


def build_rankings() -> tuple[pd.Series, dict[str, pd.Series]]:
    """Rebuild the canonical and ablated rankings from the released tables.

    The cached ``native_scored.parquet`` holds only the final top-1,000
    ligands, so the full per-ligand native summary is read instead: all 4,997
    successfully native-scored ligands, which is the population the ablations
    are defined over.
    """
    stage3 = pd.read_csv(RESULTS / "screening_full_1M_floor1000.csv")
    native = pd.read_csv(RESULTS / "screening_full_1M_floor1000_native_scored_top5000.csv")
    native = native[native.native_weighted_coverage_pct.notna()].sort_values('final_rank')
    final = pd.read_csv(RESULTS / "top_1000_glp1_mimetics_full_1M_floor1000_native_final.csv")

    canonical = pd.Series(
        np.arange(1, len(final) + 1),
        index=final.sort_values("final_rank")["zinc_id"].to_numpy(),
    )

    stage3 = stage3.assign(ligand_id=stage3["zinc_id"])
    ns = native.assign(ligand_id=native["zinc_id"])
    ns = ns.sort_values("final_rank", kind='stable')
    ns = ns.drop_duplicates("ligand_id")
    rank_col = "stage3_screen_rank"

    ablations = {
        "ablate_remove_stage1": _ranked(stage3, "ligand_id", "pair_hash_overlap_pct", False),
        "ablate_remove_stage2": _ranked(stage3, "ligand_id", "hotspot_weighted_pct", False),
        "ablate_remove_stage3": _ranked(
            ns.sort_values("cascade_score_pct", ascending=False).drop_duplicates("ligand_id"),
            "ligand_id", "cascade_score_pct", False),
        "ablate_native_top_stage3_only": _ranked(
            ns[ns[rank_col]<=5000], "ligand_id", "final_rank", True),
        "ablate_alternative_tiebreak": _ranked(ns, "ligand_id", "native_fit_rmsd_angstrom", True),
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
        tau_1000 = kendalltau(canonical, cand.reindex(canonical.index).fillna(bottom)).correlation

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
            "tau_1000": tau_1000,
            "rho_shared": rho_shared,
            "n_common": len(support),
            "tau_common": tau_common,
            "rho_common": rho_common,
            "J10": jac(10), "J50": jac(50), "J100": jac(100), "J500": jac(500),
        })
    return pd.DataFrame(rows)


def latex(df: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[tbp]",
        r"\centering",
        r"\caption{\textbf{Ranking perturbations against the final 1,000-ligand ranking.}",
        r"These analyses reorder the score tables; they are not end-to-end gate-removal",
        r"experiments. All comparisons use the same population, so $\tau$ is comparable across",
        r"rows: $\tau_{1000}$ correlates the two orderings over the top 1,000, and",
        r"$\tau_\mathrm{common}$ repeats the comparison over every ligand either ranking",
        r"reaches, with unranked ligands tied below all ranked ones. $J_k$ is the Jaccard",
        r"overlap of the top-$k$ sets. The Stage-3 subset contains only successfully native-scored",
        r"ligands whose Stage-3 rank is at most 5,000; it does not estimate scores for",
        r"unobserved members of a different native pool.}",
        r"\label{tab:ablation}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular}{@{}lrrrrr@{}}",
        r"\toprule",
        r"Ranking perturbation & $n_{\rm shared}$ & $\tau_{1000}$ & $\tau_\mathrm{common}$"
        r" & $J_{10}$ & $J_{100}$ \\",
        r"\midrule",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"{r.ablation} & {int(r.n_shared):,} & {r.tau_1000:+.3f} & {r.tau_common:+.3f}"
            f" & {r.J10:.3f} & {r.J100:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
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
