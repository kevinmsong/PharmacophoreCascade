"""
Table-level ablation and sensitivity analysis.

All ablations reorder the existing result tables without re-running the cascade.
Each ablation produces an alternate ranked list which is then compared to the
canonical final_ranked order via Kendall tau, Spearman r, and Jaccard top-k.

Entry point: run_ablations(tables, ablation_config) -> AblationResult
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .metrics import jaccard_topk, kendall_tau, spearman_r

logger = logging.getLogger(__name__)

RANK_COLS = ["ligand_id", "rank", "score", "method"]


@dataclass
class AblationResult:
    ablation_df: pd.DataFrame        # per-ablation rank-stability metrics
    sensitivity_df: pd.DataFrame     # parameter sweep results
    rank_lists: Dict[str, pd.DataFrame] = field(default_factory=dict)  # method -> ranked df


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _ranked_df(df: pd.DataFrame, id_col: str, score_col: str,
               ascending: bool, method: str) -> pd.DataFrame:
    out = df[[id_col, score_col]].copy().rename(columns={id_col: "ligand_id", score_col: "score"})
    out = out.sort_values("score", ascending=ascending, kind="stable").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["method"] = method
    return out[RANK_COLS]


def _compare_to_reference(
    candidate: pd.DataFrame,
    reference: pd.DataFrame,
    top_k_values: Sequence[int],
    method: str,
) -> dict:
    """Kendall tau, Spearman r, Jaccard top-k vs reference ranking."""
    ref = reference.set_index("ligand_id")["rank"]
    cand = candidate.set_index("ligand_id")["rank"]
    shared = ref.index.intersection(cand.index)

    row: dict = {
        "method": method,
        "n_shared": len(shared),
        "kendall_tau": kendall_tau(ref.reindex(shared).values, cand.reindex(shared).values),
        "spearman_r": spearman_r(ref.reindex(shared).values, cand.reindex(shared).values),
    }

    ref_ids  = reference.sort_values("rank")["ligand_id"].tolist()
    cand_ids = candidate.sort_values("rank")["ligand_id"].tolist()
    for k in top_k_values:
        row[f"jaccard_top{k}"] = jaccard_topk(ref_ids, cand_ids, k)

    return row


# ---------------------------------------------------------------------------
# Named ablations
# ---------------------------------------------------------------------------

def ablate_remove_stage1(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Remove Stage-1 gate: cascade_score = pair_hash_overlap_pct only."""
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]
    return _ranked_df(df, "ligand_id", "pair_hash_overlap_pct",
                      ascending=False, method="ablate_remove_stage1")


def ablate_remove_stage2(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Remove Stage-2 gate: cascade_score = hotspot_weighted_pct only."""
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]
    return _ranked_df(df, "ligand_id", "hotspot_weighted_pct",
                      ascending=False, method="ablate_remove_stage2")


def ablate_remove_stage3(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Remove Stage-3: rank the native pool by cascade_score instead of native coverage."""
    ns = tables["native_scored"].copy()
    # Deduplicate to one row per ligand
    ns = ns.sort_values("cascade_score_pct", ascending=False).drop_duplicates("ligand_id")
    return _ranked_df(ns, "ligand_id", "cascade_score_pct",
                      ascending=False, method="ablate_remove_stage3")


def ablate_native_top_stage3_only(
    tables: Dict[str, pd.DataFrame],
    n: int = 5000,
) -> pd.DataFrame:
    """Native pool = top-N by stage3_screen_rank only (ignore breadth/native-support quotas)."""
    ns = tables["native_scored"].copy()
    # Use source_input_rank as stage3 rank within the native-scored pool
    rank_col = "source_input_rank" if "source_input_rank" in ns.columns else "shortlist_rank"
    top_n = ns.nsmallest(n, rank_col)
    top_n = (
        top_n.sort_values("native_weighted_coverage_pct", ascending=False)
        .drop_duplicates("ligand_id")
    )
    return _ranked_df(top_n, "ligand_id", "native_weighted_coverage_pct",
                      ascending=False, method="ablate_native_top_stage3_only")


def ablate_alternative_tiebreak(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Alternative final ranking: RMSD primary (ascending), then native coverage (descending)."""
    ns = tables["native_scored"].copy()
    ns = (
        ns.sort_values(
            ["fit_rmsd_angstrom", "native_weighted_coverage_pct"],
            ascending=[True, False],
            kind="stable",
        )
        .drop_duplicates("ligand_id")
        .reset_index(drop=True)
    )
    ns["rank"] = np.arange(1, len(ns) + 1)
    ns["score"] = -ns["fit_rmsd_angstrom"]
    ns["method"] = "ablate_alternative_tiebreak"
    return ns[["ligand_id", "rank", "score", "method"]]


# ---------------------------------------------------------------------------
# Sensitivity: cascade weight sweep
# ---------------------------------------------------------------------------

def sweep_cascade_weights(
    tables: Dict[str, pd.DataFrame],
    weight_pairs: List[Tuple[float, float]],
    top_k_values: Sequence[int],
) -> pd.DataFrame:
    """Vary (hotspot_weight, pair_hash_weight) and measure rank stability."""
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]

    # Baseline: current published weights (0.4, 0.6)
    baseline_w = (0.4, 0.6)
    df["baseline_score"] = (
        baseline_w[0] * df["hotspot_weighted_pct"]
        + baseline_w[1] * df["pair_hash_overlap_pct"]
    )
    baseline_ranked = df.sort_values("baseline_score", ascending=False).reset_index(drop=True)
    baseline_ranked["rank"] = np.arange(1, len(baseline_ranked) + 1)
    baseline_ranked = baseline_ranked[["ligand_id", "rank"]].rename(columns={"rank": "baseline_rank"})

    records = []
    for w_h, w_p in weight_pairs:
        df["alt_score"] = w_h * df["hotspot_weighted_pct"] + w_p * df["pair_hash_overlap_pct"]
        alt_sorted = df.sort_values("alt_score", ascending=False).reset_index(drop=True)
        alt_sorted["alt_rank"] = np.arange(1, len(alt_sorted) + 1)
        merged = baseline_ranked.merge(alt_sorted[["ligand_id", "alt_rank"]], on="ligand_id")

        row: dict = {
            "hotspot_weight": w_h,
            "pair_hash_weight": w_p,
            "kendall_tau": kendall_tau(merged["baseline_rank"].values, merged["alt_rank"].values),
            "spearman_r": spearman_r(merged["baseline_rank"].values, merged["alt_rank"].values),
        }

        base_ids = baseline_ranked.sort_values("baseline_rank")["ligand_id"].tolist()
        alt_ids  = alt_sorted.sort_values("alt_rank")["ligand_id"].tolist()
        for k in top_k_values:
            row[f"jaccard_top{k}"] = jaccard_topk(base_ids, alt_ids, k)

        records.append(row)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Sensitivity: shortlist fraction sweep
# ---------------------------------------------------------------------------

def sweep_shortlist_fraction(
    tables: Dict[str, pd.DataFrame],
    fractions: List[float],
    top_k_values: Sequence[int],
) -> pd.DataFrame:
    """Vary shortlist fraction and measure how many top-k ligands are retained
    compared to the baseline 5% shortlist.
    """
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]
    n_total_stage0 = df["shortlist_rank"].max()  # approx: shortlist_rank = rank within shortlist

    # Baseline: ligands with shortlist_rank <= actual shortlist size
    # We infer shortlist size from the max shortlist_rank present
    actual_shortlist = df[df["shortlist_rank"].notna()]["ligand_id"].tolist()
    baseline_top_k = {k: set(actual_shortlist[:k]) for k in top_k_values}

    records = []
    for frac in fractions:
        # Re-cut: take top frac% by cascade_score_pct from the full stage3 table
        n_keep = max(1, int(round(len(df) * frac)))
        alt_shortlist = (
            df.sort_values("cascade_score_pct", ascending=False)
            .head(n_keep)["ligand_id"]
            .tolist()
        )

        row: dict = {"shortlist_fraction": frac, "n_shortlisted": n_keep}
        for k in top_k_values:
            alt_top = set(alt_shortlist[:k])
            row[f"jaccard_vs_baseline_top{k}"] = (
                len(alt_top & baseline_top_k[k]) / len(alt_top | baseline_top_k[k])
                if (alt_top | baseline_top_k[k]) else float("nan")
            )
        records.append(row)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_ablations(
    tables: Dict[str, pd.DataFrame],
    ablation_config: dict,
) -> AblationResult:
    """Run all named ablations and sensitivity sweeps.

    Parameters
    ----------
    tables:
        Loaded tables dict from loader.load_tables().
    ablation_config:
        Parsed ablation.yaml content.
    """
    top_k = ablation_config.get("top_k_comparison_values", [10, 50, 100, 500])
    required = ablation_config.get("required_ablations", [])

    # Reference: full_cascade ranking (use final_ranked table)
    ref_df = tables["final_ranked"].copy()
    ref_df["ligand_id"] = ref_df["zinc_id"]
    ref_df["rank"] = ref_df["final_rank"].astype(int)
    ref_df["score"] = -ref_df["final_rank"].astype(float)
    ref_df["method"] = "reference_full_cascade"

    ablation_fns = {
        "remove_stage1":          ablate_remove_stage1,
        "remove_stage2":          ablate_remove_stage2,
        "remove_stage3":          ablate_remove_stage3,
        "alternative_tiebreak":   ablate_alternative_tiebreak,
    }

    ablation_records = []
    rank_lists: Dict[str, pd.DataFrame] = {"reference_full_cascade": ref_df[RANK_COLS]}

    for ablation_name in required:
        if ablation_name == "native_top_stage3_only":
            n = ablation_config.get("native_top_stage3_n", 5000)
            alt = ablate_native_top_stage3_only(tables, n=n)
        elif ablation_name in ablation_fns:
            alt = ablation_fns[ablation_name](tables)
        else:
            logger.warning("Unknown ablation '%s', skipping", ablation_name)
            continue

        row = _compare_to_reference(alt, ref_df[RANK_COLS], top_k, ablation_name)
        ablation_records.append(row)
        rank_lists[ablation_name] = alt
        logger.info("Ablation %-40s  τ=%.3f  ρ=%.3f  J10=%.3f",
                    ablation_name, row["kendall_tau"], row["spearman_r"],
                    row.get("jaccard_top10", float("nan")))

    ablation_df = pd.DataFrame(ablation_records) if ablation_records else pd.DataFrame()

    # Sensitivity sweeps
    weight_pairs = ablation_config.get("cascade_weight_sweep", [])
    weight_pairs = [tuple(w) for w in weight_pairs]
    fractions    = ablation_config.get("shortlist_fraction_sweep", [0.05])

    sensitivity_parts = []

    if weight_pairs and "stage3" in tables:
        logger.info("Running cascade weight sensitivity sweep (%d configs) ...", len(weight_pairs))
        wdf = sweep_cascade_weights(tables, weight_pairs, top_k)
        wdf["sweep_type"] = "cascade_weights"
        sensitivity_parts.append(wdf)

    if fractions and "stage3" in tables:
        logger.info("Running shortlist fraction sensitivity sweep (%d configs) ...", len(fractions))
        sdf = sweep_shortlist_fraction(tables, fractions, top_k)
        sdf["sweep_type"] = "shortlist_fraction"
        sensitivity_parts.append(sdf)

    sensitivity_df = (
        pd.concat(sensitivity_parts, ignore_index=True)
        if sensitivity_parts
        else pd.DataFrame()
    )

    return AblationResult(
        ablation_df=ablation_df,
        sensitivity_df=sensitivity_df,
        rank_lists=rank_lists,
    )
