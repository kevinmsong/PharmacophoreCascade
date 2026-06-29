"""
Baseline ranking methods for benchmarking the full cascade.

Each function returns a pd.DataFrame with columns:
    ligand_id : str
    rank      : int  (1 = best)
    score     : float
    method    : str  (human-readable label)

All methods operate on the pre-loaded tables dict from loader.load_tables().
The universe for comparison is the Stage-3 shortlist (47,689 ligands) unless
the method is specific to the native-scored pool (4,997 ligands).
"""
from __future__ import annotations

import logging
from typing import Dict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Standard output columns
RANK_COLS = ["ligand_id", "rank", "score", "method"]


def _add_rank(df: pd.DataFrame, score_col: str, ascending: bool, method: str) -> pd.DataFrame:
    """Sort by score_col and return a RANK_COLS DataFrame."""
    out = df[["ligand_id", score_col]].copy()
    out = out.sort_values(score_col, ascending=ascending, kind="stable").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out = out.rename(columns={score_col: "score"})
    out["method"] = method
    return out[RANK_COLS]


# ---------------------------------------------------------------------------
# 1. Full cascade (as produced — lexicographic native-first ordering)
# ---------------------------------------------------------------------------

def rank_full_cascade(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Use the final_ranked table exactly as the pipeline produced it.

    This is the reference method.  The universe is the final 1,000 ligands.
    """
    df = tables["final_ranked"].copy()
    df["ligand_id"] = df["zinc_id"]
    df["score"] = -df["final_rank"].astype(float)   # lower rank = better score
    df["rank"] = df["final_rank"].astype(int)
    df["method"] = "full_cascade"
    return df[RANK_COLS]


# ---------------------------------------------------------------------------
# 2. Cascade score only (hotspot + pair-hash blend, no Stage-3 / native)
# ---------------------------------------------------------------------------

def rank_by_cascade_score(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Rank the Stage-3 shortlist by cascade_score_pct alone.

    cascade_score_pct = 0.4 * hotspot_weighted_pct + 0.6 * pair_hash_overlap_pct
    This corresponds to stopping after Stage 2.
    """
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]
    return _add_rank(df, "cascade_score_pct", ascending=False, method="cascade_score_only")


# ---------------------------------------------------------------------------
# 3. Stage-3 geometry only (no native scoring)
# ---------------------------------------------------------------------------

def rank_by_stage3_only(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Rank Stage-3 shortlist by Stage-3 geometry coverage (weighted_coverage_pct).

    weighted_coverage_pct in the stage3 audit is the anchor-weighted geometry score.
    Ties broken by fit_rmsd_angstrom (ascending).
    """
    df = tables["stage3"].copy()
    df["ligand_id"] = df["zinc_id"]

    df = df.sort_values(
        ["weighted_coverage_pct", "fit_rmsd_angstrom"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)

    df["rank"] = np.arange(1, len(df) + 1)
    df["score"] = df["weighted_coverage_pct"]
    df["method"] = "stage3_only"
    return df[RANK_COLS]


# ---------------------------------------------------------------------------
# 4. Native scoring only (from the native-scored pool)
# ---------------------------------------------------------------------------

def rank_by_native_only(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Rank the native-scored pool by native_weighted_coverage_pct alone.

    One row per unique ligand (best microstate).
    """
    df = tables["native_scored"].copy()

    # Deduplicate: keep best microstate per ligand
    df = (
        df.sort_values(
            ["native_weighted_coverage_pct", "fit_rmsd_angstrom"],
            ascending=[False, True],
            kind="stable",
        )
        .drop_duplicates(subset=["ligand_id"], keep="first")
        .reset_index(drop=True)
    )

    df["rank"] = np.arange(1, len(df) + 1)
    df["score"] = df["native_weighted_coverage_pct"]
    df["method"] = "native_only"
    return df[RANK_COLS]


# ---------------------------------------------------------------------------
# 5. Random baseline (shuffle Stage-3 shortlist)
# ---------------------------------------------------------------------------

def rank_random(tables: Dict[str, pd.DataFrame], seed: int = 42) -> pd.DataFrame:
    """Uniformly random ranking of the Stage-3 shortlist."""
    df = tables["stage3"][["zinc_id"]].copy()
    df["ligand_id"] = df["zinc_id"]
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(len(df))
    df["rank"] = shuffled + 1
    df["score"] = (len(df) - shuffled).astype(float)
    df = df.sort_values("rank").reset_index(drop=True)
    df["method"] = f"random_seed{seed}"
    return df[RANK_COLS]


# ---------------------------------------------------------------------------
# 6. Property-matched random baseline
# ---------------------------------------------------------------------------

def rank_property_matched_random(
    tables: Dict[str, pd.DataFrame],
    seed: int = 42,
    mw_bins: list | None = None,
    logp_bins: list | None = None,
) -> pd.DataFrame:
    """Random ranking stratified by MW × LogP bins.

    Within each property bin, ligands are shuffled independently, then bins
    are concatenated so rank is determined by within-bin shuffle.  This tests
    whether any observed enrichment could arise from property bias alone.
    """
    if mw_bins is None:
        mw_bins = [200, 300, 400, 500, 700]
    if logp_bins is None:
        logp_bins = [-2, 1, 3, 5, 8]

    df = tables["stage3"][["zinc_id", "mw", "logp"]].copy()
    df["ligand_id"] = df["zinc_id"]

    df["mw_bin"] = pd.cut(df["mw"], bins=mw_bins, labels=False).fillna(-1).astype(int)
    df["logp_bin"] = pd.cut(df["logp"], bins=logp_bins, labels=False).fillna(-1).astype(int)
    df["bin_key"] = df["mw_bin"].astype(str) + "_" + df["logp_bin"].astype(str)

    rng = np.random.default_rng(seed)
    shuffled_parts = []
    for _, grp in df.groupby("bin_key", sort=False):
        g = grp.copy()
        g = g.iloc[rng.permutation(len(g))].reset_index(drop=True)
        shuffled_parts.append(g)

    out = pd.concat(shuffled_parts, ignore_index=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["score"] = (len(out) - out["rank"]).astype(float)
    out["method"] = f"property_matched_random_seed{seed}"
    return out[RANK_COLS]


# ---------------------------------------------------------------------------
# Registry: all methods that operate on Stage-3 universe
# ---------------------------------------------------------------------------

def build_all_rankings(
    tables: Dict[str, pd.DataFrame],
    random_seed: int = 42,
    mw_bins: list | None = None,
    logp_bins: list | None = None,
) -> Dict[str, pd.DataFrame]:
    """Build all baseline rankings and return as a dict keyed by method name.

    Always produces:
        full_cascade, cascade_score_only, stage3_only,
        native_only, random, property_matched_random

    The universe differs by method:
        - full_cascade / native_only : native-scored pool (~5k)
        - others                     : Stage-3 shortlist (~47k)
    """
    rankings: Dict[str, pd.DataFrame] = {}

    rankings["full_cascade"] = rank_full_cascade(tables)
    logger.info("Built full_cascade ranking: %d ligands", len(rankings["full_cascade"]))

    rankings["cascade_score_only"] = rank_by_cascade_score(tables)
    logger.info("Built cascade_score_only ranking: %d ligands", len(rankings["cascade_score_only"]))

    rankings["stage3_only"] = rank_by_stage3_only(tables)
    logger.info("Built stage3_only ranking: %d ligands", len(rankings["stage3_only"]))

    rankings["native_only"] = rank_by_native_only(tables)
    logger.info("Built native_only ranking: %d ligands", len(rankings["native_only"]))

    rankings["random"] = rank_random(tables, seed=random_seed)
    logger.info("Built random ranking: %d ligands", len(rankings["random"]))

    rankings["property_matched_random"] = rank_property_matched_random(
        tables, seed=random_seed, mw_bins=mw_bins, logp_bins=logp_bins
    )
    logger.info(
        "Built property_matched_random ranking: %d ligands",
        len(rankings["property_matched_random"]),
    )

    return rankings
