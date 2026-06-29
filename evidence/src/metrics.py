"""
Ranking and enrichment metrics for benchmarking.

All enrichment functions accept:
    ranked_ids  : sequence of ligand IDs in rank order (rank 1 first)
    active_ids  : set/collection of known-active ligand IDs

When active_ids is empty or None, enrichment metrics return NaN with a warning
instead of crashing.  Rank-correlation metrics never need active_ids.
"""
from __future__ import annotations

import logging
import math
import warnings
from typing import Collection, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _check_actives(active_ids: Collection[str], metric_name: str) -> Optional[set]:
    if not active_ids:
        warnings.warn(
            f"{metric_name}: no active IDs supplied — returning NaN",
            UserWarning,
            stacklevel=3,
        )
        return None
    return set(active_ids)


def _binary_labels(ranked_ids: Sequence[str], active_set: set) -> np.ndarray:
    """1 = active, 0 = inactive, in rank order."""
    return np.array([1 if lid in active_set else 0 for lid in ranked_ids], dtype=float)


# ---------------------------------------------------------------------------
# Enrichment metrics
# ---------------------------------------------------------------------------

def enrichment_factor(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
    fraction: float,
) -> float:
    """Enrichment factor at a given top fraction of the ranked list.

    EF = (actives_in_top_fraction / total_actives) / fraction
    """
    active_set = _check_actives(active_ids, "enrichment_factor")
    if active_set is None:
        return float("nan")

    n_total = len(ranked_ids)
    n_top = max(1, int(math.ceil(n_total * fraction)))
    top_ids = set(ranked_ids[:n_top])

    n_actives_total = len(active_set)
    if n_actives_total == 0:
        return float("nan")

    n_actives_top = len(top_ids & active_set)
    return (n_actives_top / n_actives_total) / fraction


def scaffold_recovery(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
    scaffold_series: pd.Series,
    k: int,
) -> float:
    """Fraction of active scaffolds represented in the top-k ranked molecules."""
    active_set = _check_actives(active_ids, "scaffold_recovery")
    if active_set is None:
        return float("nan")

    active_scaffolds = set(
        scaffold_series.reindex(list(active_set)).dropna().astype(str).tolist()
    )
    if not active_scaffolds:
        return float("nan")

    top_scaffolds = set(
        scaffold_series.reindex(list(ranked_ids[:k])).dropna().astype(str).tolist()
    )
    return len(active_scaffolds & top_scaffolds) / len(active_scaffolds)


def top_k_recovery(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
    k: int,
) -> float:
    """Fraction of known actives recovered in the top-k ranked molecules."""
    active_set = _check_actives(active_ids, "top_k_recovery")
    if active_set is None:
        return float("nan")
    if not active_set:
        return float("nan")
    top_k = set(ranked_ids[:k])
    return len(top_k & active_set) / len(active_set)


def roc_auc(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
) -> float:
    """Area under the ROC curve (trapezoidal rule)."""
    from sklearn.metrics import roc_auc_score  # type: ignore

    active_set = _check_actives(active_ids, "roc_auc")
    if active_set is None:
        return float("nan")

    labels = _binary_labels(ranked_ids, active_set)
    if labels.sum() == 0 or labels.sum() == len(labels):
        return float("nan")

    # Scores are simply the inverse of rank index (rank 1 → highest score)
    scores = np.arange(len(ranked_ids), 0, -1, dtype=float)
    return float(roc_auc_score(labels, scores))


def pr_auc(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
) -> float:
    """Area under the precision-recall curve."""
    from sklearn.metrics import average_precision_score  # type: ignore

    active_set = _check_actives(active_ids, "pr_auc")
    if active_set is None:
        return float("nan")

    labels = _binary_labels(ranked_ids, active_set)
    if labels.sum() == 0:
        return float("nan")

    scores = np.arange(len(ranked_ids), 0, -1, dtype=float)
    return float(average_precision_score(labels, scores))


def bedroc(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
    alpha: float = 20.0,
) -> float:
    """Boltzmann-Enhanced Discrimination of Receiver Operating Characteristic.

    Implements the normalized formula from Truchon & Bayly (2007).
    alpha=20 weights early enrichment heavily (default for VS).
    BEDROC is a normalized transform of RIE and should lie on [0, 1].
    """
    active_set = _check_actives(active_ids, "bedroc")
    if active_set is None:
        return float("nan")

    labels = _binary_labels(ranked_ids, active_set)
    n = len(labels)
    n_a = int(labels.sum())
    if n_a == 0 or n_a == n:
        return float("nan")

    ra = n_a / n
    ri_sum = sum(
        math.exp(-alpha * (i + 1) / n)
        for i, lab in enumerate(labels)
        if lab == 1
    )
    rie_denom = (1 - math.exp(-alpha)) / (n * (math.exp(alpha / n) - 1))
    rie = ri_sum / (n_a * rie_denom)
    rie_max = (1 - math.exp(-alpha * ra)) / (ra * (1 - math.exp(-alpha)))
    rie_min = (1 - math.exp(alpha * ra)) / (ra * (1 - math.exp(alpha)))
    if math.isclose(rie_max, rie_min):
        return 1.0

    bedroc_score = (rie - rie_min) / (rie_max - rie_min)
    return float(min(max(bedroc_score, 0.0), 1.0))


# ---------------------------------------------------------------------------
# Diversity metrics (no actives needed)
# ---------------------------------------------------------------------------

def scaffold_diversity(
    ligand_ids: Sequence[str],
    scaffold_series: pd.Series,
) -> int:
    """Count unique Murcko scaffolds in a set of ligand IDs.

    Parameters
    ----------
    ligand_ids:
        Ordered IDs to evaluate.
    scaffold_series:
        pd.Series indexed by ligand_id with murcko_scaffold values.
    """
    present = scaffold_series.reindex(ligand_ids).dropna()
    return int(present.nunique())


# ---------------------------------------------------------------------------
# Rank correlation metrics (no actives needed)
# ---------------------------------------------------------------------------

def kendall_tau(
    rank_a: Sequence[float],
    rank_b: Sequence[float],
) -> float:
    """Kendall's tau-b between two rank sequences (same ligands, same order)."""
    if len(rank_a) != len(rank_b) or len(rank_a) < 2:
        return float("nan")
    tau, _ = stats.kendalltau(rank_a, rank_b)
    return float(tau)


def spearman_r(
    rank_a: Sequence[float],
    rank_b: Sequence[float],
) -> float:
    """Spearman rank correlation between two rank sequences."""
    if len(rank_a) != len(rank_b) or len(rank_a) < 2:
        return float("nan")
    rho, _ = stats.spearmanr(rank_a, rank_b)
    return float(rho)


def pearson_r(
    values_a: Sequence[float],
    values_b: Sequence[float],
) -> float:
    """Pearson correlation between two numeric sequences."""
    if len(values_a) != len(values_b) or len(values_a) < 2:
        return float("nan")
    r, _ = stats.pearsonr(values_a, values_b)
    return float(r)


def jaccard_topk(
    ranked_ids_a: Sequence[str],
    ranked_ids_b: Sequence[str],
    k: int,
) -> float:
    """Jaccard similarity of the top-k sets from two ranked lists."""
    top_a = set(ranked_ids_a[:k])
    top_b = set(ranked_ids_b[:k])
    union = top_a | top_b
    if not union:
        return float("nan")
    return len(top_a & top_b) / len(union)


# ---------------------------------------------------------------------------
# Convenience: compute the full enrichment metric bundle
# ---------------------------------------------------------------------------

def compute_enrichment_bundle(
    ranked_ids: Sequence[str],
    active_ids: Collection[str],
    top_k_values: Sequence[int] = (10, 25, 50, 100, 500),
    scaffold_series: pd.Series | None = None,
) -> dict:
    """Return all enrichment metrics as a flat dict.

    Gracefully returns NaN for any metric that cannot be computed.
    """
    result: dict = {}
    result["ef_1pct"]  = enrichment_factor(ranked_ids, active_ids, 0.01)
    result["ef_5pct"]  = enrichment_factor(ranked_ids, active_ids, 0.05)
    result["ef_01pct"] = enrichment_factor(ranked_ids, active_ids, 0.001)
    result["roc_auc"]  = roc_auc(ranked_ids, active_ids)
    result["pr_auc"]   = pr_auc(ranked_ids, active_ids)
    result["bedroc"]   = bedroc(ranked_ids, active_ids)
    for k in top_k_values:
        result[f"top{k}_recovery"] = top_k_recovery(ranked_ids, active_ids, k)
        if scaffold_series is not None:
            result[f"scaffold_recovery_top{k}"] = scaffold_recovery(
                ranked_ids, active_ids, scaffold_series, k
            )
    return result


def compute_correlation_bundle(
    values_a: Sequence[float],
    values_b: Sequence[float],
    label_a: str = "a",
    label_b: str = "b",
) -> dict:
    """Pearson + Spearman correlation between two score/rank sequences."""
    return {
        f"pearson_r_{label_a}_vs_{label_b}": pearson_r(values_a, values_b),
        f"spearman_r_{label_a}_vs_{label_b}": spearman_r(values_a, values_b),
    }
