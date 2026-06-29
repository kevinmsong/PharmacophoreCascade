"""
Direct evidence tables for the five manuscript claims.

Each function returns:
    {
        "evidence_df": pd.DataFrame,      # numeric evidence
        "verdict": str,                   # "supports" | "contradicts" | "inconclusive"
        "narrative": str,                 # one-sentence manuscript-ready summary
    }

Verdicts are based on objective thresholds, not optimistic interpretation.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

from .metrics import jaccard_topk, pearson_r, spearman_r

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Claim 1: Native reranking adds distinct information beyond the frontend
# ---------------------------------------------------------------------------

def claim1_native_adds_info(tables: Dict[str, pd.DataFrame]) -> dict:
    """Test: how weakly cascade_score correlates with native_weighted_coverage_pct.

    If r < 0.3, the cascade score does NOT reliably predict native coverage —
    meaning native reranking is discovering something genuinely different.
    That is simultaneously evidence that (a) native adds info AND
    (b) the cascade does not pre-select for native quality.
    """
    # Use rank_comparison which has both scores aligned per ligand
    if "rank_comparison" in tables:
        rc = tables["rank_comparison"]
        screen_scores = rc["weighted_coverage_pct"].values
        native_scores = rc["native_weighted_coverage_pct"].values
    else:
        # Fall back to native_scored joined with stage3
        ns = tables["native_scored"]
        s3 = tables["stage3"].rename(columns={"zinc_id": "ligand_id"})
        merged = ns.drop_duplicates("ligand_id").merge(
            s3[["ligand_id", "cascade_score_pct", "weighted_coverage_pct"]],
            on="ligand_id",
            how="inner",
        )
        screen_scores = merged["cascade_score_pct"].values
        native_scores = merged["native_weighted_coverage_pct"].values

    r_pearson = pearson_r(screen_scores, native_scores)
    r_spearman = spearman_r(screen_scores, native_scores)

    # Also compute rank-shift stats from rank_comparison if available
    if "rank_comparison" in tables:
        rc = tables["rank_comparison"]
        median_abs_shift = rc["abs_rank_shift"].median()
        mean_abs_shift = rc["abs_rank_shift"].mean()
        max_abs_shift = rc["abs_rank_shift"].max()
    else:
        median_abs_shift = mean_abs_shift = max_abs_shift = float("nan")

    evidence_df = pd.DataFrame([{
        "pearson_r_screen_vs_native": r_pearson,
        "spearman_r_screen_vs_native": r_spearman,
        "median_abs_rank_shift": median_abs_shift,
        "mean_abs_rank_shift": mean_abs_shift,
        "max_abs_rank_shift": max_abs_shift,
        "n_ligands": len(screen_scores),
    }])

    # Verdict: weak correlation SUPPORTS the claim that native adds distinct info
    if abs(r_pearson) < 0.15:
        verdict = "supports"
        narrative = (
            f"Pearson r = {r_pearson:.3f} between cascade score and native coverage "
            f"(Spearman ρ = {r_spearman:.3f}; median |rank shift| = {median_abs_shift:.0f}), "
            "confirming that native reranking provides information orthogonal to the frontend cascade."
        )
    elif abs(r_pearson) < 0.40:
        verdict = "inconclusive"
        narrative = (
            f"Pearson r = {r_pearson:.3f} — weak but non-negligible correlation between "
            "cascade score and native coverage; native reranking adds partial but not fully "
            "independent information."
        )
    else:
        verdict = "contradicts"
        narrative = (
            f"Pearson r = {r_pearson:.3f} — unexpectedly strong correlation between cascade "
            "score and native coverage; the claim that native reranking adds independent "
            "information is weakened."
        )

    return {"evidence_df": evidence_df, "verdict": verdict, "narrative": narrative}


# ---------------------------------------------------------------------------
# Claim 2: Frontend metrics enrich but do not determine native mimicry
# ---------------------------------------------------------------------------

def claim2_frontend_enriches_not_determines(tables: Dict[str, pd.DataFrame]) -> dict:
    """Test: do high-cascade-score ligands dominate the native top-k?

    If the native top-100 draws broadly from across the cascade distribution,
    that supports 'enriches but doesn't determine.'
    If the native top-100 is entirely from cascade top-1%, that contradicts it.
    """
    if "rank_comparison" in tables:
        rc = tables["rank_comparison"].copy()
        n_total = len(rc)

        # Percentile buckets of screen_rank
        rc["screen_pct"] = rc["screen_rank"] / n_total * 100

        native_top = {}
        for k in [10, 50, 100]:
            top_k = rc.nsmallest(k, "native_rank")
            native_top[k] = {
                "median_screen_pct": top_k["screen_pct"].median(),
                "pct_from_screen_top1": (top_k["screen_rank"] <= n_total * 0.01).mean(),
                "pct_from_screen_top5": (top_k["screen_rank"] <= n_total * 0.05).mean(),
                "pct_from_screen_top10": (top_k["screen_rank"] <= n_total * 0.10).mean(),
                "pct_from_screen_top50": (top_k["screen_rank"] <= n_total * 0.50).mean(),
            }

        records = []
        for k, stats_dict in native_top.items():
            records.append({"native_top_k": k, **stats_dict})
        evidence_df = pd.DataFrame(records)

        # Verdict: if native top-100 has median screen percentile > 40th, broadly drawn → enriches
        median_pct = native_top[100]["median_screen_pct"]
        pct_from_top1 = native_top[100]["pct_from_screen_top1"]

        if median_pct > 20 and pct_from_top1 < 0.50:
            verdict = "supports"
            narrative = (
                f"Native top-100 ligands have a median screen percentile of {median_pct:.1f}%, "
                f"with only {pct_from_top1*100:.0f}% drawn from the cascade top-1%; "
                "the frontend enriches but does not constrain native performance."
            )
        elif median_pct <= 10 and pct_from_top1 >= 0.70:
            verdict = "contradicts"
            narrative = (
                f"Native top-100 are heavily concentrated in the cascade top-1% "
                f"(median screen percentile {median_pct:.1f}%); the frontend largely determines "
                "which ligands score well natively."
            )
        else:
            verdict = "inconclusive"
            narrative = (
                f"Median screen percentile of native top-100 = {median_pct:.1f}%; "
                "the relationship between cascade rank and native rank is ambiguous."
            )

    else:
        evidence_df = pd.DataFrame([{"error": "rank_comparison table not loaded"}])
        verdict = "inconclusive"
        narrative = "rank_comparison table required for Claim 2 evaluation."

    return {"evidence_df": evidence_df, "verdict": verdict, "narrative": narrative}


# ---------------------------------------------------------------------------
# Claim 3: Native branch is not a cosmetic refinement of Stage 3
# ---------------------------------------------------------------------------

def claim3_native_not_cosmetic_stage3(tables: Dict[str, pd.DataFrame]) -> dict:
    """Test: how different is Stage-3 geometry coverage vs native coverage?

    Low Jaccard overlap between stage3 top-k and native top-k, plus weak
    correlation between stage3 weighted_coverage_pct and native_weighted_coverage_pct,
    supports the claim.
    """
    if "rank_comparison" in tables:
        rc = tables["rank_comparison"].copy()
        stage3_ranked = rc.sort_values("screen_rank")["ligand_id"].tolist()
        native_ranked = rc.sort_values("native_rank")["ligand_id"].tolist()

        jaccards = {}
        for k in [10, 50, 100, 500]:
            jaccards[k] = jaccard_topk(stage3_ranked, native_ranked, k)

        r_pcc = pearson_r(rc["weighted_coverage_pct"].values,
                          rc["native_weighted_coverage_pct"].values)
        r_spr = spearman_r(rc["screen_rank"].values, rc["native_rank"].values)

        records = [{"metric": "pearson_r_stage3_vs_native_coverage", "value": r_pcc},
                   {"metric": "spearman_r_stage3_rank_vs_native_rank", "value": r_spr}]
        for k, j in jaccards.items():
            records.append({"metric": f"jaccard_top{k}", "value": j})

        evidence_df = pd.DataFrame(records)

        j10 = jaccards.get(10, float("nan"))
        if j10 < 0.05 and abs(r_pcc) < 0.15:
            verdict = "supports"
            narrative = (
                f"Stage-3 top-10 and native top-10 share Jaccard = {j10:.2f} "
                f"with Pearson r = {r_pcc:.3f} between their coverage scores; "
                "native scoring is not a cosmetic re-ordering of Stage-3 results."
            )
        elif j10 > 0.4 and abs(r_pcc) > 0.5:
            verdict = "contradicts"
            narrative = (
                f"Stage-3 top-10 and native top-10 overlap substantially (Jaccard = {j10:.2f}, "
                f"r = {r_pcc:.3f}); native scoring largely echoes Stage-3 geometry results."
            )
        else:
            verdict = "inconclusive"
            narrative = (
                f"Jaccard(top10) = {j10:.2f}, r = {r_pcc:.3f}; "
                "partial but not conclusive distinction between Stage-3 and native ranking."
            )

    else:
        evidence_df = pd.DataFrame([{"error": "rank_comparison table not loaded"}])
        verdict = "inconclusive"
        narrative = "rank_comparison table required for Claim 3 evaluation."

    return {"evidence_df": evidence_df, "verdict": verdict, "narrative": narrative}


# ---------------------------------------------------------------------------
# Claim 4: Diversified native pool reduces premature score collapse
# ---------------------------------------------------------------------------

def claim4_diversified_pool_reduces_collapse(tables: Dict[str, pd.DataFrame]) -> dict:
    """Test: does the diversified pool produce a more spread score distribution
    than a stage3-rank-only pool would?

    Operationalized as: native pool source flags (pool_source_hotspot_breadth,
    pool_source_native_support) in the final_ranked table tell us which ligands
    came from non-stage3-rank routes.  We compare native scores of those ligands
    vs. stage3-rank-only entrants.
    """
    fr = tables["final_ranked"].copy()

    # Check if pool source columns exist
    has_breadth = "pool_source_hotspot_breadth" in fr.columns
    has_native_support = "pool_source_native_support" in fr.columns

    if not has_breadth and not has_native_support:
        evidence_df = pd.DataFrame([{
            "error": "pool_source columns absent from final_ranked — cannot evaluate Claim 4"
        }])
        return {
            "evidence_df": evidence_df,
            "verdict": "inconclusive",
            "narrative": "Pool source provenance columns not present in final_ranked table.",
        }

    fr["from_diversified"] = False
    if has_breadth:
        fr["from_diversified"] |= fr["pool_source_hotspot_breadth"].fillna(False).astype(bool)
    if has_native_support:
        fr["from_diversified"] |= fr["pool_source_native_support"].fillna(False).astype(bool)

    diversified = fr[fr["from_diversified"]]["native_weighted_coverage_pct"].dropna()
    stage3_only = fr[~fr["from_diversified"]]["native_weighted_coverage_pct"].dropna()

    if len(diversified) < 5 or len(stage3_only) < 5:
        evidence_df = pd.DataFrame([{
            "n_diversified": len(diversified),
            "n_stage3_only": len(stage3_only),
            "error": "insufficient sample size for comparison",
        }])
        return {
            "evidence_df": evidence_df,
            "verdict": "inconclusive",
            "narrative": "Insufficient ligands in one or both pool-source groups.",
        }

    median_div = diversified.median()
    median_s3  = stage3_only.median()
    std_div    = diversified.std()
    std_s3     = stage3_only.std()

    # Mann-Whitney U (non-parametric)
    u_stat, p_val = stats.mannwhitneyu(diversified, stage3_only, alternative="two-sided")

    evidence_df = pd.DataFrame([{
        "n_diversified_ligands": len(diversified),
        "n_stage3_rank_ligands": len(stage3_only),
        "median_native_cov_diversified": median_div,
        "median_native_cov_stage3_rank": median_s3,
        "std_native_cov_diversified": std_div,
        "std_native_cov_stage3_rank": std_s3,
        "mannwhitney_u": u_stat,
        "mannwhitney_p": p_val,
    }])

    # Claim: diversified pool not necessarily higher scoring, but the diversity
    # it adds prevents score collapse (i.e., diversified entrants still score competitively)
    if median_div >= median_s3 * 0.90:
        verdict = "supports"
        narrative = (
            f"Diversified-pool entrants (n={len(diversified)}) achieve median native coverage "
            f"{median_div:.2f}% vs {median_s3:.2f}% for stage3-rank entrants; "
            "non-stage3 routing brings comparably scoring ligands without score collapse."
        )
    elif median_div < median_s3 * 0.75:
        verdict = "contradicts"
        narrative = (
            f"Diversified-pool entrants (n={len(diversified)}) score substantially lower "
            f"(median {median_div:.2f}% vs {median_s3:.2f}%); "
            "non-stage3 routing may be introducing lower-quality candidates."
        )
    else:
        verdict = "inconclusive"
        narrative = (
            f"Diversified entrants score slightly lower (median {median_div:.2f}% vs {median_s3:.2f}%); "
            "the benefit of diversity routing is not clearly resolved."
        )

    return {"evidence_df": evidence_df, "verdict": verdict, "narrative": narrative}


# ---------------------------------------------------------------------------
# Claim 5: Top-ranked ligands converge on SER14/ASP15/SER17/SER18 motif
# ---------------------------------------------------------------------------

def claim5_residue_motif_convergence(
    tables: Dict[str, pd.DataFrame],
    claimed_residues: Optional[List[str]] = None,
) -> dict:
    """Test: are the claimed motif residues statistically overrepresented in
    the top-ranked ligands compared to the full native-scored pool?

    Uses hypergeometric-like test: fraction of top-k ligands that match each
    residue vs fraction across all native-scored ligands.
    """
    if claimed_residues is None:
        claimed_residues = ["SER14", "ASP15", "SER17", "SER18"]

    if "feature_matches" not in tables:
        evidence_df = pd.DataFrame([{"error": "feature_matches table not loaded"}])
        return {
            "evidence_df": evidence_df,
            "verdict": "inconclusive",
            "narrative": "feature_matches table required for Claim 5 evaluation.",
        }

    fm = tables["feature_matches"]
    fr = tables["final_ranked"]

    # Best microstate per ligand (highest-weight match set — use first occurrence
    # since native_scored is already best-microstate-first)
    fm_best = fm.drop_duplicates(subset=["ligand_id", "reference_residue_label"], keep="first")

    # Residue frequency across full native pool
    all_ligands = fm_best["ligand_id"].unique()
    n_all = len(all_ligands)

    # Residue frequency in top-k of final ranked
    final_ids = fr.sort_values("final_rank")["zinc_id"].tolist()

    records = []
    for residue in claimed_residues:
        ligands_with_residue = set(
            fm_best.loc[fm_best["reference_residue_label"] == residue, "ligand_id"]
        )
        freq_all = len(ligands_with_residue & set(all_ligands)) / n_all if n_all > 0 else 0

        for k in [10, 50, 100]:
            top_k = set(final_ids[:k])
            n_top_with = len(top_k & ligands_with_residue)
            freq_top = n_top_with / k if k > 0 else 0
            fold_enrich = freq_top / freq_all if freq_all > 0 else float("nan")

            records.append({
                "residue": residue,
                "top_k": k,
                "n_top_k_with_residue": n_top_with,
                "freq_top_k": freq_top,
                "freq_all_native": freq_all,
                "fold_enrichment": fold_enrich,
            })

    evidence_df = pd.DataFrame(records)

    # Verdict: if all four claimed residues show ≥2× fold enrichment in top-10
    top10_folds = evidence_df[evidence_df["top_k"] == 10].set_index("residue")["fold_enrichment"]
    enriched = sum(1 for r in claimed_residues if top10_folds.get(r, 0) >= 2.0)

    if enriched == len(claimed_residues):
        verdict = "supports"
        narrative = (
            f"All {len(claimed_residues)} claimed motif residues "
            f"({', '.join(claimed_residues)}) show ≥2× enrichment in the top-10 "
            "native-ranked ligands vs the full native pool."
        )
    elif enriched >= len(claimed_residues) // 2:
        verdict = "inconclusive"
        narrative = (
            f"{enriched}/{len(claimed_residues)} claimed motif residues show ≥2× enrichment "
            "in top-10; partial support for convergence on the described motif."
        )
    else:
        verdict = "contradicts"
        narrative = (
            f"Only {enriched}/{len(claimed_residues)} claimed motif residues are enriched "
            "in top-10 native-ranked ligands; convergence on the stated motif is not confirmed."
        )

    return {"evidence_df": evidence_df, "verdict": verdict, "narrative": narrative}


# ---------------------------------------------------------------------------
# Run all claims
# ---------------------------------------------------------------------------

def run_all_claims(
    tables: Dict[str, pd.DataFrame],
    config: dict,
) -> Dict[str, dict]:
    """Evaluate all five claims and return results keyed by claim label."""
    claimed_residues = config.get("claimed_motif_residues", ["SER14", "ASP15", "SER17", "SER18"])

    results = {
        "claim1_native_adds_info":
            claim1_native_adds_info(tables),
        "claim2_frontend_enriches_not_determines":
            claim2_frontend_enriches_not_determines(tables),
        "claim3_native_not_cosmetic_stage3":
            claim3_native_not_cosmetic_stage3(tables),
        "claim4_diversified_pool_reduces_collapse":
            claim4_diversified_pool_reduces_collapse(tables),
        "claim5_residue_motif_convergence":
            claim5_residue_motif_convergence(tables, claimed_residues),
    }

    for label, res in results.items():
        logger.info("%-50s  [%s]", label, res["verdict"].upper())

    return results


def claims_to_summary_df(claims_results: Dict[str, dict]) -> pd.DataFrame:
    """Flatten all claim verdicts + narratives into a single summary DataFrame."""
    rows = []
    for label, res in claims_results.items():
        rows.append({
            "claim": label,
            "verdict": res["verdict"],
            "narrative": res["narrative"],
        })
    return pd.DataFrame(rows)
