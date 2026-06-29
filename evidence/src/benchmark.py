"""
Benchmark runner: compare all ranking methods on rank-correlation and
optional enrichment metrics.

Entry point: run_benchmark(tables, config) -> BenchmarkResult
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Collection, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .baselines import build_all_rankings
from .benchmark_runner import ExternalBenchmarkArtifacts, run_external_benchmark
from .metrics import (
    compute_enrichment_bundle,
    jaccard_topk,
    kendall_tau,
    scaffold_diversity,
    spearman_r,
)

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResult:
    """All benchmark outputs in one place."""
    rankings: Dict[str, pd.DataFrame]               # method -> ranked DataFrame
    pairwise_df: pd.DataFrame                        # jaccard + rank-corr between method pairs
    summary_df: pd.DataFrame                         # per-method enrichment + diversity metrics
    scaffold_series: Optional[pd.Series] = None      # ligand_id -> murcko_scaffold
    delta_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    library_df: Optional[pd.DataFrame] = None


# ---------------------------------------------------------------------------
# Pairwise rank correlation + Jaccard
# ---------------------------------------------------------------------------

def _pairwise_metrics(
    rankings: Dict[str, pd.DataFrame],
    top_k_values: Sequence[int],
    common_universe: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Compute Kendall tau, Spearman r, and Jaccard top-k between every method pair."""
    methods = list(rankings.keys())
    records = []

    for i, m_a in enumerate(methods):
        for m_b in methods[i + 1:]:
            df_a = rankings[m_a].set_index("ligand_id")["rank"]
            df_b = rankings[m_b].set_index("ligand_id")["rank"]

            # Align on shared ligands
            shared = df_a.index.intersection(df_b.index)
            if len(shared) < 5:
                logger.warning(
                    "Fewer than 5 shared ligands between %s and %s (%d)",
                    m_a, m_b, len(shared),
                )

            ra = df_a.reindex(shared).values
            rb = df_b.reindex(shared).values

            row: dict = {
                "method_a": m_a,
                "method_b": m_b,
                "n_shared": len(shared),
                "kendall_tau": kendall_tau(ra, rb),
                "spearman_r": spearman_r(ra, rb),
            }

            ids_a = rankings[m_a].sort_values("rank")["ligand_id"].tolist()
            ids_b = rankings[m_b].sort_values("rank")["ligand_id"].tolist()

            for k in top_k_values:
                row[f"jaccard_top{k}"] = jaccard_topk(ids_a, ids_b, k)

            records.append(row)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Per-method summary (enrichment + diversity)
# ---------------------------------------------------------------------------

def _method_summary(
    rankings: Dict[str, pd.DataFrame],
    active_ids: Collection[str],
    top_k_values: Sequence[int],
    scaffold_series: Optional[pd.Series],
) -> pd.DataFrame:
    records = []
    for method, df in rankings.items():
        ranked_ids = df.sort_values("rank")["ligand_id"].tolist()

        row: dict = {"method": method, "n_ranked": len(ranked_ids)}

        # Enrichment metrics (NaN if no actives)
        bundle = compute_enrichment_bundle(
            ranked_ids,
            active_ids,
            top_k_values,
            scaffold_series=scaffold_series,
        )
        row.update(bundle)

        # Scaffold diversity at top-k
        if scaffold_series is not None:
            for k in top_k_values:
                top_ids = ranked_ids[:k]
                row[f"scaffold_diversity_top{k}"] = scaffold_diversity(top_ids, scaffold_series)

        records.append(row)

    return pd.DataFrame(records)


def _holm_adjust(pvals: Sequence[float]) -> list[float]:
    """Holm-Bonferroni step-down adjustment of a list of p-values."""
    arr = np.asarray(pvals, dtype=float)
    n = len(arr)
    if n == 0:
        return []
    order = np.argsort(arr)
    adjusted = np.empty(n, dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (n - rank) * arr[idx]
        running = max(running, val)
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


def _stratified_bootstrap_summary(
    rankings: Dict[str, pd.DataFrame],
    library_df: pd.DataFrame,
    top_k_values: Sequence[int],
    scaffold_series: Optional[pd.Series],
    iterations: int,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (summary_with_ci, paired_deltas_vs_full_cascade).

    Resampling is performed at the active-cluster level: each known active is
    resampled together with its matched decoys (keyed by ``matched_active_ligand_id``),
    which respects the active-decoy matching structure and yields confidence
    intervals dominated by the small number of actives. Paired two-sided
    bootstrap p-values (with Holm adjustment across methods within each metric)
    quantify differences between each method and the full cascade.
    """
    metric_frames: dict[str, list[dict[str, float]]] = {method: [] for method in rankings}
    delta_records: list[dict[str, float | str]] = []

    active_df = library_df.loc[library_df["label"].astype(str) == "active"].copy()
    decoy_df = library_df.loc[library_df["label"].astype(str) != "active"].copy()
    if active_df.empty or decoy_df.empty:
        summary = _method_summary(
            rankings,
            active_df["ligand_id"].astype(str).tolist(),
            top_k_values,
            scaffold_series,
        )
        return summary, pd.DataFrame()

    # Build active-cluster structure: each active plus its matched decoys.
    active_ids_all = [str(a) for a in active_df["ligand_id"].tolist()]
    clusters: dict[str, list[dict]] = {aid: [] for aid in active_ids_all}
    for _, row in active_df.iterrows():
        clusters[str(row["ligand_id"])].append(
            {"ligand_id": str(row["ligand_id"]), "label": "active",
             "murcko_scaffold": row.get("murcko_scaffold")}
        )
    grouped = "matched_active_ligand_id" in decoy_df.columns and decoy_df["matched_active_ligand_id"].notna().any()
    if grouped:
        for offset, (_, row) in enumerate(decoy_df.iterrows()):
            key = str(row.get("matched_active_ligand_id") or "")
            if key not in clusters:
                key = active_ids_all[offset % len(active_ids_all)]
            clusters[key].append(
                {"ligand_id": str(row["ligand_id"]), "label": "decoy",
                 "murcko_scaffold": row.get("murcko_scaffold")}
            )
    cluster_keys = list(clusters.keys())

    rng = np.random.default_rng(random_seed)
    rank_maps = {
        method: df.set_index("ligand_id")["rank"].to_dict()
        for method, df in rankings.items()
    }

    for iteration in range(iterations):
        if grouped:
            chosen = rng.choice(cluster_keys, size=len(cluster_keys), replace=True)
            sample_rows = []
            for ci, ckey in enumerate(chosen):
                for mi, member in enumerate(clusters[ckey]):
                    sample_rows.append(
                        {
                            "sample_id": f"{member['ligand_id']}__b{iteration}_c{ci}_{mi}",
                            "ligand_id": member["ligand_id"],
                            "label": member["label"],
                            "murcko_scaffold": member["murcko_scaffold"],
                        }
                    )
            sample_df = pd.DataFrame(sample_rows)
        else:
            active_sample = active_df.sample(
                n=len(active_df), replace=True,
                random_state=int(rng.integers(0, 2**31 - 1)),
            ).reset_index(drop=True)
            decoy_sample = decoy_df.sample(
                n=len(decoy_df), replace=True,
                random_state=int(rng.integers(0, 2**31 - 1)),
            ).reset_index(drop=True)
            sample_rows = []
            for frame in (active_sample, decoy_sample):
                for idx, row in frame.iterrows():
                    sample_rows.append(
                        {
                            "sample_id": f"{row['ligand_id']}__b{iteration}_{idx}",
                            "ligand_id": row["ligand_id"],
                            "label": row["label"],
                            "murcko_scaffold": row.get("murcko_scaffold"),
                        }
                    )
            sample_df = pd.DataFrame(sample_rows)
        active_ids = sample_df.loc[sample_df["label"] == "active", "sample_id"].tolist()
        sample_scaffold_series = sample_df.set_index("sample_id")["murcko_scaffold"]

        bootstrap_metrics: dict[str, dict] = {}
        for method, rank_map in rank_maps.items():
            ordered_ids = (
                sample_df.assign(
                    _rank=sample_df["ligand_id"].map(rank_map).fillna(np.inf),
                    _sample=sample_df["sample_id"],
                )
                .sort_values(["_rank", "_sample"], kind="stable")
                ["sample_id"]
                .tolist()
            )
            bootstrap_metrics[method] = compute_enrichment_bundle(
                ordered_ids,
                active_ids,
                top_k_values,
                scaffold_series=sample_scaffold_series,
            )
            bootstrap_metrics[method]["n_ranked"] = len(ordered_ids)
            metric_frames[method].append(bootstrap_metrics[method])

        if "full_cascade" in bootstrap_metrics:
            reference = bootstrap_metrics["full_cascade"]
            for method, bundle in bootstrap_metrics.items():
                if method == "full_cascade":
                    continue
                for metric_name, value in bundle.items():
                    if metric_name == "n_ranked":
                        continue
                    ref_value = reference.get(metric_name)
                    if pd.isna(value) or pd.isna(ref_value):
                        delta = np.nan
                    else:
                        delta = float(value) - float(ref_value)
                    delta_records.append(
                        {
                            "method": method,
                            "metric": metric_name,
                            "delta": delta,
                        }
                    )

    point_summary = _method_summary(
        rankings,
        active_df["ligand_id"].astype(str).tolist(),
        top_k_values,
        scaffold_series,
    ).set_index("method")

    enriched_rows = []
    for method, point_row in point_summary.iterrows():
        row = {"method": method}
        bootstrap_df = pd.DataFrame(metric_frames[method])
        for column, value in point_row.items():
            row[column] = value
            if column == "n_ranked":
                continue
            draws = pd.to_numeric(bootstrap_df.get(column, pd.Series(dtype=float)), errors="coerce").dropna()
            if draws.empty:
                row[f"{column}_ci_low"] = np.nan
                row[f"{column}_ci_high"] = np.nan
            else:
                row[f"{column}_ci_low"] = float(np.quantile(draws, 0.025))
                row[f"{column}_ci_high"] = float(np.quantile(draws, 0.975))
        enriched_rows.append(row)

    delta_df = pd.DataFrame(delta_records)
    if not delta_df.empty:
        agg_rows = []
        for (method, metric), grp in delta_df.groupby(["method", "metric"], dropna=False):
            draws = pd.to_numeric(grp["delta"], errors="coerce").dropna()
            if draws.empty:
                continue
            n = len(draws)
            count_ge = int((draws >= 0).sum())
            count_le = int((draws <= 0).sum())
            # two-sided Monte-Carlo p-value with the standard +1 correction
            p_value = min(1.0, (2 * min(count_ge, count_le) + 1) / (n + 1))
            agg_rows.append(
                {
                    "method": method,
                    "metric": metric,
                    "delta_mean": float(draws.mean()),
                    "delta_ci_low": float(np.quantile(draws, 0.025)),
                    "delta_ci_high": float(np.quantile(draws, 0.975)),
                    "n_bootstrap": n,
                    "p_value": float(p_value),
                }
            )
        delta_df = pd.DataFrame(agg_rows)
        if not delta_df.empty:
            delta_df["p_value_holm"] = np.nan
            for metric, grp in delta_df.groupby("metric"):
                adjusted = _holm_adjust(grp["p_value"].tolist())
                for idx, adj in zip(grp.index, adjusted):
                    delta_df.loc[idx, "p_value_holm"] = adj
            delta_df = delta_df.sort_values(["metric", "method"]).reset_index(drop=True)

    return pd.DataFrame(enriched_rows), delta_df


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_benchmark(
    tables: Dict[str, pd.DataFrame],
    config: dict,
) -> BenchmarkResult:
    """Run the full benchmarking suite.

    Parameters
    ----------
    tables:
        Dict from loader.load_tables().
    config:
        Parsed benchmark.yaml content.

    Returns
    -------
    BenchmarkResult with rankings, pairwise_df, summary_df.
    """
    random_seed = config.get("random_seed", 42)
    top_k_values: List[int] = config.get("top_k_values", [10, 25, 50, 100, 500])
    project_root = Path(__file__).resolve().parents[2]

    bench_cfg = config.get("benchmark", {})
    if bench_cfg.get("library_csv"):
        artifacts: ExternalBenchmarkArtifacts = run_external_benchmark(config, project_root)
        library_df = artifacts.library_df.copy()
        scaffold_series = library_df.set_index("ligand_id")["murcko_scaffold"]
        bootstrap_iterations = int(bench_cfg.get("bootstrap_iterations", 500))
        bootstrap_seed = int(bench_cfg.get("bootstrap_seed", random_seed))
        pairwise_df = _pairwise_metrics(artifacts.rankings, top_k_values)
        summary_df, delta_df = _stratified_bootstrap_summary(
            artifacts.rankings,
            library_df,
            top_k_values,
            scaffold_series,
            iterations=bootstrap_iterations,
            random_seed=bootstrap_seed,
        )
        return BenchmarkResult(
            rankings=artifacts.rankings,
            pairwise_df=pairwise_df,
            summary_df=summary_df,
            scaffold_series=scaffold_series,
            delta_df=delta_df,
            library_df=library_df,
        )

    prop_cfg = config.get("property_bins", {})
    mw_bins  = prop_cfg.get("mw_bins", [200, 300, 400, 500, 700])
    logp_bins = prop_cfg.get("logp_bins", [-2, 1, 3, 5, 8])

    # Load optional actives
    actives_path = config.get("data", {}).get("actives_csv")
    active_ids: List[str] = []
    if actives_path:
        import pandas as _pd
        from pathlib import Path as _Path
        ap = _Path(actives_path)
        if ap.exists():
            act_df = _pd.read_csv(ap)
            id_col = next(
                (c for c in ["zinc_id", "ligand_id", "id", "ZINC_ID"] if c in act_df.columns),
                act_df.columns[0],
            )
            active_ids = act_df[id_col].astype(str).tolist()
            logger.info("Loaded %d known actives from %s", len(active_ids), actives_path)
        else:
            logger.warning("actives_csv path not found: %s", actives_path)

    # Build scaffold series from native_scored if available (has murcko_scaffold col)
    scaffold_series: Optional[pd.Series] = None
    if "native_scored" in tables:
        ns = tables["native_scored"]
        if "murcko_scaffold" in ns.columns:
            # Best microstate per ligand
            best = ns.sort_values("native_weighted_coverage_pct", ascending=False)
            best = best.drop_duplicates(subset=["ligand_id"], keep="first")
            scaffold_series = best.set_index("ligand_id")["murcko_scaffold"]

    # Build all baseline rankings
    rankings = build_all_rankings(
        tables,
        random_seed=random_seed,
        mw_bins=mw_bins,
        logp_bins=logp_bins,
    )

    # Pairwise comparison
    logger.info("Computing pairwise rank metrics ...")
    pairwise_df = _pairwise_metrics(rankings, top_k_values)

    # Per-method summary
    logger.info("Computing per-method enrichment / diversity summary ...")
    summary_df = _method_summary(rankings, active_ids, top_k_values, scaffold_series)

    return BenchmarkResult(
        rankings=rankings,
        pairwise_df=pairwise_df,
        summary_df=summary_df,
        scaffold_series=scaffold_series,
    )
