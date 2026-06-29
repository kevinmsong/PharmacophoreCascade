"""
Run the external active-vs-decoy benchmark on an arbitrary benchmark CSV.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Dict

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import run_optimized_1M_topological_hashed_screening as screening


@dataclass
class ExternalBenchmarkArtifacts:
    library_df: pd.DataFrame
    evaluation_df: pd.DataFrame
    method_tables: Dict[str, pd.DataFrame]
    rankings: Dict[str, pd.DataFrame]
    output_dir: Path


def _ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column not in out.columns:
            out[column] = np.nan
    return out


def _resolve(path_str: str | None, root: Path) -> Path | None:
    if not path_str:
        return None
    path = Path(path_str)
    if not path.is_absolute():
        path = root / path
    return path


def load_benchmark_library(config: dict, project_root: Path) -> pd.DataFrame:
    bench_cfg = config.get("benchmark", {})
    library_path = _resolve(bench_cfg.get("library_csv"), project_root)
    if library_path is None or not library_path.exists():
        raise FileNotFoundError(
            f"Benchmark library CSV not found: {library_path}"
        )
    library_df = pd.read_csv(library_path)
    required = {"ligand_id", "smiles", "label", "murcko_scaffold"}
    missing = required.difference(library_df.columns)
    if missing:
        raise ValueError(f"Benchmark library missing required columns: {sorted(missing)}")
    if "in_domain" not in library_df.columns:
        library_df["in_domain"] = 1
    return library_df.drop_duplicates("ligand_id").reset_index(drop=True)


def _status_rank(status: str) -> int:
    order = {
        "scored": 0,
        "stage3_scored": 0,
        "aligned": 0,
        "active_state_docked": 0,
        "state_preference_docked": 0,
        "native_failed": 1,
        "stage3_failed": 2,
        "alignment_failed": 2,
        "chemistry_filtered": 3,
        "native_chemistry_excluded": 3,
        "hotspot_filtered": 4,
        "property_filtered": 5,
        "invalid": 6,
        "missing_docking": 7,
    }
    return order.get(str(status), 99)


def _ranking_from_scores(
    library_df: pd.DataFrame,
    method_df: pd.DataFrame,
    *,
    method_name: str,
    score_col: str,
    status_col: str,
    ascending: bool = False,
) -> pd.DataFrame:
    merged = library_df[["ligand_id", "label"]].merge(
        method_df[["ligand_id", score_col, status_col]],
        on="ligand_id",
        how="left",
    )
    merged["method"] = method_name
    merged["status"] = merged[status_col].fillna("missing")
    merged["score"] = pd.to_numeric(merged[score_col], errors="coerce")
    merged["success"] = merged["status"].map(_status_rank).eq(0)
    merged["_status_rank"] = merged["status"].map(_status_rank)
    merged = merged.sort_values(
        ["_status_rank", "score", "ligand_id"],
        ascending=[True, ascending, True],
        na_position="last",
        kind="stable",
    ).reset_index(drop=True)
    merged["rank"] = np.arange(1, len(merged) + 1, dtype=int)
    return merged[["ligand_id", "rank", "score", "method", "status", "label"]]


def _build_full_cascade_args(config: dict, library_size: int) -> SimpleNamespace:
    bench_cfg = config.get("benchmark", {})
    return SimpleNamespace(
        top_hotspots=int(bench_cfg.get("top_hotspots", 25)),
        pair_features=int(bench_cfg.get("pair_features", 24)),
        rerank_query_features=int(bench_cfg.get("rerank_query_features", 28)),
        pair_hash_mode=str(bench_cfg.get("pair_hash_mode", "precision_5bin")),
        native_support_max_residues=int(bench_cfg.get("native_support_max_residues", 8)),
        hotspot_min_exact=int(bench_cfg.get("hotspot_min_exact", 3)),
        hotspot_min_groups=int(bench_cfg.get("hotspot_min_groups", 2)),
        chemistry_gate_mode=str(bench_cfg.get("chemistry_gate_mode", "strict")),
        workers=int(bench_cfg.get("workers", 1)),
        rerank_conformers=int(bench_cfg.get("rerank_conformers", 16)),
        rerank_tolerance=float(bench_cfg.get("rerank_tolerance", 2.75)),
        rerank_pair_tolerance=float(bench_cfg.get("rerank_pair_tolerance", 2.75)),
        rerank_score_mode=str(bench_cfg.get("rerank_score_mode", "stage3_only")),
        native_selection_top_k=int(bench_cfg.get("native_selection_top_k", library_size)),
        native_candidate_pool_k=int(bench_cfg.get("native_candidate_pool_k", library_size)),
        native_final_top_k=int(bench_cfg.get("native_final_top_k", library_size)),
        native_rerank_max_per_scaffold=int(bench_cfg.get("native_rerank_max_per_scaffold", library_size)),
        native_rerank_pair_tolerance=float(bench_cfg.get("native_rerank_pair_tolerance", 3.0)),
        final_rank_mode=str(bench_cfg.get("final_rank_mode", "native_first")),
    )


def _build_pseudo_stage3_rows(library_df: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for _, row in library_df.iterrows():
        rows.append(
            {
                "zinc_id": row["ligand_id"],
                "smiles": row["smiles"],
                "mw": row["mw"],
                "logp": row["logp"],
                "hbd": row["hbd"],
                "hba": row["hba"],
                "rotatable_bonds": row["rotatable_bonds"],
                "lipinski": 1,
                "hotspot_weighted_pct": 0.0,
                "hotspot_bits_matched": 0,
                "hotspot_exact_matches": 0,
                "hotspot_compatible_matches": 0,
                "hotspot_exact_residue_count": 0,
                "hotspot_group_count": 0,
                "hotspot_required_groups_pass": 0,
                "pair_hash_overlap_pct": 0.0,
                "pair_hash_recall_pct": 0.0,
                "pair_hash_precision_pct": 0.0,
                "cascade_score_pct": 0.0,
                "shortlist_rank": int(row["benchmark_input_rank"]),
            }
        )
    return rows


def _evaluate_benchmark_library(
    library_df: pd.DataFrame,
    config: dict,
    project_root: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    args = _build_full_cascade_args(config, len(library_df))
    pharmacophore_path = _resolve(
        config.get("benchmark", {}).get("pharmacophore_path", "maps/pharmacophore_rigorous.json"),
        project_root,
    )
    pharmacophore = screening.load_pharmacophore(str(pharmacophore_path))
    views = screening.precompute_pharmacophore_views(
        pharmacophore,
        top_hotspots=args.top_hotspots,
        pair_features=args.pair_features,
        rerank_query_features=args.rerank_query_features,
        pair_hash_mode=args.pair_hash_mode,
        native_support_max_residues=args.native_support_max_residues,
    )

    eval_rows = []
    for _, row in library_df.iterrows():
        record, status, timings = screening.evaluate_topological_candidate_with_timings(
            row["smiles"],
            row["ligand_id"],
            views,
            hotspot_min_exact=args.hotspot_min_exact,
            hotspot_min_groups=args.hotspot_min_groups,
            pair_hash_mode=args.pair_hash_mode,
            chemistry_gate_mode=args.chemistry_gate_mode,
        )
        payload = {
            "ligand_id": row["ligand_id"],
            "smiles": row["smiles"],
            "benchmark_label": row["label"],
            "benchmark_input_rank": int(row["benchmark_input_rank"]),
            "topology_status": status,
            "stage0_property_gate_sec": float(timings["stage0_property_gate_sec"]),
            "stage1_hotspot_scoring_sec": float(timings["stage1_hotspot_scoring_sec"]),
            "stage2_pair_hash_scoring_sec": float(timings["stage2_pair_hash_scoring_sec"]),
        }
        if record is not None:
            payload.update(record)
        eval_rows.append(payload)
    return pd.DataFrame(eval_rows), views


def _rerank_stage3(
    candidate_df: pd.DataFrame,
    *,
    views: dict[str, object],
    workers: int,
    rerank_conformers: int,
    rerank_tolerance: float,
    rerank_pair_tolerance: float,
    rerank_score_mode: str,
) -> pd.DataFrame:
    if candidate_df.empty:
        return pd.DataFrame(columns=["ligand_id", "stage3_score", "stage3_status"])

    rows = candidate_df.to_dict(orient="records")
    reranked = screening.rerank_shortlist_parallel(
        rows,
        views,
        workers=workers,
        rerank_conformers=rerank_conformers,
        tolerance=rerank_tolerance,
        pair_tolerance=rerank_pair_tolerance,
        rerank_score_mode=rerank_score_mode,
        report_interval=max(100, len(rows)),
    )
    stage3_df = pd.DataFrame(reranked, columns=screening.final_output_columns())
    if stage3_df.empty:
        return pd.DataFrame(columns=["ligand_id", "stage3_score", "stage3_status"])
    stage3_df = screening.sort_stage3_final_hits(stage3_df)
    stage3_df["ligand_id"] = stage3_df["zinc_id"]
    stage3_df["stage3_score"] = stage3_df["weighted_coverage_pct"]
    stage3_df["stage3_status"] = "stage3_scored"
    stage3_df["stage3_rank"] = np.arange(1, len(stage3_df) + 1, dtype=int)
    return stage3_df


def _standard_3d_baseline(
    library_df: pd.DataFrame,
    views: dict[str, object],
    args: SimpleNamespace,
) -> pd.DataFrame:
    valid_rows = []
    status_rows = []
    for _, row in library_df.iterrows():
        mol = screening.standardize_screening_molecule(row["smiles"])
        if mol is None:
            status_rows.append({"ligand_id": row["ligand_id"], "standard3d_score": np.nan, "standard3d_status": "invalid"})
            continue
        props = screening.calculate_properties(mol)
        chemistry = screening.chemistry_alerts_for_molecule(mol)
        if not screening.passes_property_gate(props):
            status_rows.append({"ligand_id": row["ligand_id"], "standard3d_score": np.nan, "standard3d_status": "property_filtered"})
            continue
        if args.chemistry_gate_mode == "strict" and int(chemistry["chemistry_flagged"]) == 1:
            status_rows.append({"ligand_id": row["ligand_id"], "standard3d_score": np.nan, "standard3d_status": "chemistry_filtered"})
            continue

        valid_rows.append(
            {
                "ligand_id": row["ligand_id"],
                "smiles": row["smiles"],
                "mw": props["mw"],
                "logp": props["logp"],
                "hbd": props["hbd"],
                "hba": props["hba"],
                "rotatable_bonds": props["rotatable_bonds"],
                "benchmark_input_rank": int(row["benchmark_input_rank"]),
            }
        )

    baseline_df = pd.DataFrame(valid_rows)
    if baseline_df.empty:
        return pd.DataFrame(status_rows)

    pseudo_rows = _build_pseudo_stage3_rows(baseline_df)
    reranked = screening.rerank_shortlist_parallel(
        pseudo_rows,
        views,
        workers=args.workers,
        rerank_conformers=args.rerank_conformers,
        tolerance=args.rerank_tolerance,
        pair_tolerance=args.rerank_pair_tolerance,
        rerank_score_mode="stage3_only",
        report_interval=max(100, len(pseudo_rows)),
    )
    scored_df = pd.DataFrame(reranked)
    if scored_df.empty:
        scored_df = pd.DataFrame(columns=["zinc_id", "weighted_coverage_pct"])
    scored_df["ligand_id"] = scored_df.get("zinc_id")
    scored_df["standard3d_score"] = scored_df.get("weighted_coverage_pct")
    scored_df["standard3d_status"] = np.where(
        scored_df["standard3d_score"].notna(), "aligned", "alignment_failed"
    )
    out = pd.concat(
        [
            scored_df[["ligand_id", "standard3d_score", "standard3d_status"]],
            pd.DataFrame(status_rows),
        ],
        ignore_index=True,
    )
    return out.drop_duplicates("ligand_id", keep="first")


def _augment_full_cascade_failures(
    library_df: pd.DataFrame,
    evaluation_df: pd.DataFrame,
    native_final_df: pd.DataFrame,
) -> pd.DataFrame:
    native_df = native_final_df.copy()
    if native_df.empty:
        native_df = pd.DataFrame(columns=["zinc_id", "native_weighted_coverage_pct"])
    native_df["ligand_id"] = native_df.get("zinc_id")
    native_df["full_cascade_score"] = pd.to_numeric(
        native_df.get("native_weighted_coverage_pct"), errors="coerce"
    )
    native_df["full_cascade_status"] = np.where(
        native_df["full_cascade_score"].notna(), "scored", "native_failed"
    )

    missing = library_df[~library_df["ligand_id"].isin(native_df["ligand_id"])].copy()
    if not missing.empty:
        missing = missing.merge(
            evaluation_df[["ligand_id", "topology_status"]],
            on="ligand_id",
            how="left",
        )
        missing["full_cascade_score"] = np.nan
        missing["full_cascade_status"] = missing["topology_status"].fillna("native_failed")
        native_df = pd.concat(
            [
                native_df,
                missing[["ligand_id", "full_cascade_score", "full_cascade_status"]],
            ],
            ignore_index=True,
        )
    return native_df.drop_duplicates("ligand_id", keep="first")


def _aggregate_active_state_docking(docking_df: pd.DataFrame) -> pd.DataFrame:
    required = {"ligand_id", "state", "structure_id", "score"}
    missing = required.difference(docking_df.columns)
    if missing:
        raise ValueError(f"Docking results missing required columns: {sorted(missing)}")

    active_df = docking_df.loc[docking_df["state"].astype(str).str.lower() == "active"].copy()
    if active_df.empty:
        return pd.DataFrame(columns=["ligand_id", "active_state_docking_score", "active_state_docking_status"])

    best_pose = (
        active_df.groupby(["ligand_id", "structure_id"], dropna=False)["score"]
        .min()
        .reset_index()
    )
    ranked = (
        best_pose.groupby("ligand_id", dropna=False)["score"]
        .median()
        .reset_index()
        .rename(columns={"score": "active_state_docking_score"})
    )
    ranked["active_state_docking_status"] = "active_state_docked"
    return ranked


def _aggregate_state_preference_docking(docking_df: pd.DataFrame) -> pd.DataFrame:
    from glp1r_state_preference.scoring import aggregate_state_preference

    ranked = aggregate_state_preference(docking_df, bootstrap_iterations=500)
    ranked = ranked.rename(
        columns={
            "state_preference_score": "state_preference_docking_score",
        }
    )
    ranked["state_preference_docking_status"] = "state_preference_docked"
    return ranked


def _native_only_table(
    library_df: pd.DataFrame,
    args: SimpleNamespace,
    output_dir: Path,
) -> pd.DataFrame:
    """Score the ENTIRE labeled library directly with the terminal native
    pharmacophore method, bypassing the Stage 0-2 gate and Stage 3.

    All molecules (not just gate survivors) are prepared and native-scored with
    the same ligand preparation and the same native reference as the cascade's
    terminal branch, so the resulting ranking isolates whether the upstream
    cascade adds value beyond the final native-scoring stage alone.
    """
    base = library_df.copy()
    base["zinc_id"] = base["ligand_id"].astype(str)
    base["benchmark_input_rank"] = np.arange(1, len(base) + 1, dtype=int)
    pseudo = pd.DataFrame(_build_pseudo_stage3_rows(base))
    pseudo = _ensure_columns(pseudo, screening.final_output_columns())
    # The diversified native pool sorts on these fields; with all molecules
    # carrying neutral values and pool caps = library size, every molecule is
    # admitted to the native pool.
    for col in (
        "native_supported_hotspot_weighted_pct",
        "hotspot_group_count",
        "hotspot_weighted_pct",
        "cascade_score_pct",
        "weighted_coverage_pct",
        "chemistry_flagged",
    ):
        if col in pseudo.columns:
            pseudo[col] = pd.to_numeric(pseudo[col], errors="coerce").fillna(0.0)
        else:
            pseudo[col] = 0.0
    # Give the Stage-3 fields a non-degenerate spread so the native rerank's
    # internal correlation/quartile diagnostics do not crash on identical values.
    # These placeholder Stage-3 scores do not affect the native-only ranking,
    # which is taken from native weighted coverage.
    n = len(pseudo)
    spread = np.linspace(20.0, 10.0, n) if n > 1 else np.array([15.0])
    pseudo["weighted_coverage_pct"] = spread
    pseudo["cascade_score_pct"] = spread
    pseudo.insert(0, "stage3_screen_rank", np.arange(1, n + 1, dtype=int))

    native_paths = {
        "native_scored": output_dir / "native_only_native_scored.csv",
        "native_final": output_dir / "native_only_final.csv",
        "native_bundle_root": output_dir / "native_only_bundle",
        "shortlist": output_dir / "native_only_shortlist.csv",
        "full": output_dir / "native_only_stage3_full.csv",
        "top100": output_dir / "native_only_top100.csv",
        "pharmacophore_plot": output_dir / "native_only_pharmacophore.png",
        "top20_plot": output_dir / "native_only_top20.png",
        "distributions_plot": output_dir / "native_only_distributions.png",
        "run_summary": output_dir / "native_only_run_summary.json",
    }
    native_final_df, _ = screening.run_native_terminal_rerank(
        stage3_df=pseudo, output_paths=native_paths, args=args
    )
    native_df = native_final_df.copy()
    if native_df.empty:
        native_df = pd.DataFrame(columns=["zinc_id", "native_weighted_coverage_pct"])
    native_df["ligand_id"] = native_df.get("zinc_id")
    native_df["native_only_score"] = pd.to_numeric(
        native_df.get("native_weighted_coverage_pct"), errors="coerce"
    )
    native_df["native_only_status"] = np.where(
        native_df["native_only_score"].notna(), "scored", "native_failed"
    )
    out = native_df[["ligand_id", "native_only_score", "native_only_status"]].copy()
    missing = library_df[~library_df["ligand_id"].isin(out["ligand_id"])].copy()
    if not missing.empty:
        missing["native_only_score"] = np.nan
        missing["native_only_status"] = "native_failed"
        out = pd.concat(
            [out, missing[["ligand_id", "native_only_score", "native_only_status"]]],
            ignore_index=True,
        )
    return out.drop_duplicates("ligand_id", keep="first")


def run_external_benchmark(config: dict, project_root: Path) -> ExternalBenchmarkArtifacts:
    library_df = load_benchmark_library(config, project_root).copy()
    library_df["benchmark_input_rank"] = np.arange(1, len(library_df) + 1, dtype=int)
    output_dir = _resolve(
        config.get("benchmark", {}).get("output_dir", "evidence/outputs/benchmark_external"),
        project_root,
    )
    assert output_dir is not None
    output_dir.mkdir(parents=True, exist_ok=True)

    args = _build_full_cascade_args(config, len(library_df))
    evaluation_df, views = _evaluate_benchmark_library(library_df, config, project_root)
    evaluation_path = output_dir / "benchmark_evaluation.csv"
    evaluation_df.to_csv(evaluation_path, index=False)

    candidate_df = evaluation_df.loc[evaluation_df["topology_status"].isin(["candidate", "warn_only_flagged"])].copy()
    if "zinc_id" not in candidate_df.columns:
        candidate_df["zinc_id"] = candidate_df["ligand_id"]
    candidate_df = _ensure_columns(candidate_df, screening.final_output_columns() + ["shortlist_rank"])
    if "shortlist_rank" not in candidate_df.columns or candidate_df["shortlist_rank"].isna().all():
        candidate_df["shortlist_rank"] = np.arange(1, len(candidate_df) + 1, dtype=int)

    stage3_df = _rerank_stage3(
        candidate_df,
        views=views,
        workers=args.workers,
        rerank_conformers=args.rerank_conformers,
        rerank_tolerance=args.rerank_tolerance,
        rerank_pair_tolerance=args.rerank_pair_tolerance,
        rerank_score_mode=args.rerank_score_mode,
    )
    stage3_table = library_df[["ligand_id"]].merge(
        stage3_df[["ligand_id", "stage3_score", "stage3_status"]],
        on="ligand_id",
        how="left",
    ).merge(
        evaluation_df[["ligand_id", "topology_status"]],
        on="ligand_id",
        how="left",
    )
    stage3_table["stage3_status"] = stage3_table["stage3_status"].fillna(stage3_table["topology_status"]).fillna("stage3_failed")

    standard3d_table = _standard_3d_baseline(library_df, views, args)

    benchmark_output_paths = {
        "native_scored": output_dir / "benchmark_native_scored.csv",
        "native_final": output_dir / "benchmark_full_cascade.csv",
        "native_bundle_root": output_dir / "benchmark_native_bundle",
        "shortlist": output_dir / "benchmark_shortlist.csv",
        "full": output_dir / "benchmark_stage3_full.csv",
        "top100": output_dir / "benchmark_stage3_top100.csv",
        "pharmacophore_plot": output_dir / "benchmark_pharmacophore.png",
        "top20_plot": output_dir / "benchmark_top20.png",
        "distributions_plot": output_dir / "benchmark_distributions.png",
        "run_summary": output_dir / "benchmark_run_summary.json",
    }
    stage3_for_native = stage3_df.copy()
    if not stage3_for_native.empty and "stage3_screen_rank" not in stage3_for_native.columns:
        stage3_for_native.insert(0, "stage3_screen_rank", np.arange(1, len(stage3_for_native) + 1, dtype=int))
    native_final_df, _ = screening.run_native_terminal_rerank(
        stage3_df=stage3_for_native,
        output_paths=benchmark_output_paths,
        args=args,
    )
    full_cascade_table = _augment_full_cascade_failures(library_df, evaluation_df, native_final_df)

    method_tables: Dict[str, pd.DataFrame] = {
        "full_cascade": full_cascade_table,
        "stage3_only": stage3_table[["ligand_id", "stage3_score", "stage3_status"]].copy(),
        "standard_3d_pharmacophore": standard3d_table[["ligand_id", "standard3d_score", "standard3d_status"]].copy(),
    }

    bench_cfg = config.get("benchmark", {})
    if bool(bench_cfg.get("include_native_only", False)):
        method_tables["native_only"] = _native_only_table(library_df, args, output_dir)

    docking_raw_path = _resolve(bench_cfg.get("docking_raw_csv"), project_root)
    if docking_raw_path and docking_raw_path.exists():
        docking_df = pd.read_csv(docking_raw_path)
        method_tables["active_state_docking"] = _aggregate_active_state_docking(docking_df)
        method_tables["state_preference_docking_sensitivity"] = _aggregate_state_preference_docking(docking_df)

    rankings: Dict[str, pd.DataFrame] = {
        "full_cascade": _ranking_from_scores(
            library_df,
            method_tables["full_cascade"],
            method_name="full_cascade",
            score_col="full_cascade_score",
            status_col="full_cascade_status",
            ascending=False,
        ),
        "stage3_only": _ranking_from_scores(
            library_df,
            method_tables["stage3_only"],
            method_name="stage3_only",
            score_col="stage3_score",
            status_col="stage3_status",
            ascending=False,
        ),
        "standard_3d_pharmacophore": _ranking_from_scores(
            library_df,
            method_tables["standard_3d_pharmacophore"],
            method_name="standard_3d_pharmacophore",
            score_col="standard3d_score",
            status_col="standard3d_status",
            ascending=False,
        ),
    }

    if "native_only" in method_tables:
        rankings["native_only"] = _ranking_from_scores(
            library_df,
            method_tables["native_only"],
            method_name="native_only",
            score_col="native_only_score",
            status_col="native_only_status",
            ascending=False,
        )

    if "active_state_docking" in method_tables:
        rankings["active_state_docking"] = _ranking_from_scores(
            library_df,
            method_tables["active_state_docking"],
            method_name="active_state_docking",
            score_col="active_state_docking_score",
            status_col="active_state_docking_status",
            ascending=True,
        )
    if "state_preference_docking_sensitivity" in method_tables:
        rankings["state_preference_docking_sensitivity"] = _ranking_from_scores(
            library_df,
            method_tables["state_preference_docking_sensitivity"],
            method_name="state_preference_docking_sensitivity",
            score_col="state_preference_docking_score",
            status_col="state_preference_docking_status",
            ascending=False,
        )

    for method_name, ranking_df in rankings.items():
        ranking_df.to_csv(output_dir / f"{method_name}_ranking.csv", index=False)

    return ExternalBenchmarkArtifacts(
        library_df=library_df,
        evaluation_df=evaluation_df,
        method_tables=method_tables,
        rankings=rankings,
        output_dir=output_dir,
    )
