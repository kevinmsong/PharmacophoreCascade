"""
Load and validate all pipeline result tables.

Usage:
    tables = load_tables("evidence/configs/benchmark.yaml")
    stage3_df = tables["stage3"]
    final_df  = tables["final_ranked"]
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import yaml

from .schemas import (
    FEATURE_MATCHES,
    FINAL_RANKED,
    HOTSPOT_COMPAT,
    NATIVE_SCORED,
    RANK_COMPARISON,
    REFERENCE_FEATURES,
    STAGE3_AUDIT,
    validate,
)

logger = logging.getLogger(__name__)

# Keys returned in the tables dict
TABLE_KEYS = [
    "stage3",
    "final_ranked",
    "native_scored",
    "feature_matches",
    "reference_features",
    "rank_comparison",
    "hotspot_compat",
]


def _resolve(path_str: Optional[str], project_root: Path) -> Optional[Path]:
    """Resolve a path relative to the project root, or None if not given."""
    if not path_str:
        return None
    p = Path(path_str)
    if not p.is_absolute():
        p = project_root / p
    return p


def _load_csv(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV with informative error on failure."""
    if not path.exists():
        raise FileNotFoundError(
            f"[loader] Required table '{label}' not found at: {path}"
        )
    logger.info("Loading %s from %s", label, path)
    df = pd.read_csv(path, low_memory=False)
    logger.info("  -> %d rows, %d columns", len(df), len(df.columns))
    return df


def _cache_path(output_dir: Path, key: str) -> Path:
    return output_dir / "cache" / f"{key}.parquet"


def _maybe_load_cache(cache_file: Path, source_file: Path) -> Optional[pd.DataFrame]:
    """Return cached parquet if it is newer than the source CSV."""
    if cache_file.exists() and cache_file.stat().st_mtime >= source_file.stat().st_mtime:
        logger.info("Loading %s from cache", cache_file.stem)
        return pd.read_parquet(cache_file)
    return None


def _save_cache(df: pd.DataFrame, cache_file: Path) -> None:
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_file, index=False)


def load_tables(
    config_path: str | Path,
    use_cache: bool = True,
) -> Dict[str, pd.DataFrame]:
    """Load and validate all result tables described in benchmark.yaml.

    Parameters
    ----------
    config_path:
        Path to benchmark.yaml (absolute or relative to cwd).
    use_cache:
        If True, load from parquet cache when source CSV is unchanged.

    Returns
    -------
    dict with keys:
        stage3, final_ranked, native_scored, feature_matches,
        reference_features, rank_comparison, hotspot_compat (optional).
    """
    config_path = Path(config_path).resolve()
    # evidence/configs/benchmark.yaml → parent=evidence/configs → parent=evidence → parent=project root
    project_root = config_path.parent.parent.parent

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    data_cfg = cfg.get("data", {})
    output_dir = project_root / cfg.get("output_dir", "evidence/outputs")

    paths = {
        "stage3":             _resolve(data_cfg.get("stage3_audit"), project_root),
        "final_ranked":       _resolve(data_cfg.get("final_ranked"), project_root),
        "native_scored":      _resolve(data_cfg.get("native_scored"), project_root),
        "feature_matches":    _resolve(data_cfg.get("feature_matches"), project_root),
        "reference_features": _resolve(data_cfg.get("reference_features"), project_root),
        "rank_comparison":    _resolve(data_cfg.get("rank_comparison"), project_root),
        "hotspot_compat":     _resolve(data_cfg.get("hotspot_compat"), project_root),
    }

    schemas_map = {
        "stage3":             STAGE3_AUDIT,
        "final_ranked":       FINAL_RANKED,
        "native_scored":      NATIVE_SCORED,
        "feature_matches":    FEATURE_MATCHES,
        "reference_features": REFERENCE_FEATURES,
        "rank_comparison":    RANK_COMPARISON,
        "hotspot_compat":     HOTSPOT_COMPAT,
    }

    tables: Dict[str, pd.DataFrame] = {}

    for key, path in paths.items():
        if path is None:
            logger.warning("No path configured for table '%s', skipping", key)
            continue

        cache_file = _cache_path(output_dir, key)

        if use_cache:
            cached = _maybe_load_cache(cache_file, path)
            if cached is not None:
                tables[key] = validate(cached, schemas_map[key])
                continue

        df = _load_csv(path, key)
        df = validate(df, schemas_map[key])

        if use_cache:
            _save_cache(df, cache_file)

        tables[key] = df

    _log_summary(tables)
    return tables


def _log_summary(tables: Dict[str, pd.DataFrame]) -> None:
    logger.info("Loaded tables:")
    for key, df in tables.items():
        logger.info("  %-22s  %7d rows", key, len(df))
