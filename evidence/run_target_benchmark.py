#!/usr/bin/env python3
"""
Run the active-vs-decoy cascade benchmark for an arbitrary target or ablation.

This is a thin, reusable driver around ``src.benchmark.run_benchmark`` for the
external-benchmark path. It supports:

* the GLP-1R no-manual-curation ablation (automated pharmacophore + empty
  required-group gate), and
* additional peptide-receptor case studies (new library + automated interface
  pharmacophore + per-target native-reference complex).

It writes ``benchmark_summary.csv`` (point metrics + bootstrap CIs),
``benchmark_deltas.csv`` (paired deltas vs the full cascade) and
``benchmark_pairwise.csv`` to ``--output-dir``.

Usage
-----
    python evidence/run_target_benchmark.py \
        --benchmark-config evidence/configs/benchmark_glp1r_automated.yaml \
        --required-groups "" \
        --output-dir evidence/outputs/benchmark_glp1r_automated
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-config", required=True, help="Benchmark YAML (external-benchmark schema).")
    ap.add_argument(
        "--native-rerank-config",
        default=None,
        help="Study YAML (path relative to repo root) providing the native-branch "
        "reference complex/chains. Defaults to the GLP-1R study config.",
    )
    ap.add_argument(
        "--required-groups",
        default=None,
        help='Value for CASCADIA_HOTSPOT_REQUIRED_GROUPS (";"-separated). Pass "" to '
        "disable the required-group gate for automated/generic pharmacophores.",
    )
    ap.add_argument("--output-dir", default=None, help="Where to write summary CSVs.")
    args = ap.parse_args()

    # Must be set BEFORE the screening engine is imported (read at module load).
    if args.required_groups is not None:
        os.environ["CASCADIA_HOTSPOT_REQUIRED_GROUPS"] = args.required_groups

    import run_optimized_1M_topological_hashed_screening as screening

    if args.native_rerank_config:
        screening.DEFAULT_NATIVE_RERANK_CONFIG = Path(args.native_rerank_config)

    from src.benchmark import run_benchmark

    cfg_path = Path(args.benchmark_config)
    with cfg_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    print(f"[run_target_benchmark] config={cfg_path}")
    print(f"[run_target_benchmark] pharmacophore={config.get('benchmark', {}).get('pharmacophore_path')}")
    print(f"[run_target_benchmark] required_groups={tuple(screening.HOTSPOT_REQUIRED_GROUPS)}")
    print(f"[run_target_benchmark] native_rerank_config={screening.DEFAULT_NATIVE_RERANK_CONFIG}")

    # The external-benchmark path does not consult the loaded `tables` argument.
    bench = run_benchmark({}, config)

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        rel = config.get("benchmark", {}).get("output_dir", "evidence/outputs/benchmark_target")
        output_dir = (ROOT / rel).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not bench.summary_df.empty:
        bench.summary_df.to_csv(output_dir / "benchmark_summary.csv", index=False)
    if not bench.delta_df.empty:
        bench.delta_df.to_csv(output_dir / "benchmark_deltas.csv", index=False)
    if not bench.pairwise_df.empty:
        bench.pairwise_df.to_csv(output_dir / "benchmark_pairwise.csv", index=False)

    print(f"[run_target_benchmark] wrote summary to {output_dir}")
    cols = [c for c in ["method", "roc_auc", "pr_auc", "ef_1pct", "bedroc", "top10_recovery"] if c in bench.summary_df.columns]
    if cols:
        print(bench.summary_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
