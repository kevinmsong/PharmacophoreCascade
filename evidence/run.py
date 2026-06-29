#!/usr/bin/env python3
"""
GLP-1R Evidence Package — single entry point.

Usage
-----
    python evidence/run.py benchmark   # baseline comparison + rank correlations
    python evidence/run.py ablate      # ablations + sensitivity sweeps
    python evidence/run.py claims      # direct claim evidence tables
    python evidence/run.py report      # emit CSVs + Markdown + LaTeX + PDF
    python evidence/run.py all         # run all four in order

All outputs land in evidence/outputs/ (configurable via --output-dir).
All paths are resolved relative to the project root (parent of evidence/).
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

# Resolve project root so relative imports work regardless of cwd
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import yaml


def _setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )


def _load_yaml(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def _get_output_dir(config: dict, override: str | None) -> Path:
    if override:
        return Path(override).resolve()
    rel = config.get("output_dir", "evidence/outputs")
    return (_ROOT / rel).resolve()


# ---------------------------------------------------------------------------
# Sub-commands
# ---------------------------------------------------------------------------

def cmd_benchmark(args: argparse.Namespace) -> None:
    from src.loader import load_tables
    from src.benchmark import run_benchmark

    cfg_path = _HERE / "configs" / "benchmark.yaml"
    config = _load_yaml(cfg_path)
    output_dir = _get_output_dir(config, args.output_dir)

    logging.getLogger().info("Loading tables ...")
    tables = load_tables(cfg_path, use_cache=not args.no_cache)

    logging.getLogger().info("Running benchmark ...")
    bench = run_benchmark(tables, config)

    # Save pairwise + summary immediately (report step handles the rest)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not bench.summary_df.empty:
        bench.summary_df.to_csv(output_dir / "benchmark_summary.csv", index=False)
        logging.getLogger().info("Wrote benchmark_summary.csv")
    if not bench.delta_df.empty:
        bench.delta_df.to_csv(output_dir / "benchmark_deltas.csv", index=False)
        logging.getLogger().info("Wrote benchmark_deltas.csv")
    if not bench.pairwise_df.empty:
        bench.pairwise_df.to_csv(output_dir / "benchmark_pairwise.csv", index=False)
        logging.getLogger().info("Wrote benchmark_pairwise.csv")

    # Store bench + tables in a tmp file for the report step
    import pickle
    (output_dir / "cache").mkdir(exist_ok=True)
    with open(output_dir / "cache" / "_bench.pkl", "wb") as f:
        pickle.dump((bench, tables), f)

    logging.getLogger().info("Benchmark complete.")


def cmd_ablate(args: argparse.Namespace) -> None:
    from src.loader import load_tables
    from src.ablation import run_ablations

    bench_cfg_path = _HERE / "configs" / "benchmark.yaml"
    abl_cfg_path   = _HERE / "configs" / "ablation.yaml"
    bench_config = _load_yaml(bench_cfg_path)
    abl_config   = _load_yaml(abl_cfg_path)
    output_dir = _get_output_dir(bench_config, args.output_dir)

    tables = load_tables(bench_cfg_path, use_cache=not args.no_cache)

    logging.getLogger().info("Running ablations ...")
    ablation = run_ablations(tables, abl_config)

    output_dir.mkdir(parents=True, exist_ok=True)
    if not ablation.ablation_df.empty:
        ablation.ablation_df.to_csv(output_dir / "ablation_results.csv", index=False)
    if not ablation.sensitivity_df.empty:
        ablation.sensitivity_df.to_csv(output_dir / "sensitivity_results.csv", index=False)

    import pickle
    (output_dir / "cache").mkdir(exist_ok=True)
    with open(output_dir / "cache" / "_ablation.pkl", "wb") as f:
        pickle.dump(ablation, f)

    logging.getLogger().info("Ablation complete.")


def cmd_claims(args: argparse.Namespace) -> None:
    from src.loader import load_tables
    from src.claims import run_all_claims

    cfg_path = _HERE / "configs" / "benchmark.yaml"
    config = _load_yaml(cfg_path)
    output_dir = _get_output_dir(config, args.output_dir)

    tables = load_tables(cfg_path, use_cache=not args.no_cache)

    logging.getLogger().info("Evaluating manuscript claims ...")
    claims_results = run_all_claims(tables, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    from src.claims import claims_to_summary_df
    claims_to_summary_df(claims_results).to_csv(output_dir / "claims_summary.csv", index=False)

    for label, res in claims_results.items():
        edf = res["evidence_df"]
        if not edf.empty:
            edf.to_csv(output_dir / f"{label}_evidence.csv", index=False)

    import pickle
    (output_dir / "cache").mkdir(exist_ok=True)
    with open(output_dir / "cache" / "_claims.pkl", "wb") as f:
        pickle.dump(claims_results, f)

    logging.getLogger().info("Claims evaluation complete.")


def cmd_report(args: argparse.Namespace) -> None:
    import pickle
    from src.report import emit_report
    from src.benchmark import BenchmarkResult
    from src.ablation import AblationResult

    cfg_path = _HERE / "configs" / "benchmark.yaml"
    config = _load_yaml(cfg_path)
    output_dir = _get_output_dir(config, args.output_dir)
    cache_dir = output_dir / "cache"

    # Load pickled intermediates (run sub-commands first if missing)
    def _load_pkl(name, friendly):
        pkl = cache_dir / name
        if not pkl.exists():
            raise FileNotFoundError(
                f"{friendly} cache not found at {pkl}. "
                f"Run `python evidence/run.py {friendly.split()[0].lower()}` first."
            )
        with open(pkl, "rb") as f:
            return pickle.load(f)

    bench_data  = _load_pkl("_bench.pkl", "benchmark")
    bench, tables = bench_data

    ablation    = _load_pkl("_ablation.pkl", "ablation")
    claims_results = _load_pkl("_claims.pkl", "claims")

    logging.getLogger().info("Emitting report ...")
    emit_report(bench, ablation, claims_results, tables, output_dir)
    logging.getLogger().info("Report complete. Outputs in: %s", output_dir)


def cmd_all(args: argparse.Namespace) -> None:
    """Run all sub-commands in order."""
    for fn in [cmd_benchmark, cmd_ablate, cmd_claims, cmd_report]:
        t0 = time.time()
        fn(args)
        logging.getLogger().info("  -> %.1fs\n", time.time() - t0)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="GLP-1R Evidence Package",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "command",
        choices=["benchmark", "ablate", "claims", "report", "all"],
        help="Sub-command to run.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Override output directory (default: evidence/outputs).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable parquet cache; always reload from CSV.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )

    args = parser.parse_args()
    _setup_logging(args.log_level)

    dispatch = {
        "benchmark": cmd_benchmark,
        "ablate":    cmd_ablate,
        "claims":    cmd_claims,
        "report":    cmd_report,
        "all":       cmd_all,
    }
    dispatch[args.command](args)


if __name__ == "__main__":
    main()
