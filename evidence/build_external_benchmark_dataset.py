#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from src.benchmark_library_builder import (
    build_benchmark_library,
    curate_in_domain_glp1r_actives,
    load_zinc_candidate_pool,
    match_decoys_from_pool,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build the external GLP-1R benchmark library.")
    parser.add_argument("--output-dir", default=str(HERE / "data"))
    parser.add_argument("--scaffold-limit", type=int, default=10)
    parser.add_argument("--decoys-per-active", type=int, default=30)
    parser.add_argument("--max-zinc-candidates", type=int, default=250000)
    parser.add_argument("--zinc-dir", default=str(ROOT / "tmp" / "zinc"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    actives_df, exclusions_df = curate_in_domain_glp1r_actives(
        scaffold_limit=args.scaffold_limit
    )
    zinc_dir = Path(args.zinc_dir)
    zinc_paths = sorted(zinc_dir.glob("*.smi.gz"))
    if not zinc_paths:
        raise FileNotFoundError(f"No ZINC tranche files found under {zinc_dir}")

    candidate_pool_df = load_zinc_candidate_pool(
        zinc_paths,
        max_candidates=args.max_zinc_candidates,
    )
    decoys_df = match_decoys_from_pool(
        actives_df=actives_df,
        candidate_pool_df=candidate_pool_df,
        decoys_per_active=args.decoys_per_active,
    )
    benchmark_df = build_benchmark_library(actives_df, decoys_df)

    benchmark_path = output_dir / "glp1r_external_benchmark_library.csv"
    exclusion_path = output_dir / "glp1r_external_benchmark_exclusions.csv"
    benchmark_df.to_csv(benchmark_path, index=False)
    exclusions_df.to_csv(exclusion_path, index=False)

    print(f"Wrote benchmark library: {benchmark_path}")
    print(f"Wrote exclusions table: {exclusion_path}")
    print(
        f"Actives: {int((benchmark_df['label'] == 'active').sum())} | "
        f"Decoys: {int((benchmark_df['label'] == 'decoy').sum())}"
    )


if __name__ == "__main__":
    main()
