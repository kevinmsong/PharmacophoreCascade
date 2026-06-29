#!/usr/bin/env python3
"""
Build N independent matched-decoy replicate libraries per system.

Reuses each system's already-curated actives (from the published benchmark
library) and the shared local ZINC pool, drawing a fresh, equally
property-matched decoy set per random seed via
``match_decoys_from_pool(random_state=seed)``. Outputs
``evidence/data/replicates/<system>_decoyset<seed>.csv`` for the
multiple-decoy-set robustness analysis (reviewer point 5).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.benchmark_library_builder import (
    build_benchmark_library,
    load_zinc_candidate_pool,
    match_decoys_from_pool,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE / "data"

SYSTEMS = {
    "glp1r": DATA / "glp1r_external_benchmark_library.csv",
    "ghsr": DATA / "ghsr_external_benchmark_library.csv",
    "ntsr1": DATA / "ntsr1_external_benchmark_library.csv",
    "mdm2": DATA / "mdm2_external_benchmark_library.csv",
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--systems", nargs="+", default=list(SYSTEMS), choices=list(SYSTEMS))
    p.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    p.add_argument("--decoys-per-active", type=int, default=30)
    p.add_argument("--max-zinc-candidates", type=int, default=150000)
    p.add_argument("--zinc-dir", default=str(ROOT / "tmp" / "zinc"))
    p.add_argument("--output-dir", default=str(DATA / "replicates"))
    return p


def main() -> None:
    args = build_parser().parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    zinc_paths = sorted(Path(args.zinc_dir).glob("*.smi.gz"))
    if not zinc_paths:
        raise FileNotFoundError(f"No ZINC tranche files under {args.zinc_dir}")
    print(f"Loading shared ZINC pool (<= {args.max_zinc_candidates}) ...")
    pool_df = load_zinc_candidate_pool(zinc_paths, max_candidates=args.max_zinc_candidates)
    print(f"  candidates: {len(pool_df)}")

    for system in args.systems:
        lib = pd.read_csv(SYSTEMS[system])
        actives = lib.loc[lib["label"].astype(str) == "active"].copy().reset_index(drop=True)
        print(f"\n=== {system}: {len(actives)} actives ===")
        for seed in args.seeds:
            decoys = match_decoys_from_pool(
                actives_df=actives,
                candidate_pool_df=pool_df,
                decoys_per_active=args.decoys_per_active,
                random_state=int(seed),
            )
            library = build_benchmark_library(actives, decoys)
            path = out_dir / f"{system}_decoyset{seed}.csv"
            library.to_csv(path, index=False)
            n_act = int((library["label"] == "active").sum())
            n_dec = int((library["label"] == "decoy").sum())
            print(f"  seed {seed}: actives={n_act} decoys={n_dec} -> {path.name}")


if __name__ == "__main__":
    main()
