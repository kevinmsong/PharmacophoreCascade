#!/usr/bin/env python3
"""
Build retrospective active-vs-decoy benchmark libraries for additional
peptide-receptor / protein-protein-interface case studies.

Reuses the (now target-agnostic) curation + decoy-matching machinery in
``src.benchmark_library_builder``. Actives are curated from ChEMBL for each
target; decoys are property-matched, scaffold-distinct, ECFP4-dissimilar
molecules drawn once from the local ZINC screening pool and shared across
targets. Outputs ``evidence/data/<target>_external_benchmark_library.csv`` and
``_exclusions.csv`` plus a small ``additional_benchmarks_summary.json``.

The ZINC pool load dominates runtime; it is performed once and reused.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.benchmark_library_builder import (
    build_benchmark_library,
    curate_in_domain_actives,
    load_zinc_candidate_pool,
    match_decoys_from_pool,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# Per-target curation specs. Peptide-GPCR agonists use functional EC50 with an
# agonist assay-term requirement; the MDM2-p53 protein-protein interface uses
# binding/inhibition potencies (IC50/Ki) with no agonist requirement.
TARGET_SPECS: dict[str, dict] = {
    "ghsr": {
        "name": "Ghrelin receptor (GHSR)",
        "chembl_target": "CHEMBL4616",
        "standard_types": ("EC50",),
        "default_potency_type": "EC50",
        "require_assay_terms": True,
        "assay_terms": ("agonist",),
        "excluded_terms": ("antagonist", "inverse agonist", "allosteric", "pam", "nam"),
        "potency_threshold_nM": 15000.0,
        "scaffold_limit": 10,
    },
    "ntsr1": {
        "name": "Neurotensin receptor 1 (NTSR1)",
        "chembl_target": "CHEMBL4123",
        "standard_types": ("EC50",),
        "default_potency_type": "EC50",
        "require_assay_terms": True,
        "assay_terms": ("agonist",),
        "excluded_terms": ("antagonist", "inverse agonist", "allosteric", "pam", "nam"),
        "potency_threshold_nM": 15000.0,
        "scaffold_limit": 10,
    },
    "mdm2": {
        "name": "MDM2-p53 interface",
        "chembl_target": "CHEMBL5023",
        "standard_types": ("IC50", "Ki"),
        "default_potency_type": "IC50",
        "require_assay_terms": False,
        "assay_terms": (),
        "excluded_terms": (),
        "potency_threshold_nM": 10000.0,
        "scaffold_limit": 10,
    },
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", nargs="+", default=list(TARGET_SPECS), choices=list(TARGET_SPECS))
    p.add_argument("--output-dir", default=str(HERE / "data"))
    p.add_argument("--decoys-per-active", type=int, default=30)
    p.add_argument("--max-zinc-candidates", type=int, default=200000)
    p.add_argument("--zinc-dir", default=str(ROOT / "tmp" / "zinc"))
    return p


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) Curate actives per target first (cheap, network), so we can skip the
    #    expensive ZINC load if nothing curated.
    actives: dict[str, "object"] = {}
    exclusions: dict[str, "object"] = {}
    for key in args.targets:
        spec = TARGET_SPECS[key]
        print(f"\n=== Curating actives: {spec['name']} ({spec['chembl_target']}) ===")
        act_df, exc_df = curate_in_domain_actives(
            spec["chembl_target"],
            standard_types=spec["standard_types"],
            default_potency_type=spec["default_potency_type"],
            require_assay_terms=spec["require_assay_terms"],
            assay_terms=spec["assay_terms"],
            excluded_terms=spec["excluded_terms"],
            scaffold_limit=spec["scaffold_limit"],
            potency_threshold_nM=spec["potency_threshold_nM"],
        )
        actives[key] = act_df
        exclusions[key] = exc_df
        print(f"  in-domain actives: {len(act_df)}")

    # 2) Load the ZINC decoy pool once, reused across all targets.
    zinc_dir = Path(args.zinc_dir)
    zinc_paths = sorted(zinc_dir.glob("*.smi.gz"))
    if not zinc_paths:
        raise FileNotFoundError(f"No ZINC tranche files under {zinc_dir}")
    print(f"\n=== Loading shared ZINC decoy pool (<= {args.max_zinc_candidates}) ===")
    candidate_pool_df = load_zinc_candidate_pool(zinc_paths, max_candidates=args.max_zinc_candidates)
    print(f"  ZINC in-domain candidates loaded: {len(candidate_pool_df)}")

    summary: dict[str, dict] = {}
    for key in args.targets:
        spec = TARGET_SPECS[key]
        act_df = actives[key]
        if act_df.empty:
            print(f"\n[skip] {spec['name']}: no in-domain actives.")
            summary[key] = {"name": spec["name"], "actives": 0, "decoys": 0, "status": "no_actives"}
            continue
        print(f"\n=== Matching decoys: {spec['name']} ({len(act_df)} actives) ===")
        decoys_df = match_decoys_from_pool(
            actives_df=act_df,
            candidate_pool_df=candidate_pool_df,
            decoys_per_active=args.decoys_per_active,
        )
        library_df = build_benchmark_library(act_df, decoys_df)
        lib_path = output_dir / f"{key}_external_benchmark_library.csv"
        exc_path = output_dir / f"{key}_external_benchmark_exclusions.csv"
        library_df.to_csv(lib_path, index=False)
        exclusions[key].to_csv(exc_path, index=False)
        n_act = int((library_df["label"] == "active").sum())
        n_dec = int((library_df["label"] == "decoy").sum())
        print(f"  wrote {lib_path.name}: actives={n_act}, decoys={n_dec}")
        summary[key] = {
            "name": spec["name"],
            "chembl_target": spec["chembl_target"],
            "actives": n_act,
            "decoys": n_dec,
            "library_csv": str(lib_path),
            "status": "ok",
        }

    (output_dir / "additional_benchmarks_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print("\n=== Summary ===")
    for key, info in summary.items():
        print(f"  {key:6s}: actives={info['actives']:3d}  decoys={info['decoys']:4d}  ({info['status']})")


if __name__ == "__main__":
    main()
