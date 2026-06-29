#!/usr/bin/env python3
"""
Export a machine-readable per-molecule record for each benchmark system.

For every benchmarked molecule, joins the curated library (canonical SMILES,
label, ChEMBL molecule/assay identifiers, property-matching distance, maximum
active Tanimoto, Murcko scaffold) with the per-method ranks and scores from the
headline benchmark run. Output: evidence/data/machine_readable/<system>_benchmark_scored.csv
(reviewer point 16; P5 transparency).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE / "data"
OUT = DATA / "machine_readable"

SYSTEMS = {
    "glp1r": (DATA / "glp1r_external_benchmark_library.csv", ROOT / "evidence/outputs/benchmark_glp1r_full/benchmark_external"),
    "ghsr": (DATA / "ghsr_external_benchmark_library.csv", ROOT / "evidence/outputs/benchmark_ghsr_full/benchmark_external"),
    "ntsr1": (DATA / "ntsr1_external_benchmark_library.csv", ROOT / "evidence/outputs/benchmark_ntsr1_full/benchmark_external"),
    "mdm2": (DATA / "mdm2_external_benchmark_library.csv", ROOT / "evidence/outputs/benchmark_mdm2_full/benchmark_external"),
}
METHODS = ["full_cascade", "native_only", "stage3_only", "standard_3d_pharmacophore"]
KEEP = ["ligand_id", "label", "canonical_smiles", "chembl_id", "pref_name",
        "potency_value", "potency_type", "assay_description", "murcko_scaffold",
        "matched_active_ligand_id", "property_distance", "max_active_tanimoto",
        "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for system, (lib_path, bench_dir) in SYSTEMS.items():
        if not lib_path.exists():
            continue
        lib = pd.read_csv(lib_path)
        cols = [c for c in KEEP if c in lib.columns]
        rec = lib[cols].copy()
        n_methods = 0
        for method in METHODS:
            rk = bench_dir / f"{method}_ranking.csv"
            if not rk.exists():
                continue
            r = pd.read_csv(rk)[["ligand_id", "rank", "score", "status"]].rename(
                columns={"rank": f"{method}_rank", "score": f"{method}_score", "status": f"{method}_status"}
            )
            rec = rec.merge(r, on="ligand_id", how="left")
            n_methods += 1
        out_path = OUT / f"{system}_benchmark_scored.csv"
        rec.to_csv(out_path, index=False)
        print(f"{system}: {len(rec)} molecules, {n_methods} methods -> {out_path.name}")


if __name__ == "__main__":
    main()
