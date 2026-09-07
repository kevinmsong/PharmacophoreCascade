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
import json

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
        manifest=ROOT/f'evidence/outputs/alert_disabled_revision/benchmark_{system}_full.json'
        state=json.loads(manifest.read_text())
        assert state['status']=='completed' and state['settings']['chemistry_gate_mode']=='warn_only'
        lib = pd.read_csv(lib_path)
        assert lib.ligand_id.is_unique
        cols = [c for c in KEEP if c in lib.columns]
        rec = lib[cols].copy()
        n_methods = 0
        for method in METHODS:
            rk = bench_dir / f"{method}_ranking.csv"
            r = pd.read_csv(rk)[["ligand_id", "rank", "score", "status"]].rename(
                columns={"rank": f"{method}_rank", "score": f"{method}_score", "status": f"{method}_status"}
            )
            rec = rec.merge(r, on="ligand_id", how="left",validate='one_to_one')
            ties=pd.read_csv(bench_dir.parent/'benchmark_tie_groups.csv')
            ties=ties[ties.method==method][['ligand_id','evaluation_tie_group','evaluation_rank']].rename(
                columns={'evaluation_tie_group':f'{method}_tie_group','evaluation_rank':f'{method}_evaluation_rank'})
            rec=rec.merge(ties,on='ligand_id',how='left',validate='one_to_one')
            assert rec[f'{method}_evaluation_rank'].notna().all()
            assert rec[f'{method}_rank'].notna().all(),f'Incomplete ranks: {system}/{method}'
            n_methods += 1
        evaluation=pd.read_csv(bench_dir/'benchmark_evaluation.csv')
        gate_columns=[c for c in ['ligand_id','topology_status','chemistry_flagged','pains_alert','reactive_flags'] if c in evaluation]
        rec=rec.merge(evaluation[gate_columns],on='ligand_id',how='left',validate='one_to_one')
        rec['chemistry_gate_mode']='warn_only'
        rec['native_pains_filter']=False
        out_path = OUT / f"{system}_benchmark_scored.csv"
        rec.to_csv(out_path, index=False)
        print(f"{system}: {len(rec)} molecules, {n_methods} methods -> {out_path.name}")


if __name__ == "__main__":
    main()
