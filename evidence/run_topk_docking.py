#!/usr/bin/env python3
"""
Redock the current top-ranked GLP-1R final-ranked ligands against active and
inactive receptor structures with AutoDock Vina (reviewer points 9-11).

Reuses the prepared receptor PDBQT files and the per-structure docking boxes
from the prior docking study, generates 3D conformers + Meeko PDBQT for the
current top-N final-ranked ligands, runs Vina with fully documented parameters,
and aggregates per-ligand best affinity by receptor state. Active-state
preference is reported only for ligands docked successfully against both states.
"""
from __future__ import annotations

import subprocess
import json
import hashlib
import time
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "GLP1_top_ligand_analysis" / "src"))
from glp1r_state_preference import docking as dk  # noqa: E402

RECEPTORS = ROOT / "GLP1_top_ligand_analysis/derived_715_topological/receptors"
OUT = ROOT / "evidence/outputs/docking_top10"
SEEDS = [11, 29, 47]
EXHAUSTIVENESS = 16
NUM_MODES = 9


def main() -> None:
    started=time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    top = pd.read_csv(ROOT / "evidence/outputs/revision/top10_for_docking.csv")
    boxes = pd.read_csv(ROOT / "evidence/outputs/revision/docking_boxes.csv")
    assert len(top)==10 and top.ligand_id.is_unique
    for b in boxes.itertuples():
        if not (RECEPTORS/f'{b.conformer_id}.pdbqt').exists():
            raise FileNotFoundError(f'Missing fixed receptor: {b.conformer_id}')

    # 1) 3D conformers -> SDF (title = ligand_id) -> Meeko PDBQT
    sdf = OUT / "top10_ligands.sdf"
    w = Chem.SDWriter(str(sdf))
    prepared = 0
    prep_records=[]
    for _, r in top.iterrows():
        m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
        if m is None:
            prep_records.append({'ligand_id':r['ligand_id'],'status':'invalid_smiles'})
            continue
        m = Chem.AddHs(m)
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260313
        if AllChem.EmbedMolecule(m, params) != 0:
            prep_records.append({'ligand_id':r['ligand_id'],'status':'embedding_failed'})
            continue
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=400)
        except Exception:
            pass
        m.SetProp("_Name", str(r["ligand_id"]))
        w.write(m)
        prepared += 1
        prep_records.append({'ligand_id':r['ligand_id'],'status':'conformer_prepared'})
    w.close()
    print(f"prepared {prepared} ligand conformers")
    lig_dir = OUT / "ligands_pdbqt"
    dk.convert_sdf_to_pdbqt_dir(sdf, lig_dir)
    for record in prep_records:
        if record['status']=='conformer_prepared':
            record['status']='prepared' if (lig_dir/f"{record['ligand_id']}.pdbqt").exists() else 'pdbqt_preparation_failed'
    pd.DataFrame(prep_records).to_csv(OUT/'ligand_preparation_status.csv',index=False)

    vina = dk.ensure_binary("vina")
    rows = []
    n_jobs = len(boxes) * len(top) * len(SEEDS)
    done = 0
    for _, b in boxes.iterrows():
        rec = RECEPTORS / f"{b['conformer_id']}.pdbqt"
        for _, r in top.iterrows():
            lig = lig_dir / f"{r['ligand_id']}.pdbqt"
            for seed in SEEDS:
                if not lig.exists():
                    rows.append({'structure_id':b['structure_id'],'state':b['state'],
                        'ligand_id':r['ligand_id'],'final_rank':int(r['analysis_rank']),
                        'seed':seed,'best_affinity':np.nan,'status':'ligand_preparation_failed','returncode':None})
                    done+=1
                    continue
                jobdir = OUT / "runs" / f"{b['structure_id']}__{r['ligand_id']}__seed{seed}"
                jobdir.mkdir(parents=True, exist_ok=True)
                conf = jobdir / "vina.conf"
                outp = jobdir / "out.pdbqt"
                logp = jobdir / "vina.log"
                dk.write_vina_config(
                    conf, rec, lig,
                    [float(b["center_x"]), float(b["center_y"]), float(b["center_z"])],
                    [float(b["size_x"]), float(b["size_y"]), float(b["size_z"])],
                    EXHAUSTIVENESS, NUM_MODES, int(seed), outp, logp,
                )
                res = subprocess.run([vina, "--config", str(conf), "--cpu", "8"],
                                     capture_output=True, text=True)
                logp.write_text(res.stdout + "\n" + res.stderr, encoding="utf-8")
                poses = dk.parse_vina_log(logp)
                best = min((p["affinity_kcal_mol"] for p in poses), default=np.nan)
                rows.append({"structure_id": b["structure_id"], "state": b["state"],
                             "ligand_id": r["ligand_id"], "final_rank": int(r["analysis_rank"]),
                             "seed": int(seed), "best_affinity": best,
                             "status":"scored" if np.isfinite(best) and res.returncode==0 else "docking_failed",
                             "returncode":res.returncode})
                done += 1
                if done % 10 == 0:
                    print(f"  {done}/{n_jobs} jobs", flush=True)
    df = pd.DataFrame(rows)
    assert len(df)==n_jobs
    df.to_csv(OUT / "docking_results.csv", index=False)
    if not (df.status=='scored').any():
        raise RuntimeError('No successful docking jobs; inspect preparation and Vina logs')

    # 2) aggregate per ligand: best over active / inactive structures
    agg = []
    for lig, g in df.groupby("ligand_id"):
        rank = int(g["final_rank"].iloc[0])
        ba = g.loc[(g["state"] == "active")&(g.status=='scored'), "best_affinity"].min()
        bi = g.loc[(g["state"] == "inactive")&(g.status=='scored'), "best_affinity"].min()
        agg.append({"ligand_id": lig, "final_rank": rank, "best_active": ba, "best_inactive": bi,
                    "active_pref": (ba - bi) if (pd.notna(ba) and pd.notna(bi)) else np.nan})
    adf = pd.DataFrame(agg).sort_values("final_rank")
    adf.to_csv(ROOT / "evidence/outputs/revision/docking_top10_summary.csv", index=False)
    print("\n=== per-ligand summary ===")
    print(adf.round(2).to_string(index=False))
    both = adf.dropna(subset=["active_pref"])
    inputs=[ROOT/'evidence/outputs/revision/top10_for_docking.csv',ROOT/'evidence/outputs/revision/docking_boxes.csv']
    inputs += [RECEPTORS/f'{b.conformer_id}.pdbqt' for b in boxes.itertuples()]
    manifest={'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'vina_version':subprocess.run([vina,'--version'],capture_output=True,text=True,check=True).stdout.strip(),
        'ligand_ids':top.ligand_id.tolist(),'seeds':SEEDS,'exhaustiveness':EXHAUSTIVENESS,
        'num_modes':NUM_MODES,'etkdg_seed':20260313,'mmff_max_iters':400,
        'requested_jobs':n_jobs,'job_status_counts':df.status.value_counts().to_dict(),
        'ligands_with_both_states':len(both),'cached_ligand_scores_used':False,'elapsed_seconds':time.perf_counter()-started}
    (OUT/'rerun_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(f"\nligands docked in both states: {len(both)}/{len(adf)}")
    if len(both):
        print(f"median best active affinity: {both['best_active'].median():.2f}; "
              f"median active-state preference: {both['active_pref'].median():.2f}; "
              f"prefer active: {(both['active_pref']<0).sum()}/{len(both)}")


if __name__ == "__main__":
    main()
