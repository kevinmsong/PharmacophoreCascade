"""Paired, freshly scored production benchmarks for a 1,000-ligand shortlist floor.

Only input identities/SMILES/labels are read from archived production evaluations.
All Stage 0-3 computations and native preparation/scoring are run afresh. Stage 0-2
and the Stage-3 union are shared between the paired policies; downstream selection
and native scoring are executed separately when shortlisted IDs differ.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"evidence/outputs/absolute_floor"
sys.path.insert(0,str(ROOT))
SYSTEMS=("ghsr","ntsr1","glp1r","mdm2")


def init_target_worker(views,args,caps,groups):
    import run_optimized_1M_topological_hashed_screening as s
    s.HOTSPOT_REQUIRED_GROUPS=tuple(groups)
    s._init_prescreen(views,args,caps)


def evaluate_target(supplier,n,views,args,caps,groups):
    import run_optimized_1M_topological_hashed_screening as s
    # Pass empty target groups explicitly: Windows may drop empty env values
    # when spawning a child interpreter.
    with Pool(args.prescreen_workers,initializer=init_target_worker,
              initargs=(views,args,caps,groups)) as pool:
        yield from pool.imap(s._prescreen_one,s.iter_supplier(supplier,n),chunksize=128)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_system(system,workers):
    out=OUT/system
    out.mkdir(parents=True,exist_ok=True)
    os.environ["CASCADIA_HOTSPOT_REQUIRED_GROUPS"]=("ECD anchoring;Upper TMD activation pocket" if system=="glp1r" else "")
    import run_optimized_1M_topological_hashed_screening as s
    s.HOTSPOT_REQUIRED_GROUPS=("ECD anchoring","Upper TMD activation pocket") if system=="glp1r" else ()
    cfg_path=ROOT/f"evidence/configs/benchmark_{system}_full.yaml"
    cfg=yaml.safe_load(cfg_path.read_text())
    bench=cfg["benchmark"]
    native_cfg=ROOT/("GLP1_top_ligand_analysis/configs/study_top1000_full_1M_native_tuned_no_docking.yaml" if system=="glp1r" else f"GLP1_top_ligand_analysis/configs/study_{system}_benchmark.yaml")
    s.DEFAULT_NATIVE_RERANK_CONFIG=native_cfg
    source=ROOT/f"evidence/data/production_frozen/{system}.csv"
    library=pd.read_csv(source)
    assert library.ligand_id.is_unique
    library.to_csv(out/"input_library.csv",index=False)
    supplier=out/"input_library.smi.gz"
    with gzip.open(supplier,"wt",encoding="utf-8") as f:
        for row in library.itertuples():
            f.write(f"{row.smiles} {row.ligand_id}\n")
    pharm=s.load_pharmacophore(ROOT/bench["pharmacophore_path"])
    # Preserve the historical benchmark's fixed per-type caps.
    s.TYPE_CAPS=s._BASE_TYPE_CAPS.copy()
    # Archived production benchmarks used 0.25/0.75, unlike the headline's
    # 0.4/0.6; pin their actual coefficients to isolate the floor intervention.
    s.CASCADE_STAGE_WEIGHTS={"hotspot":0.25,"pair_hash":0.75}
    views=s.precompute_pharmacophore_views(pharm,top_hotspots=25,pair_features=24,
        rerank_query_features=28,pair_hash_mode="precision_5bin",
        native_support_max_residues=int(bench.get("native_support_max_residues",0)))
    args=SimpleNamespace(prescreen_workers=workers,cascade_hotspot_weight=0.25,
        hotspot_min_exact=3,hotspot_min_groups=2,pair_hash_mode="precision_5bin",chemistry_gate_mode="warn_only")
    manifest={"status":"running","system":system,"started_utc":datetime.now(timezone.utc).isoformat(),
        "scientific_settings":{"shortlist_fraction":0.05,"shortlist_floor":1000,"shortlist_basis":"stage012_candidates",
        "cascade_weights":s.CASCADE_STAGE_WEIGHTS,"type_caps":s.TYPE_CAPS,"required_groups":s.HOTSPOT_REQUIRED_GROUPS,
        "stage3_conformers":16,"native_config":str(native_cfg.relative_to(ROOT)),"native_max_per_scaffold":8},
        "workers":workers,"native_workers":workers,"chemistry_gate_mode":"warn_only","native_pains_filter":False,
        "source_sha256":digest(source),"config_sha256":digest(cfg_path),
        "engine_sha256":digest(ROOT/"run_optimized_1M_topological_hashed_screening.py"),
        "native_config_sha256":digest(native_cfg),"pharmacophore_sha256":digest(ROOT/bench["pharmacophore_path"]),
        "cached_scores_used":False}
    (out/"provenance.json").write_text(json.dumps(manifest,indent=2))
    started=time.perf_counter()
    records=[]
    candidates=[]
    evaluations=evaluate_target(supplier,len(library),views,args,s.TYPE_CAPS,s.HOTSPOT_REQUIRED_GROUPS)
    for idx,(src,(record,status,timing)) in enumerate(zip(library.to_dict("records"),evaluations),1):
        row=dict(src,status=status,cascade_score_pct=0.0,gate=("candidate" if record is not None else "stage1_fail" if status=="hotspot_filtered" else "stage0_fail"))
        row.update(record or {})
        records.append(row)
        if record is not None:
            candidates.append(record)
        if idx%5000==0:
            print(f"{system} fresh Stage 0-2: {idx}/{len(library)}",flush=True)
    ev=pd.DataFrame(records)
    ev.to_csv(out/"stage012_evaluation.csv",index=False)
    stage012_sec=time.perf_counter()-started
    if (ev.status=='chemistry_filtered').any():
        raise RuntimeError('Structural-alert exclusion occurred despite the revised non-excluding setting')
    # Use the complete production tie-break contract; retain actual feature metadata.
    candidates.sort(key=lambda r:(-r['cascade_score_pct'],-r['hotspot_group_count'],
        -r['hotspot_weighted_pct'],-r['pair_hash_overlap_pct'],str(r['zinc_id'])))
    n_base=s.shortlist_count(len(candidates),5,0)
    n_floor=s.shortlist_count(len(candidates),5,1000)
    union=candidates[:n_floor]
    for i,r in enumerate(union,1): r['shortlist_rank']=i
    pd.DataFrame(union).to_csv(out/"stage3_union_input.csv",index=False)
    t=time.perf_counter()
    rows=s.rerank_shortlist_parallel(union,views,workers=workers,rerank_conformers=16,
        tolerance=2.75,pair_tolerance=2.75,rerank_score_mode="stage3_only",report_interval=100)
    union_df=s.sort_stage3_final_hits(pd.DataFrame(rows,columns=s.final_output_columns()))
    union_df.to_csv(out/"stage3_union_scored.csv",index=False)
    union_stage3_sec=time.perf_counter()-t
    summaries=[]
    all_survival=[]
    policy_outputs={}
    for policy,k in [("percentage",n_base),("floor1000",n_floor)]:
        folder=out/policy
        folder.mkdir(exist_ok=True)
        ids=set(r['zinc_id'] for r in candidates[:k])
        shortlist=pd.DataFrame(candidates[:k])
        shortlist.to_csv(folder/"shortlist.csv",index=False)
        stage3=union_df[union_df.zinc_id.isin(ids)].copy().reset_index(drop=True)
        stage3.insert(0,"stage3_screen_rank",np.arange(1,len(stage3)+1))
        stage3.to_csv(folder/"stage3_scored.csv",index=False)
        paths=s.resolve_output_paths(folder,"production")
        prod_args=SimpleNamespace(native_candidate_pool_k=20000,native_selection_top_k=5000,
            native_final_top_k=1000,native_rerank_max_per_scaffold=8,
            native_rerank_pair_tolerance=3.0,final_rank_mode="native_first",native_workers=workers)
        if policy=="floor1000" and n_floor==n_base:
            # Same IDs, ranks, and inputs: execute once and explicitly report shared results.
            native_final,native_summary,selected,scored,native_sec=policy_outputs['percentage']
            shared=True
        else:
            t=time.perf_counter()
            native_final,native_summary=s.run_native_terminal_rerank(stage3_df=stage3,output_paths=paths,args=prod_args)
            native_sec=time.perf_counter()-t
            selected=pd.read_csv(paths['native_bundle_root']/"input/native_selected_input.csv")
            scored=pd.read_csv(paths['native_scored'])
            shared=False
        policy_outputs[policy]=(native_final,native_summary,selected,scored,native_sec)
        selected.to_csv(folder/"native_selected.csv",index=False)
        scored.to_csv(folder/"native_scored.csv",index=False)
        native_final.to_csv(folder/"final_ranked.csv",index=False)
        success=set(scored.loc[scored.native_weighted_coverage_pct.notna(),'zinc_id'])
        stage_sets=[('input',set(library.ligand_id)),('stage0_pass',set(ev.loc[ev.gate.isin(['candidate','stage1_fail']),'ligand_id'])),
            ('stage012_pass',set(r['zinc_id'] for r in candidates)),('shortlist',ids),('stage3_success',set(stage3.zinc_id)),
            ('native_selected',set(selected.zinc_id)),('native_success',success),('final_ranked',set(native_final.zinc_id))]
        for stage,stage_ids in stage_sets:
            counts=library.loc[library.ligand_id.isin(stage_ids),'label'].value_counts()
            all_survival.append(dict(system=system,policy=policy,stage=stage,total=len(stage_ids),
                actives=int(counts.get('active',0)),decoys=int(counts.get('decoy',0)),background=int(counts.get('background',0))))
        molecule=library.copy()
        for stage,stage_ids in stage_sets: molecule[stage]=molecule.ligand_id.isin(stage_ids)
        molecule.to_csv(folder/"per_molecule_survival.csv",index=False)
        summaries.append(dict(system=system,policy=policy,n_input=len(library),n_candidates=len(candidates),n_shortlist=k,
            stage012_sec=stage012_sec,union_stage3_sec=union_stage3_sec,native_wall_sec=native_sec,
            native_results_shared_with_percentage=shared,native_summary=native_summary))
    pd.DataFrame(all_survival).to_csv(out/"stage_survival.csv",index=False)
    (out/"policy_summary.json").write_text(json.dumps(summaries,indent=2))
    manifest.update(status="completed",completed_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.perf_counter()-started)
    (out/"provenance.json").write_text(json.dumps(manifest,indent=2))
    print(pd.DataFrame(all_survival).to_string(index=False),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--system',choices=SYSTEMS)
    ap.add_argument('--workers',type=int,default=3)
    a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if a.system:
        run_system(a.system,a.workers)
    else:
        for system in SYSTEMS:
            with open(OUT/f'{system}.log','w',encoding='utf-8') as log:
                subprocess.run([sys.executable,'-u',__file__,'--system',system,'--workers',str(a.workers)],
                    cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)


if __name__=='__main__': main()
