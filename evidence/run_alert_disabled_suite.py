"""Resumable, fresh rerun of every analysis experiment used by the ACS revision.

Two independent lanes permit the full benchmarks and production pairs to execute
beside the headline run. Completed historical folders are moved out of the active
outputs before a job starts. Resume trusts only this revision's completed manifest.
"""
import argparse
import copy
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime,timezone
from pathlib import Path
import yaml

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'tmp/alert_disabled_revision'
META=ROOT/'evidence/outputs/alert_disabled_revision'
SPECS={
    'glp1r':('benchmark_glp1r_full.yaml',None,None),
    'ghsr':('benchmark_ghsr_full.yaml','study_ghsr_benchmark.yaml',''),
    'ntsr1':('benchmark_ntsr1_full.yaml','study_ntsr1_benchmark.yaml',''),
    'mdm2':('benchmark_mdm2_full.yaml','study_mdm2_benchmark.yaml',''),
}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def archive(path):
    path=path.resolve()
    assert path.is_relative_to(ROOT.resolve()) and path!=ROOT.resolve()
    if not path.exists():return
    destination=(WORK/'prior_outputs'/path.relative_to(ROOT)).resolve()
    assert destination.is_relative_to(WORK.resolve())
    if destination.exists(): raise RuntimeError(f'Archive already exists; inspect interrupted run: {path}')
    destination.parent.mkdir(parents=True,exist_ok=True)
    path.rename(destination)


def job(name,command,outputs,inputs,settings):
    META.mkdir(parents=True,exist_ok=True)
    manifest=META/f'{name}.json'
    if manifest.exists():
        old=json.loads(manifest.read_text())
        if old['status']=='completed':
            print(f'Already completed: {name}',flush=True);return
        raise RuntimeError(f'Inspect unfinished job before restart: {manifest}')
    for path in outputs:archive(path)
    record={'job':name,'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),
            'command':command,'settings':settings,'cached_scores_used':False,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [
                ROOT/'run_optimized_1M_topological_hashed_screening.py',
                ROOT/'evidence/src/benchmark_runner.py',ROOT/'evidence/parallel_native.py',Path(__file__)]}}
    manifest.write_text(json.dumps(record,indent=2))
    started=time.perf_counter()
    print(f'Starting {name}',flush=True)
    with open(META/f'{name}.log','w',encoding='utf-8') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    record.update(status='completed' if result.returncode==0 else 'failed',
                  returncode=result.returncode,elapsed_seconds=time.perf_counter()-started,
                  completed_utc=datetime.now(timezone.utc).isoformat())
    manifest.write_text(json.dumps(record,indent=2))
    if result.returncode:raise RuntimeError(f'{name} failed; see its log')
    print(f'Completed {name} in {record["elapsed_seconds"]/60:.1f} min',flush=True)


def benchmark(system,workers,variant=None,seed=None):
    base,native,groups=SPECS[system]
    cfg=yaml.safe_load((ROOT/'evidence/configs'/base).read_text())
    name=f'benchmark_{system}_full'
    if variant=='automated':
        cfg=yaml.safe_load((ROOT/'evidence/configs/benchmark_glp1r_automated.yaml').read_text())
        groups='';name='benchmark_glp1r_automated'
    if variant=='7ki0':native='study_glp1r_7ki0.yaml';name='benchmark_glp1r_7ki0'
    out=ROOT/'evidence/outputs'/name
    if seed is not None:
        name=f'{system}_seed{seed}'
        out=ROOT/'evidence/outputs/decoy_replicates'/name
        cfg['benchmark']['library_csv']=f'evidence/data/replicates/{system}_decoyset{seed}.csv'
        cfg['benchmark']['include_native_only']=False
        cfg['benchmark']['bootstrap_iterations']=0
    settings=cfg['benchmark']
    settings.update(chemistry_gate_mode='warn_only',workers=workers,native_workers=workers,
                    cascade_hotspot_weight=.25,feature_cap_mode='fixed',
                    output_dir=str((out/'benchmark_external').relative_to(ROOT)))
    cfg_folder=META/'configs';cfg_folder.mkdir(parents=True,exist_ok=True)
    cfg_path=cfg_folder/f'{name}.yaml'
    cfg_path.write_text(yaml.safe_dump(cfg,sort_keys=False))
    cmd=[sys.executable,'-u','evidence/run_target_benchmark.py','--benchmark-config',str(cfg_path),'--output-dir',str(out)]
    native_path=ROOT/'GLP1_top_ligand_analysis/configs'/(native or 'study_top1000_full_1M_native_tuned_no_docking.yaml')
    cmd+=['--native-rerank-config',str(native_path)]
    if groups is not None:cmd+=['--required-groups',groups]
    job(name,cmd,[out],[cfg_path,ROOT/settings['library_csv'],ROOT/settings['pharmacophore_path'],native_path],
        {'chemistry_gate_mode':'warn_only','native_pains_filter':False,'workers':workers,
         'cascade_hotspot_weight':.25,'feature_caps':'fixed','required_groups':groups,
         'variant':variant,'decoy_seed':seed})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lane',choices=['benchmarks','production'],required=True)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    META.mkdir(parents=True,exist_ok=True)
    if args.lane=='production':
        for system in ['ghsr','ntsr1','glp1r','mdm2']:
            cmd=[sys.executable,'-u','evidence/run_absolute_floor_benchmarks.py','--system',system,'--workers',str(args.workers)]
            job(f'production_{system}',cmd,[ROOT/f'evidence/outputs/absolute_floor/{system}'],
                [ROOT/f'evidence/data/production_frozen/{system}.csv',ROOT/f'evidence/configs/benchmark_{system}_full.yaml'],
                {'chemistry_gate_mode':'warn_only','native_pains_filter':False,'workers':args.workers,'policies':['5% only','5% + floor 1000']})
    else:
        for system in SPECS:benchmark(system,args.workers)
        for variant in ['automated','7ki0']:benchmark('glp1r',args.workers,variant=variant)
        for system in SPECS:
            for seed in range(1,6):benchmark(system,args.workers,seed=seed)


if __name__=='__main__':main()
