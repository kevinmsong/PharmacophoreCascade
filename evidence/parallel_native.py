"""Process-parallel execution of unchanged, deterministic per-ligand native routines.

Each worker calls the original preparation/scoring functions. Ordered chunks are
joined, and best-microstate selection is repeated globally. No conformer count,
random seed, feature matching, score, or selection rule is changed.
"""
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from functools import partial
import math
import csv
import re
from pathlib import Path
import shutil
import pandas as pd


def modules():
    import run_optimized_1M_topological_hashed_screening as engine
    ligands=engine.glp1_analysis_import('glp1r_state_preference.ligands')
    native=engine.glp1_analysis_import('glp1r_state_preference.native_pharmacophore')
    from rdkit import Chem
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps | Chem.PropertyPickleOptions.CoordsAsDouble)
    return ligands,native


def prepare_worker(kwargs):
    ligands,_=modules()
    result=ligands.prepare_library(**kwargs)
    return result,kwargs['output_sdf']


def prepare_parallel(workers,**kwargs):
    ligands,_=modules()
    with open(kwargs['input_csv'],newline='',encoding='utf-8') as handle:
        reader=csv.reader(handle);header=next(reader);source=list(reader)
    if kwargs.get('row_limit') is not None: source=source[:kwargs['row_limit']]
    if not source:
        return ligands._serial_prepare_library(**kwargs)
    folder=Path(kwargs['output_table']).parent/'parallel_prepare'
    folder.mkdir(parents=True,exist_ok=True)
    jobs=[]
    size=max(1,min(64,math.ceil(len(source)/workers)))
    for i,start in enumerate(range(0,len(source),size)):
        chunk=folder/f'chunk_{i:04d}';chunk.mkdir(exist_ok=True)
        input_csv=chunk/'input.csv'
        # Preserve original numeric text exactly; a pandas read/write/read
        # round trip can otherwise perturb metadata in the last decimal place.
        with open(input_csv,'w',newline='',encoding='utf-8') as handle:
            writer=csv.writer(handle);writer.writerow(header);writer.writerows(source[start:start+size])
        job=dict(kwargs,input_csv=input_csv,output_table=chunk/'prepared.csv',
            output_failures=chunk/'failures.csv',output_sdf=chunk/'prepared.sdf',progress_every=None)
        job.pop('row_limit',None);jobs.append(job)
    prepared=[];failures=[];sdfs=[];stats=Counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,(result,sdf) in enumerate(pool.map(prepare_worker,jobs),1):
            prepared.append(result.records);failures.append(result.failures);sdfs.append(sdf)
            stats.update(result.stats)
            print(f'[progress] Parallel native preparation: {i}/{len(jobs)} chunks',flush=True)
    records=pd.concat(prepared,ignore_index=True)
    failed=pd.concat(failures,ignore_index=True)
    records.to_csv(kwargs['output_table'],index=False);failed.to_csv(kwargs['output_failures'],index=False)
    with open(kwargs['output_sdf'],'wb') as out:
        offset=0
        for sdf in sdfs:
            data=Path(sdf).read_bytes()
            # SDWriter numbers property headers within each file. Restore the
            # global record numbers so concatenation is byte-identical to serial.
            data=re.sub(rb'(>  <[^>]+>  \()(\d+)(\)[ \t]*\r?\n)',
                lambda m:m[1]+str(int(m[2])+offset).encode()+m[3],data)
            out.write(data);offset+=data.count(b'$$$$')
    return ligands.LigandPreparationResult(records=records,failures=failed,stats=dict(stats))


def score_worker(kwargs):
    _,native=modules()
    return native.score_prepared_native_library(**kwargs)


def score_parallel(workers,**kwargs):
    _,native=modules()
    from rdkit import Chem
    frame=kwargs['prepared_df']
    molecules=native.molecule_lookup(kwargs['prepared_sdf'])
    folder=Path(kwargs['prepared_sdf']).parent/'parallel_scoring'
    folder.mkdir(parents=True,exist_ok=True)
    size=max(1,min(128,math.ceil(len(frame)/workers)))
    jobs=[]
    for i,start in enumerate(range(0,len(frame),size)):
        sub=frame.iloc[start:start+size].copy()
        sdf=folder/f'chunk_{i:04d}.sdf'
        writer=Chem.SDWriter(str(sdf))
        for mid in sub.microstate_id.astype(str):
            if mid in molecules: writer.write(molecules[mid])
        writer.close()
        jobs.append(dict(kwargs,prepared_df=sub,prepared_sdf=sdf,progress_every=None))
    scored=[];matches=[];failures=[];mapping={};updated={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,result in enumerate(pool.map(score_worker,jobs),1):
            scored.append(result.scored_df);matches.append(result.match_df);failures.append(result.scoring_failures_df)
            mapping.update(result.mapping_cache);updated.update(result.molecules)
            print(f'[progress] Parallel native scoring: {i}/{len(jobs)} chunks',flush=True)
    scored_df=pd.concat(scored,ignore_index=True)
    order=[('source_input_rank',True),('ligand_id',True),('native_weighted_coverage_pct',False),
           ('matched_reference_features',False),('microstate_id',True)]
    order=[(k,asc) for k,asc in order if k in scored_df]
    scored_df=scored_df.sort_values([k for k,a in order],ascending=[a for k,a in order]).reset_index(drop=True)
    return native.NativeLibraryScoringResult(molecules=updated,mapping_cache=mapping,
        scored_df=scored_df,match_df=pd.concat(matches,ignore_index=True),
        best_df=native.select_best_microstates(scored_df),scoring_failures_df=pd.concat(failures,ignore_index=True))


def install(workers):
    if workers<=1: return
    ligands,native=modules()
    if not hasattr(ligands,'_serial_prepare_library'):
        ligands._serial_prepare_library=ligands.prepare_library
        native._serial_score_prepared_native_library=native.score_prepared_native_library
    ligands.prepare_library=partial(prepare_parallel,workers)
    native.score_prepared_native_library=partial(score_parallel,workers)
