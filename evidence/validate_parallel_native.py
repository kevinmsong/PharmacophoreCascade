"""Numerically validate chunked native preparation/scoring against serial routines."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from evidence import parallel_native as parallel


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-dir',type=Path,default=ROOT/'evidence/outputs/alert_disabled_revision/parallel_validation')
    args=ap.parse_args()
    ligands,native=parallel.modules()
    source=ROOT/'tmp/chart_lock_smoke_results/screening_floor_smoke_native_terminal_bundle/input/native_selected_input.csv'
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    opts=dict(input_csv=source,id_column='zinc_id',smiles_column='smiles',assay_ph=7.4,
              max_tautomers_per_microstate=8,max_microstates_per_ligand=2)
    serial=ligands.prepare_library(**opts,output_table=out/'serial.csv',output_failures=out/'serial_fail.csv',output_sdf=out/'serial.sdf')
    chunked=parallel.prepare_parallel(2,**opts,output_table=out/'parallel.csv',output_failures=out/'parallel_fail.csv',output_sdf=out/'parallel.sdf')
    pd.testing.assert_frame_equal(serial.records,chunked.records,check_dtype=False,atol=1e-12,rtol=1e-12)
    pd.testing.assert_frame_equal(serial.failures,chunked.failures)
    assert (out/'serial.sdf').read_bytes()==(out/'parallel.sdf').read_bytes(),'Prepared SDF differs'
    ref=native.build_reference_features(ROOT/'GLP1_top_ligand_analysis/6X18_GLP1_GLP1R.pdb','R','P',4.5)
    opts=dict(prepared_df=serial.records,prepared_sdf=out/'serial.sdf',reference_features=ref,pair_tolerance=3.0,embed_conformers=4,progress_every=None)
    a=native.score_prepared_native_library(**opts)
    b=parallel.score_parallel(2,**opts)
    for field in ['scored_df','best_df','match_df','scoring_failures_df']:
        aa=getattr(a,field).reset_index(drop=True);bb=getattr(b,field).reset_index(drop=True)
        pd.testing.assert_frame_equal(aa,bb,check_dtype=False,atol=1e-12,rtol=1e-12)
    assert a.mapping_cache==b.mapping_cache
    for key in a.mapping_cache:
        for ca,cb in zip(a.molecules[key].GetConformers(),b.molecules[key].GetConformers()):
            np.testing.assert_allclose(ca.GetPositions(),cb.GetPositions(),atol=1e-12,rtol=0)
    (out/'validation.json').write_text(json.dumps({'status':'passed','ligands':int(serial.records.ligand_id.nunique()),
        'prepared_microstates':len(serial.records),'prepared_sdf_byte_identical':True,'numeric_tolerance':1e-12,
        'structural_alert_exclusions':False,'cached_scores_used':False,
        'source_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'parallel_source_sha256':hashlib.sha256((ROOT/'evidence/parallel_native.py').read_bytes()).hexdigest()},indent=2))
    print(f'PASS: {len(serial.records)} prepared microstates, identical SDF, scores/mappings/conformers agree to 1e-12',flush=True)


if __name__=='__main__': main()
