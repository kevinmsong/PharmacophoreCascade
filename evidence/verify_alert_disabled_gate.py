"""Check the revised Stage-0 gate independently on all four labeled libraries."""
import sys
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'reference_implementation'))
import run_optimized_1M_topological_hashed_screening as engine
from stage_scoring_reference import stage0


def main():
    records=[]
    for system in ['glp1r','ghsr','ntsr1','mdm2']:
        source=pd.read_csv(ROOT/f'evidence/data/{system}_external_benchmark_library.csv')
        for row in source.itertuples():
            ref=stage0(row.smiles)
            mol=engine.standardize_screening_molecule(row.smiles)
            passed=mol is not None and engine.passes_property_gate(engine.calculate_properties(mol))
            assert ref['pass']==passed,f'Stage-0 disagreement: {system}/{row.ligand_id}'
            alert=engine.chemistry_alerts_for_molecule(mol) if mol is not None else {}
            records.append({'system':system,'ligand_id':row.ligand_id,'label':row.label,
                'property_gate_pass':passed,'chemistry_flagged':int(alert.get('chemistry_flagged',0)),
                'pains_alert':alert.get('pains_alert',''),'reactive_flags':alert.get('reactive_flags',''),
                'stage0_pass':passed,'independent_reference_agrees':True})
    frame=pd.DataFrame(records)
    out=ROOT/'evidence/outputs/alert_disabled_revision'
    out.mkdir(parents=True,exist_ok=True)
    frame.to_csv(out/'independent_stage0_check.csv',index=False)
    active=frame[frame.label=='active']
    print(active.groupby('system')[['stage0_pass','chemistry_flagged']].sum().to_string())
    print(f'Independent gate agreement: {len(frame)}/{len(frame)} labeled molecules')


if __name__=='__main__':main()
