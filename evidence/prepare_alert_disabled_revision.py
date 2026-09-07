"""Freeze input libraries and switch current study configs to non-excluding alerts."""
from pathlib import Path
import hashlib
import json
import shutil
from datetime import datetime,timezone
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'tmp/alert_disabled_revision/before'


def main():
    ARCHIVE.mkdir(parents=True,exist_ok=True)
    files=list((ROOT/'evidence/configs').glob('*.yaml'))
    files+=list((ROOT/'GLP1_top_ligand_analysis/configs').glob('*.yaml'))
    files+=list((ROOT/'ACS_Omega_resubmission').glob('*.tex'))
    for path in files:
        target=ARCHIVE/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists(): raise RuntimeError(f'Backup exists: {target}')
        shutil.copy2(path,target)
    frozen=ROOT/'evidence/data/production_frozen'
    frozen.mkdir(exist_ok=True)
    for system in ['glp1r','ghsr','ntsr1','mdm2']:
        source=ROOT/f'evidence/outputs/production_{system}/stage012_evaluation.csv'
        frame=pd.read_csv(source,usecols=['ligand_id','smiles','label'])
        frame.to_csv(frozen/f'{system}.csv',index=False)
    changes=[]
    for path in files:
        if path.suffix!='.yaml': continue
        text=path.read_text(encoding='utf-8')
        updated=text.replace('chemistry_gate_mode: strict','chemistry_gate_mode: warn_only')
        updated=updated.replace('pains_filter: true','pains_filter: false')
        if updated!=text:
            path.write_text(updated,encoding='utf-8')
            changes.append(str(path.relative_to(ROOT)))
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),
        'user_direction':'Disable structural-alert exclusions; retain standardization and property limits. Rerun all analyses and replace current deliverables.',
        'stage0_mode':'warn_only','native_pains_filter':False,
        'shortlist_policy':'5% with minimum 1000, capped by eligible candidates',
        'updated_configs':changes,
        'frozen_input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in frozen.glob('*.csv')}}
    (ROOT/'tmp/alert_disabled_revision/revision_settings.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
