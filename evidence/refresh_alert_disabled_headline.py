"""Promote a completed fresh headline and redock its actual new top 10."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd
import yaml
from run_alert_disabled_suite import archive

ROOT=Path(__file__).resolve().parents[1]
HEAD=ROOT/'results/absolute_floor_1000'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--wait',action='store_true');args=ap.parse_args()
    while True:
        status=json.loads((HEAD/'rerun_provenance.json').read_text())['status']
        if status=='completed':break
        if status!='running' or not args.wait:raise RuntimeError(f'Headline status: {status}')
        time.sleep(30)
    manifest=META/'headline_promotion.json'
    if manifest.exists():raise RuntimeError('Headline already promoted; inspect current docking manifest before repeating')
    aliases={
        'screening_full_1M_floor1000.csv':'screening_full_1M_topological_hashed.csv',
        'screening_full_1M_floor1000_shortlist.csv':'screening_full_1M_topological_hashed_shortlist.csv',
        'screening_full_1M_floor1000_native_scored_top5000.csv':'screening_full_1M_topological_hashed_native_scored_top5000.csv',
        'top_1000_glp1_mimetics_full_1M_floor1000_native_final.csv':'top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv',
        'screening_full_1M_floor1000_run_summary.json':'screening_full_1M_topological_hashed_run_summary.json'}
    records={}
    for src,dst in aliases.items():
        source=HEAD/src;target=ROOT/'results'/dst
        assert (ROOT/'tmp/alert_disabled_revision/prior_headline'/dst).exists(),'Missing prior-result backup'
        shutil.copy2(source,target)
        records[dst]=hashlib.sha256(target.read_bytes()).hexdigest()
    # Analysis consumers read the actual new bundle directly; no old top-1000
    # diagnostic study is allowed to masquerade as the new full native cohort.
    bundle='results/absolute_floor_1000/screening_full_1M_floor1000_native_terminal_bundle/analysis/'
    for config in (ROOT/'evidence/configs').glob('benchmark*.yaml'):
        value=yaml.safe_load(config.read_text())
        data=value.get('data',{})
        for key,file in [('native_scored','microstate_native_mapping.csv'),('feature_matches','microstate_native_feature_matches.csv'),
            ('reference_features','native_reference_features.csv'),('rank_comparison','native_rank_comparison.csv'),
            ('hotspot_compat','ligand_hotspot_compatibility.csv')]:
            if key in data:data[key]=bundle+file
        config.write_text(yaml.safe_dump(value,sort_keys=False))
    final=pd.read_csv(HEAD/'top_1000_glp1_mimetics_full_1M_floor1000_native_final.csv').sort_values('final_rank').head(10)
    smiles='canonical_smiles' if 'canonical_smiles' in final else 'smiles'
    top=pd.DataFrame({'analysis_rank':final.final_rank,'ligand_id':final.zinc_id,'canonical_smiles':final[smiles]})
    assert len(top)==10 and top.ligand_id.is_unique
    revision=ROOT/'evidence/outputs/revision'
    for path in [revision/'top10_for_docking.csv',revision/'docking_top10_summary.csv',ROOT/'evidence/outputs/docking_top10']:
        archive(path)
    top.to_csv(revision/'top10_for_docking.csv',index=False)
    manifest.write_text(json.dumps({'status':'completed','canonical_alias_sha256':records,
        'headline_source':str(HEAD.relative_to(ROOT)),'new_top10':top.ligand_id.tolist()},indent=2))
    print('Fresh headline promoted; starting all 150 new top-10 docking jobs',flush=True)
    with open(META/'fresh_top10_docking.log','w',encoding='utf-8') as log:
        subprocess.run([sys.executable,'-u','evidence/run_topk_docking.py'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Fresh top-10 docking completed',flush=True)


if __name__=='__main__':main()
