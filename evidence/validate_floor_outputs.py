"""Independent membership and count audit of completed paired floor outputs."""
import argparse
import json
import math
from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'evidence/outputs/absolute_floor'
STAGES=['input','stage0_pass','stage012_pass','shortlist','stage3_success',
        'native_selected','native_success','final_ranked']


def audit(system):
    folder=DATA/system
    provenance=json.loads((folder/'provenance.json').read_text())
    assert provenance['status']=='completed',f'{system} incomplete'
    library=pd.read_csv(folder/'input_library.csv')
    assert library.ligand_id.is_unique
    universe=set(library.ligand_id)
    table=pd.read_csv(folder/'stage_survival.csv')
    evaluation=pd.read_csv(folder/'stage012_evaluation.csv')
    eligible=set(evaluation.loc[evaluation.gate=='candidate','ligand_id'])
    assert set(evaluation.ligand_id)==universe
    source_names={'shortlist':'shortlist.csv','stage3_success':'stage3_scored.csv',
                  'native_selected':'native_selected.csv','native_success':'native_scored.csv',
                  'final_ranked':'final_ranked.csv'}
    policies={}
    for policy in ['percentage','floor1000']:
        target=folder/policy
        per=pd.read_csv(target/'per_molecule_survival.csv')
        assert per.ligand_id.is_unique and set(per.ligand_id)==universe
        assert_frame_equal(library.sort_values('ligand_id').reset_index(drop=True),
                           per[library.columns].sort_values('ligand_id').reset_index(drop=True))
        prior=universe
        counts=table[table.policy==policy].set_index('stage')
        for stage in STAGES:
            assert set(per[stage].unique()).issubset({True,False})
            ids=set(per.loc[per[stage],'ligand_id'])
            assert ids<=prior,f'{system}/{policy}/{stage} contains a molecule lost earlier'
            if stage=='stage012_pass': assert ids==eligible
            if stage in source_names:
                source=pd.read_csv(target/source_names[stage])
                assert source.zinc_id.is_unique
                if stage=='native_success': source=source[source.native_weighted_coverage_pct.notna()]
                assert ids==set(source.zinc_id),f'{system}/{policy}/{stage} membership mismatch'
            measured=per.loc[per[stage],'label'].value_counts()
            assert int(counts.loc[stage,'total'])==len(ids)
            for label,col in [('active','actives'),('decoy','decoys'),('background','background')]:
                assert int(counts.loc[stage,col])==int(measured.get(label,0))
            prior=ids
        shortlist=pd.read_csv(target/'shortlist.csv')
        requested=max(math.ceil(.05*len(eligible)),1000 if policy=='floor1000' else 1)
        assert len(shortlist)==min(len(eligible),requested)
        policies[policy]=shortlist.zinc_id.tolist()
    a,b=policies['percentage'],policies['floor1000']
    assert a==b[:len(a)],f'{system} floor is not an extension of the same ranking'
    shared=a==b
    if shared:
        for name in ['stage3_scored.csv','native_selected.csv','native_scored.csv','final_ranked.csv']:
            assert_frame_equal(pd.read_csv(folder/'percentage'/name),pd.read_csv(folder/'floor1000'/name))
    result={'system':system,'input_molecules':len(library),'shortlist_percentage':len(a),
            'shortlist_floor':len(b),'identical_policy_outputs':shared,
            'all_stage_memberships_and_label_counts_verified':True}
    (folder/'independent_output_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--systems',nargs='+',default=['glp1r','ghsr','ntsr1','mdm2'])
    args=parser.parse_args()
    for system in args.systems: audit(system)


if __name__=='__main__': main()
