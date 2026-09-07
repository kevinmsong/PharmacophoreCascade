"""Collect only completed revision outputs into one manuscript data object."""
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
META=ROOT/'evidence/outputs/alert_disabled_revision'
HEAD=ROOT/'results/absolute_floor_1000'


def json_safe(value):
    if isinstance(value,dict):return {key:json_safe(item) for key,item in value.items()}
    if isinstance(value,list):return [json_safe(item) for item in value]
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def completed(name):
    value=json.loads((META/(name+'.json')).read_text())
    if value['status']!='completed':raise RuntimeError(f'Incomplete experiment: {name}')
    return value


def check_benchmark_alert_settings(name, folder):
    execution=completed(name)
    assert execution['settings']['chemistry_gate_mode']=='warn_only'
    assert execution['settings']['native_pains_filter'] is False
    evaluation=pd.read_csv(folder/'benchmark_external/benchmark_evaluation.csv')
    assert not evaluation.topology_status.eq('chemistry_filtered').any(), name
    native=[]
    for bundle in ['benchmark_native_bundle','native_only_bundle']:
        path=folder/'benchmark_external'/bundle/'manifests/analysis_manifest.json'
        if not path.exists():
            continue
        manifest=json.loads(path.read_text(encoding='utf-8'))
        assert manifest['pains_filter_enabled'] is False, str(path)
        assert manifest['counts']['chemistry_excluded_selected_ligands']==0, str(path)
        native.append(str(path.relative_to(ROOT)))
    assert native, f'No actual native execution manifest for {name}'
    return dict(job=name, evaluated_molecules=len(evaluation),
                property_exclusions=int(evaluation.topology_status.eq('property_filtered').sum()),
                structural_alert_exclusions=0, native_execution_manifests=native)


def main():
    benchmarks={}
    alert_audit=[]
    for name in ['glp1r_full','ghsr_full','ntsr1_full','mdm2_full','glp1r_automated','glp1r_7ki0']:
        completed('benchmark_'+name)
        alert_audit.append(check_benchmark_alert_settings('benchmark_'+name,ROOT/f'evidence/outputs/benchmark_{name}'))
        assert json.loads((ROOT/f'evidence/outputs/benchmark_{name}/tie_aware_evaluation.json').read_text())['status']=='completed'
        table=pd.read_csv(ROOT/f'evidence/outputs/benchmark_{name}/benchmark_summary.csv')
        benchmarks[name]=table.set_index('method').to_dict(orient='index')
    replicates=[]
    for system in ['glp1r','ghsr','ntsr1','mdm2']:
        completed('production_'+system)
        for seed in range(1,6):
            completed(f'{system}_seed{seed}')
            alert_audit.append(check_benchmark_alert_settings(f'{system}_seed{seed}',ROOT/f'evidence/outputs/decoy_replicates/{system}_seed{seed}'))
            assert json.loads((ROOT/f'evidence/outputs/decoy_replicates/{system}_seed{seed}/tie_aware_evaluation.json').read_text())['status']=='completed'
            table=pd.read_csv(ROOT/f'evidence/outputs/decoy_replicates/{system}_seed{seed}/benchmark_summary.csv')
            for row in table.to_dict('records'):
                replicates.append(dict(system=system,seed=seed,**row))
    provenance=json.loads((HEAD/'rerun_provenance.json').read_text())
    if provenance['status']!='completed':raise RuntimeError('Headline not complete')
    head=json.loads((HEAD/'screening_full_1M_floor1000_run_summary.json').read_text())
    if provenance['chemistry_gate_mode']!='warn_only' or head['native_rerank']['pains_filter_enabled']:
        raise RuntimeError('Headline did not use the authorized alert settings')
    assert head['native_rerank']['counts']['chemistry_excluded_selected_ligands']==0
    (META/'alert_setting_execution_audit.json').write_text(json.dumps(dict(
        status='passed', benchmark_runs_checked=len(alert_audit), benchmarks=alert_audit,
        headline_native_structural_alert_exclusions=0,
        production_checks='Per-system independent_output_audit.json files under evidence/outputs/absolute_floor'),
        indent=2),encoding='utf-8')
    native=pd.read_csv(HEAD/'screening_full_1M_floor1000_native_scored_top5000.csv')
    native=native[native.native_weighted_coverage_pct.notna()].copy().sort_values('final_rank')
    native['native_cohort_rank']=np.arange(1,len(native)+1)
    native['screen_cohort_rank']=native.stage3_screen_rank.rank(method='first').astype(int)
    native['cohort_rank_shift']=native.screen_cohort_rank-native.native_cohort_rank
    top=native.head(1000)
    cov=native.native_weighted_coverage_pct
    upstream={}
    for column in ['screen_weighted_coverage_pct','pair_hash_overlap_pct','pair_hash_recall_pct',
                   'hotspot_compatible_matches','cascade_score_pct']:
        upstream[column]={'pearson_r':float(native[column].corr(cov)),
                          'spearman_rho':float(native[column].corr(cov,method='spearman'))}
    diagnostics={'n':len(native),'coverage_median':float(cov.median()),
        'coverage_q25':float(cov.quantile(.25)),'coverage_q75':float(cov.quantile(.75)),
        'coverage_min':float(cov.min()),'coverage_max':float(cov.max()),
        'rank_spearman':float(native.screen_cohort_rank.corr(native.native_cohort_rank,method='spearman')),
        'top1000_median_absolute_cohort_rank_shift':float(top.cohort_rank_shift.abs().median()),
        'top10_overlap':len(set(native.head(10).zinc_id)&set(native.nsmallest(10,'screen_cohort_rank').zinc_id)),
        'upstream_correlations':upstream}
    native.to_csv(META/'headline_native_cohort_diagnostics.csv',index=False)
    result={'head':head,'benchmarks':benchmarks,'replicates':replicates,'diagnostics':diagnostics,
        'benchmark_evaluation':'Equal status and score remain tied; EF, BEDROC and top-k recovery average over within-tie permutations. Identifier-ordered estimates are superseded.',
        'top10':native.head(10).to_dict('records'),
        'reference_verification':json.loads((META/'reference_verification/verification_summary.json').read_text()),
        'floor':json.loads((ROOT/'evidence/outputs/absolute_floor/analysis_summary.json').read_text())}
    (META/'manuscript_results.json').write_text(json.dumps(json_safe(result),indent=2,allow_nan=False))
    print('Collected complete, current results for all manuscript analyses')


if __name__=='__main__':main()
