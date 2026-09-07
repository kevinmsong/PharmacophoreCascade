"""Replace legacy analysis deliverables with current, provenance-checked outputs."""
import json
import argparse
import pickle
import shutil
import subprocess
from pathlib import Path
import pandas as pd
import yaml
from src.loader import load_tables
from src.benchmark import BenchmarkResult
from src.ablation import AblationResult,sweep_cascade_weights
from src.claims import run_all_claims
from src.report import emit_report

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/outputs'
META=OUT/'alert_disabled_revision'
HEAD=ROOT/'results/absolute_floor_1000'


def quote(path):
    return "'"+str(path).replace("'","''")+"'"


def archive(path):
    if path.is_symlink() or path.is_junction():
        raise RuntimeError(f'Refusing to archive through an existing results alias: {path}')
    path=path.resolve()
    assert path.is_relative_to(ROOT.resolve()) and path!=ROOT.resolve()
    if not path.exists():return
    base=(ROOT/'tmp/alert_disabled_revision/postprocess_prior').resolve()
    target=(base/path.relative_to(ROOT)).resolve()
    assert target.is_relative_to(base)
    if target.exists():raise RuntimeError(f'Postprocess archive exists; inspect before retry: {path}')
    target.parent.mkdir(parents=True,exist_ok=True)
    path.rename(target)


def main():
    data=json.loads((META/'manuscript_results.json').read_text())
    # Root-level reports and caches are prior results; archive them before
    # overwriting their established deliverable paths.
    for old in list(OUT.iterdir()):
        if old.is_file() or old.name=='cache':archive(old)
    old_bundle=ROOT/'results/screening_full_1M_topological_hashed_native_terminal_bundle'
    new_bundle=HEAD/'screening_full_1M_floor1000_native_terminal_bundle'
    assert new_bundle.is_dir()
    archive(old_bundle)
    subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',
        'New-Item -ItemType Junction -Path '+quote(old_bundle)+' -Target '+quote(new_bundle)+' | Out-Null'],check=True)
    cfg=yaml.safe_load((ROOT/'evidence/configs/benchmark.yaml').read_text())
    tables=load_tables(ROOT/'evidence/configs/benchmark.yaml',use_cache=False)
    cohort=pd.read_csv(META/'headline_native_cohort_diagnostics.csv')
    tables['best_native']=cohort
    rc=pd.DataFrame({'ligand_id':cohort.zinc_id,'weighted_coverage_pct':cohort.screen_weighted_coverage_pct,
        'native_weighted_coverage_pct':cohort.native_weighted_coverage_pct,'screen_rank':cohort.screen_cohort_rank,
        'native_rank':cohort.native_cohort_rank,'rank_shift':cohort.cohort_rank_shift,
        'abs_rank_shift':cohort.cohort_rank_shift.abs()})
    tables['rank_comparison']=rc
    rc.to_csv(META/'current_native_rank_comparison.csv',index=False)
    benchdir=OUT/'benchmark_glp1r_full'
    rankings={m:pd.read_csv(benchdir/'benchmark_external'/f'{m}_ranking.csv')
              for m in data['benchmarks']['glp1r_full']}
    library=pd.read_csv(ROOT/'evidence/data/glp1r_external_benchmark_library.csv')
    bench=BenchmarkResult(rankings=rankings,pairwise_df=pd.read_csv(benchdir/'benchmark_pairwise.csv'),
        summary_df=pd.read_csv(benchdir/'benchmark_summary.csv'),delta_df=pd.read_csv(benchdir/'benchmark_deltas.csv'),
        scaffold_series=library.set_index('ligand_id').murcko_scaffold,library_df=library)
    ab=pd.read_csv(OUT/'ablation_common/ablation_common_support.csv').rename(
        columns={'ablation':'method','tau_1000':'kendall_tau','rho_shared':'spearman_r',
                 'J10':'jaccard_top10','J50':'jaccard_top50','J100':'jaccard_top100','J500':'jaccard_top500'})
    acfg=yaml.safe_load((ROOT/'evidence/configs/ablation.yaml').read_text())
    sensitivity=sweep_cascade_weights(tables,acfg['cascade_weight_sweep'],[10,50,100,500])
    sensitivity['sweep_type']='cascade_weights'
    sensitivity['analysis_scope']='Reordering observed successful Stage-3 table; no upstream rescoring or new native pool'
    ablation=AblationResult(ablation_df=ab,sensitivity_df=sensitivity)
    claims=run_all_claims(tables,cfg)
    emit_report(bench,ablation,claims,tables,OUT)
    from PIL import Image
    for stem in ['benchmark_plots','ablation_sensitivity']:
        assert (OUT/(stem+'.pdf')).is_file(),f'Missing current report figure: {stem}'
        with Image.open(OUT/(stem+'.png')) as raster:
            assert min(raster.info.get('dpi',(0,0)))>=599,stem
    cache=OUT/'cache';cache.mkdir(exist_ok=True)
    for name,value in [('_bench.pkl',(bench,tables)),('_ablation.pkl',ablation),('_claims.pkl',claims)]:
        with open(cache/name,'wb') as f:pickle.dump(value,f)
    for key,table in tables.items():table.to_parquet(cache/(key+'.parquet'),index=False)
    finish_outputs(data,cohort)


def finish_outputs(data,cohort):
    report=OUT/'benchmark_report.md'
    report.write_text('# Current alert-disabled revision\n\nStructural alerts are recorded without exclusion; property limits remain. '
        'The headline shortlist has a capped minimum of 1,000. Retrospective statistics below use the freshly '
        'executed GLP-1R benchmark; internal ranking and motif diagnostics use the new headline cohort. '
        'Ranking perturbations reorder observed score tables and do not imply gate-removal reruns.\n\n'+report.read_text(encoding='utf-8'),encoding='utf-8')
    aliases=[]
    for system in ['glp1r','ghsr','ntsr1','mdm2']:
        current=OUT/'absolute_floor'/system
        target=OUT/f'production_{system}'
        archive(target);target.mkdir()
        for src,dst in [(current/'stage012_evaluation.csv','stage012_evaluation.csv'),
            (current/'floor1000/native_scored.csv','prod_native_scored.csv'),
            (current/'floor1000/final_ranked.csv','prod_native_final.csv'),
            (current/'floor1000/per_molecule_survival.csv','per_molecule_survival.csv')]:shutil.copy2(src,target/dst)
        stages=pd.read_csv(current/'stage_survival.csv');stages=stages[stages.policy=='floor1000'].copy()
        bystage=stages.set_index('stage');actives=int(bystage.loc['input','actives'])
        stages['active_survival_pct']=100*stages.actives/actives
        stages.to_csv(target/'stage_survival.csv',index=False)
        summary={'background_size':30000,'shortlist_fraction':.05,'shortlist_floor':1000,
            'chemistry_gate_mode':'warn_only','native_pains_filter':False,
            'n_candidates':int(bystage.loc['stage012_pass','total']),
            'n_shortlist':int(bystage.loc['shortlist','total']),
            'actives_in_final':int(bystage.loc['final_ranked','actives']),
            'stage_survival':stages.to_dict('records'),'current_source':str(current.relative_to(ROOT))}
        (target/'production_summary.json').write_text(json.dumps(summary,indent=2))
        aliases.append({'path':str(target.relative_to(ROOT)),'source':str(current.relative_to(ROOT)),'policy':'floor1000'})
        # Older small-benchmark folders refer explicitly to the new expanded
        # full benchmark rather than retaining obsolete numerical results.
        legacy=OUT/f'benchmark_{system}'
        if legacy.exists():
            archive(legacy);legacy.mkdir()
            (legacy/'CURRENT_RESULTS.md').write_text(f'Current results: `../benchmark_{system}_full/`.\n\n'
                'The current expanded benchmark and all original 10-active decoy replicates were rerun with '
                'structural-alert exclusions disabled. See `../decoy_replicates/` for the replicate universes.\n')
    efficiency=OUT/'efficiency'
    if efficiency.exists():archive(efficiency)
    efficiency.mkdir()
    pairs=pd.read_csv(OUT/'absolute_floor/policy_comparison.csv')
    pairs.to_csv(efficiency/'executed_policy_comparison.csv',index=False)
    (efficiency/'README.md').write_text('Current efficiency/retention results are actual paired executions in '
        '`../absolute_floor/`. Inferred native counts from the earlier shortlist sweep are superseded. '
        'The percentage denominator is the Stage-1/2 candidate pool in these production benchmarks.\n')
    for legacy_name in ['benchmark_external','_smoke_glp1r_full']:
        legacy=OUT/legacy_name
        if legacy.exists():archive(legacy)
    # Keep the original root benchmark path as a current-results alias.
    subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',
        'New-Item -ItemType Junction -Path '+quote(OUT/'benchmark_external')+' -Target '+
        quote(OUT/'benchmark_glp1r_full/benchmark_external')+' | Out-Null'],check=True)
    revision=OUT/'revision'
    for name in ['ablation_curated_vs_automated.csv','docking_orthogonal_summary.csv','rank_reconciliation.json',
                 'revision_numbers.json','tableS3_down_reconciled.csv','tableS3_up_reconciled.csv']:
        archive(revision/name)
    pd.DataFrame([{'metric':key,'curated':data['benchmarks']['glp1r_full']['full_cascade'][key],
        'automated_no_curation':data['benchmarks']['glp1r_automated']['full_cascade'][key]}
        for key in ['roc_auc','pr_auc','ef_1pct','bedroc','top10_recovery']]).to_csv(
            revision/'ablation_curated_vs_automated.csv',index=False)
    docking=pd.read_csv(revision/'docking_top10_summary.csv')
    docking['overall_best']=docking[['best_active','best_inactive']].min(axis=1)
    docking.to_csv(revision/'docking_orthogonal_summary.csv',index=False)
    cohort.nlargest(5,'cohort_rank_shift').to_csv(revision/'tableS3_up_reconciled.csv',index=False)
    cohort.nsmallest(5,'cohort_rank_shift').to_csv(revision/'tableS3_down_reconciled.csv',index=False)
    (revision/'rank_reconciliation.json').write_text(json.dumps(data['diagnostics'],indent=2))
    (revision/'revision_numbers.json').write_text(json.dumps(data,indent=2,allow_nan=False))
    # The native-only comparator bypasses Stage 3. Its inherited correlation
    # plotting shim used placeholder Stage-3 values; remove those diagnostics
    # from current deliverables while preserving genuine native scores.
    for name in data['benchmarks']:
        bundle=OUT/f'benchmark_{name}/benchmark_external/native_only_bundle'
        if not bundle.exists():continue
        reports=bundle/'reports'
        if reports.exists():archive(reports)
        reports.mkdir()
        (reports/'README.md').write_text('Native-only scoring bypasses Stages 0–3. Stage-3/native correlation '
            'plots and rank-shift diagnostics are not applicable. Only actual native scores and method-ranking '
            'files should be analyzed. The benchmark driver used placeholder Stage-3 columns for an inherited '
            'output interface; those columns are not measured screen scores.\n',encoding='utf-8')
        manifest_path=bundle/'manifests/analysis_manifest.json'
        manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
        archive(manifest_path)
        manifest['diagnostics']={'status':'not_applicable','always_on':False,'diagnostics':[],
            'reason':'Native-only bypasses Stage 3; correlation diagnostics based on adapter placeholders are excluded.',
            'outputs':{}}
        manifest['input_adapter_note']='Stage-0–3 fields inherited from the native-only adapter are execution placeholders, not measured screen scores. Genuine native scores and fitted geometry are retained.'
        manifest_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
        for filename in ['native_correlation_metrics.csv','native_correlation_quartiles.csv','native_rank_comparison.csv']:
            path=bundle/'analysis'/filename
            if path.exists():archive(path)
    (META/'output_refresh.json').write_text(json.dumps({'status':'completed','production_aliases':aliases,
        'root_report_source':'fresh GLP-1R full benchmark plus fresh headline diagnostics',
        'native_only_stage3_diagnostics':'not applicable; excluded from current outputs'},indent=2))
    print('Refreshed root reports, caches, production aliases, and diagnostic scope')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume-after-root-report',action='store_true',
        help='Resume the inspected interruption after fresh reports/caches were written, before report prefix and aliases.')
    args=parser.parse_args()
    if args.resume_after_root_report:
        alias=ROOT/'results/screening_full_1M_topological_hashed_native_terminal_bundle'
        assert alias.is_junction() and alias.resolve()==(HEAD/'screening_full_1M_floor1000_native_terminal_bundle').resolve()
        assert all((OUT/'cache'/name).exists() for name in ['_bench.pkl','_ablation.pkl','_claims.pkl'])
        assert not (OUT/'benchmark_report.md').read_text(encoding='utf-8').startswith('# Current alert-disabled revision')
        data=json.loads((META/'manuscript_results.json').read_text(encoding='utf-8'))
        cohort=pd.read_csv(META/'headline_native_cohort_diagnostics.csv')
        finish_outputs(data,cohort)
    else:
        main()
