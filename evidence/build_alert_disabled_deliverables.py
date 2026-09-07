"""Wait for fresh experiments, then rebuild their dependent revision artifacts."""
import argparse
import json
import subprocess
import sys
import time
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
META=ROOT/'evidence/outputs/alert_disabled_revision'


def dependencies():
    names=['benchmark_'+x for x in ['glp1r_full','ghsr_full','ntsr1_full','mdm2_full','glp1r_automated','glp1r_7ki0']]
    names+=['production_'+s for s in ['glp1r','ghsr','ntsr1','mdm2']]
    names += [f'{s}_seed{i}' for s in ['glp1r','ghsr','ntsr1','mdm2'] for i in range(1,6)]
    for name in names:
        path=META/(name+'.json')
        if not path.exists():return False
        state=json.loads(path.read_text())
        if state['status']=='failed':raise RuntimeError(f'Experiment failed: {name}')
        if state['status']!='completed':return False
    promotion=META/'headline_promotion.json'
    docking=ROOT/'evidence/outputs/docking_top10/rerun_manifest.json'
    if not promotion.exists() or not docking.exists():return False
    fresh=json.loads(promotion.read_text());dock=json.loads(docking.read_text())
    if fresh['new_top10']!=dock['ligand_ids']:raise RuntimeError('Docking top 10 differs from fresh headline')
    return True


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--wait',action='store_true');ap.add_argument('--resume',action='store_true')
    args=ap.parse_args()
    while not dependencies():
        if not args.wait:raise RuntimeError('Fresh experiments are incomplete')
        time.sleep(30)
    phases=[
        ('tie_aware_evaluation','evidence/reevaluate_tied_benchmarks.py'),
        ('validate_floor','evidence/validate_floor_outputs.py'),
        ('floor_analysis','evidence/analyze_absolute_floor.py'),
        ('attrition','evidence/analyze_alert_disabled_attrition.py'),
        ('replicate_aggregate','evidence/aggregate_replicates.py'),
        ('machine_readable','evidence/export_machine_readable.py'),
        ('equivalence','evidence/analyze_equivalence.py'),
        ('ranking_perturbations','evidence/analyze_ablation_common_support.py'),
        ('collect_results','evidence/collect_alert_disabled_results.py'),
        ('publication_tables','evidence/write_alert_disabled_tables.py'),
        ('manuscript','evidence/write_alert_disabled_manuscript.py'),
        ('supporting_information','evidence/write_alert_disabled_si.py'),
        ('reviewer_response','evidence/write_alert_disabled_response.py'),
        ('root_outputs','evidence/refresh_alert_disabled_outputs.py'),
        ('independent_metrics','reproduce/reproduce_benchmarks.py'),
        ('main_figures','evidence/make_ieee_figures.py','--target','acs'),
        ('headline_diagnostic_figures','evidence/refresh_headline_diagnostic_figures.py'),
        ('si_figures','evidence/make_ieee_appendix_figures.py','--target','acs'),
        ('graphical_abstract','evidence/make_floor_toc_graphic.py'),
        ('manuscript_numbers','evidence/check_manuscript_numbers.py','--dir','ACS_Omega_resubmission'),
        ('reviewer_coverage','evidence/check_reviewer_coverage.py','--dir','ACS_Omega_resubmission'),
        ('figure_accessibility','evidence/check_figure_accessibility.py','--dir','ACS_Omega_resubmission'),
    ]
    statepath=META/'deliverable_build.json'
    states=json.loads(statepath.read_text()) if statepath.exists() else {}
    if states and not args.resume:raise RuntimeError('Build state exists; use --resume after inspecting it')
    for name,*command in phases:
        if states.get(name,{}).get('status')=='completed':continue
        print(f'Building {name}',flush=True)
        started=time.perf_counter()
        with open(META/f'build_{name}.log','w',encoding='utf-8') as log:
            result=subprocess.run([sys.executable,'-u',*command],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        states[name]={'status':'completed' if result.returncode==0 else 'failed',
            'returncode':result.returncode,'command':command,'elapsed_seconds':time.perf_counter()-started,
            'finished_utc':datetime.now(timezone.utc).isoformat()}
        statepath.write_text(json.dumps(states,indent=2))
        if result.returncode:raise RuntimeError(f'Build failed: {name}; inspect build_{name}.log')
    print('Current data, manuscript sources, figures, and numerical checks ready for document rendering and visual QA',flush=True)


if __name__=='__main__':main()
