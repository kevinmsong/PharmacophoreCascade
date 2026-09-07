"""Numerically compare independently recomputed metrics with fresh benchmark summaries."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
META=ROOT/'evidence/outputs/alert_disabled_revision'


def main():
    spec=importlib.util.spec_from_file_location('independent_metrics',ROOT/'reproduce/reproduce_benchmarks.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[]
    mapping={'ROC-AUC':'roc_auc','PR-AUC':'pr_auc','EF1%':'ef_1pct','EF5%':'ef_5pct',
             'BEDROC':'bedroc','top10_rec':'top10_recovery','top25_rec':'top25_recovery'}
    for system in module.SYSTEMS:
        state=json.loads((META/f'benchmark_{system}_full.json').read_text())
        assert state['status']=='completed'
        raw=pd.read_csv(module.DATA/f'{system}_benchmark_scored.csv')
        actives=int((raw.label=='active').sum())
        for method in module.METHODS:
            ranks=raw[f'{method}_rank'].to_numpy()
            np.testing.assert_array_equal(np.sort(ranks),np.arange(1,len(raw)+1))
        independent=module.run_system(system)
        emitted=pd.read_csv(ROOT/f'evidence/outputs/benchmark_{system}_full/benchmark_summary.csv').set_index('method')
        for method in module.METHODS:
            for alternate,key in mapping.items():
                expected=float(emitted.loc[method,key])*(actives if alternate.startswith('top') else 1)
                actual=float(independent.loc[method,alternate])
                error=abs(actual-expected)
                assert error<1e-12,f'{system}/{method}/{key}: {actual} != {expected}'
                rows.append({'system':system,'method':method,'metric':key,'independent':actual,
                    'emitted':expected,'absolute_difference':error})
    pd.DataFrame(rows).to_csv(META/'independent_metric_comparison.csv',index=False)
    (META/'independent_metric_validation.json').write_text(json.dumps({'status':'passed','comparisons':len(rows),
        'maximum_absolute_error':max(row['absolute_difference'] for row in rows),'tolerance':1e-12},indent=2))
    print(f'All {len(rows)} independently recomputed metric comparisons agree within 1e-12')


if __name__=='__main__':main()
