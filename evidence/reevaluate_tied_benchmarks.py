"""Replace identifier-ordered enrichment statistics with tie-aware statistics.

Uses freshly executed molecular scores; never rescoring molecules or changing
their screening/selection outcomes. Original evaluation summaries are archived.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from tie_aware_metrics import clustered_counts, groups_from_ranking

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / 'evidence/outputs/alert_disabled_revision'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def jobs():
    return ['benchmark_' + name for name in ['glp1r_full', 'ghsr_full', 'ntsr1_full', 'mdm2_full', 'glp1r_automated', 'glp1r_7ki0']] + [
        f'{system}_seed{seed}' for system in ['glp1r', 'ghsr', 'ntsr1', 'mdm2'] for seed in range(1, 6)]


def process(job):
    execution = json.loads((META / (job + '.json')).read_text())
    assert execution['status'] == 'completed'
    cfg = yaml.safe_load((META / 'configs' / (job + '.yaml')).read_text())['benchmark']
    folder = ROOT / cfg['output_dir']
    target = folder.parent
    marker = target / 'tie_aware_evaluation.json'
    library = ROOT / cfg['library_csv']
    rankings = sorted(folder.glob('*_ranking.csv'))
    signature = {str(path.relative_to(ROOT)): sha(path) for path in
                 [Path(__file__), Path(__file__).with_name('tie_aware_metrics.py'), library, *rankings]}
    if marker.exists():
        prior=json.loads(marker.read_text())
        if prior.get('source_sha256') == signature and prior.get('summary_sha256') == sha(target/'benchmark_summary.csv'):
            print('Already tie-aware: ' + job, flush=True)
            return
    backup = META / 'identifier_order_evaluation' / job
    backup.mkdir(parents=True, exist_ok=True)
    for name in ['benchmark_summary.csv', 'benchmark_deltas.csv']:
        path = target / name
        if path.exists() and not (backup / name).exists():
            shutil.copy2(path, backup / name)
    old = pd.read_csv(backup / 'benchmark_summary.csv').set_index('method')
    lib = pd.read_csv(library)
    count = int(cfg.get('bootstrap_iterations', 5000))
    seed = int(cfg.get('bootstrap_seed', 42))
    rows, audit, boots, comparisons = [], [], {}, []
    for path in rankings:
        method = path.name.removesuffix('_ranking.csv')
        frame = pd.read_csv(path)
        assert frame.ligand_id.is_unique and set(frame.ligand_id) == set(lib.ligand_id)
        point, boot = clustered_counts(frame, lib, draws=count, seed=seed)
        boots[method] = boot
        row = {'method': method, 'n_ranked': len(frame), **point}
        for metric, value in point.items():
            if boot:
                row[metric + '_ci_low'], row[metric + '_ci_high'] = np.quantile(boot[metric], [.025, .975])
            else:
                row[metric + '_ci_low'] = row[metric + '_ci_high'] = np.nan
            comparisons.append(dict(job=job, method=method, metric=metric,
                                    identifier_ordered=float(old.loc[method, metric]), tie_aware=value,
                                    correction='tie handling plus actual-cutoff EF normalization' if metric.startswith('ef_') else 'tie handling'))
        rows.append(row)
        ordered, group = groups_from_ranking(frame)
        ordered['evaluation_tie_group'] = group
        ordered['evaluation_rank'] = ordered.groupby('evaluation_tie_group')['rank'].transform('mean')
        audit.extend(ordered[['ligand_id', 'method', 'evaluation_tie_group', 'evaluation_rank']].to_dict('records'))
        print(f'{job} / {method}: tie-aware ROC-AUC {point["roc_auc"]:.6f}', flush=True)
    pd.DataFrame(rows).to_csv(target / 'benchmark_summary.csv', index=False)
    pd.DataFrame(audit).to_csv(target / 'benchmark_tie_groups.csv', index=False)
    pd.DataFrame(comparisons).to_csv(target / 'identifier_tie_sensitivity.csv', index=False)
    deltas = []
    for method, boot in boots.items():
        if method == 'full_cascade' or not boot:
            continue
        for metric, draw in boot.items():
            delta = draw - boots['full_cascade'][metric]
            p = min(1., (2 * min((delta >= 0).sum(), (delta <= 0).sum()) + 1) / (len(delta) + 1))
            deltas.append(dict(method=method, metric=metric, delta_mean=float(delta.mean()),
                delta_ci_low=float(np.quantile(delta, .025)), delta_ci_high=float(np.quantile(delta, .975)),
                n_bootstrap=len(delta), p_value=float(p)))
    delta = pd.DataFrame(deltas, columns=['method', 'metric', 'delta_mean', 'delta_ci_low', 'delta_ci_high', 'n_bootstrap', 'p_value'])
    delta['p_value_holm'] = np.nan
    for _, group in delta.groupby('metric'):
        order = group.sort_values('p_value').index
        adjusted = np.maximum.accumulate(delta.loc[order, 'p_value'].to_numpy() * np.arange(len(order), 0, -1))
        delta.loc[order, 'p_value_holm'] = np.minimum(adjusted, 1.)
    delta.to_csv(target / 'benchmark_deltas.csv', index=False)
    np.savez_compressed(target / 'tie_aware_bootstrap.npz', **{
        method + '__' + metric: values for method, boot in boots.items() for metric, values in boot.items()})
    marker.write_text(json.dumps(dict(status='completed', source_sha256=signature,
        summary_sha256=sha(target/'benchmark_summary.csv'), bootstrap_iterations=count,
        seed=seed, completed_utc=datetime.now(timezone.utc).isoformat(),
        rule='Equal status and score are tied; missing scores tie within failure status. ROC-AUC and average precision preserve ties. EF, BEDROC and top-k recovery average over all within-tie permutations.',
        enrichment_factor_normalization='Expected active prevalence within the actual ceil(fraction*N) selected positions divided by the full-library active prevalence; the old nominal-fraction denominator is superseded.',
        unchanged='Molecular scores, original total-order ranking files, candidate selections, and production floor retention.',
        original_summary_archive=str(backup.relative_to(ROOT)),
        scope='Primary enrichment metrics and paired bootstrap. Pairwise original-list correlations remain descriptive total-order diagnostics; scaffold enrichment based on identifier tie breaks is omitted.'), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job', choices=jobs())
    args = parser.parse_args()
    for job in [args.job] if args.job else jobs():
        process(job)


if __name__ == '__main__':
    main()
