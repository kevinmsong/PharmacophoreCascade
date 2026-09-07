"""Validate tie handling against sklearn and exhaustive within-tie permutations."""
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from tie_aware_metrics import evaluate, metrics_from_counts, clustered_counts


def case(blocks):
    records = []
    for group, labels in enumerate(blocks):
        for label in labels:
            records.append(dict(ligand_id=f'ID{len(records)}', rank=len(records)+1,
                status='scored', score=-group, label='active' if label else 'decoy'))
    frame = pd.DataFrame(records)
    actual = evaluate(frame)
    y = frame.label.eq('active').to_numpy()
    np.testing.assert_allclose(actual['roc_auc'], roc_auc_score(y, frame.score), atol=1e-14)
    np.testing.assert_allclose(actual['pr_auc'], average_precision_score(y, frame.score), atol=1e-14)
    samples = []
    for arrangements in itertools.product(*[sorted(set(itertools.permutations(block))) for block in blocks]):
        labels = np.asarray([x for block in arrangements for x in block])
        n = len(labels); positives = labels.sum(); alpha = 20
        exp_sum = np.exp(-alpha * (np.flatnonzero(labels) + 1) / n).sum()
        ra = positives/n
        rie = exp_sum / (positives * (1-np.exp(-alpha))/(n*np.expm1(alpha/n)))
        upper = (1-np.exp(-alpha*ra))/(ra*(1-np.exp(-alpha)))
        lower = (1-np.exp(alpha*ra))/(ra*(1-np.exp(alpha)))
        row = {'bedroc': (rie-lower)/(upper-lower), 'roc_auc': roc_auc_score(labels, -np.arange(n))}
        for name, frac in [('ef_1pct', .01), ('ef_5pct', .05), ('ef_01pct', .001)]:
            cutoff=max(1, int(np.ceil(n*frac)))
            row[name] = labels[:cutoff].sum()/positives*n/cutoff
        for k in [10,25,50]: row[f'top{k}_recovery'] = labels[:k].sum()/positives
        samples.append(row)
    expected = pd.DataFrame(samples).mean()
    for metric, value in expected.items(): np.testing.assert_allclose(actual[metric], value, atol=1e-13)
    renamed = frame.copy()
    renamed['ligand_id'] = ['CHEMBL' + str(i) if label else 'DECOY' + str(i) for i,label in enumerate(y)]
    renamed = renamed.sample(frac=1, random_state=9)
    assert evaluate(renamed) == actual
    return len(samples)


def main():
    counts = [case(blocks) for blocks in [
        [[1,0],[1,1,0],[0]], [[1,0,1,0,1,0]],
        [[1],[0],[0],[0],[0],[1],[0],[0],[0],[1,0,1]],
        [[1],[0],[1],[0]],
    ]]
    # Integer multiplicities must equal explicitly repeated observations.
    frame = pd.DataFrame(dict(ligand_id=['a','d1','b','d2'], rank=[1,2,3,4],
        status=['scored']*4, score=[1,1,0,0], label=['active','decoy','active','decoy']))
    lib = frame.assign(matched_active_ligand_id=['a','a','b','b'])
    _, boot = clustered_counts(frame, lib, draws=20, seed=42)
    selected = np.random.default_rng(42).choice(2,size=(20,2),replace=True)
    for i, draw in enumerate(selected):
        expanded = pd.concat([frame.iloc[2*x:2*x+2] for x in draw]).sort_values('rank').copy()
        expanded['rank'] = np.arange(1,len(expanded)+1)
        expected = evaluate(expanded)
        for metric in expected: np.testing.assert_allclose(boot[metric][i],expected[metric],atol=1e-13)
    out = Path(__file__).resolve().parents[1]/'evidence/outputs/alert_disabled_revision/tie_metric_validation.json'
    out.write_text(json.dumps(dict(status='passed', exhaustive_orderings=counts,
        bootstrap_expansion_checks=20, sklearn_roc_and_average_precision=True,
        identifier_and_row_order_invariance=True),indent=2))
    print('Tie-aware metrics pass sklearn, exhaustive permutations, identifier invariance and explicit bootstrap expansion checks.')


if __name__ == '__main__': main()
