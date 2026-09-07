"""Identifier-invariant benchmark metrics, including expected retrieval within ties.

Score groups retain the declared method-status ordering. Molecules with equal
status and equal numerical score (including missing-score groups) are tied.
No label or ligand identifier resolves a tie.
"""
import numpy as np
import pandas as pd


def groups_from_ranking(frame):
    frame = frame.sort_values('rank').copy()
    keys = [(str(s), None if pd.isna(v) else float(v)) for s, v in zip(frame.status, frame.score)]
    mapping = {}
    group = []
    for key in keys:
        if key not in mapping:
            mapping[key] = len(mapping)
        group.append(mapping[key])
    group = np.asarray(group)
    assert np.all(np.diff(group) >= 0), 'Tie groups must be contiguous in declared status/score order'
    return frame, group


def metrics_from_counts(totals, actives, alpha=20.0):
    """Evaluate one or many count arrays, with best score group first."""
    n = np.atleast_2d(np.asarray(totals, dtype=float))
    a = np.atleast_2d(np.asarray(actives, dtype=float))
    assert n.shape == a.shape and np.all(a <= n) and np.all(a >= 0)
    cn, ca = n.cumsum(axis=1), a.cumsum(axis=1)
    count, positives = cn[:, -1], ca[:, -1]
    negatives = count - positives
    assert np.all(positives > 0) and np.all(negatives > 0)
    previous = cn - n
    probability = np.divide(a, n, out=np.zeros_like(a), where=n > 0)
    precision = np.divide(ca, cn, out=np.zeros_like(a), where=cn > 0)
    decoys = n - a
    roc = (a * (negatives[:, None] - decoys.cumsum(axis=1) + .5 * decoys)).sum(axis=1) / (positives * negatives)
    out = {'roc_auc': roc, 'pr_auc': (a * precision).sum(axis=1) / positives}
    def recovered(k):
        admitted = np.minimum(np.maximum(np.asarray(k)[:, None] - previous, 0), n)
        return (admitted * probability).sum(axis=1) / positives
    for name, frac in [('ef_1pct', .01), ('ef_5pct', .05), ('ef_01pct', .001)]:
        cutoff = np.minimum(count, np.maximum(1, np.ceil(count * frac)))
        out[name] = recovered(cutoff) * count / cutoff
    for k in [10, 25, 50]:
        out[f'top{k}_recovery'] = recovered(np.minimum(count, k))
    step = alpha / count
    weight_sum = np.exp(-step[:, None] * (previous + 1)) * (-np.expm1(-step[:, None] * n)) / (-np.expm1(-step[:, None]))
    exponential_sum = (weight_sum * probability).sum(axis=1)
    ra = positives / count
    rie = exponential_sum / (positives * (-np.expm1(-alpha)) / (count * np.expm1(step)))
    high = (-np.expm1(-alpha * ra)) / (ra * (-np.expm1(-alpha)))
    low = (-np.expm1(alpha * ra)) / (ra * (-np.expm1(alpha)))
    out['bedroc'] = np.clip((rie - low) / (high - low), 0, 1)
    return out


def evaluate(frame):
    ordered, group = groups_from_ranking(frame)
    totals = np.bincount(group)
    actives = np.bincount(group, weights=ordered.label.eq('active').to_numpy())
    return {k: float(v[0]) for k, v in metrics_from_counts(totals, actives).items()}


def clustered_counts(frame, library, draws=5000, seed=42):
    """Point metrics and paired cluster bootstrap, without expanding sampled rows."""
    ordered, group = groups_from_ranking(frame)
    lib = library.set_index('ligand_id').loc[ordered.ligand_id]
    active_ids = library.loc[library.label.eq('active'), 'ligand_id'].tolist()
    cluster_map = {name: i for i, name in enumerate(active_ids)}
    owners = [lid if label == 'active' else owner for lid, label, owner in
              zip(lib.index, lib.label, lib.matched_active_ligand_id)]
    assert all(owner in cluster_map for owner in owners), 'Each decoy needs its matched active'
    cluster = np.asarray([cluster_map[owner] for owner in owners])
    total = np.zeros((len(active_ids), int(group.max()) + 1))
    active = np.zeros_like(total)
    np.add.at(total, (cluster, group), 1)
    np.add.at(active, (cluster, group), lib.label.eq('active').astype(int).to_numpy())
    point = {k: float(v[0]) for k, v in metrics_from_counts(total.sum(axis=0), active.sum(axis=0)).items()}
    if draws == 0:
        return point, {}
    chosen = np.random.default_rng(seed).choice(len(active_ids), size=(draws, len(active_ids)), replace=True)
    multiplicity = np.zeros((draws, len(active_ids)))
    np.add.at(multiplicity, (np.repeat(np.arange(draws), len(active_ids)), chosen.ravel()), 1)
    boot = metrics_from_counts(multiplicity @ total, multiplicity @ active)
    return point, boot
