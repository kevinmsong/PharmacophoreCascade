# Reproducibility package

This directory lets anyone recompute every retrospective enrichment number reported
in the manuscript *A staged topology-to-native pharmacophore cascade for scalable
native rescoring at peptide-receptor interfaces* **without** the production screening
engine, directly from the released per-molecule benchmark scores.

## What is released here

The conclusions rest on the retrospective benchmarks, and the full benchmark evidence
is public in this repository:

| Artifact | Location |
|---|---|
| Per-molecule scored tables (label, per-method rank/score/status, descriptors, ChEMBL/decoy provenance) | `evidence/data/machine_readable/{glp1r,ghsr,ntsr1,mdm2}_benchmark_scored.csv` |
| Benchmark configurations (gates, conformers, tolerances, seeds) | `evidence/configs/benchmark_*.yaml` |
| Pre-computed benchmark summaries with bootstrap CIs | `evidence/outputs/benchmark_*/benchmark_summary.csv` |
| Decoy-replicate summaries (5 sets/system) | `evidence/outputs/decoy_replicates/` |
| Receptor + native pharmacophores | `maps/` |
| Pinned environment | `requirements.txt` |

Each scored table has one row per benchmarked molecule (10 actives + 300 matched
decoys = 310 per system) with the `label` column and, for each method
(`full_cascade`, `native_only`, `stage3_only`, `standard_3d_pharmacophore`), a
`*_rank`, `*_score`, and `*_status` column. Molecules that failed an intermediate
stage carry a non-`scored` status and are ranked below scored molecules, exactly as
in the manuscript.

## Reproduce the reported metrics

```bash
pip install -r requirements.txt
python reproduce/reproduce_benchmarks.py            # all four systems
python reproduce/reproduce_benchmarks.py --system glp1r
```

This recomputes, per method and system, ROC-AUC, PR-AUC, EF1%, EF5%, BEDROC
(alpha = 20), and top-k active recovery. The output matches the manuscript, e.g.
for GLP-1R the full cascade gives ROC-AUC = 0.800, PR-AUC = 0.465, EF1% = 30,
BEDROC = 0.546, 5/10 actives in the top 10; the native-only baseline matches it
(ROC-AUC = 0.809); and across systems the native-only baseline equals or exceeds
the full cascade, confirming that the terminal native scoring supplies the
enrichment while the staged cascade supplies scalability.

### Metric definitions (match the Methods section)

- **ROC-AUC / PR-AUC** — global ranking quality (scikit-learn); PR-AUC is more
  sensitive to the low active prevalence.
- **EF k%** — fraction of actives recovered in the top `ceil(k% x N)` ranked
  molecules, divided by `k/100`.
- **BEDROC (alpha = 20)** — early-recognition-weighted score (Truchon & Bayly,
  *J. Chem. Inf. Model.* 2007, 47, 488-508).
- **top-k recovery** — number of actives ranked within the top k positions.

Grouped-bootstrap 95% confidence intervals and paired p-values (each active
resampled with its matched decoys) are provided pre-computed in the
`benchmark_*/benchmark_summary.csv` tables.

## What is not released

The production screening engine that executes the million-compound scan is under
active development for separate applications and is therefore not open-sourced; it
is available from the corresponding author under a reasonable-use agreement. It is
not required to reproduce any reported result: the Methods specify the exact
operation, parameters, and closed-form scoring functions of every stage, and the
released per-molecule scores above are sufficient to recompute all enrichment
metrics with the script in this directory.
