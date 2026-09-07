# NTSR1 | 5% only | 20 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.692` (95% bootstrap CI `0.407` to `0.891`)
- Spearman `rho = 0.698` (95% bootstrap CI `0.381` to `0.881`)
- Kendall `tau = 0.585` (95% bootstrap CI `0.322` to `0.775`)
- Linear slope: `4.46` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows a strong positive association between screening weighted coverage and native peptide-contact 3D mimicry.
- Highest native 3D coverage ligand: `CHEMBL3315213` (92.9% native coverage at 17.209% screen weighted coverage)
- Highest screen weighted coverage ligand: `CHEMBL3315213` (17.209% screen weighted coverage but 92.9% native coverage)
- Median absolute rank shift between the two scoring systems: `2.5` positions
- Top-10 overlap between screen ranking and native ranking: `7` ligands

## Quartile Summary

- `Q1 10.978-12.531%`: `n=5`; native median `58.3%`; native range `47.2%` to `82.7%`
- `Q2 12.814-14.353%`: `n=5`; native median `82.7%`; native range `58.3%` to `82.7%`
- `Q3 14.411-15.927%`: `n=5`; native median `82.7%`; native range `72.3%` to `82.7%`
- `Q4 15.994-17.209%`: `n=5`; native median `82.7%`; native range `82.7%` to `92.9%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `CHEMBL1312502` (screen rank `18` to native rank `2`)
- Strongest downward reranking under native 3D scoring: `CHEMBL463811` (screen rank `5` to native rank `14`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `8.63` native-coverage percentage points
- Linear-fit MAE: `6.49` native-coverage percentage points
- Residual standard deviation: `8.86`
- Spearman correlation between fitted values and absolute residuals: `-0.644`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
