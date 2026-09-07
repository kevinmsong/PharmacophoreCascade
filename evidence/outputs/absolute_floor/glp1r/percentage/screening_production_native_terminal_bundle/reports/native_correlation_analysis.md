# peptide-contactR | 5% only | 647 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.225` (95% bootstrap CI `0.134` to `0.310`)
- Spearman `rho = 0.163` (95% bootstrap CI `0.087` to `0.238`)
- Kendall `tau = 0.125` (95% bootstrap CI `0.067` to `0.182`)
- Linear slope: `0.32` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows a weak positive association between screening weighted coverage and native peptide-contact 3D mimicry.
- Highest native 3D coverage ligand: `CHEMBL6027640` (27.9% native coverage at 15.295% screen weighted coverage)
- Highest screen weighted coverage ligand: `CHEMBL6027640` (15.295% screen weighted coverage but 27.9% native coverage)
- Median absolute rank shift between the two scoring systems: `189.0` positions
- Top-10 overlap between screen ranking and native ranking: `1` ligands

## Quartile Summary

- `Q1 6.545-10.250%`: `n=162`; native median `22.3%`; native range `13.4%` to `22.5%`
- `Q2 10.253-10.986%`: `n=162`; native median `22.3%`; native range `16.4%` to `27.2%`
- `Q3 10.989-11.579%`: `n=161`; native median `22.3%`; native range `16.4%` to `27.2%`
- `Q4 11.586-15.295%`: `n=162`; native median `22.3%`; native range `16.4%` to `27.9%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `ZINCf50000009Fc3` (screen rank `632` to native rank `80`)
- Strongest downward reranking under native 3D scoring: `ZINCf50000009CVi` (screen rank `11` to native rank `623`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `1.80` native-coverage percentage points
- Linear-fit MAE: `1.25` native-coverage percentage points
- Residual standard deviation: `1.80`
- Spearman correlation between fitted values and absolute residuals: `-0.529`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
