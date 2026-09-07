# MDM2-p53 | 5% only | 1,152 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.188` (95% bootstrap CI `0.127` to `0.248`)
- Spearman `rho = 0.148` (95% bootstrap CI `0.089` to `0.205`)
- Kendall `tau = 0.102` (95% bootstrap CI `0.061` to `0.144`)
- Linear slope: `1.03` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows a weak positive association between screening weighted coverage and native peptide-contact 3D mimicry.
- Highest native 3D coverage ligand: `CHEMBL6006846` (87.8% native coverage at 15.874% screen weighted coverage)
- Highest screen weighted coverage ligand: `ZINCk700000gshd2` (18.411% screen weighted coverage but 56.0% native coverage)
- Median absolute rank shift between the two scoring systems: `307.0` positions
- Top-10 overlap between screen ranking and native ranking: `0` ligands

## Quartile Summary

- `Q1 8.956-13.993%`: `n=288`; native median `47.7%`; native range `29.0%` to `72.0%`
- `Q2 13.998-14.953%`: `n=288`; native median `48.8%`; native range `33.4%` to `74.7%`
- `Q3 14.954-15.635%`: `n=288`; native median `49.2%`; native range `34.6%` to `74.7%`
- `Q4 15.649-18.411%`: `n=288`; native median `56.0%`; native range `29.0%` to `87.8%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `CHEMBL3217780` (screen rank `1103` to native rank `6`)
- Strongest downward reranking under native 3D scoring: `ZINCj5000002eWzt` (screen rank `63` to native rank `1150`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `7.59` native-coverage percentage points
- Linear-fit MAE: `6.42` native-coverage percentage points
- Residual standard deviation: `7.59`
- Spearman correlation between fitted values and absolute residuals: `0.056`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
