# peptide-contactR | 5% + floor | 4,997 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.049` (95% bootstrap CI `0.023` to `0.074`)
- Spearman `rho = 0.055` (95% bootstrap CI `0.027` to `0.081`)
- Kendall `tau = 0.039` (95% bootstrap CI `0.019` to `0.057`)
- Linear slope: `0.14` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows little evidence of a positive association between screening weighted coverage and native peptide-contact 3D mimicry.
- Highest native 3D coverage ligand: `ZINCk700000Gfz1H` (30.0% native coverage at 12.787% screen weighted coverage)
- Highest screen weighted coverage ligand: `ZINCk700000BmVa0` (16.920% screen weighted coverage but 24.8% native coverage)
- Median absolute rank shift between the two scoring systems: `1457.0` positions
- Top-10 overlap between screen ranking and native ranking: `0` ligands

## Quartile Summary

- `Q1 12.046-12.539%`: `n=1250`; native median `24.6%`; native range `17.8%` to `30.0%`
- `Q2 12.539-12.889%`: `n=1249`; native median `24.2%`; native range `16.9%` to `30.0%`
- `Q3 12.890-13.387%`: `n=1249`; native median `24.2%`; native range `15.1%` to `29.1%`
- `Q4 13.388-16.920%`: `n=1249`; native median `24.6%`; native range `16.4%` to `29.1%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `ZINCj700000JVICU` (screen rank `4955` to native rank `6`)
- Strongest downward reranking under native 3D scoring: `ZINCk8000004kfyb` (screen rank `98` to native rank `4978`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `1.89` native-coverage percentage points
- Linear-fit MAE: `1.52` native-coverage percentage points
- Residual standard deviation: `1.89`
- Spearman correlation between fitted values and absolute residuals: `-0.036`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
