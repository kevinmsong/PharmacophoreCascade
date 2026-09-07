# GHSR | 5% + floor | 271 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.505` (95% bootstrap CI `0.400` to `0.597`)
- Spearman `rho = 0.466` (95% bootstrap CI `0.357` to `0.565`)
- Kendall `tau = 0.350` (95% bootstrap CI `0.266` to `0.430`)
- Linear slope: `2.73` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows a moderate positive association between screening weighted coverage and native peptide-contact 3D mimicry.
- Highest native 3D coverage ligand: `CHEMBL130229` (92.6% native coverage at 14.536% screen weighted coverage)
- Highest screen weighted coverage ligand: `CHEMBL4164040` (17.728% screen weighted coverage but 81.8% native coverage)
- Median absolute rank shift between the two scoring systems: `50.0` positions
- Top-10 overlap between screen ranking and native ranking: `2` ligands

## Quartile Summary

- `Q1 6.396-11.270%`: `n=68`; native median `63.7%`; native range `27.0%` to `77.4%`
- `Q2 11.278-12.234%`: `n=68`; native median `64.7%`; native range `50.4%` to `79.8%`
- `Q3 12.240-13.360%`: `n=67`; native median `64.2%`; native range `51.5%` to `79.3%`
- `Q4 13.364-17.728%`: `n=68`; native median `77.4%`; native range `51.5%` to `92.6%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `ZINCf5000000a5ui` (screen rank `256` to native rank `32`)
- Strongest downward reranking under native 3D scoring: `ZINCj600000JU75f` (screen rank `12` to native rank `220`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `8.00` native-coverage percentage points
- Linear-fit MAE: `6.10` native-coverage percentage points
- Residual standard deviation: `8.01`
- Spearman correlation between fitted values and absolute residuals: `0.065`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
