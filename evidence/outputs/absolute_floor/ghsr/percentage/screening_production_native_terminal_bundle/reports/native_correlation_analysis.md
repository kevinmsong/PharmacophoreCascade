# GHSR | 5% only | 14 native-scored ligands Native Correlation Analysis

## Summary

- Pearson `r = 0.449` (95% bootstrap CI `0.048` to `0.738`)
- Spearman `rho = 0.284` (95% bootstrap CI `-0.303` to `0.733`)
- Kendall `tau = 0.211` (95% bootstrap CI `-0.208` to `0.602`)
- Linear slope: `2.18` native-coverage percentage points per 1 percentage point of screen weighted coverage
- Interpretation: The ligand-level comparison shows a moderate positive association between screening weighted coverage and native peptide-contact 3D mimicry. At least one bootstrap confidence interval includes zero.
- Highest native 3D coverage ligand: `CHEMBL130229` (92.6% native coverage at 14.536% screen weighted coverage)
- Highest screen weighted coverage ligand: `CHEMBL4164040` (17.728% screen weighted coverage but 81.8% native coverage)
- Median absolute rank shift between the two scoring systems: `3.0` positions
- Top-10 overlap between screen ranking and native ranking: `8` ligands

## Quartile Summary

- `Q1 10.548-12.486%`: `n=4`; native median `70.8%`; native range `64.2%` to `79.8%`
- `Q2 12.916-14.536%`: `n=3`; native median `82.4%`; native range `51.5%` to `92.6%`
- `Q3 14.778-15.119%`: `n=3`; native median `79.3%`; native range `79.3%` to `79.3%`
- `Q4 15.190-17.728%`: `n=4`; native median `79.3%`; native range `63.1%` to `81.8%`

## Rank Reordering

- Strongest upward reranking under native 3D scoring: `CHEMBL4171914` (screen rank `12` to native rank `4`)
- Strongest downward reranking under native 3D scoring: `ZINCj600000JU75f` (screen rank `4` to native rank `13`)
- Interpretation: the native 3D score meaningfully reshuffles the screen-derived ordering, so screen `weighted_coverage_pct` and peptide-interface 3D mimicry should be treated as related but non-interchangeable prioritization signals.

## Residual Diagnostics

- Linear-fit RMSE: `9.13` native-coverage percentage points
- Linear-fit MAE: `6.85` native-coverage percentage points
- Residual standard deviation: `9.47`
- Spearman correlation between fitted values and absolute residuals: `-0.569`
- Residual diagnostics are used here as a simple check on model adequacy; they do not upgrade the analysis from descriptive association to a causal or mechanistic model.

## Figures

![Correlation scatter](native_correlation_scatter.png)

![Quartile summary](native_correlation_quartiles.png)

![Rank reranking](native_correlation_rank_reranking.png)

![Residual diagnostics](native_correlation_residual_diagnostics.png)

See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.

## Scope of the diagnostic ranks

These descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.
