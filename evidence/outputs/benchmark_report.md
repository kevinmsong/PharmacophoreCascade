# Current alert-disabled revision

Structural alerts are recorded without exclusion; property limits remain. The headline shortlist has a capped minimum of 1,000. Retrospective statistics below use the freshly executed GLP-1R benchmark; internal ranking and motif diagnostics use the new headline cohort. Ranking perturbations reorder observed score tables and do not imply gate-removal reruns.

# GLP-1R Evidence Package: Benchmark Report

## 1. Baseline Method Comparison

| method                    |   n_ranked |   roc_auc |   pr_auc |   ef_1pct |   ef_5pct |   ef_01pct |   top10_recovery |   top25_recovery |   top50_recovery |   bedroc |   roc_auc_ci_low |   roc_auc_ci_high |   pr_auc_ci_low |   pr_auc_ci_high |   ef_1pct_ci_low |   ef_1pct_ci_high |   ef_5pct_ci_low |   ef_5pct_ci_high |   ef_01pct_ci_low |   ef_01pct_ci_high |   top10_recovery_ci_low |   top10_recovery_ci_high |   top25_recovery_ci_low |   top25_recovery_ci_high |   top50_recovery_ci_low |   top50_recovery_ci_high |   bedroc_ci_low |   bedroc_ci_high |
|:--------------------------|-----------:|----------:|---------:|----------:|----------:|-----------:|-----------------:|-----------------:|-----------------:|---------:|-----------------:|------------------:|----------------:|-----------------:|-----------------:|------------------:|-----------------:|------------------:|------------------:|-------------------:|------------------------:|-------------------------:|------------------------:|-------------------------:|------------------------:|-------------------------:|----------------:|-----------------:|
| full_cascade              |        310 |     0.753 |    0.363 |    23.250 |     7.926 |     31.000 |            0.401 |            0.421 |            0.453 |    0.480 |            0.592 |             0.889 |           0.152 |            0.633 |            8.004 |            31.000 |            2.574 |            13.562 |            31.000 |             31.000 |                   0.117 |                    0.700 |                   0.158 |                    0.707 |                   0.200 |                    0.733 |           0.177 |            0.743 |
| native_only               |        310 |     0.758 |    0.361 |    23.250 |     7.893 |     31.000 |            0.401 |            0.417 |            0.443 |    0.477 |            0.693 |             0.866 |           0.150 |            0.626 |            7.959 |            31.000 |            2.505 |            13.562 |            31.000 |             31.000 |                   0.115 |                    0.700 |                   0.151 |                    0.706 |                   0.200 |                    0.728 |           0.172 |            0.742 |
| stage3_only               |        310 |     0.727 |    0.284 |    15.500 |     5.812 |     31.000 |            0.300 |            0.300 |            0.400 |    0.357 |            0.545 |             0.887 |           0.072 |            0.558 |            0.000 |            31.000 |            1.938 |            11.625 |             0.000 |             31.000 |                   0.000 |                    0.500 |                   0.100 |                    0.600 |                   0.200 |                    0.800 |           0.114 |            0.632 |
| standard_3d_pharmacophore |        310 |     0.717 |    0.259 |    15.500 |     3.875 |     31.000 |            0.200 |            0.300 |            0.400 |    0.300 |            0.596 |             0.848 |           0.057 |            0.530 |            0.000 |            31.000 |            0.000 |             9.688 |             0.000 |             31.000 |                   0.000 |                    0.500 |                   0.000 |                    0.600 |                   0.100 |                    0.700 |           0.041 |            0.581 |



## 2. Pairwise Rank Correlation

| method_a                  | method_b                  |   n_shared |   kendall_tau |   spearman_r |   jaccard_top10 |   jaccard_top25 |   jaccard_top50 |
|:--------------------------|:--------------------------|-----------:|--------------:|-------------:|----------------:|----------------:|----------------:|
| full_cascade              | stage3_only               |        310 |         0.574 |        0.731 |           0.176 |           0.136 |           0.235 |
| full_cascade              | standard_3d_pharmacophore |        310 |         0.151 |        0.206 |           0.111 |           0.111 |           0.163 |
| full_cascade              | native_only               |        310 |         0.667 |        0.689 |           1.000 |           0.852 |           0.695 |
| stage3_only               | standard_3d_pharmacophore |        310 |         0.576 |        0.538 |           0.333 |           0.316 |           0.493 |
| stage3_only               | native_only               |        310 |         0.242 |        0.346 |           0.176 |           0.111 |           0.190 |
| standard_3d_pharmacophore | native_only               |        310 |         0.266 |        0.381 |           0.111 |           0.087 |           0.136 |



## 3. Paired Bootstrap Deltas vs Full Cascade

| method                    | metric         |   delta_mean |   delta_ci_low |   delta_ci_high |   n_bootstrap |   p_value |   p_value_holm |
|:--------------------------|:---------------|-------------:|---------------:|----------------:|--------------:|----------:|---------------:|
| native_only               | roc_auc        |        0.007 |         -0.065 |           0.119 |          5000 |     0.957 |          1.000 |
| native_only               | pr_auc         |       -0.003 |         -0.013 |           0.008 |          5000 |     0.559 |          0.726 |
| native_only               | ef_1pct        |       -0.006 |         -0.090 |           0.000 |          5000 |     1.000 |          1.000 |
| native_only               | ef_5pct        |       -0.045 |         -0.246 |           0.000 |          5000 |     0.752 |          0.752 |
| native_only               | ef_01pct       |       -0.001 |          0.000 |           0.000 |          5000 |     1.000 |          1.000 |
| native_only               | top10_recovery |       -0.001 |         -0.006 |           0.000 |          5000 |     1.000 |          1.000 |
| native_only               | top25_recovery |       -0.005 |         -0.022 |           0.000 |          5000 |     0.682 |          0.869 |
| native_only               | top50_recovery |       -0.012 |         -0.050 |           0.000 |          5000 |     0.682 |          1.000 |
| native_only               | bedroc         |       -0.004 |         -0.018 |          -0.000 |          5000 |     0.007 |          0.020 |
| stage3_only               | roc_auc        |       -0.026 |         -0.085 |           0.024 |          5000 |     0.368 |          1.000 |
| stage3_only               | pr_auc         |       -0.091 |         -0.311 |           0.063 |          5000 |     0.363 |          0.726 |
| stage3_only               | ef_1pct        |       -7.405 |        -23.250 |           6.200 |          5000 |     0.587 |          1.000 |
| stage3_only               | ef_5pct        |       -2.023 |         -6.004 |           0.000 |          5000 |     0.251 |          0.501 |
| stage3_only               | ef_01pct       |       -2.961 |        -31.000 |           0.000 |          5000 |     1.000 |          1.000 |
| stage3_only               | top10_recovery |       -0.139 |         -0.400 |           0.000 |          5000 |     0.192 |          0.384 |
| stage3_only               | top25_recovery |       -0.083 |         -0.300 |           0.132 |          5000 |     0.435 |          0.869 |
| stage3_only               | top50_recovery |       -0.001 |         -0.300 |           0.233 |          5000 |     1.000 |          1.000 |
| stage3_only               | bedroc         |       -0.115 |         -0.328 |           0.018 |          5000 |     0.197 |          0.197 |
| standard_3d_pharmacophore | roc_auc        |       -0.035 |         -0.145 |           0.088 |          5000 |     0.516 |          1.000 |
| standard_3d_pharmacophore | pr_auc         |       -0.123 |         -0.359 |           0.047 |          5000 |     0.204 |          0.612 |
| standard_3d_pharmacophore | ef_1pct        |       -7.909 |        -23.357 |           6.200 |          5000 |     0.576 |          1.000 |
| standard_3d_pharmacophore | ef_5pct        |       -3.472 |         -7.836 |           0.000 |          5000 |     0.074 |          0.221 |
| standard_3d_pharmacophore | ef_01pct       |       -2.999 |        -31.000 |           0.000 |          5000 |     1.000 |          1.000 |
| standard_3d_pharmacophore | top10_recovery |       -0.197 |         -0.409 |           0.000 |          5000 |     0.090 |          0.269 |
| standard_3d_pharmacophore | top25_recovery |       -0.136 |         -0.329 |           0.000 |          5000 |     0.151 |          0.454 |
| standard_3d_pharmacophore | top50_recovery |       -0.087 |         -0.300 |           0.100 |          5000 |     0.567 |          1.000 |
| standard_3d_pharmacophore | bedroc         |       -0.173 |         -0.397 |          -0.001 |          5000 |     0.047 |          0.095 |



## 4. Ablation Results

| method                                |   n_shared |   tau_shared |   kendall_tau |   spearman_r |   n_common |   tau_common |   rho_common |   jaccard_top10 |   jaccard_top50 |   jaccard_top100 |   jaccard_top500 |
|:--------------------------------------|-----------:|-------------:|--------------:|-------------:|-----------:|-------------:|-------------:|----------------:|----------------:|-----------------:|-----------------:|
| Hotspot-score ordering                |       1000 |       -0.023 |        -0.023 |       -0.034 |      49757 |        0.057 |        0.070 |           0.000 |           0.000 |            0.000 |            0.000 |
| Pair-overlap ordering                 |       1000 |        0.107 |         0.107 |        0.160 |      49757 |        0.066 |        0.081 |           0.000 |           0.010 |            0.010 |            0.032 |
| Cascade score only within native pool |       1000 |        0.055 |         0.055 |        0.081 |      49757 |        0.420 |        0.433 |           0.000 |           0.020 |            0.015 |            0.098 |
| Observed Stage-3 top-5,000 subset     |        424 |        1.000 |         0.230 |        1.000 |      49757 |        0.282 |        0.285 |           0.111 |           0.299 |            0.274 |            0.269 |
| Native-fit RMSD ordering              |       1000 |        0.681 |         0.681 |        0.667 |      49757 |        0.417 |        0.427 |           0.000 |           0.000 |            0.000 |            0.020 |



## 5. Sensitivity Sweeps

|   hotspot_weight |   pair_hash_weight |   kendall_tau |   spearman_r |   jaccard_top10 |   jaccard_top50 |   jaccard_top100 |   jaccard_top500 | sweep_type      | analysis_scope                                                                         |
|-----------------:|-------------------:|--------------:|-------------:|----------------:|----------------:|-----------------:|-----------------:|:----------------|:---------------------------------------------------------------------------------------|
|            0.600 |              0.400 |         0.715 |        0.880 |           0.538 |           0.471 |            0.429 |            0.490 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |
|            0.500 |              0.500 |         0.834 |        0.958 |           0.667 |           0.639 |            0.653 |            0.681 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |
|            0.400 |              0.600 |         1.000 |        1.000 |           1.000 |           1.000 |            1.000 |            1.000 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |
|            0.200 |              0.800 |         0.684 |        0.856 |           0.250 |           0.408 |            0.449 |            0.499 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |
|            1.000 |              0.000 |         0.386 |        0.549 |           0.176 |           0.333 |            0.227 |            0.167 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |
|            0.000 |              1.000 |         0.488 |        0.653 |           0.111 |           0.176 |            0.205 |            0.302 | cascade_weights | Reordering observed successful Stage-3 table; no upstream rescoring or new native pool |



## 6. Manuscript Claim Evidence


| Claim | Verdict | Summary |
|---|---|---|

| claim1_native_adds_info | [SUPPORTS] | Pearson r = 0.049 between cascade score and native coverage (Spearman ρ = 0.055; median \|rank shift\| = 1440), confirming that native reranking provides information orthogonal to the frontend cascade. |

| claim2_frontend_enriches_not_determines | [INCONCLUSIVE] | Native top-100 ligands have a median screen percentile of 46.9%, with only 1% drawn from the cascade top-1%; native ordering draws from a broad range of observed screen ranks. These are within-cohort rank diagnostics; no native scores for excluded molecules were measured, so frontend enrichment is not established by this analysis. |

| claim3_native_not_cosmetic_stage3 | [SUPPORTS] | Stage-3 top-10 and native top-10 share Jaccard = 0.00 with Pearson r = 0.049 between their coverage scores; native scoring is not a cosmetic re-ordering of Stage-3 results. |

| claim4_diversified_pool_reduces_collapse | [INCONCLUSIVE] | Insufficient ligands in one or both pool-source groups. |

| claim5_residue_motif_convergence | [CONTRADICTS] | Only 0/4 claimed motif residues are enriched in top-10 native-ranked ligands; convergence on the stated motif is not confirmed. |


### Claim Detail Tables


#### claim1_native_adds_info

|   pearson_r_screen_vs_native |   spearman_r_screen_vs_native |   median_abs_rank_shift |   mean_abs_rank_shift |   max_abs_rank_shift |   n_ligands |
|-----------------------------:|------------------------------:|------------------------:|----------------------:|---------------------:|------------:|
|                       0.0491 |                        0.0550 |               1440.0000 |             1623.1887 |            4939.0000 |   4997.0000 |



#### claim2_frontend_enriches_not_determines

|   native_top_k |   median_screen_pct |   pct_from_screen_top1 |   pct_from_screen_top5 |   pct_from_screen_top10 |   pct_from_screen_top50 |
|---------------:|--------------------:|-----------------------:|-----------------------:|------------------------:|------------------------:|
|        10.0000 |             67.1303 |                 0.0000 |                 0.1000 |                  0.1000 |                  0.2000 |
|        50.0000 |             44.4166 |                 0.0200 |                 0.0600 |                  0.0800 |                  0.5600 |
|       100.0000 |             46.8881 |                 0.0100 |                 0.0400 |                  0.0700 |                  0.5300 |



#### claim3_native_not_cosmetic_stage3

| metric                                |   value |
|:--------------------------------------|--------:|
| pearson_r_stage3_vs_native_coverage   |  0.0491 |
| spearman_r_stage3_rank_vs_native_rank |  0.0567 |
| jaccard_top10                         |  0.0000 |
| jaccard_top50                         |  0.0101 |
| jaccard_top100                        |  0.0101 |
| jaccard_top500                        |  0.0515 |



#### claim4_diversified_pool_reduces_collapse

|   n_diversified |   n_stage3_only | error                                   |
|----------------:|----------------:|:----------------------------------------|
|            1000 |               0 | insufficient sample size for comparison |



#### claim5_residue_motif_convergence

| residue   |   top_k |   n_top_k_with_residue |   freq_top_k |   freq_all_native |   fold_enrichment |
|:----------|--------:|-----------------------:|-------------:|------------------:|------------------:|
| SER14     |      10 |                     10 |       1.0000 |            0.9756 |            1.0250 |
| SER14     |      50 |                     49 |       0.9800 |            0.9756 |            1.0045 |
| SER14     |     100 |                     99 |       0.9900 |            0.9756 |            1.0148 |
| ASP15     |      10 |                     10 |       1.0000 |            0.7030 |            1.4224 |
| ASP15     |      50 |                     48 |       0.9600 |            0.7030 |            1.3655 |
| ASP15     |     100 |                     98 |       0.9800 |            0.7030 |            1.3940 |
| SER17     |      10 |                     10 |       1.0000 |            0.9286 |            1.0769 |
| SER17     |      50 |                     29 |       0.5800 |            0.9286 |            0.6246 |
| SER17     |     100 |                     79 |       0.7900 |            0.9286 |            0.8508 |
| SER18     |      10 |                     10 |       1.0000 |            0.8895 |            1.1242 |
| SER18     |      50 |                     44 |       0.8800 |            0.8895 |            0.9893 |
| SER18     |     100 |                     94 |       0.9400 |            0.8895 |            1.0567 |

