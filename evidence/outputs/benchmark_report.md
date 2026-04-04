# GLP-1R Evidence Package: Benchmark Report

## 1. Baseline Method Comparison

| method                    |   n_ranked |   ef_1pct |   ef_1pct_ci_low |   ef_1pct_ci_high |   ef_5pct |   ef_5pct_ci_low |   ef_5pct_ci_high |   ef_01pct |   ef_01pct_ci_low |   ef_01pct_ci_high |   roc_auc |   roc_auc_ci_low |   roc_auc_ci_high |   pr_auc |   pr_auc_ci_low |   pr_auc_ci_high |   bedroc |   bedroc_ci_low |   bedroc_ci_high |   top10_recovery |   top10_recovery_ci_low |   top10_recovery_ci_high |   scaffold_recovery_top10 |   scaffold_recovery_top10_ci_low |   scaffold_recovery_top10_ci_high |   top25_recovery |   top25_recovery_ci_low |   top25_recovery_ci_high |   scaffold_recovery_top25 |   scaffold_recovery_top25_ci_low |   scaffold_recovery_top25_ci_high |   top50_recovery |   top50_recovery_ci_low |   top50_recovery_ci_high |   scaffold_recovery_top50 |   scaffold_recovery_top50_ci_low |   scaffold_recovery_top50_ci_high |   scaffold_diversity_top10 |   scaffold_diversity_top10_ci_low |   scaffold_diversity_top10_ci_high |   scaffold_diversity_top25 |   scaffold_diversity_top25_ci_low |   scaffold_diversity_top25_ci_high |   scaffold_diversity_top50 |   scaffold_diversity_top50_ci_low |   scaffold_diversity_top50_ci_high |
|:--------------------------|-----------:|----------:|-----------------:|------------------:|----------:|-----------------:|------------------:|-----------:|------------------:|-------------------:|----------:|-----------------:|------------------:|---------:|----------------:|-----------------:|---------:|----------------:|-----------------:|-----------------:|------------------------:|-------------------------:|--------------------------:|---------------------------------:|----------------------------------:|-----------------:|------------------------:|-------------------------:|--------------------------:|---------------------------------:|----------------------------------:|-----------------:|------------------------:|-------------------------:|--------------------------:|---------------------------------:|----------------------------------:|---------------------------:|----------------------------------:|-----------------------------------:|---------------------------:|----------------------------------:|-----------------------------------:|---------------------------:|----------------------------------:|-----------------------------------:|
| full_cascade              |    310.000 |    30.000 |           10.000 |            40.000 |    10.000 |            4.000 |            16.000 |    100.000 |           100.000 |            100.000 |     0.800 |            0.636 |             0.940 |    0.465 |           0.177 |            0.789 |    0.546 |           0.229 |            0.834 |            0.500 |                   0.200 |                    0.700 |                     0.500 |                            0.250 |                             0.667 |            0.500 |                   0.200 |                    0.800 |                     0.500 |                            0.286 |                             0.800 |            0.500 |                   0.200 |                    0.800 |                     0.500 |                            0.286 |                             0.800 |                      7.000 |                               nan |                                nan |                     18.000 |                               nan |                                nan |                     32.000 |                               nan |                                nan |
| stage3_only               |    310.000 |    20.000 |            0.000 |            40.000 |     6.000 |            2.000 |            12.000 |    100.000 |             0.000 |            100.000 |     0.742 |            0.569 |             0.898 |    0.286 |           0.061 |            0.613 |    0.358 |           0.095 |            0.686 |            0.300 |                   0.000 |                    0.600 |                     0.300 |                            0.000 |                             0.500 |            0.300 |                   0.100 |                    0.700 |                     0.300 |                            0.143 |                             0.600 |            0.400 |                   0.148 |                    0.800 |                     0.400 |                            0.200 |                             0.714 |                     10.000 |                               nan |                                nan |                     19.000 |                               nan |                                nan |                     37.000 |                               nan |                                nan |
| standard_3d_pharmacophore |    310.000 |    20.000 |            0.000 |            40.000 |     4.000 |            0.000 |            10.000 |    100.000 |             0.000 |            100.000 |     0.727 |            0.570 |             0.878 |    0.261 |           0.047 |            0.573 |    0.304 |           0.047 |            0.627 |            0.200 |                   0.000 |                    0.500 |                     0.200 |                            0.000 |                             0.400 |            0.300 |                   0.100 |                    0.600 |                     0.300 |                            0.125 |                             0.500 |            0.400 |                   0.100 |                    0.700 |                     0.400 |                            0.167 |                             0.667 |                      9.000 |                               nan |                                nan |                     24.000 |                               nan |                                nan |                     40.000 |                               nan |                                nan |



## 2. Pairwise Rank Correlation

| method_a     | method_b                  |   n_shared |   kendall_tau |   spearman_r |   jaccard_top10 |   jaccard_top25 |   jaccard_top50 |
|:-------------|:--------------------------|-----------:|--------------:|-------------:|----------------:|----------------:|----------------:|
| full_cascade | stage3_only               |        310 |         0.585 |        0.741 |           0.176 |           0.136 |           0.250 |
| full_cascade | standard_3d_pharmacophore |        310 |         0.169 |        0.229 |           0.111 |           0.111 |           0.163 |
| stage3_only  | standard_3d_pharmacophore |        310 |         0.584 |        0.548 |           0.429 |           0.351 |           0.515 |



## 3. Paired Bootstrap Deltas vs Full Cascade

| method                    | metric                  |   delta_mean |   delta_ci_low |   delta_ci_high |
|:--------------------------|:------------------------|-------------:|---------------:|----------------:|
| stage3_only               | bedroc                  |       -0.182 |         -0.436 |          -0.000 |
| standard_3d_pharmacophore | bedroc                  |       -0.239 |         -0.512 |          -0.017 |
| stage3_only               | ef_01pct                |       -8.800 |       -100.000 |           0.000 |
| standard_3d_pharmacophore | ef_01pct                |       -8.800 |       -100.000 |           0.000 |
| stage3_only               | ef_1pct                 |      -11.900 |        -30.000 |           0.000 |
| standard_3d_pharmacophore | ef_1pct                 |      -12.120 |        -30.000 |           0.000 |
| stage3_only               | ef_5pct                 |       -3.908 |        -10.000 |           0.000 |
| standard_3d_pharmacophore | ef_5pct                 |       -5.608 |        -12.000 |           0.000 |
| stage3_only               | pr_auc                  |       -0.177 |         -0.430 |          -0.011 |
| standard_3d_pharmacophore | pr_auc                  |       -0.205 |         -0.463 |          -0.026 |
| stage3_only               | roc_auc                 |       -0.057 |         -0.137 |           0.017 |
| standard_3d_pharmacophore | roc_auc                 |       -0.069 |         -0.209 |           0.091 |
| stage3_only               | scaffold_recovery_top10 |       -0.199 |         -0.400 |           0.000 |
| standard_3d_pharmacophore | scaffold_recovery_top10 |       -0.264 |         -0.500 |           0.000 |
| stage3_only               | scaffold_recovery_top25 |       -0.170 |         -0.368 |           0.000 |
| standard_3d_pharmacophore | scaffold_recovery_top25 |       -0.202 |         -0.400 |           0.000 |
| stage3_only               | scaffold_recovery_top50 |       -0.039 |         -0.200 |           0.167 |
| standard_3d_pharmacophore | scaffold_recovery_top50 |       -0.109 |         -0.286 |           0.000 |
| stage3_only               | top10_recovery          |       -0.198 |         -0.500 |           0.000 |
| standard_3d_pharmacophore | top10_recovery          |       -0.260 |         -0.500 |           0.000 |
| stage3_only               | top25_recovery          |       -0.173 |         -0.400 |           0.000 |
| standard_3d_pharmacophore | top25_recovery          |       -0.204 |         -0.500 |           0.000 |
| stage3_only               | top50_recovery          |       -0.041 |         -0.300 |           0.200 |
| standard_3d_pharmacophore | top50_recovery          |       -0.109 |         -0.300 |           0.000 |



## 4. Ablation Results

| method                 |   n_shared |   kendall_tau |   spearman_r |   jaccard_top10 |   jaccard_top50 |   jaccard_top100 |   jaccard_top500 |
|:-----------------------|-----------:|--------------:|-------------:|----------------:|----------------:|-----------------:|-----------------:|
| remove_stage1          |       1000 |         0.118 |        0.176 |           0.000 |           0.010 |            0.010 |            0.038 |
| remove_stage2          |       1000 |        -0.020 |       -0.030 |           0.000 |           0.000 |            0.000 |            0.000 |
| remove_stage3          |         80 |         0.018 |        0.031 |           0.000 |           0.010 |            0.015 |            0.026 |
| native_top_stage3_only |         52 |         0.431 |        0.608 |           0.053 |           0.031 |            0.026 |            0.028 |
| alternative_tiebreak   |         80 |         0.022 |        0.026 |           0.000 |           0.010 |            0.005 |            0.016 |



## 5. Sensitivity Sweeps

|   hotspot_weight |   pair_hash_weight |   kendall_tau |   spearman_r |   jaccard_top10 |   jaccard_top50 |   jaccard_top100 |   jaccard_top500 | sweep_type         |   shortlist_fraction |   n_shortlisted |   jaccard_vs_baseline_top10 |   jaccard_vs_baseline_top50 |   jaccard_vs_baseline_top100 |   jaccard_vs_baseline_top500 |
|-----------------:|-------------------:|--------------:|-------------:|----------------:|----------------:|-----------------:|-----------------:|:-------------------|---------------------:|----------------:|----------------------------:|----------------------------:|-----------------------------:|-----------------------------:|
|            0.600 |              0.400 |         0.717 |        0.881 |           0.538 |           0.493 |            0.429 |            0.488 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|            0.500 |              0.500 |         0.835 |        0.958 |           0.667 |           0.667 |            0.639 |            0.678 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|            0.400 |              0.600 |         1.000 |        1.000 |           1.000 |           1.000 |            1.000 |            1.000 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|            0.200 |              0.800 |         0.685 |        0.857 |           0.250 |           0.449 |            0.471 |            0.504 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|            1.000 |              0.000 |         0.388 |        0.552 |           0.176 |           0.333 |            0.227 |            0.167 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|            0.000 |              1.000 |         0.489 |        0.654 |           0.111 |           0.220 |            0.220 |            0.304 | cascade_weights    |              nan     |         nan     |                     nan     |                     nan     |                      nan     |                      nan     |
|          nan     |            nan     |       nan     |      nan     |         nan     |         nan     |          nan     |          nan     | shortlist_fraction |                0.010 |         477.000 |                       0.000 |                       0.000 |                        0.000 |                        0.007 |
|          nan     |            nan     |       nan     |      nan     |         nan     |         nan     |          nan     |          nan     | shortlist_fraction |                0.025 |        1192.000 |                       0.000 |                       0.000 |                        0.000 |                        0.007 |
|          nan     |            nan     |       nan     |      nan     |         nan     |         nan     |          nan     |          nan     | shortlist_fraction |                0.050 |        2384.000 |                       0.000 |                       0.000 |                        0.000 |                        0.007 |
|          nan     |            nan     |       nan     |      nan     |         nan     |         nan     |          nan     |          nan     | shortlist_fraction |                0.100 |        4769.000 |                       0.000 |                       0.000 |                        0.000 |                        0.007 |



## 6. Manuscript Claim Evidence


| Claim | Verdict | Summary |
|---|---|---|

| claim1_native_adds_info | [SUPPORTS] | Pearson r = 0.038 between cascade score and native coverage (Spearman ρ = 0.000; median \|rank shift\| = 307), confirming that native reranking provides information orthogonal to the frontend cascade. |

| claim2_frontend_enriches_not_determines | [SUPPORTS] | Native top-100 ligands have a median screen percentile of 47.2%, with only 2% drawn from the cascade top-1%; the frontend enriches but does not constrain native performance. |

| claim3_native_not_cosmetic_stage3 | [SUPPORTS] | Stage-3 top-10 and native top-10 share Jaccard = 0.00 with Pearson r = 0.038 between their coverage scores; native scoring is not a cosmetic re-ordering of Stage-3 results. |

| claim4_diversified_pool_reduces_collapse | [INCONCLUSIVE] | Insufficient ligands in one or both pool-source groups. |

| claim5_residue_motif_convergence | [CONTRADICTS] | Only 0/4 claimed motif residues are enriched in top-10 native-ranked ligands; convergence on the stated motif is not confirmed. |


### Claim Detail Tables


#### claim1_native_adds_info

|   pearson_r_screen_vs_native |   spearman_r_screen_vs_native |   median_abs_rank_shift |   mean_abs_rank_shift |   max_abs_rank_shift |   n_ligands |
|-----------------------------:|------------------------------:|------------------------:|----------------------:|---------------------:|------------:|
|                       0.0375 |                        0.0002 |                307.0000 |              339.0880 |             962.0000 |   1000.0000 |



#### claim2_frontend_enriches_not_determines

|   native_top_k |   median_screen_pct |   pct_from_screen_top1 |   pct_from_screen_top5 |   pct_from_screen_top10 |   pct_from_screen_top50 |
|---------------:|--------------------:|-----------------------:|-----------------------:|------------------------:|------------------------:|
|        10.0000 |             47.7000 |                 0.0000 |                 0.0000 |                  0.0000 |                  0.5000 |
|        50.0000 |             44.9500 |                 0.0400 |                 0.0600 |                  0.1400 |                  0.5600 |
|       100.0000 |             47.2500 |                 0.0200 |                 0.0700 |                  0.1400 |                  0.5400 |



#### claim3_native_not_cosmetic_stage3

| metric                                |   value |
|:--------------------------------------|--------:|
| pearson_r_stage3_vs_native_coverage   |  0.0375 |
| spearman_r_stage3_rank_vs_native_rank | -0.0051 |
| jaccard_top10                         |  0.0000 |
| jaccard_top50                         |  0.0309 |
| jaccard_top100                        |  0.0753 |
| jaccard_top500                        |  0.3333 |



#### claim4_diversified_pool_reduces_collapse

|   n_diversified |   n_stage3_only | error                                   |
|----------------:|----------------:|:----------------------------------------|
|            1000 |               0 | insufficient sample size for comparison |



#### claim5_residue_motif_convergence

| residue   |   top_k |   n_top_k_with_residue |   freq_top_k |   freq_all_native |   fold_enrichment |
|:----------|--------:|-----------------------:|-------------:|------------------:|------------------:|
| SER14     |      10 |                      1 |       0.1000 |            0.9990 |            0.1001 |
| SER14     |      50 |                      4 |       0.0800 |            0.9990 |            0.0801 |
| SER14     |     100 |                      8 |       0.0800 |            0.9990 |            0.0801 |
| ASP15     |      10 |                      1 |       0.1000 |            0.9370 |            0.1067 |
| ASP15     |      50 |                      4 |       0.0800 |            0.9370 |            0.0854 |
| ASP15     |     100 |                      8 |       0.0800 |            0.9370 |            0.0854 |
| SER17     |      10 |                      1 |       0.1000 |            0.9970 |            0.1003 |
| SER17     |      50 |                      4 |       0.0800 |            0.9970 |            0.0802 |
| SER17     |     100 |                      8 |       0.0800 |            0.9970 |            0.0802 |
| SER18     |      10 |                      1 |       0.1000 |            0.9960 |            0.1004 |
| SER18     |      50 |                      4 |       0.0800 |            0.9960 |            0.0803 |
| SER18     |     100 |                      8 |       0.0800 |            0.9960 |            0.0803 |

