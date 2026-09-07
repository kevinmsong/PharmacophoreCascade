# Analysis report

Numbers here are the ones the manuscript reports, regenerated from
`evidence/outputs/`. Enrichment preserves score ties: molecules with equal status and
equal score form a tie group that no ligand identifier resolves, and EF, BEDROC, and
top-k recovery average over the possible within-tie orders.

Structural alerts do not exclude (`chemistry_gate_mode=warn_only`, native
`pains_filter=false`); PAINS and reactive-group matches are annotations. Molecular
standardization and the MW, LogP, HBD, and HBA limits are active. The shortlist keeps
the larger of the 5% count or 1,000 molecules, capped by eligible candidates.

## Million-compound GLP-1R screen

1,000,000 ZINC inputs. Stage 0 admits 997,590; 990,191 pass the hotspot gate; 49,880
enter 3D scoring; 49,757 receive Stage-3 scores. The native branch scores 4,997 ligands
successfully and emits 1,000 final-ranked candidates. End-to-end wall time is 12,336.7 s
(3.43 h) including the native branch once, on 12 workers. The 1,000-molecule minimum is
inactive at this library scale, since 5% of Stage-0 passers is already 49,880.

## Final active retention by shortlist rule

| System | Rule | Shortlist | Native scored | Final actives | Input actives |
|:---|:---|---:|---:|---:|---:|
| GLP-1R | 5% only | 1,464 | 647 | 8 | 10 |
| GLP-1R | 5% + min. 1,000 | 1,464 | 647 | 8 | 10 |
| GHSR | 5% only | 14 | 14 | 9 | 50 |
| GHSR | 5% + min. 1,000 | 277 | 271 | 39 | 50 |
| NTSR1 | 5% only | 20 | 20 | 18 | 50 |
| NTSR1 | 5% + min. 1,000 | 398 | 343 | 35 | 50 |
| MDM2-p53 | 5% only | 1,566 | 1,152 | 43 | 50 |
| MDM2-p53 | 5% + min. 1,000 | 1,566 | 1,152 | 43 | 50 |

Each pair shares one evaluation of Stages 0-2 and the larger Stage-3 union. Differing
shortlists get separate native preparation and scoring; identical shortlists share one
execution. Counts include measured preparation, selection, and scoring failures, rather
than treating shortlist admission as a successful native score. The extra native-branch
cost of the minimum is 2.4 to 8.8 min for GHSR and 3.3 to 8.9 min for NTSR1.

## Retrospective metrics

### GLP-1R

|                           |   roc_auc |   pr_auc |   ef_1pct |   bedroc |   top10_recovery |
|:--------------------------|----------:|---------:|----------:|---------:|-----------------:|
| full_cascade              |     0.753 |    0.363 |    23.250 |    0.480 |            0.401 |
| native_only               |     0.758 |    0.361 |    23.250 |    0.477 |            0.401 |
| stage3_only               |     0.727 |    0.284 |    15.500 |    0.357 |            0.300 |
| standard_3d_pharmacophore |     0.717 |    0.259 |    15.500 |    0.300 |            0.200 |

### GHSR

|                           |   roc_auc |   pr_auc |   ef_1pct |   bedroc |   top10_recovery |
|:--------------------------|----------:|---------:|----------:|---------:|-----------------:|
| full_cascade              |     0.867 |    0.610 |    31.000 |    0.740 |            0.200 |
| native_only               |     0.926 |    0.680 |    31.000 |    0.748 |            0.200 |
| stage3_only               |     0.852 |    0.408 |    19.375 |    0.589 |            0.120 |
| standard_3d_pharmacophore |     0.794 |    0.231 |    15.500 |    0.364 |            0.100 |

### NTSR1

|                           |   roc_auc |   pr_auc |   ef_1pct |   bedroc |   top10_recovery |
|:--------------------------|----------:|---------:|----------:|---------:|-----------------:|
| full_cascade              |     0.785 |    0.390 |    25.188 |    0.561 |            0.164 |
| native_only               |     0.690 |    0.127 |     8.221 |    0.312 |            0.059 |
| stage3_only               |     0.778 |    0.362 |    19.375 |    0.525 |            0.160 |
| standard_3d_pharmacophore |     0.560 |    0.142 |    13.563 |    0.250 |            0.100 |

### MDM2-p53

|                           |   roc_auc |   pr_auc |   ef_1pct |   bedroc |   top10_recovery |
|:--------------------------|----------:|---------:|----------:|---------:|-----------------:|
| full_cascade              |     0.945 |    0.705 |    31.000 |    0.820 |            0.200 |
| native_only               |     0.945 |    0.705 |    31.000 |    0.820 |            0.200 |
| stage3_only               |     0.339 |    0.023 |     0.000 |    0.028 |            0.000 |
| standard_3d_pharmacophore |     0.339 |    0.023 |     0.000 |    0.028 |            0.000 |

Native-only scoring exceeds the single-pass 3D pharmacophore on all four systems and
matches or exceeds the full cascade on three. Equivalence against a prespecified
+/-0.05 ROC-AUC margin is established only for MDM2-p53.

## Native diagnostics

Across the 4,997 successfully native-scored ligands: median coverage 24.65%
(IQR 23.25 to 24.87, range 15.09 to 29.97). Stage-3 and native coverage correlate at
Pearson r = 0.049, Spearman rho = 0.055. Ordering the same set by Stage-3 rank and by
final native rank gives rho = 0.057, their top-10 sets share 0 molecules, and the median
absolute within-cohort rank shift among the final top 1,000 is 1,940 positions. The
terminal stage, not the upstream ordering, determines the output.

Ranks are computed within the common native-scored set. Coverage-only ranks in native
diagnostic bundles are labeled separately and may differ because of their tie rules.

## Docking of the top 10

| ligand_id        |   final_rank |   best_active |   best_inactive |   active_pref |
|:-----------------|-------------:|--------------:|----------------:|--------------:|
| ZINCk700000Gfz1H |            1 |         -7.18 |           -5.71 |         -1.47 |
| ZINCk700001bQgxo |            2 |         -7.34 |           -6.22 |         -1.12 |
| ZINCk500000BlHKL |            3 |         -7.07 |           -5.88 |         -1.19 |
| ZINCj600000JWMq0 |            4 |         -8.45 |           -7.07 |         -1.39 |
| ZINCk500000gdMQW |            5 |         -7.05 |           -6.13 |         -0.92 |
| ZINCk600000EmAGZ |            6 |         -8.56 |           -7.25 |         -1.31 |
| ZINCk5000003XjkE |            7 |         -7.21 |           -6.12 |         -1.08 |
| ZINCk800000gvPBH |            8 |         -6.83 |           -6.40 |         -0.44 |
| ZINCk600000Ja5bY |            9 |         -8.70 |           -7.24 |         -1.46 |
| ZINCk600000y8tEo |           10 |         -7.42 |           -6.28 |         -1.14 |

All 10 favor the active state; median difference -1.16 kcal/mol (IQR -1.37 to -1.09),
exact signed-rank W = 0, two-sided p = 0.0020. The active-state ensemble has three
structures and the inactive-state ensemble two, so best-per-state sampling is unequal
and the paired test cannot isolate a receptor-state effect.

These computational measurements prioritize candidates for testing. They do not
demonstrate binding, state selectivity, or agonism.
