# Topological-to-3D Native Pharmacophore Cascade

Data and workflow documentation accompanying:

> **Topological-to-3D Native Pharmacophore Cascade for Ultra-High-Throughput Ligand Discovery at Peptide-Receptor Interfaces**
>
> Kevin Song, Lei Ye, Jianyi Zhang
>
> Department of Biomedical Engineering, The University of Alabama at Birmingham
>
> *Journal of Chemical Information and Modeling* (2026)

---

## Overview

This repository contains the data files, configuration, structural inputs, benchmark evidence, and result tables produced by the Topological-to-3D Native Pharmacophore Cascade applied to GLP-1R small-molecule mimetic discovery. The cascade screens ~1 million ZINC compounds through four sequential stages plus a terminal native branch:

1. **Stage 0** -- Physicochemical property and chemistry gating
2. **Stage 1** -- Receptor-informed hotspot bitmask compatibility scoring
3. **Stage 2** -- Typed pair-hash relational filtering
4. **Stage 3** -- Geometry-aware 3D conformer reranking
5. **Native Branch** -- Diversified pool assembly, scaffold-aware selection, and terminal peptide-contact reranking against the active GLP-1/GLP1R complex

The final output is a ranked set of 1,000 peptide-mimicking small molecules scored for native peptide-contact mimicry.

### Pipeline Architecture

```
 1M ZINC compounds
      |
 Stage 0: Property + chemistry gate
      |
 Stage 1: Hotspot bitmask scoring
      |
 Stage 2: Pair-hash prescreen
      |
 Shortlist (top 5% by cascade score)
      |
 Stage 3: 3D conformer rerank
      |
      +---> Stage-3 audit outputs (CSV, plots)
      |
 Diversified native pool (20,000 ligands)
      |
 Scaffold-aware selection (5,000 ligands)
      |
 Native scoring vs. 6X18 peptide-contact pharmacophore
      |
 Final native top-1,000 ranking
```

The full workflow diagram is in [`virtual_screening_pipeline.mmd`](virtual_screening_pipeline.mmd) (Mermaid format).

---

## Repository Contents

### Structural Inputs

| Path | Description |
|------|-------------|
| `structures/6X18_GLP1_GLP1R.pdb` | Active GLP-1/GLP1R cryo-EM complex (PDB 6X18) |
| `structures/AF-P43220-F1-model_v6_GLP1R.pdb` | AlphaFold GLP1R model |
| `structures/complex.pdb` | Processed complex used for pharmacophore derivation |
| `pocket/mimic_pocket.pdb` | Extracted binding pocket |

### Pharmacophore and Interface Data

| Path | Description |
|------|-------------|
| `maps/pharmacophore_rigorous.json` | 102-feature receptor-side pharmacophore with curated contact annotations |
| `maps/glp1r_interface.json` | GLP1R interface residue map (77 residues) |

### Configuration

| Path | Description |
|------|-------------|
| `config/pipeline.yaml` | Pipeline configuration |
| `evidence/configs/benchmark.yaml` | Retrospective benchmark configuration |
| `evidence/configs/ablation.yaml` | Ablation study configuration |

### Primary Results

| Path | Description |
|------|-------------|
| `results/top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv` | **Primary output:** final native-ranked top-1,000 ligands |
| `results/screening_full_1M_topological_hashed_native_scored_top5000.csv` | Full native-scored 5,000-ligand table |
| `results/top_100_glp1_mimetics_full_1M_topological_hashed.csv` | Stage-3 audit top-100 |
| `results/screening_full_1M_topological_hashed_run_summary.json` | Run parameters and summary statistics |

### Result Figures

| Path | Description |
|------|-------------|
| `results/pharmacophore_3d_full_1M_topological_hashed.png` | 3D pharmacophore visualization |
| `results/top_20_glp1_mimetics_full_1M_topological_hashed.png` | Top-20 ligand structures |
| `results/property_distributions_full_1M_topological_hashed.png` | Property distribution plots |
| `publication/publication_figures/figure1_pipeline.png` | Pipeline architecture (Figure 1) |
| `publication/publication_figures/figure2_pharmacophore.png` | Receptor pharmacophore (Figure 2) |
| `publication/publication_figures/figure3_top_20_structures.png` | Top-20 structures (Figure 3) |
| `publication/publication_figures/figure4_top_native_ligand_structural_overlap.png` | Structural overlap analysis (Figure 4) |
| `publication/publication_figures/figure5_correlation_grid.png` | Metric correlation grid (Figure 5) |
| `publication/publication_figures/figure6_native_weighted_coverage_profile.png` | Coverage profile (Figure 6) |
| `publication/publication_figures/figure7_top_native_peptide_mimicry_report.png` | Peptide mimicry report (Figure 7) |

### Native Terminal Analysis

Summary analysis files from the native terminal reranking stage are in `results/screening_full_1M_topological_hashed_native_terminal_bundle/`:

- `analysis/` -- Correlation metrics, quartiles, rank comparison, hotspot frequencies, reference features, receptor contacts
- `reports/` -- Diagnostic plots (correlation scatter, rank reranking, quartiles, residuals) and figure captions
- `manifests/analysis_manifest.json` -- Machine-readable output inventory

### Benchmark and Evidence Data

| Path | Description |
|------|-------------|
| `evidence/data/glp1r_external_benchmark_library.csv` | 310-molecule benchmark library (10 actives + 300 matched decoys) |
| `evidence/data/glp1r_external_benchmark_exclusions.csv` | Curated ChEMBL exclusion log |
| `evidence/outputs/benchmark_summary.csv` | Benchmark summary metrics (ROC-AUC, PR-AUC, EF, BEDROC) |
| `evidence/outputs/benchmark_plots.pdf` | Benchmark performance plots |
| `evidence/outputs/benchmark_report.md` | Benchmark narrative report |
| `evidence/outputs/ablation_results.csv` | Ablation study results |
| `evidence/outputs/ablation_sensitivity.pdf` | Ablation sensitivity plots |
| `evidence/outputs/claim*_evidence.csv` | Per-claim supporting evidence tables |
| `evidence/outputs/claims_summary.csv` | Claims evidence summary |
| `evidence/outputs/*_table.tex` | LaTeX table fragments used in the manuscript |
| `evidence/outputs/benchmark_external/` | External benchmark ranking CSVs per method |

### GLP1R Ligand Analysis

| Path | Description |
|------|-------------|
| `GLP1_top_ligand_analysis/configs/` | Study configuration files |
| `GLP1_top_ligand_analysis/references/` | Reference PDB structures (5VEW, 6LN2, 7KI0, 7LCJ) |
| `GLP1_top_ligand_analysis/reports_715_topological/` | QC and screening comparison reports |

### Other Files

| Path | Description |
|------|-------------|
| `virtual_screening_pipeline.mmd` | Full pipeline flowchart (Mermaid) |
| `ANALYSIS_REPORT.md` | Detailed workflow report with stage summaries |
| `requirements.txt` | Python dependencies |

---

## Key Results Summary

- **Library:** 1,000,000 ZINC compounds (tranches H17--H20)
- **Stage-3 shortlist:** 47,812 ligands (top 5% by cascade score)
- **Native pool:** 20,000 diversified candidates (12k by Stage-3 rank, 4k by hotspot breadth, 4k by native-supported hotspot score)
- **Native scored:** 4,997 ligands after scaffold-aware selection and preparation
- **Final output:** 1,000 native-ranked peptide-mimicking small molecules

### Retrospective Benchmark (310 molecules: 10 actives + 300 decoys)

| Method | ROC-AUC | PR-AUC | EF1% | EF5% | BEDROC |
|--------|---------|--------|------|------|--------|
| Full Cascade | 0.800 | 0.465 | 30 | 10 | 0.546 |
| Stage-3 Only | 0.742 | 0.286 | 20 | 6 | 0.358 |
| Standard 3D Pharmacophore | 0.727 | 0.261 | 20 | 4 | 0.304 |

---

## Final Ranking Criteria

Ligands in the primary output are sorted by:

1. `native_weighted_coverage_pct` (descending)
2. `native_matched_reference_features` (descending)
3. `native_fit_rmsd_angstrom` (ascending)
4. `native_mean_pair_distance_error_angstrom` (ascending)
5. `stage3_screen_rank` (ascending)
6. `zinc_id` (ascending)

---

## Note

This repository provides data, configuration, and results for reproducibility and review. Source code for the pipeline implementation is not included in this public release.

---

## License

All rights reserved. This repository is provided for peer review and academic reference in connection with the accompanying manuscript submission.
