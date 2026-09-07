# Reproduce the current revision

The current per-molecule records contain 310 GLP-1R molecules (10 actives, 300 decoys) and 1,550 molecules for each other full benchmark (50 actives, 1,500 decoys). The separate 20 decoy-replicate runs use the original 10-active, 300-decoy sets. Do not conflate those populations.

Benchmark enrichment preserves score ties; EF, BEDROC, and top-k recovery average over within-tie orderings. Identifier-ordered enrichment estimates are superseded (Table S15). This evaluation correction does not alter molecular scores or observed production retention.

Structural-alert exclusions are disabled: `chemistry_gate_mode=warn_only` and native `pains_filter=false`. PAINS/reactive matches remain annotations. Molecular standardization and the MW, LogP, HBD, and HBA limits remain active.

The revised shortlist keeps the larger of the 5% count or 1,000 molecules, capped by eligible candidates. The headline run preserves its Stage-0-pass percentage denominator and 0.4/0.6 hotspot/pair weights. Production-policy pairs preserve their original Stage-1/2 candidate denominator and 0.25/0.75 weights. Typed-feature caps are fixed at 2, 2, 4, 6, 6, 6. Retrospective method comparisons are uncapped; the native-only comparator bypasses Stages 0–3 by definition.

## Recompute measurements from the fresh score tables

Run from the project root with the environment in `requirements.txt`:

```bash
python evidence/reevaluate_tied_benchmarks.py
python evidence/export_machine_readable.py
python reproduce/reproduce_benchmarks.py
python evidence/verify_current_metrics.py
python evidence/validate_tie_aware_metrics.py
python evidence/validate_floor_outputs.py
python evidence/analyze_absolute_floor.py
python evidence/analyze_alert_disabled_attrition.py
python evidence/aggregate_replicates.py
python evidence/analyze_equivalence.py
python evidence/analyze_ablation_common_support.py
python evidence/collect_alert_disabled_results.py
python evidence/make_ieee_figures.py --target acs
python evidence/refresh_headline_diagnostic_figures.py
python evidence/make_ieee_appendix_figures.py --target acs
python evidence/make_floor_toc_graphic.py
python evidence/check_manuscript_numbers.py --dir ACS_Omega_resubmission
python evidence/check_reviewer_coverage.py --dir ACS_Omega_resubmission
python evidence/check_figure_accessibility.py --dir ACS_Omega_resubmission
```

The independent reference check and its exact input paths are documented in `reference_implementation/README.md`. Its current verification exports quantify residual numerical differences rather than asserting exact identity. Stage 3 and native conformer/alignment scoring are not independently reimplemented; their fresh outputs and serial/parallel validation are provided.

ROC-AUC and average precision preserve equal-status/equal-score ties, including missing scores within a failure status. EF, BEDROC, and top-k retrieval are expected values averaged over within-tie orderings. Expected active counts may be fractional. An audit found that ligand-ID tie breaks could favor ChEMBL actives; those identifier-ordered enrichment estimates are superseded. Original ranks remain available for traceability, and evaluation tie groups are released explicitly. EF k% uses the fraction of actives recovered among the top ceil(k% × N), divided by the actual selected fraction ceil(k% × N)/N. This corrects the earlier nominal-fraction denominator, which overstated EF when the percentage cutoff rounded upward. BEDROC uses alpha=20. Grouped bootstrap resamples each active together with its matched decoys (5,000 iterations, seed 42). Equivalence requires the whole 90% paired-difference interval within ±0.05 ROC-AUC. Docking uses an exact signed-rank sign-flip test with average ranks for ties, zero differences excluded, and only successful paired state results.

## Execute the molecular experiments again

The exact executed commands, input/source hashes, settings, completion records, and effective worker configurations are in `evidence/outputs/alert_disabled_revision/` and `results/absolute_floor_1000/rerun_provenance.json`. The run drivers create fresh scores and prepared conformers and deliberately refuse to overwrite an already recorded experiment without inspecting its state. The current build driver supports `--resume` for completed analysis phases; inspect a failed phase before retrying a mutation that archives outputs.

## Document rendering

Build all four LaTeX documents first to update figure/table references. Then run `python evidence/make_acs_docx.py` and render the Word files with `evidence/render_word_documents.py`. Numerical and cross-reference checks supplement visual inspection of every final page; they do not substitute for it. The current publication figures are 600-dpi PNGs with vector PDFs where applicable.

The manuscript's code-availability statement governs external distribution of the production engine. Per-molecule metric reproduction does not require that engine. This local revision operation does not publish materials to an external repository.
