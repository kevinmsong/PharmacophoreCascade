# Reproducibility package

This directory lets anyone recompute every retrospective enrichment number reported
in the manuscript *Native Peptide-Contact Pharmacophore Scoring at Library Scale: A
Staged Cascade and Its Efficiency–Retention Trade-off* **without** the production
screening engine, directly from the released per-molecule benchmark scores.

Three other directories carry the rest of the released material:

| Directory | What it does |
|---|---|
| `reference_implementation/` | A readable reimplementation of Stage 0–2 scoring, written from the manuscript's equations alone, plus a script that checks it against the engine's released scores. See its README for the measured agreement. |
| `evidence/` | The analysis scripts behind the paper's figures and tables: active attrition, the efficiency–retention sweep, TOST equivalence testing and the paired docking test, the common-support ablation recomputation, the figure scripts, and three checkers covering figure accessibility, numerical consistency, and manuscript cross-references. |
| `evidence/outputs/` | Their outputs, including the per-molecule attrition table and the shortlist sweep. |

To regenerate the analyses and figures from released data:

```bash
python evidence/analyze_active_attrition.py --check
python evidence/analyze_efficiency_retention.py
python evidence/analyze_equivalence.py
python evidence/analyze_ablation_common_support.py
python evidence/make_ieee_figures.py                 # --target acs for the ACS build
python evidence/make_ieee_appendix_figures.py        # --target acs for the ACS build
python evidence/check_figure_accessibility.py
python reference_implementation/verify_against_released_scores.py
```

The manuscript sources are checked against these outputs rather than transcribed
by hand. `evidence/check_manuscript_numbers.py` re-derives every headline
statistic from the CSVs above and confirms the value appears in the LaTeX, and
`evidence/check_response_crossrefs.py` confirms the response letter names the
figures and tables the built manuscript actually numbers that way:

```bash
python evidence/check_manuscript_numbers.py --dir <submission folder>
python evidence/check_response_crossrefs.py --dir <submission folder>
```

## What is released here

The conclusions rest on the retrospective benchmarks, and the full benchmark evidence
is public in this repository:

| Artifact | Location |
|---|---|
| Per-molecule scored tables (label, per-method rank/score/status, descriptors, ChEMBL/decoy provenance) | `evidence/data/machine_readable/{glp1r,ghsr,ntsr1,mdm2}_benchmark_scored.csv` |
| Benchmark configurations (gates, conformers, tolerances, seeds) | `evidence/configs/benchmark_*.yaml` |
| Pre-computed benchmark summaries with bootstrap CIs | `evidence/outputs/benchmark_*/benchmark_summary.csv` |
| Decoy-replicate summaries (5 sets/system) | `evidence/outputs/decoy_replicates/` |
| Receptor + native pharmacophores | `maps/` |
| Pinned environment | `requirements.txt` |

Each scored table has one row per benchmarked molecule (10 actives + 300 matched
decoys = 310 per system) with the `label` column and, for each method
(`full_cascade`, `native_only`, `stage3_only`, `standard_3d_pharmacophore`), a
`*_rank`, `*_score`, and `*_status` column. Molecules that failed an intermediate
stage carry a non-`scored` status and are ranked below scored molecules, exactly as
in the manuscript.

## Reproduce the reported metrics

```bash
pip install -r requirements.txt
python reproduce/reproduce_benchmarks.py            # all four systems
python reproduce/reproduce_benchmarks.py --system glp1r
```

This recomputes, per method and system, ROC-AUC, PR-AUC, EF1%, EF5%, BEDROC
(alpha = 20), and top-k active recovery. The output matches the manuscript, e.g.
for GLP-1R the full cascade gives ROC-AUC = 0.800, PR-AUC = 0.465, EF1% = 30,
BEDROC = 0.546, 5/10 actives in the top 10; the native-only baseline matches it
(ROC-AUC = 0.809); and across systems the native-only baseline equals or exceeds
the full cascade, confirming that the terminal native scoring supplies the
enrichment while the staged cascade supplies scalability.

### Metric definitions (match the Methods section)

- **ROC-AUC / PR-AUC** — global ranking quality (scikit-learn); PR-AUC is more
  sensitive to the low active prevalence.
- **EF k%** — fraction of actives recovered in the top `ceil(k% x N)` ranked
  molecules, divided by `k/100`.
- **BEDROC (alpha = 20)** — early-recognition-weighted score (Truchon & Bayly,
  *J. Chem. Inf. Model.* 2007, 47, 488-508).
- **top-k recovery** — number of actives ranked within the top k positions.

Grouped-bootstrap 95% confidence intervals and paired p-values (each active
resampled with its matched decoys) are provided pre-computed in the
`benchmark_*/benchmark_summary.csv` tables. `evidence/analyze_equivalence.py`
recomputes the bootstrap from scratch and adds the TOST equivalence tests and the
paired Wilcoxon test on the docking deltas.

## What is not released

The production screening engine that executes the million-compound scan is under
active development for separate applications and is therefore not open-sourced; it
is available from the corresponding author under a reasonable-use agreement. It is
not required to reproduce any reported result.

Two things make that claim checkable rather than asserted. First, the released
per-molecule scores are sufficient to recompute every enrichment metric with the
script in this directory. Second, `reference_implementation/` reimplements the
Stage 0–2 scoring equations independently of the engine and verifies them against
its released output; it reproduces the Stage-0 and Stage-1 gate outcomes for 100% of
the GLP-1R benchmark and the continuous scores to within 0.03–0.05 percentage points
on average. Writing it surfaced two errors in earlier descriptions of the method — a
rotatable-bond bound in the Stage-0 equation that the code does not apply, and a
receptor-feature selection rule described as "top N by weight" when it is not — both
of which are corrected in the current manuscript.

Stage 3 and the terminal native branch are not reimplemented, because both depend on
RDKit conformer embedding rather than on a closed-form expression. Their per-molecule
outputs are released instead.
