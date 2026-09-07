# Native peptide-contact pharmacophore screening

Data, analysis code, and an independent reference implementation for
**High-throughput native peptide-contact pharmacophore scoring at library scale: a staged cascade and its efficiency-retention trade-off**,
by Kevin Song, John Zhang, Lei Ye, and Jianyi Zhang.

The paper asks which stage of a staged virtual-screening cascade actually supplies its
retrieval. The answer is the terminal one: a pharmacophore built from the bound peptide's
own receptor contacts beats a conventional single-pass 3D pharmacophore on all four
systems tested, and matches or exceeds the full cascade on three of them. The upstream
stages supply tractability, at a cost in active retention that shortlist depth controls.

## Headline screen

One million ZINC compounds (tranches H17-H20) screened against GLP-1R:

| Stage | Molecules |
|---|---|
| Input library | 1,000,000 |
| Stage 0, standardization and property limits | 997,590 |
| Stage 1, receptor hotspot gate | 990,191 |
| Cascade-score shortlist (5%, minimum 1,000) | 49,880 |
| Successful Stage-3 scores | 49,757 |
| Native candidate pool, then scaffold-capped selection | 20,000 to 5,000 |
| Successful native ligand scores | 4,997 |
| Final ranked output | 1,000 |

End-to-end wall time 12,336.7 s (3.43 h) on an eight-core, 16-thread workstation with 12
workers, including the native branch once. Three-dimensional work runs on 5% of the
library and terminal peptide-contact scoring on 0.5%.

## Retrospective enrichment

ROC-AUC on matched labeled universes (10 actives with 300 decoys for GLP-1R; 50 with
1,500 for the others). Enrichment preserves score ties: molecules with equal status and
equal score are tied, no identifier resolves a tie, and EF, BEDROC, and top-k recovery
average over the possible within-tie orders.

| System | Full cascade | Native-only | Stage-3 only | Single-pass 3D |
|---|---|---|---|---|
| GLP-1R | 0.753 | 0.758 | 0.727 | 0.717 |
| GHSR | 0.867 | 0.926 | 0.852 | 0.794 |
| NTSR1 | 0.785 | 0.690 | 0.778 | 0.560 |
| MDM2-p53 | 0.945 | 0.945 | 0.339 | 0.339 |

## Shortlist depth against active retention

Paired production runs hold every other setting fixed and change only the shortlist rule.

| System | Shortlist (5% only) | Shortlist (5% + min. 1,000) | Actives retained | Native branch |
|---|---|---|---|---|
| GLP-1R | 1,464 | 1,464 | 8/10 both | 17.6 min, shared |
| GHSR | 14 | 277 | 9/50 to 39/50 | 2.4 to 8.8 min |
| NTSR1 | 20 | 398 | 18/50 to 35/50 | 3.3 to 8.9 min |
| MDM2-p53 | 1,566 | 1,566 | 43/50 both | 41.7 min, shared |

A fixed percentage is the wrong control once an upstream gate has already reduced the
pool to a few hundred molecules. At million-compound scale the percentage term already
exceeds 1,000, so the minimum never binds. Structural alerts are recorded rather than
excluding, so 158 of 160 actives clear Stage 0; the two exceptions fail the property
limits after standardization.

## Settings

Molecular standardization and the MW, LogP, HBD, and HBA limits are active.
`chemistry_gate_mode=warn_only` and native `pains_filter=false`, so PAINS and
reactive-group matches are annotations rather than exclusions. The shortlist keeps the
larger of the 5% count or 1,000 molecules, capped by eligible candidates. The headline
run uses a Stage-0-pass percentage denominator and 0.4/0.6 hotspot/pair weights;
production-policy pairs use the Stage-1/2 candidate denominator and 0.25/0.75 weights.
Typed-feature caps are 2, 2, 4, 6, 6, 6. Retrospective method comparisons are uncapped,
and the native-only comparator bypasses Stages 0-3 by definition.

## What is here

| Artifact | Location |
|---|---|
| Independent Stage 0-2 implementation, and its verification script | `reference_implementation/` |
| Recompute every reported enrichment metric from per-molecule records | `reproduce/` |
| Analysis, benchmark, and figure source | `evidence/` |
| Machine-readable per-molecule benchmark records | `evidence/data/machine_readable/` |
| Retrospective experiments, all four systems | `evidence/outputs/benchmark_{glp1r,ghsr,ntsr1,mdm2}_full/` |
| Curation and 7KI0 reference checks | `evidence/outputs/benchmark_glp1r_automated/`, `benchmark_glp1r_7ki0/` |
| Decoy replicates | `evidence/outputs/decoy_replicates/` |
| Paired shortlist-policy inputs, scores, and survival | `evidence/outputs/absolute_floor/` |
| Top-10 docking | `evidence/outputs/docking_top10/` |
| Ranked screening outputs and native bundle | `results/` |
| Receptor and native pharmacophores, stage configs, structures | `config/`, `maps/`, `pocket/`, `structures/` |
| 600-dpi figures and the graphical abstract | `ACS_Omega_resubmission/` |

The production screening engine that executes the million-compound scan is not released;
it is available from the corresponding author under a reasonable-use agreement. Everything
needed to recompute the reported statistics from the released per-molecule data, and to
check the Stage 0-2 scoring functions against the published equations, is in this
repository.

## Reproduction

See [`reproduce/README.md`](reproduce/README.md) and
[`reference_implementation/README.md`](reference_implementation/README.md).
Install with `pip install -r requirements.txt` (Python 3.12.3, RDKit 2024.09.6,
scikit-learn 1.3.2, NumPy 1.26.4, SciPy 1.11.4, pandas 2.2.3, AutoDock Vina 1.2.7,
Meeko 0.7.1).

Measured times reflect concurrent workstation workloads and are not portable speedup
estimates. A serial/parallel check on 51 prepared microstates gave byte-identical
prepared SDFs, with scores, feature mappings, and coordinates agreeing to 1e-12.

## Interpretation

Retrospective retrieval, peptide-feature coverage, and docking scores do not establish
biological activity. Actives and decoys come from different databases, so property
matching and repeated decoy draws cannot eliminate source-related bias. Only one non-GPCR
interface is included. Ranking perturbations reorder observed score tables; they do not
infer scores for excluded molecules and are not end-to-end gate-removal experiments.
Prospective binding and signaling measurements remain necessary.
