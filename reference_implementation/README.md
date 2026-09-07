# Reference implementation of Stage 0–2 cascade scoring

This directory lets a third party check that the equations printed in the
manuscript are the equations the screen actually computes.

`stage_scoring_reference.py` is a plain, unoptimised reimplementation of
Stages 0–2 written from the published equations and parameter values alone. It
shares **no code with the production screening engine** and makes no attempt to
be fast; it exists to be read and checked, not to screen a million compounds.

## What is implemented

| Manuscript | Implemented here |
|---|---|
| Eq. (1) Stage-0 gate | `stage0()` — standardization and property envelope; PAINS/reactive matches are annotations without exclusion |
| Eq. (2) hotspot fraction `H` | `stage1()` — exact / compatible family matching at γ = 0.5, plus the three-part gate |
| Eq. (4) typed pair-hash | `ligand_pair_catalogue()`, `receptor_pair_query()`, `stage2()` — recall, precision, and the F1 overlap `O_pair` |
| Eq. (3) cascade score | `evaluate()` — `100 × (0.4 H + 0.6 O_pair)` |

Stage 3 and the terminal native branch are **not** reimplemented: both need 3D
conformer generation and alignment, which depends on RDKit's embedding rather
than on a closed-form expression. The released per-molecule Stage-3 and native
scores support reanalysis of those stages directly.

## Run it

```bash
# Score the bundled examples (actives, decoys, and each kind of gate failure)
python reference_implementation/stage_scoring_reference.py --examples

# Score arbitrary molecules
python reference_implementation/stage_scoring_reference.py \
    --smiles "CC(C)(C)Nc1nc2cc(Cl)c(Cl)cc2nc1S(C)(=O)=O"

# Check the reimplementation against the engine's released scores
python reference_implementation/verify_against_released_scores.py --benchmark-evaluation evidence/outputs/benchmark_glp1r_full/benchmark_external/benchmark_evaluation.csv --audit results/absolute_floor_1000/screening_full_1M_floor1000_shortlist.csv --output-dir evidence/outputs/alert_disabled_revision/reference_verification
```

## Measured agreement with the production engine

`verify_against_released_scores.py` compares this implementation against the
released per-molecule GLP-1R benchmark table (310 molecules) and against a
random sample of the headline million-compound run's audit table.

| Quantity | Agreement |
|---|---|
| Stage-0 pass/fail | **100.0%** (310 / 310 molecules) |
| Stage-1 gate pass/fail | **100.0%** (310 / 310 Stage-0 survivors) |
| Hotspot fraction `H` | Pearson *r* = 0.994349; mean abs. difference 0.046373 percentage points |
| Pair overlap `O_pair` | Pearson *r* = 0.999537; mean abs. difference 0.007533 pp |
| Cascade score | Pearson *r* = 0.996551; mean abs. difference 0.017926 pp |
| Eq. (3) on the engine's own columns | max absolute error 1.07 × 10⁻¹⁴ pp |

The continuous check uses 500 randomly sampled molecules from the fresh
49,880-molecule headline shortlist (seed 42). Maximum absolute residuals are
2.586892 pp for H, 0.757347 pp for pair overlap, and 1.450871 pp for cascade
score. The implementation is not numerically identical for every molecule;
the full residual records are included with the verification summary.

A separate Stage-0 check (`evidence/verify_alert_disabled_gate.py`) agrees on
all 4,960 labeled molecules across four systems. All 34 alert-flagged actives
are admitted. Of 160 actives, 158 pass Stage 0; two NTSR1 compounds fail the
retained property limits. An optional `exclude_alerts=True` argument exists
only for explicit historical comparisons; the default is non-excluding.

## Two things this check surfaced

Reimplementing from the published description exposed two places where the
manuscript's earlier wording did not match the code. Both are corrected in the
current manuscript, and the corrections are what allow the agreement above:

1. **Rotatable-bond count is not a Stage-0 gate condition.** The screen
   computes and reports it, but `passes_property_gate` tests only MW, LogP,
   HBD, and HBA. Eq. (1) previously listed a rotatable-bond bound.
2. **Receptor features are not simply the top *N* by weight.** They are
   collapsed to one per (residue, complementary family), ordered by curated
   priority, and extended with up to eight native-supported residues, with
   backfill disabled. The revised scoring also pins the original fixed
   typed-feature caps (2, 2, 4, 6, 6, 6), rather than later size-scaled defaults.

## Scope

This is a verification aid, not the production system. The engine that executes
the million-compound scan is under active development for separate applications
and is not released; it is available from the corresponding author under a
reasonable-use agreement. The retrospective enrichment metrics are reproducible from the accompanying
per-molecule scores with `reproduce/reproduce_benchmarks.py`, and the scoring
functions themselves are checkable here.

## Requirements

`rdkit`, `numpy`, `pandas` — all pinned in the repository's `requirements.txt`.
