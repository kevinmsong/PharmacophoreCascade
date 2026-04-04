# GLP-1 / GLP1R Small Molecule Mimic Discovery

## Chart-Locked Final Workflow Report

**Generated:** 2026-03-13 05:51:16
**Workflow name:** `Topological-to-3D Native Pharmacophore Cascade`
**Primary workflow:** `run_optimized_1M_topological_hashed_screening.py`
**Primary native output:** `results\top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv`
**Stage-3 audit table:** `results\screening_full_1M_topological_hashed.csv`
**Native bundle manifest:** `results\screening_full_1M_topological_hashed_native_terminal_bundle\manifests\analysis_manifest.json`
**Run summary:** `results\screening_full_1M_topological_hashed_run_summary.json`

---

## 1. Project Overview

This workflow follows the finalized chart contract for the Topological-to-3D Native Pharmacophore Cascade:
a topological prescreen, a stage-3 geometry-aware rerank, and a native
GLP-1/GLP1R terminal rerank that produces the final top-1000 ligand table.
Stage-3 tables and plots remain required audit artifacts, but they are no
longer the primary endpoint.

The default cascade is:

1. Stage 0 property + chemistry gate
2. Stage 1 physicochemical hotspot bitmask
3. Stage 2 typed pair-hash prescreen
4. Stage 3 geometry-aware rerank on the top 5% shortlist
5. Diversified 20,000-ligand native pool with 12k / 4k / 4k source quotas
6. Scaffold-aware native selection of up to 5,000 ligands under a Murcko cap
7. Native reference-feature scoring with best conformer per microstate and best state per ligand
8. Final native-first ranking to the primary top-1000 output

---

## 2. Structural Context

- Interface residues from `maps/glp1r_interface.json`: 77
- GLP1R chain(s): R
- Residue list preview: THR29, VAL30, SER31, LEU32, TRP33, THR35, VAL36, TRP39, ARG43, PHE66

---

## 3. Curated Receptor Pharmacophore

- Total pharmacophore features: 102
- Curated feature count: 21
- Curated receptor residues: 17
- Native-supported residues: 8
- Curated weight bonus: +2.0

### Curated Contact Groups

### ECD anchoring

- GLP-1 residues: `Phe28, Ile29, Leu32, Val33`
- GLP1R residues: `LEU32, TRP39, ASP67, TYR69, ARG121, LEU123, GLU128`

### Upper TMD activation pocket

- GLP-1 residues: `His7, Glu9, Thr13, Ser14, Ser17, Ser18`
- GLP1R residues: `TYR145, ARG190, LYS197, TRP297, THR298, ARG299, LEU388, SER392`

### ECL1 support

- GLP-1 residues: `Trp31`
- GLP1R residues: `GLN211, HIS212`


### Top Weighted Features

| Rank | Type | Residue | Weight | Curated Group | GLP-1 Residues |
|-----:|:-----|:--------|-------:|:--------------|:---------------|
| 1 | negative | GLU128 OE1 | 6.35 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |
| 2 | positive | ARG299 NH1 | 5.91 | Upper TMD activation pocket | His7, Glu9, Thr13, Ser14, Ser17, Ser18 |
| 3 | negative | ASP67 OD1 | 5.85 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |
| 4 | positive | ARG121 NH1 | 5.84 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |
| 5 | positive | LYS197 NZ | 5.76 | Upper TMD activation pocket | His7, Glu9, Thr13, Ser14, Ser17, Ser18 |
| 6 | positive | ARG190 NH1 | 5.63 | Upper TMD activation pocket | His7, Glu9, Thr13, Ser14, Ser17, Ser18 |
| 7 | aromatic | HIS212 CD2 | 5.29 | ECL1 support | Trp31 |
| 8 | aromatic | TYR145 CD1 | 5.23 | Upper TMD activation pocket | His7, Glu9, Thr13, Ser14, Ser17, Ser18 |
| 9 | aromatic | TYR145 CZ | 5.23 | Upper TMD activation pocket | His7, Glu9, Thr13, Ser14, Ser17, Ser18 |
| 10 | aromatic | TYR69 CD1 | 5.13 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |
| 11 | aromatic | TYR69 CZ | 5.13 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |
| 12 | aromatic | TRP39 CD1 | 5.12 | ECD anchoring | Phe28, Ile29, Leu32, Val33 |

---

## 4. Query Composition

### Stage 1 Hotspot Query

Selection metadata not available.

### Stage 2 Pair Query

Selection metadata not available.

### Stage 3 Geometry Query

Selection metadata not available.

---

## 5. Primary Native Results

Primary native-final results not yet available. Run the chart-locked workflow first.


---

## 6. Stage-3 Audit Outputs

Stage-3 audit outputs are not yet available.


---

## 7. Diagnostics

- Diagnostics status: `unavailable`
- Always-on contract: `True`
- Reason: `missing_diagnostics_metadata`
- No diagnostics files were recorded for this run.

---

## 8. Pipeline Timing

Run timing summary not yet available.


---

## 9. File Inventory

### Primary Native Outputs

- `results\top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv`
- `results\screening_full_1M_topological_hashed_native_scored_top5000.csv`
- `results\screening_full_1M_topological_hashed_native_terminal_bundle\manifests\analysis_manifest.json`

### Stage-3 Audit Outputs

- `results\screening_full_1M_topological_hashed.csv`
- `results\top_100_glp1_mimetics_full_1M_topological_hashed.csv`
- `results/pharmacophore_3d_full_1M_topological_hashed.png`
- `results/top_20_glp1_mimetics_full_1M_topological_hashed.png`
- `results/property_distributions_full_1M_topological_hashed.png`
- `results\screening_full_1M_topological_hashed_run_summary.json`

### Workflow Source

- `virtual_screening_pipeline.mmd`

*Report generated from the chart-locked Topological-to-3D Native Pharmacophore Cascade workflow.*
