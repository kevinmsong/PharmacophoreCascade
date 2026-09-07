"""Refresh project/reproduction prose and the current quantitative analysis report."""
import json
import shutil
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
META=ROOT/'evidence/outputs/alert_disabled_revision'
SUB=ROOT/'ACS_Omega_resubmission'


def replace(path,text):
    backup=ROOT/'tmp/alert_disabled_revision/documentation_before'/path.relative_to(ROOT)
    if path.exists() and not backup.exists():
        backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,backup)
    path.write_text(text,encoding='utf-8')


def main():
    data=json.loads((META/'manuscript_results.json').read_text())
    h=data['head'];n=h['native_rerank']['counts'];c=h['counts']
    settings='''Benchmark enrichment preserves score ties; EF, BEDROC, and top-k recovery average over within-tie orderings. Identifier-ordered enrichment estimates are superseded (Table S15). This evaluation correction does not alter molecular scores or observed production retention.

Structural-alert exclusions are disabled: `chemistry_gate_mode=warn_only` and native `pains_filter=false`. PAINS/reactive matches remain annotations. Molecular standardization and the MW, LogP, HBD, and HBA limits remain active.

The revised shortlist keeps the larger of the 5% count or 1,000 molecules, capped by eligible candidates. The headline run preserves its Stage-0-pass percentage denominator and 0.4/0.6 hotspot/pair weights. Production-policy pairs preserve their original Stage-1/2 candidate denominator and 0.25/0.75 weights. Typed-feature caps are fixed at 2, 2, 4, 6, 6, 6. Retrospective method comparisons are uncapped; the native-only comparator bypasses Stages 0–3 by definition.
'''
    headline=f'''The fresh screen evaluates 1,000,000 original ZINC inputs. Stage 0 admits {c['property_pass']:,}; {c['hotspot_pass']:,} pass the hotspot gate; {c['shortlist_size']:,} enter 3D scoring; and {c['final_hits']:,} receive successful Stage-3 scores. The native branch successfully scores {n['successful_best_ligands']:,} ligands and emits {n['final_rows']:,} final-ranked candidates. End-to-end wall time is {h['timings']['total_pipeline_sec']/3600:.2f} hours, including the native branch once. The 1,000-molecule floor is nonbinding at this library scale. Comparison with the archived strict-filter run does not isolate the floor, because alert handling also changed.
'''
    pairs=pd.DataFrame(data['floor']['benchmarks'])
    tab=pairs[['system_label','policy','shortlist_n','native_success_n','final_ranked_actives','n_actives']].rename(
        columns={'system_label':'System','policy':'Policy','shortlist_n':'Shortlist','native_success_n':'Native scored',
                 'final_ranked_actives':'Final actives','n_actives':'Input actives'}).to_markdown(index=False)
    overview='''# Native peptide-contact pharmacophore screening: current ACS revision

Materials for **High-throughput native peptide-contact pharmacophore scoring at library scale: a staged cascade and its efficiency-retention trade-off**, by Kevin Song, John Zhang, Lei Ye, and Jianyi Zhang.

The current submission is in [`ACS_Omega_resubmission/`](ACS_Omega_resubmission/). The 2026-09-06 revision recomputes the headline screen, four expanded retrospective comparisons, curation and 7KI0 reference checks, 20 original decoy-replicate experiments, four paired production-policy experiments, and docking of the new final top 10. Current completion manifests distinguish executed results from historical material.

## Settings used

'''+settings+'''
## Fresh headline result

'''+headline+'''
## Measured production-policy pairs

'''+tab+'''

Each pair shares fresh Stages 0–2 and the larger Stage-3 union. Different shortlists receive separate native preparation/scoring; identical shortlists explicitly share one fresh native execution. Counts include actual preparation, selection, and scoring failures, rather than treating shortlist admission as a successful native score.

## Current artifacts

| Artifact | Location |
|---|---|
| Clean and highlighted manuscript, SI, reviewer response | `ACS_Omega_resubmission/` (Word, PDF, LaTeX) |
| 600-dpi figures and graphical abstract | `ACS_Omega_resubmission/` (PNG and vector PDF where applicable) |
| Fresh headline tables and native bundle | `results/absolute_floor_1000/` |
| Established headline filename aliases | `results/screening_full_1M_topological_hashed*` and `results/top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv` |
| Full retrospective experiments | `evidence/outputs/benchmark_{glp1r,ghsr,ntsr1,mdm2}_full/` |
| Curation and reference checks | `evidence/outputs/benchmark_glp1r_automated/`, `benchmark_glp1r_7ki0/` |
| Original 10-active decoy replicates | `evidence/outputs/decoy_replicates/` |
| Paired production inputs, scores, survival | `evidence/outputs/absolute_floor/` |
| Machine-readable benchmark records | `evidence/data/machine_readable/` |
| New top-10 docking jobs | `evidence/outputs/docking_top10/` |
| Complete current quantitative manuscript object | `evidence/outputs/alert_disabled_revision/manuscript_results.json` |
| Commands, source/input hashes, execution manifests, validation | `evidence/outputs/alert_disabled_revision/` |
| Independent Stage-0–2 implementation | `reference_implementation/` |
| Recompute tie-aware enrichment from per-molecule records | `reproduce/reproduce_benchmarks.py` |

The established native-bundle and root GLP-1R benchmark paths are local junction aliases to the fresh outputs. Historical outputs are archived under `tmp/alert_disabled_revision/` and are not the current scientific results.

The earlier `publication/`, `ACS_Omega_submission/`, and `IEEE_Access_submission/` folders are retained as historical manuscript drafts. Their `SUPERSEDED.md` files point to the current ACS resubmission; use `ACS_Omega_resubmission/` for the updated submission documents.

## Reproduction

See [`reproduce/README.md`](reproduce/README.md). Fresh rerun drivers are `evidence/run_absolute_floor_headline.py` and `evidence/run_alert_disabled_suite.py`; their completed manifests record the exact commands. `refresh_alert_disabled_headline.py` promotes the new headline and executes the new docking study. `build_alert_disabled_deliverables.py` first replaces identifier-ordered enrichment with tie-aware evaluation and then rebuilds dependent analyses after verifying completion. Run `reevaluate_tied_benchmarks.py` after molecular benchmark execution and before consuming its summaries; the original molecular runner emits legacy identifier-ordered statistics as an intermediate output. These revision drivers deliberately refuse blind overwriting of an already recorded run; use their manifests to inspect or resume a partial build.

The headline uses 12 workers. Benchmark jobs initially request two workers; future jobs may use four after the production lane finishes, eight after the headline finishes and 12 after docking finishes. Each such job records its full effective configuration and actual worker count in `benchmark_execution.json`. Future native invocations may increase workers as the headline and docking finish; `native_worker_events` records each branch allocation. This changes execution parallelism, not molecular parameters or per-ligand seeds. The fresh serial/parallel check covers 51 prepared states from five ligands, with identical prepared SDF content and scores/mappings/coordinates agreeing to 1e-12. Measured times reflect concurrent workstation workloads and are not portable speedup estimates.

## Interpretation

Retrospective retrieval, peptide-feature coverage, and docking scores do not establish biological activity. Source-related active/decoy bias, limited active sets and receptor sampling, and missing prospective binding/signaling validation remain limitations. Ranking perturbations reorder observed score tables; they do not infer scores for excluded molecules or substitute for end-to-end gate-removal experiments. Native-only Stage-3 diagnostic placeholders are not measurements and are excluded from the current diagnostic figures.
'''
    replace(ROOT/'README.md',overview)
    for legacy in ['publication','ACS_Omega_submission','IEEE_Access_submission']:
        if (ROOT/legacy).exists():
            replace(ROOT/legacy/'SUPERSEDED.md',
                '# Historical manuscript draft\n\n'
                'This folder preserves an earlier submission or dated draft. Its results and figures are superseded. '
                'The current manuscript, Supporting Information, reviewer response, abstract, and figure assets are in '
                '[ACS_Omega_resubmission](../ACS_Omega_resubmission/README.md). '
                'The [project README](../README.md) identifies the fresh scientific outputs and execution records.\n')
    report='# Current analysis report: non-excluding structural alerts\n\n'+settings+'\n'+headline+'\n## Paired final retention\n\n'+tab+'\n\n## Retrospective metrics\n\n'
    for key,label in [('glp1r_full','GLP-1R'),('ghsr_full','GHSR'),('ntsr1_full','NTSR1'),('mdm2_full','MDM2–p53')]:
        frame=pd.DataFrame.from_dict(data['benchmarks'][key],orient='index')[['roc_auc','pr_auc','ef_1pct','bedroc','top10_recovery']]
        report+='### '+label+'\n\n'+frame.to_markdown(floatfmt='.3f')+'\n\n'
    report+='## Native diagnostics\n\n```json\n'+json.dumps(data['diagnostics'],indent=2)+'\n```\n\n'
    report+='The diagnostic ranks use the same successfully native-scored cohort. Coverage-only ranks in native diagnostic bundles are separately labeled and may differ because of their tie rules.\n\n'
    report+='## Docking of the new top 10\n\n'+pd.read_csv(ROOT/'evidence/outputs/revision/docking_top10_summary.csv').to_markdown(index=False,floatfmt='.2f')+'\n\n'
    report+='These computational measurements prioritize candidates for testing; they do not demonstrate binding, state selectivity, or agonism.\n'
    replace(ROOT/'ANALYSIS_REPORT.md',report)
    reproduction='''# Reproduce the current revision

The current per-molecule records contain 310 GLP-1R molecules (10 actives, 300 decoys) and 1,550 molecules for each other full benchmark (50 actives, 1,500 decoys). The separate 20 decoy-replicate runs use the original 10-active, 300-decoy sets. Do not conflate those populations.

'''+settings+'''
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
'''
    replace(ROOT/'reproduce/README.md',reproduction)
    replace(SUB/'README.md','''# Current ACS Omega revision deliverables

This folder contains the clean manuscript (`main`), highlighted manuscript (`main_highlighted`), Supporting Information, and response to reviewers in Word, PDF, and LaTeX forms. Raster figures and the graphical abstract are exported at 600 dpi; vector PDF counterparts are supplied where applicable.

'''+settings+'\n'+headline+'\nThe current data and execution records are described in the [project README](../README.md). The final artifact inventory records paths and SHA-256 hashes.\n')
    print('Updated project, analysis, reproduction, and submission documentation')


if __name__=='__main__':main()
