"""Refresh author replies while preserving all reviewer quotations verbatim."""
import json
import re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUB=ROOT/'ACS_Omega_resubmission'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def main():
    d=json.loads((META/'manuscript_results.json').read_text());b=d['benchmarks']
    p=pd.DataFrame(d['floor']['benchmarks']).set_index(['system','policy'])
    metric=lambda s,m,k:float(b[s][m][k])
    def pair(s):return f"{int(p.loc[(s,'percentage'),'final_ranked_actives'])}/50 to {int(p.loc[(s,'floor1000'),'final_ranked_actives'])}/50"
    eq=pd.read_csv(ROOT/'evidence/outputs/equivalence/equivalence_tost.csv')
    eq=eq[eq.metric=='roc_auc']
    dock=json.loads((ROOT/'evidence/outputs/equivalence/docking_wilcoxon.json').read_text())
    h=d['head'];nc=h['native_rerank']['counts'];ref=d['reference_verification']['continuous']
    replies={}
    replies['1.1']=(rf'''We have separated terminal scoring performance from the upstream selection mechanism. The new comparisons use the same labeled universe within each system and include native-only, full-cascade, Stage-3-only, and conventional single-pass 3D methods. All runs were recomputed from molecular inputs with structural-alert exclusions disabled while retaining property limits. The full numerical comparison is in Table~2 and Figures~2 and 3. We do not infer that each cascade gate improves retrieval merely because the terminal scoring method is useful.

The revised headline screen was also executed from the original million-compound inputs with the minimum-shortlist rule. It shortlists {h['counts']['shortlist_size']:,} molecules, successfully native-scores {nc['successful_best_ligands']:,} ligands, and writes {nc['final_rows']:,} final candidates in {h['timings']['total_pipeline_sec']/3600:.2f}~h. This measured total includes the native branch once. We removed the earlier extrapolation to exhaustive million-ligand native scoring because that comparison was not measured.''',
        r'The abstract, Results, Discussion, Table~2, and Table~S1 now use the completed revised experiments and measured runtime.')
    replies['1.2']=(rf'''We now test the shortlist intervention directly rather than recommend an unexecuted setting. The rule retains the larger of 5\% and 1,000 molecules, capped by the eligible pool. Fresh paired production runs hold non-excluding structural alerts and all scoring settings fixed while changing only the shortlist minimum. The floor changes final active retention from {pair('ghsr')} for GHSR and from {pair('ntsr1')} for NTSR1. Shortlist admission, successful native scoring, and final ranking are counted separately.

The floor is nonbinding in the million-compound run because its percentage term exceeds 1,000. The newly computed headline cohort differs from the archived strict-filter run because alert handling also changed; Table~S12 explicitly does not interpret that archive comparison as a floor-only experiment. Table~S14 reports measured additional native-branch work for the matched policy pairs.''',
        r'Figure~5, Table~5, Figure~S5, and Tables~S12--S14 report the executed intervention; Tables~4, S10, and S11 account for actual active survival.')
    reps=pd.DataFrame(d['replicates'])
    reptext=' '.join(f"{label}: {reps[(reps.system==s)&(reps.method=='full_cascade')].roc_auc.mean():.3f} $\\pm$ {reps[(reps.system==s)&(reps.method=='full_cascade')].roc_auc.std(ddof=1):.3f}." for s,label in [('glp1r','GLP-1R'),('ghsr','GHSR'),('ntsr1','NTSR1'),('mdm2','MDM2--p53')])
    replies['1.3']=(r'''We retain the distinction between retrospective enrichment and prospective activity. Every labeled benchmark and all 20 decoy-replicate experiments were rescored with the revised alert settings. The expanded comparisons contain 50 actives for each additional target; GLP-1R remains limited to 10. The five decoy replicates per system use the original 10-active sets. Full-cascade ROC-AUC mean $\pm$ sample standard deviation across those replicates is: '''+reptext+r'''

These repeated draws test variation in decoy selection but do not remove a possible systematic difference between ChEMBL actives and ZINC decoys. We state that limitation directly. Per-molecule labels, molecular structures, curation identifiers, and method ranks/scores are provided for independent recalculation. Experimental binding and signaling validation remain outstanding.''',
        r'Table~S8, Methods, the decoy-robustness Results subsection, and Discussion distinguish expanded benchmarks, decoy replicates, source bias, and prospective validation.')
    replies['2.1']=(rf'''We recomputed stage survival per active, including the outcomes after shortlisting. Structural-alert exclusions are now disabled in both the screening engine and native preparation, while property limits remain. An independent implementation agrees with all 4,960 labeled Stage-0 decisions. Of 160 actives, 158 pass; all 34 alert-flagged actives are admitted, and two NTSR1 actives still fail property limits.

The floor is then evaluated as a separate intervention under those same alert settings. Final active retention changes from {pair('ghsr')} for GHSR and from {pair('ntsr1')} for NTSR1. Every active absent from the final ranking has a recorded first limiting stage, including hotspot gating, shortlist selection, preparation, native-pool/scaffold selection, and final ranking. The table no longer treats admission to native processing as a successful native score.''',
        r'Figure~4 and Tables~4, S10, and S11 report individual causes and measured survival; Figure~5 and Table~5 compare the executed shortlist policies.')
    eqtext='\n\n'.join(rf'{r.system}: full-minus-native ROC-AUC difference {r.observed_delta:+.3f}, 90\% interval [{r.ci90_low:+.3f}, {r.ci90_high:+.3f}]; '+('equivalence established.' if r.equivalent else 'equivalence not established.') for r in eq.itertuples())
    replies['2.2']=(r'''We no longer use a non-significant difference as evidence of equivalence. The paired grouped-bootstrap analysis was recomputed from the new method scores with equal-score ties retained, with a prespecified ROC-AUC equivalence margin of $\pm0.05$. Equivalence requires the entire 90\% interval inside that margin. The results are:

'''+eqtext+r'''

The small GLP-1R active set limits precision. We report the observed differences, interval widths, and equivalence outcome together.''',
        r'Table~3 and the equivalence Results subsection report the revised intervals and explicit decision rule; Methods describes matched-active-group resampling.')
    replies['2.3']=(rf'''The earlier strict-filter interpretation has been replaced by the new measured result. At MDM2--p53, full-cascade ROC-AUC is {metric('mdm2_full','full_cascade','roc_auc'):.3f}, native-only is {metric('mdm2_full','native_only','roc_auc'):.3f}, and conventional single-pass 3D is {metric('mdm2_full','standard_3d_pharmacophore','roc_auc'):.3f}. The corresponding full-cascade and native-only BEDROC values are {metric('mdm2_full','full_cascade','bedroc'):.3f} and {metric('mdm2_full','native_only','bedroc'):.3f}.

We separate this measured target result from a general claim about interface geometry. A single protein-protein system is insufficient to establish a flat-interface boundary of applicability. We also retain the distinction between overall discrimination and early retrieval, rather than describing either as a substitute for the other.''',
        r'Table~2, Figure~3, the cross-system Results subsection, and Discussion use the fresh MDM2--p53 data and remove the unsupported geometry generalization.')
    replies['2.4']=(rf'''The new final top 10 were prepared and redocked against the same five receptor structures, with three seeds per structure and exhaustiveness 16. Successful paired active/inactive results were obtained for {dock['n']}/10 ligands; {dock['n_favoring_active']} favor the active state. The median active-minus-inactive difference is {dock['median_delta']:.2f}~kcal/mol (interquartile range {dock['iqr_low']:.2f} to {dock['iqr_high']:.2f}). The exact signed-rank permutation test gives $W={dock['wilcoxon_W']:g}$ and two-sided $p={dock['p_two_sided']:.4f}$.

The test enumerates sign flips, uses average ranks for tied absolute differences, and excludes zero differences. Preparation and docking failures remain in the job accounting. The best-score comparison samples three active-state structures and two inactive-state structures; these unequal opportunities mean the paired test cannot isolate a receptor-state effect. We describe this as an orthogonal computational score comparison, not evidence of binding, state selectivity, or agonism.''',
        r'Figure~6, Table~S9, the docking Results subsection, and Methods report the new ligands, measured differences, exact test, and limits of interpretation.')
    replies['2.5']=(rf'''The released independent implementation was updated to the authorized non-excluding alert setting and the fixed typed-feature caps used by the screening runs. It agrees on Stage-0 decisions for all 4,960 labeled benchmark molecules; on the 310-molecule GLP-1R benchmark it also reproduces Stage-1 decisions. On 500 randomly sampled freshly shortlisted headline molecules, Pearson correlations are {ref['engine_H']['pearson_r']:.3f} for hotspot fraction, {ref['engine_O']['pearson_r']:.3f} for pair overlap, and {ref['engine_cascade']['pearson_r']:.3f} for cascade score. Mean absolute differences are {ref['engine_H']['mean_abs_difference']:.3f}, {ref['engine_O']['mean_abs_difference']:.3f}, and {ref['engine_cascade']['mean_abs_difference']:.3f} percentage points, respectively.

The independent implementation is not numerically identical for every molecule; complete residuals are provided. Emitted engine scores satisfy the documented weighted blend to numerical precision. Serial/parallel native preparation and scoring were separately checked on 51 freshly prepared states from five ligands, preserving prepared SDF content and numerical scores/mappings. The revision package includes inputs, configurations, commands, source hashes, per-molecule outputs, and completion manifests.''',
        r'The independent-reference Results subsection, Methods, reference-implementation README, and machine-readable verification exports document agreement, residual differences, and execution provenance.')
    replies['2.6a']=(r'''The workflow diagram explains the sequence in plain language and now shows retained property limits, non-excluding alert annotations, and the capped shortlist minimum. The glossary defines the scoring terms. Figure captions identify the population and quantity plotted, and the graphical abstract has been refactored around the settings actually tested and the measured retention comparison. Figures use 600-dpi raster exports and vector versions where applicable, with colorblind-safe colors and redundant markers, hatching, or labels.''',
        r'Scheme~1, Table~1, the graphical abstract, figure captions, and the abstract describe the revised workflow and measured outcomes.')
    replies['2.6b']=(r'''The ranking comparison now uses a common support across every row, assigning unranked molecules tied ranks below the ranked set. We separately report Kendall correlation over the canonical top 1,000, using the same missing-rank rule, and the shared-set size. These computations use the complete freshly scored tables.

We also corrected the experiment description: these are table-level ranking perturbations, not independent end-to-end runs with a gate removed. The Stage-3 top-5,000 subset is restricted to molecules actually native-scored within that original rank threshold; scores for unobserved molecules are not inferred. Native-fit RMSD ordering is identified as its own ordering criterion rather than mislabeled a secondary tie-break.''',
        r'Table~7 and the ranking-perturbation Results subsection state the support, missing-rank treatment, observed-subset limit, and what the analyses can establish.')
    replies['2.6c']=(rf'''The reference-structure sensitivity run was repeated with the independent 7KI0 native reference under the revised alert settings. Full-cascade ROC-AUC changes from {metric('glp1r_full','full_cascade','roc_auc'):.3f} to {metric('glp1r_7ki0','full_cascade','roc_auc'):.3f}; native-only ROC-AUC changes from {metric('glp1r_full','native_only','roc_auc'):.3f} to {metric('glp1r_7ki0','native_only','roc_auc'):.3f}. Full-cascade BEDROC changes from {metric('glp1r_full','full_cascade','bedroc'):.3f} to {metric('glp1r_7ki0','full_cascade','bedroc'):.3f}, and expected top-10 active recovery under score ties changes from {10*metric('glp1r_full','full_cascade','top10_recovery'):.2f} to {10*metric('glp1r_7ki0','full_cascade','top10_recovery'):.2f} of 10.

We report early retrieval alongside global discrimination. This two-reference comparison does not establish robustness to the full conformational ensemble, which remains a limitation.''',
        r'The reference-structure Results subsection, Methods, and Discussion report the fresh 7KI0 comparison and its limited coverage of receptor conformations.')
    baseline=ROOT/'tmp/alert_disabled_revision/before/ACS_Omega_resubmission/response_to_reviewers.tex'
    text=baseline.read_text(encoding='utf-8')
    quotes=re.findall(r'\\begin\{reviewerquote\}.*?\\end\{reviewerquote\}',text,re.S)
    for key,(answer,changes) in replies.items():
        start=text.index(r'\comment{Comment '+key)
        response=text.index(r'\response',start)
        ends=[m.start() for m in re.finditer(r'\\(?:comment\{|section\*\{)',text[response:])]
        end=response+min(ends) if ends else text.index(r'\end{document}',response)
        text=text[:response]+r'\response '+answer+'\n\n'+r'\changes '+changes+'\n\n'+text[end:]
    start=text.index('Dear Editor,');end=text.index(r'\section*{Reviewer 1}',start)
    intro=r'''Dear Editor,

Thank you for the detailed reviews. We have revised the manuscript to distinguish terminal native peptide-contact scoring from the selection effects of the staged cascade, and have recomputed the experiments under a consistent structural-alert policy.

\section*{Summary of principal changes}
\begin{itemize}
\item PAINS/reactive alerts are recorded without exclusion; molecular-property limits and standardization remain active. All headline, retrospective, curation, reference-structure, decoy-replicate, and production experiments were rerun from molecular inputs.
\item The capped 1,000-molecule minimum on the 5\% shortlist was tested through fresh paired production runs and a complete million-compound screen. Counts distinguish shortlist admission, successful native scoring, and final ranking.
\item The abstract, quantitative Results, Discussion, tables, figures, graphical abstract, and Supporting Information now use these results. End-to-end runtime includes the native branch once, and unmeasured speedup claims have been removed.
\item An evaluation audit identified ligand-identifier bias when equal scores were ordered by ID. All primary enrichment estimates, intervals, paired tests, and replicate summaries now preserve score ties; EF, BEDROC, and top-$k$ recovery average over within-tie orderings (Table~S15). EF also uses the actual rounded cutoff fraction in its normalization. Production survival counts are unaffected.
\item Equivalence testing, docking of the new top 10, ranking perturbations, independent verification, and per-active attrition were recomputed. Limitations of the active sets, decoy source, receptor sampling, and computational validation remain explicit.
\end{itemize}

'''
    text=text[:start]+intro+text[end:]
    start=text.index(r'\response',text.index(r'\section*{Reviewer 1}'));end=text.index(r'\comment{Comment 1.1}',start)
    text=text[:start]+r'\response We have revised the contribution around the distinction between terminal scoring and upstream selection, and tested the selection changes directly. The responses below identify the newly executed comparisons and their limits.'+'\n\n'+text[end:]
    start=text.index(r'\section*{Editorial and style points}');end=text.index(r'\end{document}',start)
    text=text[:start]+r'''\section*{Editorial and presentation changes}
All quantitative passages and graphical summaries have been refreshed from the new outputs. Figures are exported at 600 dpi with colorblind-safe palettes and vector counterparts where applicable. Tables use consistent numerical precision, clear headings, and repeated headers for multipage material. Manuscript references and the final Word/PDF layouts are checked against the regenerated figure and table order.

Sincerely,\\
Kevin Song, John Zhang, Lei Ye, and Jianyi Zhang

'''+text[end:]
    assert re.findall(r'\\begin\{reviewerquote\}.*?\\end\{reviewerquote\}',text,re.S)==quotes,'Reviewer quotations changed'
    (SUB/'response_to_reviewers.tex').write_text(text,encoding='utf-8')
    print('Updated all 11 detailed replies; reviewer quotations preserved verbatim')


if __name__=='__main__':main()
