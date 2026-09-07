"""Rebuild quantitative manuscript passages from the completed revision data."""
import json
import re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUB=ROOT/'ACS_Omega_resubmission'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def paragraph(text):
    return '\\new{'+text.strip()+'}\n\n'


def heading(title):
    return '\\subsection{\\new{'+title+'}}\n'


def figure(file,label,title,caption):
    return ('\\begin{figure}[tbp]\n\\centering\n'
        '\\includegraphics[width=\\textwidth]{'+file+'.pdf}\n'
        '\\caption{\\new{\\textbf{'+title+'} '+caption+'}}\n'
        '\\label{'+label+'}\n\\end{figure}\n\n')


def table(name):
    return '\\input{'+name+'.tex}\n\n'


def main():
    data=json.loads((META/'manuscript_results.json').read_text())
    head=data['head'];counts=head['counts'];native=head['native_rerank']['counts']
    bench=data['benchmarks'];diag=data['diagnostics'];ref=data['reference_verification']
    policies=pd.DataFrame(data['floor']['benchmarks']).set_index(['system','policy'])
    eq=pd.read_csv(ROOT/'evidence/outputs/equivalence/equivalence_tost.csv')
    eq=eq[eq.metric=='roc_auc']
    docking=json.loads((ROOT/'evidence/outputs/equivalence/docking_wilcoxon.json').read_text())
    dock=pd.read_csv(ROOT/'evidence/outputs/revision/docking_top10_summary.csv')
    ablation=pd.read_csv(ROOT/'evidence/outputs/ablation_common/ablation_common_support.csv')
    reps=pd.DataFrame(data['replicates'])
    def metric(system,method,key):return float(bench[system][method][key])
    def policy(system,arm,col):return int(policies.loc[(system,arm),col])
    def retention(system,arm):
        return f"{policy(system,arm,'final_ranked_actives')}/{policy(system,arm,'n_actives')}"
    hours=head['timings']['total_pipeline_sec']/3600
    keys=['glp1r_full','ghsr_full','ntsr1_full','mdm2_full']
    gains=[metric(s,'native_only','roc_auc')-metric(s,'standard_3d_pharmacophore','roc_auc') for s in keys]
    outperform=sum(x>0 for x in gains)
    native_ge=sum(metric(s,'native_only','roc_auc')>=metric(s,'full_cascade','roc_auc')-1e-12 for s in keys)
    ga,gb=retention('ghsr','percentage'),retention('ghsr','floor1000')
    na,nb=retention('ntsr1','percentage'),retention('ntsr1','floor1000')
    for system in ['glp1r','mdm2']:
        assert policy(system,'percentage','shortlist_n')==policy(system,'floor1000','shortlist_n')>=1000
    abstract=paragraph(rf"""Peptide-contact pharmacophores offer a way to rank small molecules against
distributed peptide-binding interfaces, but upstream filters can discard known actives.
We recompute a staged screening workflow across GLP-1R, GHSR, NTSR1, and MDM2--p53,
retaining molecular-property limits while recording PAINS/reactive alerts without
exclusion. Tie-aware native-only scoring exceeds a conventional single-pass 3D pharmacophore
on {outperform}/4 systems, with ROC-AUC differences from {min(gains):.3f} to {max(gains):.3f},
and matches or exceeds the full cascade on {native_ge}/4 systems. A fresh
million-compound GLP-1R run shortlists {counts['shortlist_size']:,} molecules and produces
{native['successful_best_ligands']:,} successful native ligand scores in {hours:.2f}~h.
The revised shortlist keeps the larger of 5\% or 1,000 molecules, capped by eligible
candidates; the floor is nonbinding at million-compound scale. In paired production
benchmarks, it changes final active retention from {ga} to {gb} for GHSR and from
{na} to {nb} for NTSR1. Across all four active sets, non-excluding alerts allow
158/160 molecules through Stage 0; two still fail the retained property limits.
The floor addresses an additional shortlist loss but cannot reverse upstream
hotspot exclusions or downstream selection failures. Complete stage-survival
records, repeated decoy benchmarks, reference-structure and curation checks, and
an independent implementation of the prescreening stages accompany the results. These computational
rankings prioritize candidates for testing and do not establish biological activity.""")
    results=heading('Million-scale screening with non-excluding structural alerts')
    results+=paragraph(rf"""The complete GLP-1R screen was rerun from the original H17--H20 ZINC
inputs, without reusing molecular scores or prepared conformers.\cite{{Irwin2020}}
Of 1,000,000 inputs, {counts['property_pass']:,} pass Stage 0 and
{counts['hotspot_pass']:,} pass the hotspot gate. The 5\% rule, with a minimum of
1,000 eligible molecules, gives a {counts['shortlist_size']:,}-molecule shortlist.
Because the percentage term exceeds the minimum, the floor does not change
shortlist size at this scale. Disabling structural-alert exclusions does change
the screened cohort, so comparison with the earlier strict-filter run is not
interpreted as a floor-only effect (Table~S12).""")
    results+=paragraph(rf"""Fresh Stage-3 reranking yields {counts['final_hits']:,} successful ligand scores.
The diversified native branch assembles {native['candidate_pool_size']:,} candidates,
selects {native['selected_ligands']:,} with an eight-per-scaffold cap,\cite{{Bemis1996}}
prepares {native['prepared_microstates']:,} microstates, and successfully scores
{native['scored_microstates']:,} states from {native['successful_best_ligands']:,} ligands.
The final output contains {native['final_rows']:,} ranked candidates. Measured
end-to-end wall time is {head['timings']['total_pipeline_sec']:,.1f}~s
({hours:.2f}~h), using 12 workers for prescreening, Stage 3, and the native branch
on an eight-core, 16-thread workstation. This total includes native preparation
and scoring once. No full-library native-only runtime was measured, so we do not
assert an extrapolated speedup (Table~S1).""")
    results+=heading('Retrospective enrichment and the native-only comparison')
    results+=paragraph(rf"""The four methods rank the same labeled universe within each target:
10 in-domain actives and 300 matched decoys for GLP-1R, and 50 actives with
1,500 decoys for each other system. These uncapped comparisons measure scoring
and upstream-gate effects; the production shortlist is tested separately.
For GLP-1R, full-cascade ROC-AUC is {metric('glp1r_full','full_cascade','roc_auc'):.3f},
compared with {metric('glp1r_full','native_only','roc_auc'):.3f} for native-only,
{metric('glp1r_full','stage3_only','roc_auc'):.3f} for Stage-3-only, and
{metric('glp1r_full','standard_3d_pharmacophore','roc_auc'):.3f} for single-pass 3D
pharmacophore scoring. The full cascade has expected recovery of
{10*metric('glp1r_full','full_cascade','top10_recovery'):.2f} of 10 actives in the top 10
and has BEDROC {metric('glp1r_full','full_cascade','bedroc'):.3f}
(\figref{{fig:benchmark}}; \tabref{{tab:crosssystem}}).""")
    results+=paragraph(r'''All enrichment statistics preserve score ties. An evaluation audit found that the previous ligand-ID tie break could favor ChEMBL actives over differently named decoys, including within unscored groups. We therefore recomputed point estimates, grouped-bootstrap intervals, paired comparisons, and replicate summaries without using identifiers to resolve equal scores. Table~S15 quantifies the resulting differences. Expected top-$k$ recovery averages over possible within-tie orders and can be fractional; production-policy survival counts remain observed integers.''')
    results+=figure('fig3_glp1r_benchmark','fig:benchmark','Retrospective GLP-1R benchmark.',
        r'Four methods evaluated on the same 310 labeled molecules under the revised alert settings. '
        r'Bars show point estimates and intervals show 95\% matched-group bootstrap intervals. '
        r'Color and hatching distinguish methods; larger values indicate better retrieval.')
    results+=table('crosssystem_table')
    results+=heading('Equivalence testing across systems')
    for row in eq.itertuples():
        interpretation='meets the equivalence criterion' if row.equivalent else 'does not establish equivalence'
        results+=paragraph(rf"""For {row.system}, the observed full-cascade minus native-only
ROC-AUC difference is {row.observed_delta:+.3f}, with a 90\% grouped-bootstrap
interval of [{row.ci90_low:+.3f}, {row.ci90_high:+.3f}]. This {interpretation}
at the prespecified margin of $\pm0.05$ ROC-AUC.""")
    results+=paragraph(r"""Equivalence requires the whole 90\% interval to fall inside that margin;
a non-significant difference alone is insufficient. The small GLP-1R active set
limits precision, and observed differences and interval widths are reported
together rather than collapsed into a claim of parity (\tabref{tab:tost}).""")
    results+=table('tost_table')
    results+=heading('Enrichment across all four systems')
    results+=paragraph(rf"""Native-only scoring exceeds single-pass 3D on {outperform} of four systems,
with ROC-AUC differences ranging from {min(gains):.3f} to {max(gains):.3f}.
It matches or exceeds the full cascade on {native_ge} systems. These comparisons
separate terminal peptide-contact scoring from the upstream cascade and are
recomputed under non-excluding alerts; conclusions from the earlier strict-filter
configuration are not carried over (\figref{{fig:crosssystem}}).""")
    results+=paragraph(rf"""At MDM2--p53, full-cascade and native-only ROC-AUC are
{metric('mdm2_full','full_cascade','roc_auc'):.3f} and
{metric('mdm2_full','native_only','roc_auc'):.3f}, respectively; their BEDROC values
are {metric('mdm2_full','full_cascade','bedroc'):.3f} and
{metric('mdm2_full','native_only','bedroc'):.3f}. This single protein-protein
interface is insufficient to establish a general rule about flat interfaces.
Its measured result is reported alongside the three GPCR systems rather than
used to infer a geometry-dependent failure from the old alert-filter setting.""")
    results+=figure('fig4_cross_system','fig:crosssystem','Enrichment across all four systems.',
        r'Full cascade, native-only, and single-pass 3D evaluated on matched labeled universes. '
        r'Panels show ROC-AUC, EF1\%, and BEDROC; intervals are 95\% matched-group bootstrap intervals. '
        r'All configurations record structural alerts without exclusion.')
    results+=heading('Where actives are lost with structural-alert exclusions disabled')
    results+=paragraph(r"""An independent Stage-0 implementation agrees with the engine on all
4,960 labeled benchmark molecules. Among 160 actives, 158 pass Stage 0:
10/10 for GLP-1R, 50/50 for GHSR, 48/50 for NTSR1, and 50/50 for MDM2--p53.
The 34 actives carrying PAINS/reactive flags are annotated and admitted at
this stage. Two NTSR1 actives remain outside the property envelope after
standardization. Thus structural-alert matching is no longer a cause of
active loss; it remains visible in the released records.""")
    results+=paragraph(r"""For the revised production setting, \tabref{tab:attrition} and
\figref{fig:attrition} identify the first limiting stage for each active,
including the shortlist, 3D scoring, native selection, native preparation/scoring,
and the final ranking limit. Counts describe measured final survival, not
shortlist admission inferred to imply native scoring. Table~S10 identifies
individual losses, and Table~S11 reports stage-by-stage surviving actives.""")
    results+=table('attrition_table')
    results+=figure('fig5_attrition','fig:attrition','Fate of every in-domain active under revised production settings.',
        r'The 5\% shortlist has a minimum of 1,000 eligible molecules and structural alerts do not exclude. '
        r'Segments count the first stage of loss; retained actives are present in the final native ranking. '
        r'Hatching and labels provide information independently of hue.')
    results+=heading('The trade between efficiency and active retention')
    results+=paragraph(rf"""Fresh matched policy pairs isolate the shortlist minimum while holding the
revised alert settings and all subsequent scoring and selection rules fixed.
For GHSR, the floor expands the shortlist from
{policy('ghsr','percentage','shortlist_n'):,} to {policy('ghsr','floor1000','shortlist_n'):,}
and changes final active retention from {ga} to {gb}. For NTSR1, the corresponding
shortlists contain {policy('ntsr1','percentage','shortlist_n'):,} and
{policy('ntsr1','floor1000','shortlist_n'):,} molecules, with final retention changing
from {na} to {nb}. GLP-1R and MDM2--p53 have identical policy shortlists because
their 5\% counts already exceed 1,000 (\figref{{fig:efficiency}}; \tabref{{tab:efficiency}}).""")
    results+=paragraph(r"""The rule is a minimum on a percentage shortlist, capped by the number
available; it is not a fixed 1,000-molecule maximum. It avoids an additional severe
cut when the upstream pool is small but cannot restore hotspot-gate losses or
guarantee survival through downstream preparation and scaffold constraints.
Figure~S5 and Table~S13 give the complete stage counts. Table~S14 reports measured
native-branch wall times, including selection, preparation, scoring, and outputs;
we do not characterize the additional work as negligible.""")
    results+=figure('fig6_efficiency','fig:efficiency','The trade between efficiency and active retention under a measured floor.',
        r'(a) Actual shortlist sizes. (b) Actives present in the final native ranking. '
        r'Orange circles denote 5\% only and blue squares denote 5\% plus the capped 1,000-molecule minimum. '
        r'Coincident values are offset for visibility. Native execution is shared only for identical policy inputs.')
    results+=table('efficiency_table')
    results+=heading('Robustness to decoy selection and reference structure')
    replicate_summaries=[]
    for key,label in [('glp1r','GLP-1R'),('ghsr','GHSR'),('ntsr1','NTSR1'),('mdm2','MDM2--p53')]:
        values=reps[(reps.system==key)&(reps.method=='full_cascade')].roc_auc
        replicate_summaries.append(rf'{values.mean():.3f} $\pm$ {values.std(ddof=1):.3f} for {label} '
            rf'(range {values.min():.3f}--{values.max():.3f})')
    results+=paragraph(r'Across five independently sampled matched-decoy sets per 10-active benchmark, '
        r'full-cascade ROC-AUC (mean $\pm$ standard deviation) is '
        + '; '.join(replicate_summaries) + r' (Table~S8).')
    results+=paragraph(r"""These repeated draws assess sensitivity to decoy selection, while
all draws still share the same ZINC source and curation procedure. They cannot
exclude a systematic ChEMBL-versus-ZINC difference. The expanded 50-active
comparisons and the 10-active decoy replicates are therefore reported separately.""")
    results+=paragraph(rf"""Rebuilding the GLP-1R native reference from the independent active,
peptide-bound 7KI0 complex changes full-cascade ROC-AUC from
{metric('glp1r_full','full_cascade','roc_auc'):.3f} to
{metric('glp1r_7ki0','full_cascade','roc_auc'):.3f}, native-only ROC-AUC from
{metric('glp1r_full','native_only','roc_auc'):.3f} to
{metric('glp1r_7ki0','native_only','roc_auc'):.3f}, and full-cascade BEDROC from
{metric('glp1r_full','full_cascade','bedroc'):.3f} to
{metric('glp1r_7ki0','full_cascade','bedroc'):.3f}. Top-10 active recovery changes from
{10*metric('glp1r_full','full_cascade','top10_recovery'):.2f} to
{10*metric('glp1r_7ki0','full_cascade','top10_recovery'):.2f} of 10 (expected counts under score ties).
Reference dependence is assessed on both overall discrimination and early retrieval;
neither measurement establishes robustness to all receptor conformations.""")
    results+=heading('Internal diagnostics and ranking perturbations')
    upstream=diag['upstream_correlations']['screen_weighted_coverage_pct']
    results+=paragraph(rf"""Across {diag['n']:,} successfully native-scored headline ligands,
Stage-3 and native coverage have Pearson $r = {upstream['pearson_r']:.3f}$ and
Spearman $\rho = {upstream['spearman_rho']:.3f}$. Ordering the same cohort by
Stage-3 rank and final native rank gives $\rho = {diag['rank_spearman']:.3f}$;
their top-10 sets share {diag['top10_overlap']} molecules. Among the final top 1,000,
the median absolute within-cohort rank shift is
{diag['top1000_median_absolute_cohort_rank_shift']:,.0f} positions.
These comparisons use a common native-scored universe, rather than subtracting
ranks defined over different populations (\tabref{{tab:evidence}}; Figure~S1).""")
    results+=table('evidence_table')
    results+=paragraph(rf"""Ranking perturbations are recomputed from the fresh score tables
(\tabref{{tab:ablation}}). Correlations are reported over the canonical top 1,000
and over the common support of {int(ablation.n_common.iloc[0]):,} ligands, assigning
unranked molecules tied ranks below ranked ones. These are table-level ranking
analyses, not independent end-to-end executions of an altered screening pipeline.
Consequently they assess sensitivity of the reported ordering without claiming
to recover molecules excluded by an upstream gate.""")
    results+=table('ablation_table')
    results+=heading('Sensitivity to manual pharmacophore curation')
    results+=paragraph(rf"""The automated GLP-1R pharmacophore removes hand-assigned contact groups,
weight bonuses, native-supported expansion, and curated rescue. Its freshly
recomputed full-cascade ROC-AUC is
{metric('glp1r_automated','full_cascade','roc_auc'):.3f}, compared with
{metric('glp1r_full','full_cascade','roc_auc'):.3f} for the curated model;
BEDROC is {metric('glp1r_automated','full_cascade','bedroc'):.3f} versus
{metric('glp1r_full','full_cascade','bedroc'):.3f}. The automated and curated models
have expected recovery of {10*metric('glp1r_automated','full_cascade','top10_recovery'):.2f} and
{10*metric('glp1r_full','full_cascade','top10_recovery'):.2f} actives in their top 10,
respectively. Global discrimination and early retrieval are both reported;
one does not substitute for the other (\tabref{{tab:ablation_curation}}).""")
    results+=table('ablation_curation_table')
    results+=heading('Properties of the new native-ranked cohort')
    top=pd.DataFrame(data['top10']);first=top.iloc[0]
    results+=paragraph(rf"""The top-ranked ligand, {first.zinc_id}, has native weighted coverage
{first.native_weighted_coverage_pct:.2f}\%, matches
{int(first.native_matched_reference_features)} retained reference features,
and has native-fit RMSD {first.native_fit_rmsd_angstrom:.2f}~\angstrom\ and
mean pair-distance error {first.native_mean_pair_distance_error_angstrom:.2f}~\angstrom.
The top-10 cohort spans {top.native_weighted_coverage_pct.min():.2f}--{top.native_weighted_coverage_pct.max():.2f}\%
native coverage, {top.mw.min():.1f}--{top.mw.max():.1f}~Da, and LogP
{top.logp.min():.2f} to {top.logp.max():.2f}. Figure~S3 shows the leading native
overlay and Figure~S4 gives the top-20 structures. The broader native-scored
cohort has median coverage {diag['coverage_median']:.2f}\% and interquartile range
{diag['coverage_q25']:.2f}--{diag['coverage_q75']:.2f}\% (Figure~S2).
These scores quantify reference-feature recovery, not binding affinity or potency.""")
    results+=table('top10_table')
    results+=heading('Orthogonal docking of the new top-ranked ligands')
    valid=dock.dropna(subset=['best_active','best_inactive'])
    results+=paragraph(rf"""The newly ranked top 10 were redocked against the same three active-state
and two inactive-state GLP-1R structures using three seeds per structure and
exhaustiveness 16.\cite{{Trott2010,Eberhardt2021}} Valid results against both
states were obtained for {docking['n']}/10 ligands. Their median best active-state
affinity is {valid.best_active.median():.2f}~kcal/mol. Among the paired results,
{docking['n_favoring_active']}/{docking['n']} favor the active state, with median
active-minus-inactive difference {docking['median_delta']:.2f}~kcal/mol
(interquartile range {docking['iqr_low']:.2f} to {docking['iqr_high']:.2f}).
The two-sided exact signed-rank permutation test gives
$W = {docking['wilcoxon_W']:g}$ and $p = {docking['p_two_sided']:.4f}$
(\figref{{fig:docking}}; Table~S9).""")
    results+=paragraph(r"""This is an orthogonal computational scoring comparison. The ligand
population is selected by the fresh pharmacophore ranking, and preparation and
docking failures are retained in the output accounting. Docking scores do not
establish receptor engagement, state selectivity, or functional agonism;
those claims require experimental measurements. The active-state ensemble has
three structures and the inactive-state ensemble two, so taking the best score
per state gives unequal sampling opportunities. The paired test describes these
selected ensembles and does not isolate a receptor-state effect.""")
    results+=figure('fig7_docking','fig:docking','Redocking of the top-10 final-ranked ligands.',
        r'(a) Best affinity against active- and inactive-state structures. '
        r'(b) Paired active-minus-inactive differences, their median, and the exact signed-rank test. '
        r'Only ligands with successful results against both states contribute to the paired analysis.')
    results+=heading('Independent reference implementation')
    c=ref['continuous']
    results+=paragraph(rf"""The independent implementation reproduces Stage-0 and Stage-1 decisions
for all 310 GLP-1R benchmark molecules under the revised configuration. On a
random sample of {c['n']} freshly shortlisted headline molecules, Pearson
correlations with the engine are {c['engine_H']['pearson_r']:.3f} for hotspot
fraction, {c['engine_O']['pearson_r']:.3f} for pair overlap, and
{c['engine_cascade']['pearson_r']:.3f} for cascade score; corresponding mean
absolute differences are {c['engine_H']['mean_abs_difference']:.3f},
{c['engine_O']['mean_abs_difference']:.3f}, and
{c['engine_cascade']['mean_abs_difference']:.3f} percentage points.
The engine's emitted cascade score satisfies its stated blend to numerical
precision. The independent implementation is not numerically identical for every
molecule; full comparison records, including the residual differences, are
released. It implements Stages 0--2; conformer generation and native alignment
are checked through fresh execution and serial/parallel validation.""")
    results+='\\FloatBarrier\n\n'
    body=(SUB/'manuscript_body.tex').read_text(encoding='utf-8')
    body=re.sub(r'(?<=\\begin\{abstract\}\n).*?(?=\n\\end\{abstract\})',lambda m:abstract,body,count=1,flags=re.S)
    starts=[body.find(r'\subsection{Million-scale tractability}'),
            body.find(r'\subsection{\new{Million-scale screening with non-excluding structural alerts}}')]
    start=min(x for x in starts if x>=0)
    end=body.index(r'\section{Methods}',start)
    body=body[:start]+results+body[end:]
    # Discussion is rewritten after the measured comparisons, not inherited
    # from the strict-filter configuration.
    discussion=paragraph(rf"""The revised experiments separate scoring performance from two
configuration choices: excluding structural-alert matches and imposing a fractional
shortlist after an upstream gate. With alerts recorded rather than excluded,
native-only scoring exceeds single-pass 3D on {outperform}/4 systems and matches
or exceeds the full cascade on {native_ge}/4. The native reference and the cascade
must therefore be judged separately: a useful terminal score does not establish
that every upstream filter improves retrieval.""")
    discussion+=paragraph(rf"""The Stage-0 change admits the 34 alert-flagged actives previously
excluded by that configuration while retaining the property envelope. The floor
then addresses a different loss. In matched runs with identical non-excluding
alerts, final retention changes from {ga} to {gb} for GHSR and from {na} to {nb}
for NTSR1. At million-compound scale the percentage term already exceeds 1,000,
so the floor is inactive. The measured screen produces
{native['successful_best_ligands']:,} successful native ligand scores in
{hours:.2f}~h. This is a measured operating point, not a benchmark against
exhaustive native scoring of the entire library.""")
    discussion+=paragraph(r"""A practical configuration should preserve molecular validity and
the intended property limits, record structural alerts transparently, and avoid
reducing a small eligible pool to a handful of molecules solely because a fixed
percentage was inherited from a large-library workflow. The minimum remains
bounded by availability and does not guarantee final survival. Native-pool quotas,
scaffold caps, failed preparation, and final rank limits must remain visible in
the accounting. Retaining an alert-flagged molecule is a screening decision;
the annotation still informs subsequent assay design and compound evaluation.""")
    discussion+=paragraph(r"""Several limitations remain. The GLP-1R set contains only 10 in-domain
actives; the other systems contain 50 each. Actives and decoys originate from
different databases, so property matching and repeated decoy draws cannot
eliminate source-related bias. Only one non-GPCR interface is included, and the
curation and reference-structure checks cover a limited set of alternatives.
The conventional comparator is a single-pass 3D pharmacophore; the top-10 docking
exercise does not replace a full retrospective comparison against a docking
pipeline at matched computational cost. Neither retrospective enrichment,
pharmacophore coverage, nor favorable docking demonstrates prospective biological
activity. Binding and cellular signaling measurements remain necessary.""")
    discussion+=paragraph(r"""The contribution is therefore a reproducible scoring and selection
analysis with a tested shortlist intervention, supported by complete stage-survival
records. The revised data package allows readers to distinguish a molecule
admitted to a shortlist, one successfully native-scored, and one present in the
final ranking, and to recompute the reported enrichment and uncertainty estimates.""")
    start=body.index(r'\section{Discussion}');end=body.index(r'\section*{Associated Content}',start)
    body=body[:start]+'\\section{Discussion}\n'+discussion+body[end:]
    body=body.replace('free of charge at \\url{https://pubs.acs.org/doi/XXXXXXX}.','as part of the accompanying revision package.')
    body=body.replace('(Figures S1--S4, Tables S1--S11)','(Figures S1--S5, Tables S1--S15)')
    body=body.replace('Materials accompanying this study are publicly available at\n\\url{https://github.com/kevinmsong/PharmacophoreCascade}:',
        'The revision data package contains the following materials; the project repository is\n\\url{https://github.com/kevinmsong/PharmacophoreCascade}:')
    body=body.replace('Peptide-mediated receptor activation and protein-protein recognition are important\nbut difficult targets for small-molecule discovery.',
        'Peptide-mediated receptor activation and protein-protein recognition are important,\nbut difficult targets, for small-molecule discovery.')
    intro_start=body.find('Four contributions follow:')
    if intro_start<0:
        intro_start=body.index('We compare the full cascade with a native-only baseline')
        last='independent implementation make the scoring and selection decisions reviewable.'
    else:last='per-molecule benchmark data for all four systems.'
    intro_end=body.index(last,intro_start)+len(last)
    body=body[:intro_start]+('We compare the full cascade with a native-only baseline and a conventional 3D\n'
        'pharmacophore, test equivalence with explicit margins, and repeat the scoring\n'
        'under non-excluding structural alerts. We then evaluate a minimum on the\n'
        'percentage shortlist through fresh paired production runs and a complete\n'
        'million-compound rerun. Per-molecule records, sensitivity analyses, and an\n'
        'independent implementation make the scoring and selection decisions reviewable.')+body[intro_end:]
    body=body.replace('and \\tabref{tab:ablation} confirms\nthat pool composition is the most consequential choice tested.',
        'and the observed-subset analysis in \\tabref{tab:ablation} characterizes ranking sensitivity without inferring scores for unobserved ligands.')
    start=body.index(r'\new{Component ablations replace one pipeline element')
    end=body.index(r'\new{To assess sensitivity to the reference structure',start)
    body=body[:start]+paragraph(r'''Table-level ranking perturbations reorder the fresh score tables using hotspot score, pair overlap, cascade score within the native cohort, an observed Stage-3 top-5,000 subset, or native-fit RMSD. They do not rerun an altered gate or infer native scores for excluded molecules. Kendall correlations use the canonical top 1,000 and a fixed common support, with unranked molecules tied below ranked molecules. A six-configuration weight sweep similarly reorders the observed Stage-3 table. The separate uncurated benchmark is a full fresh run: its receptor pharmacophore is rebuilt from interface atoms with automated feature priority, $k$-means contact groups, no weight bonuses, no curated rescue, and the GLP-1R required-group condition disabled. Gate thresholds, conformer counts, tolerances, labeled inputs, and seeds remain fixed.''')+body[end:]
    body=body.replace('with a paired Wilcoxon signed-rank test.',
        'with an exact signed-rank permutation test over successful paired results. Absolute differences receive average ranks for ties; zero differences are excluded. All sign flips are enumerated for the two-sided test, and all preparation/docking failures remain in the job accounting.')
    body=body.replace('No conclusion depends on it: the\nconclusions rest on the retrospective benchmarks, recomputable from the released\nper-molecule data,',
        'The retrospective enrichment statistics are recomputable from the accompanying per-molecule data,')
    body=body.replace('which was not previously possible.',
        'while exposing the measured residual differences.')
    body=body.replace('active and scaffold recovery;', 'expected active recovery;')
    body=body.replace('The final native ranking is independent of the Stage-3 order.',
        'The final ranking prioritizes native coverage and fit, with Stage-3 rank retained as a later tie-break.')
    body=body.replace(r'\subsection{\new{Statistical analysis}}',r'\subsection{\new{Statistical analysis}}'+'\n'+paragraph(r'''For retrospective evaluation, method status determines the ordering of successful and unsuccessful outcomes, followed by the emitted score. Equal status and equal score, including missing scores within a failure status, form a tie group. ROC-AUC assigns half credit to tied active--decoy pairs, and PR-AUC is average precision evaluated at tied score thresholds. EF, BEDROC, and top-$k$ recovery are their expected values over uniformly ordered molecules within each tie group. For EF at a nominal fraction $f$, the actual cutoff is $m=\lceil fN\rceil$ and EF equals the expected active fraction among those $m$ molecules divided by the library active prevalence; the denominator is $m/N$, not an unrounded nominal fraction. Grouped bootstrap samples retain these ties and jointly resample each active with its matched decoys (5,000 draws, seed 42). Original identifier-ordered lists remain available for traceability, but their tie breaks do not enter the revised enrichment statistics. Observed production selections and final survival are evaluated separately.'''),1)
    (SUB/'manuscript_body.tex').write_text(body,encoding='utf-8')
    print('Rebuilt abstract, quantitative Results, Discussion, and current availability statements')


if __name__=='__main__':main()
