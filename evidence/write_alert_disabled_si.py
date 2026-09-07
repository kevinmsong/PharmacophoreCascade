"""Refresh SI captions, run accounting, ranking tables, and experiment descriptions."""
import json
import re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUB=ROOT/'ACS_Omega_resubmission'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def longtable(caption,label,columns,rows,alignment):
    header=' & '.join(columns)+r' \\'
    return '\n'.join([r'{\footnotesize',r'\setlength{\tabcolsep}{4pt}',
        r'\begin{longtable}{'+alignment+'}',r'\caption{'+caption+r'}\label{'+label+r'}\\',
        r'\toprule',header,r'\midrule',r'\endfirsthead',r'\toprule',header,r'\midrule',r'\endhead']+
        [' & '.join(map(str,row))+r' \\' for row in rows]+[r'\bottomrule',r'\end{longtable}','}',''])


def main():
    data=json.loads((META/'manuscript_results.json').read_text())
    h=data['head'];c=h['counts'];n=h['native_rerank']['counts'];t=h['timings'];nt=h['native_rerank']['timings']
    d=data['diagnostics'];top=pd.DataFrame(data['top10']);first=top.iloc[0]
    source=SUB/'supporting_information.tex'
    text=source.read_text(encoding='utf-8')
    text=text.replace('Figure S1 through Figure S4','Figure S1 through Figure S5').replace('Table S1 through Table S11','Table S1 through Table S15')
    text=text.replace(r'\item Figure S4. Structures of the top 20 native-ranked ligands.',
        r'\item Figure S4. Structures of the top 20 native-ranked ligands.'+'\n'+r'\item Figure S5. Active survival in the paired shortlist-policy runs.')
    text=text.replace('Every in-domain active removed before native scoring, with the stage and specific cause.',
        'Every in-domain active absent from the final ranking, with the first limiting stage and cause.')
    text=text.replace(r'\item Table S11. Stage-by-stage active survival under production constraints.',
        r'\item Table S11. Stage-by-stage active survival under production constraints.'+'\n'+
        r'\item Table S12. Fresh headline run and the archived strict-filter configuration.'+'\n'+
        r'\item Table S13. Complete molecule and active counts for paired shortlist policies.'+'\n'+
        r'\item Table S14. Measured native-branch wall times for paired shortlist policies.'+'\n'+
        r'\item Table S15. Identifier-ordered versus tie-aware benchmark evaluation.')
    captions={
        'fig:s1_upstream':rf'\textbf{{Upstream screen metrics against terminal native coverage}} across {d["n"]:,} successfully native-scored headline ligands. Hexbins show point density, dashed lines are least-squares fits, and each panel reports Pearson and Spearman correlations computed from the displayed cohort. These are internal associations among scored candidates and do not measure enrichment relative to molecules excluded upstream.',
        'fig:s2_coverage':rf'\textbf{{Distribution of native weighted coverage}} across {d["n"]:,} successfully scored ligands. (a) Histogram: median {d["coverage_median"]:.2f}\%, interquartile range {d["coverage_q25"]:.2f}--{d["coverage_q75"]:.2f}\%, and full range {d["coverage_min"]:.2f}--{d["coverage_max"]:.2f}\%. (b) Coverage in final-rank order, with the reported top 1,000 indicated.',
        'fig:s3_overlay':rf'\textbf{{Best-scoring native overlay for the top-ranked ligand}}, {first.zinc_id}, against the retained GLP-1 peptide-contact reference from PDB 6X18. Filled markers denote matched reference features and open markers denote unmatched features. Shape and color indicate family; marker area indicates weight. Its {int(first.native_matched_reference_features)} matched features give {first.native_weighted_coverage_pct:.2f}\% coverage, native-fit RMSD {first.native_fit_rmsd_angstrom:.2f}~\angstrom, and mean pair-distance error {first.native_mean_pair_distance_error_angstrom:.2f}~\angstrom.',
        'fig:s4_top20':r'\textbf{Structures of the top 20 ligands ranked by native peptide-interface mimicry}, in final-rank order from the fresh million-compound run. Each entry identifies the ligand and its native weighted coverage. Structural-alert matches remain annotations; the displayed ranks do not establish potency or assay suitability.'}
    captions['fig:s3_overlay']+=r' Gray bonds show the actual aligned heavy-atom ligand structure; black crosses mark its matched feature centroids, with dotted segments to the corresponding reference features. Coordinates are projected onto the first two principal axes of the reference. The reported fit metrics are three-dimensional; the aligned SDF independently reproduces the RMSD within its coordinate-rounding precision.'
    for label,caption in captions.items():
        pattern=r'\\caption\{.*?\}\s*(?=\\label\{'+re.escape(label)+r'\})'
        # Restrict the search to the figure carrying this label to avoid crossing floats.
        end=text.index(r'\label{'+label+'}');start=text.rfind(r'\begin{figure}',0,end)
        block=text[start:end]
        block=re.sub(r'\\caption\{.*\}\s*$',lambda _:r'\caption{'+caption+'}\n',block,flags=re.S)
        text=text[:start]+block+text[end:]
    marker=r'\section*{Supplementary Tables}'
    s5=(r'\begin{figure}[H]'+'\n'+r'\centering'+'\n'+
        r'\includegraphics[width=\textwidth]{figA5_floor_survival.pdf}'+'\n'+
        r'\caption{\textbf{Active survival across the measured shortlist-policy pairs.} Both policies record structural alerts without exclusion and retain property limits. Lines trace actual survival through gating, shortlisting, 3D scoring, native selection, successful native scoring, and the final ranking. The 1,000-molecule floor is capped by available candidates. The two curves coincide for GLP-1R and MDM2--p53.}'+'\n'+
        r'\label{fig:s5_floor}'+'\n'+r'\end{figure}'+'\n'+r'\clearpage'+'\n\n')
    text=text.replace(marker,s5+marker,1)
    rows=[['Starting library','1,000,000 ZINC molecules; H17--H20'],['Workers','12 prescreen / 12 Stage 3 / 12 native'],
        ['Stage-1 / Stage-2 query features','25 / 24'],['Shortlist rule',r'$\min(N_{\rm eligible},\max(\lceil0.05N_{\rm Stage0}\rceil,1000))$'],
        ['Hotspot / pair weight','0.4 / 0.6'],['Typed-feature caps','Fixed: 2, 2, 4, 6, 6, 6'],
        ['Stage-3 conformers / query anchors','16 / 28'],['Stage-3 distance / pair tolerance',r'2.75 / 2.75~\angstrom'],
        ['Pair-hash mode',r'precision\_5bin'],['Structural alerts','Recorded without exclusion; native PAINS filter disabled'],
        ['Native pool','12,000 Stage-3; 4,000 hotspot breadth; 4,000 native-supported hotspot'],
        ['Native selection','Up to 5,000 ligands; maximum eight per Murcko scaffold'],
        ['Native preparation','Up to two charge states; eight tautomers per charge state; four conformers per state'],
        ['Final rank mode',r'native\_first']]
    for label,key in [('Stage-0 pass','property_pass'),('Stage-1 pass','hotspot_pass'),('Shortlist','shortlist_size'),('Successful Stage-3 ligands','final_hits')]:rows.append([label,f'{c[key]:,}'])
    for label,key in [('Native candidate pool','candidate_pool_size'),('Selected native ligands','selected_ligands'),('Prepared states','prepared_microstates'),('Scored states','scored_microstates'),('Successfully native-scored ligands','successful_best_ligands'),('Final ranking','final_rows')]:rows.append([label,f'{n[key]:,}'])
    for label,key in [('Stage-0 accumulated task time','stage0_property_gate_sec'),('Stage-1 accumulated task time','stage1_hotspot_scoring_sec'),('Stage-2 accumulated task time','stage2_pair_hash_scoring_sec'),('Stage-3 wall time','stage3_rerank_sec')]:rows.append([label,f'{t[key]:,.1f}~s'])
    for label,key in [('Native selection wall time','selection_sec'),('Native preparation wall time','library_prepare_sec'),('Native scoring wall time','scoring_sec'),('Native branch total wall time','total_sec')]:rows.append([label,f'{nt[key]:,.1f}~s'])
    rows.append([r'\textbf{End-to-end wall time}',rf'\textbf{{{t["total_pipeline_sec"]:,.1f}~s ({t["total_pipeline_sec"]/3600:.2f}~h)}}'])
    table=longtable(r'\textbf{Audited run configuration, counts, and timings for the GLP-1R screen.} '
        r'Stage-0--2 task times sum per-molecule elapsed '
        r'times across workers and are not additive wall-clock phases. End-to-end wall time includes '
        r'the native branch once; concurrent benchmark jobs share the workstation. No exhaustive '
        r'million-ligand native-only runtime was measured.',
        'tab:s1_runconfig',['Field','Value'],rows,r'@{}p{.36\textwidth}p{.57\textwidth}@{}')
    (SUB/'current_run_config_table.tex').write_text(table,encoding='utf-8')
    begin=text.index(r'{\footnotesize',text.index(marker));end=text.index(r'\newpage',begin)
    text=text[:begin]+r'\input{current_run_config_table.tex}'+'\n\n'+text[end:]
    native=pd.read_csv(META/'headline_native_cohort_diagnostics.csv');rows=[]
    for group in [native.nlargest(5,'cohort_rank_shift'),native.nsmallest(5,'cohort_rank_shift')]:
        for r in group.itertuples():rows.append([r.zinc_id,f'{r.screen_cohort_rank:,}',f'{r.native_cohort_rank:,}',f'{r.cohort_rank_shift:+,}',f'{r.screen_weighted_coverage_pct:.2f}',f'{r.native_weighted_coverage_pct:.2f}'])
    table=longtable(r'\textbf{Ligands showing the largest rank changes under native reranking.} '
        r'The first five rows show the largest upward shifts and the next five the largest downward shifts. '
        r'Both ranks are defined over the same successfully native-scored cohort; positive shifts mean improvement.',
        'tab:s3_rankshifts',['Ligand ID','Screen rank','Native rank','Shift',r'Screen (\%)',r'Native (\%)'],rows,'lrrrrr')
    (SUB/'current_rank_shift_table.tex').write_text(table,encoding='utf-8')
    endlabel=text.index(r'\label{tab:s3_rankshifts}')
    begin=text.rfind(r'{\footnotesize',0,endlabel);end=text.index(r'\end{longtable}',endlabel)+len(r'\end{longtable}')
    end=text.index('}',end)+1
    text=text[:begin]+r'\input{current_rank_shift_table.tex}'+text[end:]
    text=re.sub(r'These three interfaces were chosen.*?(?=\\input\{si_new_systems_counts_table)',
        'The three additional interfaces broaden the evaluation beyond GLP-1R. Two are peptide-GPCRs and one is a protein-protein interface. One non-GPCR system cannot establish a general relationship between interface geometry and screening performance.\n\n',text,flags=re.S)
    start=text.index('To verify that the reported enrichment');end=text.index(r'\input{decoy_robustness_table',start)
    text=text[:start]+('Five original matched-decoy sets per system were rescored from molecular inputs with structural-alert exclusions disabled. Each contains 10 actives and 300 decoys; these replicate sets are distinct from the expanded 50-active benchmarks. Table~S8 reports the measured mean and sample standard deviation. All decoys share the same ZINC source, so repeated draws cannot exclude systematic source-related bias.\n\n')+text[end:]
    dock=json.loads((ROOT/'evidence/outputs/equivalence/docking_wilcoxon.json').read_text())
    start=text.index('As a structure-based check');end=text.index(r'\input{top10_docking_table',start)
    text=text[:start]+rf'''The new top 10 were prepared and redocked against three active-state (6X18, 7KI0, 7LCJ) and two inactive-state (5VEW, 6LN2) GLP-1R structures, with three seeds per structure and exhaustiveness 16. Of 10 selected ligands, {dock['n']} have successful results against both states. The median active-minus-inactive difference is {dock['median_delta']:.2f}~kcal/mol; {dock['n_favoring_active']} paired ligands favor the active state. An exact signed-rank test using exhaustive sign flips, average ranks for ties, and zero differences excluded gives $W={dock['wilcoxon_W']:g}$ and two-sided $p={dock['p_two_sided']:.4f}$. This is a computational score comparison, not evidence of agonism. The raw output accounts for every requested job, including preparation and docking failures.

'''+text[end:]
    start=text.index("Each system's actives");end=text.index(r'\input{attrition_full_table',start)
    text=text[:start]+r'''Each system's original actives, matched decoys, and 30,000 background molecules were recomputed under non-excluding structural alerts and retained property limits. Table~S10 identifies every active absent from the final native ranking, including losses after shortlist admission. Table~S11 gives actual stage survival under the 5\% plus 1,000-molecule minimum policy. All 34 alert-flagged actives pass Stage 0; only two NTSR1 actives fail the property envelope. Later losses remain attributable to the hotspot gate, shortlist, preparation, native selection, or final ranking as recorded per molecule.

'''+text[end:]
    extra=r'''
\clearpage
\section*{Fresh Headline Rerun and Paired Shortlist Policies}
The headline screen was executed again from all original million-compound inputs with structural-alert exclusions disabled and the capped 1,000-molecule minimum on its 5\% shortlist. Its percentage denominator remains the Stage-0-pass count. Table~S12 compares this result with the archived strict-filter run; because both alert handling and the floor differ, the archive comparison does not isolate the floor.

Fresh production-policy pairs isolate the floor while holding non-excluding alerts fixed. These pairs retain the original production-benchmark denominator: candidates surviving Stages 1/2. Each pair shares fresh prescreening and the larger Stage-3 union, with separate native preparation and scoring when policy inputs differ. Identical shortlists explicitly share a single fresh native execution. Tables~S13 and S14 report complete counts and measured native-branch wall times, respectively.

\input{floor_headline_table.tex}
\clearpage
\input{floor_survival_table.tex}
\clearpage
\input{floor_timing_table.tex}

\clearpage
\section*{Tie-Aware Benchmark Evaluation}
Ligand identifiers encode different naming conventions for ChEMBL actives and ZINC decoys. Ordering equal scores by those identifiers inflated some earlier rank-based enrichment estimates. All primary benchmark statistics in this revision instead retain equal-status, equal-score ties, including missing scores within failure status. ROC-AUC and average precision operate at tied score thresholds. EF, BEDROC, and top-$k$ retrieval average over every possible within-tie ordering, so expected active counts can be fractional. EF is normalized by the actual selected fraction $\lceil fN\rceil/N$, correcting the earlier nominal-fraction denominator when rounding changes the selected fraction. All grouped-bootstrap intervals and paired comparisons were recomputed with the same tie treatment and EF normalization. The molecular scores, production shortlists, and observed final survival counts are unchanged by this evaluation correction.

Table~S15 compares the superseded identifier-ordered calculation with the current primary ROC-AUC and PR-AUC. Machine-readable records provide both original ranks and evaluation tie groups; the independent reproduction script computes the current statistics from those groups. Validation checks cover exhaustive within-tie permutations, scikit-learn ROC-AUC and average precision, identifier invariance, and explicit expansion of bootstrap samples.

\input{identifier_tie_table.tex}

'''
    text=text.replace(r'\end{document}',extra+r'\end{document}')
    text=text.replace('figA4_top20_structures.png','figA4_top20_structures.pdf')
    source.write_text(text,encoding='utf-8')
    print('Updated all quantitative SI passages and added measured floor experiments')


if __name__=='__main__':main()
