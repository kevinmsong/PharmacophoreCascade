"""Generate publication tables exclusively from completed revision records."""
import json
import shutil
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUB=ROOT/'ACS_Omega_resubmission'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def write_table(name,caption,label,columns,rows,alignment=None):
    alignment=alignment or 'l'+'r'*(len(columns)-1)
    lines=[r'\begin{table}[tbp]',r'\centering',r'\caption{'+caption+'}',
           r'\label{'+label+'}',r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',
           r'\begin{tabular}{@{}'+alignment+'@{}}',r'\toprule',
           ' & '.join(columns)+r' \\',r'\midrule']
    lines+=[' & '.join(map(str,row))+r' \\' for row in rows]
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}','']
    text='\n'.join(lines)
    (SUB/(name+'.tex')).write_text(text,encoding='utf-8')
    (META/(name+'.tex')).write_text(text,encoding='utf-8')


def main():
    d=json.loads((META/'manuscript_results.json').read_text())
    b=d['benchmarks'];diag=d['diagnostics']
    systems={'glp1r':'GLP-1R','ghsr':'GHSR','ntsr1':'NTSR1','mdm2':'MDM2--p53'}
    methods={'full_cascade':'Full cascade','native_only':'Native-only',
             'stage3_only':'Stage-3 only','standard_3d_pharmacophore':'Single-pass 3D'}
    metrics=['roc_auc','pr_auc','ef_1pct','bedroc']
    rows=[]
    for system,label in systems.items():
        data=b[system+'_full']
        best={k:max(data[m][k] for m in methods) for k in metrics}
        for index,(method,methodlabel) in enumerate(methods.items()):
            cells=[]
            for k in metrics:
                value=data[method][k];cell=f'{value:.1f}' if k=='ef_1pct' else f'{value:.3f}'
                if abs(value-best[k])<1e-12:cell=r'\textbf{'+cell+'}'
                cells.append(cell)
            rows.append([label if index==0 else '',methodlabel]+cells)
    write_table('crosssystem_table',r'\textbf{Retrospective enrichment across all four systems.} '
        r'Each method ranks the same labeled universe: 10 GLP-1R actives or 50 actives per other system, '
        r'with 30 matched decoys per active. Structural alerts are recorded without exclusion. '
        r'Equal scores remain tied: EF and BEDROC average over within-tie orders. '
        r'BEDROC uses $\alpha=20$. Bold indicates the highest value within each system.',
        'tab:crosssystem',['System','Method','ROC-AUC','PR-AUC',r'EF1\%','BEDROC'],rows,'llrrrr')
    upstream=diag['upstream_correlations']['screen_weighted_coverage_pct']
    write_table('evidence_table',r'\textbf{Comparison of screen-stage and native rankings.} '
        f"The common cohort contains {diag['n']:,} successfully native-scored headline ligands. "
        r'Both ranks are computed within that same cohort; the Stage-3 ordering uses the emitted '
        r'Stage-3 screen rank, and the native ordering uses the final rank. No ranks from different '
        r'population sizes are subtracted.',
        'tab:evidence',['Diagnostic','Value'],[
        [r'Stage-3 vs native coverage: Pearson $r$ / Spearman $\rho$',f"{upstream['pearson_r']:.3f} / {upstream['spearman_rho']:.3f}"],
        [r'Stage-3 vs native cohort rank: Spearman $\rho$',f"{diag['rank_spearman']:.3f}"],
        ['Median absolute cohort-rank shift, final top 1,000',f"{diag['top1000_median_absolute_cohort_rank_shift']:,.0f}"],
        ['Shared molecules in cohort top-10 sets',f"{diag['top10_overlap']} of 10"]])
    rows=[]
    for key,label in [('roc_auc','ROC-AUC'),('pr_auc','PR-AUC'),('ef_1pct',r'EF1\%'),('bedroc',r'BEDROC ($\alpha=20$)'),('top10_recovery','Actives in top 10')]:
        vals=[b[x]['full_cascade'][key] for x in ['glp1r_full','glp1r_automated']]
        # EF carries one decimal everywhere; the other metrics carry three.
        if key=='top10_recovery':cells=[f'{x*10:.2f} of 10' for x in vals]
        elif key=='ef_1pct':cells=[f'{x:.1f}' for x in vals]
        else:cells=[f'{x:.3f}' for x in vals]
        rows.append([label]+cells)
    rows.append(['Single-pass 3D ROC-AUC']+[f"{b[x]['standard_3d_pharmacophore']['roc_auc']:.3f}" for x in ['glp1r_full','glp1r_automated']])
    write_table('ablation_curation_table',r'\textbf{Cascade performance with and without manual curation (GLP-1R).} '
        r'Both columns use the same 10 actives and 300 decoys, with structural alerts recorded rather than '
        r'excluding. The automated model removes hand-assigned contact groups, weight bonuses, native-supported '
        r'expansion, and curated rescue. Top-10 values are expected active counts under score ties.',
        'tab:ablation_curation',['Metric','Curated','Automated'],rows)
    rows=[]
    for r in pd.DataFrame(d['top10']).itertuples():
        rows.append([int(r.final_rank),r.zinc_id,f'{r.native_weighted_coverage_pct:.2f}',
            int(r.native_matched_reference_features),f'{r.native_fit_rmsd_angstrom:.2f}',
            f'{r.native_mean_pair_distance_error_angstrom:.2f}',f'{r.stage3_screen_rank:,}',f'{r.mw:.1f}',f'{r.logp:.2f}'])
    write_table('top10_table',r'\textbf{Top 10 ligands ranked by native peptide-interface mimicry.} '
        r'Coverage is the percentage of retained peptide-reference weight recovered by the best prepared state. '
        r'Feat. is the matched feature count; RMSD and pair error are in \angstrom. Stage-3 ranks refer '
        r'to the full successful Stage-3 cohort. MW is in Da. Structures appear in Figure~S4.',
        'tab:top10',['Rank','ZINC ID',r'Cov. (\%)','Feat.','RMSD','Pair err.','3D rank','MW','LogP'],rows,'rlrrrrrrr')
    reps=pd.DataFrame(d['replicates']);rows=[]
    for key,label in systems.items():
        sub=reps[(reps.system==key)&(reps.method=='full_cascade')]
        assert len(sub)==5
        rows.append([label]+[f'{sub[k].mean():.3f} $\\pm$ {sub[k].std(ddof=1):.3f}' for k in ['roc_auc','pr_auc','bedroc']])
    write_table('decoy_robustness_table',r'\textbf{Robustness of the full cascade across five independent matched-decoy sets.} '
        r'Values are mean $\pm$ sample standard deviation over five runs per system. Each replicate '
        r'contains the 10-active set and 300 matched decoys, with structural alerts recorded rather than '
        r'excluding; these replicates are distinct from the 50-active comparisons.',
        'tab:decoy_robustness',['System','ROC-AUC','PR-AUC','BEDROC'],rows)
    dock=pd.read_csv(ROOT/'evidence/outputs/revision/docking_top10_summary.csv')
    # An en dash, not "---": the manuscript sets no em dashes anywhere.
    fmt=lambda x:'--' if pd.isna(x) else f'{x:.2f}'
    rows=[[int(r.final_rank),r.ligand_id,fmt(r.best_active),fmt(r.best_inactive),fmt(r.active_pref)] for r in dock.itertuples()]
    write_table('top10_docking_table',r'\textbf{Redocking of the top-10 final-ranked GLP-1R ligands.} '
        r'Affinities and active-minus-inactive differences are in kcal/mol. Each state value is the best '
        r'successful result across its receptor structures and three seeds. Missing paired results are '
        r'shown as dashes and excluded from the paired test; all requested jobs remain in the raw output.',
        'tab:top10_docking',['Rank','Ligand ID','Active','Inactive',r'$\Delta$'],rows,'rlrrr')
    # The statistical writers now emit booktabs directly, matching every other
    # table here, so these two only need copying into the submission directory.
    for filename,folder in [('ablation_table','ablation_common'),('tost_table','equivalence')]:
        text=(ROOT/'evidence/outputs'/folder/(filename+'.tex')).read_text()
        assert r'\hline' not in text, f'{filename} still emits \\hline rules'
        (SUB/(filename+'.tex')).write_text(text,encoding='utf-8')
    rows=[]
    for system,label in systems.items():
        audit=pd.read_csv(ROOT/f'evidence/outputs/benchmark_{system}_full/identifier_tie_sensitivity.csv')
        for i,(method,methodlabel) in enumerate(methods.items()):
            a=audit[audit.method==method].set_index('metric')
            rows.append([label if i==0 else '',methodlabel]+[
                f'{a.loc[metric,col]:.3f}' for metric in ['roc_auc','pr_auc'] for col in ['identifier_ordered','tie_aware']])
    write_table('identifier_tie_table',r'\textbf{Sensitivity of benchmark enrichment to how equal scores are ordered.} '
        r'The ID columns order molecules with equal scores by ligand identifier; the tied columns, used '
        r'throughout this work, leave them tied and average over the possible within-tie orders. Molecules '
        r'with equal status and score are tied, including unscored molecules within a failure status. '
        r'Both columns evaluate the same molecular scores on the same molecules.',
        'tab:identifier_ties',['System','Method','ROC ID','ROC tied','PR ID','PR tied'],rows,'llrrrr')
    print('Wrote current main and supplementary quantitative tables')


if __name__=='__main__':main()
