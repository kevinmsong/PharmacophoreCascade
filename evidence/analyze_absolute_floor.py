"""Build publication tables and 600-dpi/vector figures from completed floor reruns."""
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'evidence/outputs/absolute_floor'
SUB=ROOT/'ACS_Omega_resubmission'
HEAD=ROOT/'results/absolute_floor_1000'
ARCHIVE=ROOT/'tmp/alert_disabled_revision/prior_headline'
SYSTEMS={'glp1r':'GLP-1R','ghsr':'GHSR','ntsr1':'NTSR1','mdm2':'MDM2–p53'}
#: LaTeX spelling of the same labels; MDM2--p53 is an en dash in the source.
TEX_SYSTEMS={'glp1r':'GLP-1R','ghsr':'GHSR','ntsr1':'NTSR1','mdm2':'MDM2--p53'}
#: "Minimum" rather than "floor" throughout: the reviewers asked for plainer terms.
POLICIES={'percentage':'5% only','floor1000':'5% + minimum 1,000'}


def load_completed():
    rows=[]
    survival=[]
    for key,label in SYSTEMS.items():
        path=DATA/key
        provenance=json.loads((path/'provenance.json').read_text())
        if provenance['status']!='completed': raise RuntimeError(f'{key} has not completed')
        stages=pd.read_csv(path/'stage_survival.csv')
        survival.append(stages)
        summaries=json.loads((path/'policy_summary.json').read_text())
        for summary in summaries:
            group=stages[stages.policy==summary['policy']].set_index('stage')
            rec={'system':key,'system_label':label,'policy':summary['policy'],
                 'n_input':int(group.loc['input','total']),'n_actives':int(group.loc['input','actives']),
                 'native_wall_seconds':summary['native_wall_sec'],
                 'stage012_seconds':summary['stage012_sec'],'union_stage3_seconds':summary['union_stage3_sec'],
                 'shared_native_results':summary['native_results_shared_with_percentage']}
            for stage in group.index:
                rec[stage+'_n']=int(group.loc[stage,'total'])
                rec[stage+'_actives']=int(group.loc[stage,'actives'])
            rec['final_retention_pct']=100*rec['final_ranked_actives']/rec['n_actives']
            rows.append(rec)
    return pd.DataFrame(rows),pd.concat(survival,ignore_index=True)


def write_table(rows):
    lines=[r'\begin{table}[tbp]',r'\centering',r'\caption{\textbf{What shortlist depth changes under production constraints.} '
        r'Each matched pair of shortlist rules shares one evaluation of Stages 0--2 and the Stage-3 union. '
        r'Native preparation and scoring are repeated when the shortlist differs. '
        r'The 1,000-molecule minimum is capped by the number of eligible candidates. '
        r'Native counts are successful ligand scores; final retention additionally reflects '
        r'scaffold selection and the final top-1,000 limit.}',r'\label{tab:efficiency}',
        r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{@{}llrrrr@{}}',r'\toprule',
        r'System & Shortlist rule & Shortlist & Native & Shortlist actives & Final actives \\',r'\midrule']
    for key,label in TEX_SYSTEMS.items():
        for r in rows[rows.system==key].to_dict('records'):
            policy=r'5\% only' if r['policy']=='percentage' else r'5\% + min.\ 1,000'
            lines.append(f"{label} & {policy} & {r['shortlist_n']:,} & {r['native_success_n']:,} & "
                f"{r['shortlist_actives']}/{r['n_actives']} & {r['final_ranked_actives']}/{r['n_actives']} \\\\")
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}','']
    (SUB/'efficiency_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    shutil.copy2(SUB/'efficiency_table.tex',DATA/'efficiency_table.tex')


def figures(rows,stages):
    fs.use_target('acs');fs.apply()
    colors=[fs.OKABE_ITO['orange'],fs.OKABE_ITO['blue']]
    markers=['o','s']
    fig,axes=plt.subplots(1,2,figsize=(6.5,3.0),gridspec_kw={'width_ratios':[1,1.12]})
    labels=list(SYSTEMS.values())
    y=np.arange(4)[::-1]
    for ax in axes:
        ax.set_yticks(y,labels);ax.grid(axis='y',visible=False)
    for i,(key,label) in enumerate(SYSTEMS.items()):
        sub=rows[rows.system==key].set_index('policy')
        counts=[sub.loc[p,'shortlist_n'] for p in POLICIES]
        rets=[sub.loc[p,'final_retention_pct'] for p in POLICIES]
        axes[0].plot(counts,[y[i],y[i]],color='#777777',lw=1,zorder=1)
        axes[1].plot(rets,[y[i],y[i]],color='#777777',lw=1,zorder=1)
        for j,p in enumerate(POLICIES):
            dy=(-.07 if j==0 else .07) if counts[0]==counts[1] else 0
            axes[0].scatter(counts[j],y[i]+dy,s=26,marker=markers[j],color=colors[j],edgecolor='white',linewidth=.4,zorder=3)
            axes[0].annotate(f'{int(counts[j]):,}',(counts[j],y[i]+dy),xytext=(0,-13 if j==0 else 7),textcoords='offset points',ha='center',fontsize=7)
            dy=(-.07 if j==0 else .07) if rets[0]==rets[1] else 0
            axes[1].scatter(rets[j],y[i]+dy,s=26,marker=markers[j],color=colors[j],edgecolor='white',linewidth=.4,zorder=3)
            axes[1].annotate(f"{int(sub.loc[p,'final_ranked_actives'])}/{int(sub.loc[p,'n_actives'])}",(rets[j],y[i]+dy),
                xytext=(0,-13 if j==0 else 7),textcoords='offset points',ha='center',fontsize=7)
    axes[0].set_xscale('log');axes[0].set_xlim(7,5000)
    axes[0].set_xlabel('Molecules shortlisted (log scale)')
    axes[1].set_xlim(0,103);axes[1].set_xlabel('Actives in final ranking (%)')
    for ax in axes: ax.set_ylim(-.55,3.5)
    for letter,ax in zip('ab',axes): fs.panel_label(ax,f'({letter})',dx=-.20,dy=1.06)
    handles=[Line2D([],[],color=c,marker=m,linestyle='none',label=POLICIES[p]) for c,m,p in zip(colors,markers,POLICIES)]
    fig.legend(handles=handles,loc='lower center',ncol=2,bbox_to_anchor=(.52,-.005),fontsize=8)
    fig.tight_layout(rect=[0,.11,1,1],w_pad=1.5)
    fs.save(fig,SUB,'fig6_efficiency')
    stage_order=['stage0_pass','stage012_pass','shortlist','stage3_success','native_selected','native_success','final_ranked']
    ticks=['Stage 0','Gate','Shortlist','3D','Selected','Scored','Final']
    fig,axes=plt.subplots(2,2,figsize=(6.5,4.7),sharex=True)
    for ax,(key,label),letter in zip(axes.flat,SYSTEMS.items(),'abcd'):
        n=int(rows[rows.system==key].n_actives.iloc[0])
        for j,p in enumerate(POLICIES):
            sub=stages[(stages.system==key)&(stages.policy==p)].set_index('stage')
            vals=sub.loc[stage_order,'actives'].to_numpy()
            ax.plot(range(7),vals,marker=markers[j],color=colors[j],linestyle='--' if j==0 else '-',label=POLICIES[p],lw=1.3)
        ax.set_title(f'{label} ({n} input actives)',fontsize=8.5)
        ax.set_ylim(0,n*1.12);ax.set_ylabel('Actives retained')
        ax.set_xticks(range(7),ticks,rotation=38,ha='right',fontsize=6.7)
        ax.grid(axis='x',visible=False);fs.panel_label(ax,f'({letter})',dx=-.16,dy=1.12)
    fig.legend(handles=handles,loc='lower center',ncol=2,fontsize=8)
    fig.tight_layout(rect=[0,.07,1,1],h_pad=1.6,w_pad=1.7)
    fs.save(fig,SUB,'figA5_floor_survival')


def headline_compare():
    provenance=json.loads((HEAD/'rerun_provenance.json').read_text())
    if provenance['status']!='completed': raise RuntimeError('Headline rerun has not completed')
    old=json.loads((ARCHIVE/'screening_full_1M_topological_hashed_run_summary.json').read_text())
    new=json.loads((HEAD/'screening_full_1M_floor1000_run_summary.json').read_text())
    old_short=pd.read_csv(ARCHIVE/'screening_full_1M_topological_hashed_shortlist.csv')
    new_short=pd.read_csv(HEAD/'screening_full_1M_floor1000_shortlist.csv')
    old_final=pd.read_csv(ARCHIVE/'top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv')
    new_final=pd.read_csv(HEAD/'top_1000_glp1_mimetics_full_1M_floor1000_native_final.csv')
    common=old_final.merge(new_final,on='zinc_id',suffixes=('_old','_new'))
    comparisons = {}
    for name, old_path, new_path, columns in [
        ('stage3', ARCHIVE/'screening_full_1M_topological_hashed.csv',
         HEAD/'screening_full_1M_floor1000.csv',
         ['weighted_coverage_pct','fit_rmsd_angstrom','mean_pair_distance_error_angstrom']),
        ('all_native', ARCHIVE/'screening_full_1M_topological_hashed_native_scored_top5000.csv',
         HEAD/'screening_full_1M_floor1000_native_scored_top5000.csv',
         ['native_weighted_coverage_pct','native_fit_rmsd_angstrom','native_mean_pair_distance_error_angstrom'])]:
        a=pd.read_csv(old_path);b=pd.read_csv(new_path)
        matched=a.merge(b,on='zinc_id',suffixes=('_old','_new'))
        comparisons[name]={'archived_n':len(a),'fresh_n':len(b),'common_n':len(matched),
            'same_order':a.zinc_id.tolist()==b.zinc_id.tolist(),
            'max_abs_deltas':{col:float((matched[col+'_old']-matched[col+'_new']).abs().max()) for col in columns}}
    result={'comparison_scope':'Archived strict-filter run versus fresh alert-disabled floor run; differences are not attributed to the floor alone.',
        'chemistry_gate_mode':new['settings'].get('chemistry_gate_mode','unknown') if 'settings' in new else provenance['chemistry_gate_mode'],
        'archived_counts':old['counts'],'rerun_counts':new['counts'],
        'rerun_native_counts':new['native_rerank']['counts'],
        'shortlist_ids_same_order_as_archive':old_short.zinc_id.tolist()==new_short.zinc_id.tolist(),
        'shortlist_overlap_with_archive':len(set(old_short.zinc_id)&set(new_short.zinc_id)),
        'final_ids_same_order_as_archive':old_final.zinc_id.tolist()==new_final.zinc_id.tolist(),
        'final_overlap_with_archive':len(common),
        'native_score_max_abs_delta_common_final':float((common.native_weighted_coverage_pct_old-common.native_weighted_coverage_pct_new).abs().max()),
        'floor_changes_current_run_shortlist_size':new['counts']['shortlist_size']!=new['counts']['percentage_only_shortlist_size'],
        'rerun_wall_seconds':provenance['elapsed_seconds'],'pipeline_wall_seconds':new['timings']['total_pipeline_sec']}
    result['full_cohort_checks']=comparisons
    (DATA/'headline_comparison.json').write_text(json.dumps(result,indent=2))
    return result


def supplementary_tables(rows,stages,headline):
    old=headline['archived_counts'];new=headline['rerun_counts']
    native=headline['rerun_native_counts']
    original=json.loads((ARCHIVE/'screening_full_1M_topological_hashed_run_summary.json').read_text())['native_rerank']['counts']
    comparisons=[('Input compounds',old['total_scanned'],new['total_scanned']),
        ('Stage-0 survivors',old['property_pass'],new['property_pass']),
        ('Stage-1/2 candidates',old['hotspot_pass'],new['hotspot_pass']),
        ('Stage-3 shortlist',old['shortlist_size'],new['shortlist_size']),
        ('Successful Stage-3 scores',old['final_hits'],new['final_hits']),
        ('Native candidate pool',original['candidate_pool_size'],native['candidate_pool_size']),
        ('Selected for native scoring',original['selected_ligands'],native['selected_ligands']),
        ('Successful native ligand scores',original['successful_best_ligands'],native['successful_best_ligands']),
        ('Final ranking',original['final_rows'],native['final_rows'])]
    # The archived screen differed from the current one in alert mode, shortlist
    # minimum, and worker count at once, so a two-column table of the two runs
    # cannot attribute any difference to a single setting. The counts stay in the
    # released CSV as provenance; the document does not present them as a result.
    pd.DataFrame(comparisons,columns=['stage','archived','floor_rerun']).to_csv(DATA/'headline_stage_counts.csv',index=False)
    names={'input':'Input','stage0_pass':'Stage 0','stage012_pass':'Stages 1/2','shortlist':'Shortlist',
           'stage3_success':'Successful 3D','native_selected':'Native selection','native_success':'Successful native','final_ranked':'Final ranking'}
    lines=[r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',
        r'\begin{longtable}{@{}p{.13\textwidth}p{.25\textwidth}>{\raggedleft\arraybackslash}p{.16\textwidth}>{\raggedleft\arraybackslash}p{.10\textwidth}>{\raggedleft\arraybackslash}p{.16\textwidth}>{\raggedleft\arraybackslash}p{.10\textwidth}@{}}',
        r'\caption{\textbf{Stage-by-stage counts for each pair of shortlist rules.} '
        r'Entries are total ligands and in-domain actives, respectively. All systems use the same '
        r'30,000-molecule background, fixed 0.25/0.75 production-benchmark weights, and the candidate '
        r'count as the percentage denominator.}\label{tab:floor-survival}\\',
        r'\toprule',r' & & \multicolumn{2}{c}{5\% only} & \multicolumn{2}{c}{5\% + min.\ 1,000} \\',
        r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}',
        r'System & Stage & Total & Actives & Total & Actives \\',r'\midrule',r'\endfirsthead',
        r'\toprule',r'System & Stage & Total & Actives & Total & Actives \\',r'\midrule',r'\endhead']
    for key,label in TEX_SYSTEMS.items():
        base=stages[(stages.system==key)&(stages.policy=='percentage')].set_index('stage')
        floor=stages[(stages.system==key)&(stages.policy=='floor1000')].set_index('stage')
        for stage,name in names.items():
            a=base.loc[stage];b=floor.loc[stage]
            lines.append(f'{label} & {name} & {int(a.total):,} & {int(a.actives)} & {int(b.total):,} & {int(b.actives)} \\\\')
        lines.append(r'\midrule')
    # The loop leaves a \midrule after the last system; a longtable must close on a
    # \bottomrule or its final row renders with no rule beneath it.
    lines[-1] = r'\bottomrule'
    lines += [r'\end{longtable}',r'\setlength{\tabcolsep}{6pt}',r'\normalsize','']
    (SUB/'floor_survival_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    lines=[r'\begin{table}[htbp]',r'\centering',r'\caption{\textbf{Measured native-branch wall times for each pair of shortlist rules.} '
        r'Values cover native selection, ligand preparation, scoring, and output generation. '
        r'Stages 0--2 and the larger Stage-3 union are shared within each pair and are excluded here. '
        r'Identical shortlists share one native execution. These runs execute alongside the million-compound '
        r'screen, so the elapsed times are workload-specific measurements rather than portable speedup estimates.}',
        r'\label{tab:floor-timing}',r'\footnotesize',r'\setlength{\tabcolsep}{6pt}',
        r'\begin{tabular}{@{}lrrl@{}}',r'\toprule',
        r'System & 5\% only (min) & 5\% + min.\ 1,000 (min) & Native execution \\',r'\midrule']
    for key,label in TEX_SYSTEMS.items():
        r=rows[rows.system==key].set_index('policy')
        shared=bool(r.loc['floor1000','shared_native_results'])
        lines.append(f"{label} & {r.loc['percentage','native_wall_seconds']/60:.2f} & {r.loc['floor1000','native_wall_seconds']/60:.2f} & "
                     +('Shared (identical input)' if shared else 'Separate per policy')+r' \\')
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}','']
    (SUB/'floor_timing_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    for stem in ['floor_survival_table','floor_timing_table']:
        shutil.copy2(SUB/(stem+'.tex'),DATA/(stem+'.tex'))


def main():
    rows,stages=load_completed()
    rows.to_csv(DATA/'policy_comparison.csv',index=False)
    stages.to_csv(DATA/'all_stage_survival.csv',index=False)
    headline=headline_compare()
    write_table(rows);figures(rows,stages)
    supplementary_tables(rows,stages,headline)
    summary={'headline':headline,'benchmarks':rows.to_dict('records')}
    (DATA/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
    print(rows[['system','policy','shortlist_n','native_success_n','final_ranked_actives','native_wall_seconds']].to_string(index=False))


if __name__=='__main__': main()
