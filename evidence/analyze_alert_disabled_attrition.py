"""Measured final-retention accounting under non-excluding alerts and a 1,000 floor."""
import json
import sys
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference_implementation'))
from stage_scoring_reference import stage0
DATA=ROOT/'evidence/outputs/absolute_floor'
OUT=ROOT/'evidence/outputs/attrition'
SUB=ROOT/'ACS_Omega_resubmission'
SYSTEMS={'glp1r':'GLP-1R','ghsr':'GHSR','ntsr1':'NTSR1','mdm2':'MDM2-p53'}
STAGES=[('stage0_pass','Stage 0','property envelope'),
        ('stage012_pass','Stage 1','hotspot gate'),
        ('shortlist','Shortlist','below shortlist cut'),
        ('stage3_success','Stage 3','3D scoring failure'),
        ('native_selected','Native selection','native pool or scaffold cap'),
        ('native_success','Native scoring','native preparation or scoring failure'),
        ('final_ranked','Final top-1,000','below final rank limit')]


def tex(text):
    # SYSTEMS keeps the hyphen because its labels are also the CSV keys the
    # figures index on; only the typeset copy takes the manuscript's en dash.
    return (str(text).replace('&',r'\&').replace('_',r'\_').replace('%',r'\%')
            .replace('#',r'\#').replace('MDM2-p53','MDM2--p53'))


def main():
    rows=[];survival=[]
    for key,label in SYSTEMS.items():
        if json.loads((DATA/key/'provenance.json').read_text())['status']!='completed':
            raise RuntimeError(f'Incomplete production rerun: {key}')
        per=pd.read_csv(DATA/key/'floor1000/per_molecule_survival.csv')
        stages=pd.read_csv(DATA/key/'stage_survival.csv')
        survival.append(stages[stages.policy=='floor1000'])
        for row in per[per.label=='active'].itertuples():
            gate=stage0(row.smiles)
            assert gate['pass']==row.stage0_pass
            stage,cause,detail='-','retained','present in the final native ranking'
            for col,name,reason in STAGES:
                if not getattr(row,col):
                    stage,cause=name,reason
                    detail={'stage0_pass':gate['reason'],
                        'stage012_pass':'required hotspot matches or contact groups not met',
                        'shortlist':'outside the 5% shortlist with a minimum of 1,000 eligible candidates',
                        'stage3_success':'no successful conformer-level Stage-3 score',
                        'native_selected':'not selected by native-pool quotas and the eight-per-scaffold cap',
                        'native_success':'no successful native score after ligand preparation',
                        'final_ranked':'native rank below the final top-1,000 cutoff'}[col]
                    if col=='stage0_pass' and gate.get('mol') is None:cause='unparseable'
                    break
            rows.append({'system':label,'ligand_id':row.ligand_id,'smiles':row.smiles,
                'policy':'5% + floor1000; alerts non-excluding','stage_lost':stage,
                'reason_class':cause,'reason_detail':detail,'structural_alerts':'; '.join(gate.get('alerts',[])),
                'stage0_pass':row.stage0_pass,'final_ranked':row.final_ranked})
    frame=pd.DataFrame(rows)
    OUT.mkdir(parents=True,exist_ok=True)
    frame.to_csv(OUT/'active_attrition_per_molecule.csv',index=False)
    summary=frame.groupby(['system','stage_lost','reason_class']).size().reset_index(name='n')
    summary.to_csv(OUT/'active_attrition_summary.csv',index=False)
    lines=[r'\begin{table}[tbp]',r'\centering',
        r'\caption{\textbf{In-domain actives removed at each stage, by cause.} '
        r'Production runs record structural alerts without excluding, retain the molecular-property limits, '
        r'and use a 5\% shortlist with a minimum of 1,000 eligible molecules. '
        r'Counts follow measured preparation, selection, and scoring outcomes. '
        r'Table~S10 identifies individual losses; Table~S11 reports surviving actives.}',
        r'\label{tab:attrition}',r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',
        r'\begin{tabular}{@{}llrrrr@{}}',r'\toprule',
        r'Stage & Cause & GLP-1R & GHSR & NTSR1 & MDM2--p53 \\',r'\midrule']
    for _,stage,_ in STAGES:
        selected=frame[frame.stage_lost==stage]
        for cause in selected.reason_class.unique():
            counts=[str(int(((selected.system==s)&(selected.reason_class==cause)).sum())) for s in SYSTEMS.values()]
            lines.append(tex(stage)+' & '+tex(cause)+' & '+' & '.join(counts)+r' \\')
    lines.append(r'\midrule')
    cells=[]
    for label in SYSTEMS.values():
        sub=frame[frame.system==label]
        cells.append(f'{int(sub.final_ranked.sum())}/{len(sub)}')
    lines += [r'Final ranking & Retained & '+' & '.join(cells)+r' \\',r'\bottomrule',r'\end{tabular}',r'\end{table}','']
    (SUB/'attrition_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    lost=frame[frame.stage_lost!='-']
    lines=[r'\footnotesize',r'\setlength{\tabcolsep}{4pt}',
        r'\begin{longtable}{@{}p{.10\textwidth}p{.21\textwidth}p{.16\textwidth}p{.42\textwidth}@{}}',
        r'\caption{\textbf{In-domain actives absent from the final ranking under production constraints.} '
        r'Structural alerts are annotations, not exclusions. Causes refer to the first failed or limiting stage.}'
        r'\label{tab:attrition_full}\\',r'\toprule',r'System & Ligand & Stage & Cause \\',r'\midrule',
        r'\endfirsthead',r'\toprule',r'System & Ligand & Stage & Cause \\',r'\midrule',r'\endhead']
    for row in lost.itertuples():
        lines.append(' & '.join(tex(v) for v in [row.system,row.ligand_id,row.stage_lost,row.reason_detail])+r' \\')
    lines += [r'\bottomrule',r'\end{longtable}',r'\normalsize']
    (SUB/'attrition_full_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    stages=pd.concat(survival)
    lines=[r'\begin{table}[htbp]',r'\centering',
        r'\caption{\textbf{Stage-by-stage active survival under production constraints.} '
        r'Each labeled library sits in a 30,000-molecule background, with structural alerts recorded '
        r'rather than excluding and a minimum of 1,000 molecules on the 5\% shortlist. '
        r'Entries count actives and include the measured native-preparation and scoring outcomes.}',
        r'\label{tab:production}',r'\footnotesize',r'\setlength{\tabcolsep}{6pt}',
        r'\begin{tabular}{@{}lrrrr@{}}',r'\toprule',
        r'Stage & GLP-1R & GHSR & NTSR1 & MDM2--p53 \\',r'\midrule']
    for column,label in [('input','Input')]+[(c,n) for c,n,_ in STAGES]:
        cells=[str(int(stages[(stages.system==key)&(stages.stage==column)].actives.iloc[0])) for key in SYSTEMS]
        lines.append(tex(label)+' & '+' & '.join(cells)+r' \\')
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
    (SUB/'production_survival_table.tex').write_text('\n'.join(lines),encoding='utf-8')
    print(frame.groupby(['system','stage_lost']).size().to_string())


if __name__=='__main__':main()
