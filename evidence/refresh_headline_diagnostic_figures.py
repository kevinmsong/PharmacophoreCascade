"""Refresh legacy headline figure paths from fresh Stage-3 results at 600 dpi."""
import hashlib
import io
import json
import shutil
from pathlib import Path

import fitz
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import figstyle as fs

ROOT=Path(__file__).resolve().parents[1]
HEAD=ROOT/'results/absolute_floor_1000'
SUB=ROOT/'ACS_Omega_resubmission'
META=ROOT/'evidence/outputs/alert_disabled_revision'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def backup(path):
    assert path.resolve().is_relative_to(ROOT.resolve())
    saved=ROOT/'tmp/alert_disabled_revision/headline_figures_before'/path.relative_to(ROOT)
    if path.exists() and not saved.exists():
        saved.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,saved)


def main():
    assert json.loads((HEAD/'rerun_provenance.json').read_text())['status']=='completed'
    source=HEAD/'screening_full_1M_floor1000.csv'
    before=sha(source)
    frame=pd.read_csv(source).sort_values('stage3_screen_rank')
    assert len(frame)==49757
    for path in HEAD.glob('*.png'):
        backup(path)
    # This is the same receptor-feature input, with its current, labeled PCA
    # projection replacing the crowded original 3D scatter. Preserve aliases.
    for ext in ['.png','.pdf']:
        shutil.copy2(SUB/('fig2_pharmacophore'+ext),HEAD/('pharmacophore_3d_full_1M_floor1000'+ext))
    fs.use_target('acs');fs.apply()
    fig,axes=plt.subplots(2,3,figsize=(6.5,4.3))
    fig.suptitle(f'Successful Stage-3 headline cohort (n = {len(frame):,})',fontsize=9,weight='bold')
    for ax,column,label,color in zip(axes[0],['cascade_score_pct','weighted_coverage_pct','mw'],
            ['Cascade score (%)','Stage-3 weighted coverage (%)','Molecular weight (Da)'],
            ['blue','orange','grey']):
        ax.hist(frame[column],bins=40,color=fs.OKABE_ITO[color],edgecolor='white',linewidth=.15)
        ax.set(xlabel=label,ylabel='Ligands')
    for ax,x,y,xlab,ylab in [(axes[1,0],'hotspot_bits_matched','pair_hash_overlap_pct','Matched hotspot features','Pair overlap (%)'),
                            (axes[1,1],'cascade_score_pct','weighted_coverage_pct','Cascade score (%)','Stage-3 coverage (%)')]:
        ax.hexbin(frame[x],frame[y],gridsize=30,mincnt=1,bins='log',cmap='Blues',linewidths=0)
        ax.set(xlabel=xlab,ylabel=ylab)
        ax.set_title('Density: darker = more ligands',fontsize=6)
    axes[1,2].hist(frame.rotatable_bonds,bins=range(int(frame.rotatable_bonds.max())+2),
                   color=fs.OKABE_ITO['blue'],edgecolor='white',linewidth=.3)
    axes[1,2].set(xlabel='Rotatable bonds',ylabel='Ligands')
    fig.tight_layout(rect=[0,0,1,.95],w_pad=1.4,h_pad=1.5)
    fs.save(fig,HEAD,'property_distributions_full_1M_floor1000')
    top=frame.head(20)
    top.to_csv(HEAD/'top_20_stage3_figure_data.csv',index=False)
    mols=[Chem.MolFromSmiles(s) for s in top.smiles]
    assert all(m is not None for m in mols)
    for mol in mols:
        Chem.rdDepictor.Compute2DCoords(mol)
    legends=[f'{int(r.stage3_screen_rank)}. {r.zinc_id}\n{r.weighted_coverage_pct:.2f}% Stage-3 coverage' for r in top.itertuples()]
    stem='top_20_glp1_mimetics_full_1M_floor1000'
    for vector in [False,True]:
        drawer=(rdMolDraw2D.MolDraw2DSVG if vector else rdMolDraw2D.MolDraw2DCairo)(3900,4250,975,850)
        opts=drawer.drawOptions();opts.legendFontSize=58;opts.legendFraction=.20
        opts.fixedFontSize=48;opts.bondLineWidth=4.5;opts.useBWAtomPalette()
        drawer.DrawMolecules(mols,legends=legends);drawer.FinishDrawing()
        drawing=drawer.GetDrawingText()
        if not vector:
            Image.open(io.BytesIO(drawing)).save(HEAD/(stem+'.png'),dpi=(600,600))
        else:
            (HEAD/(stem+'.svg')).write_text(drawing,encoding='utf-8')
            with fitz.open(stream=drawing.encode('utf-8'),filetype='svg') as svg:
                with fitz.open(stream=svg.convert_to_pdf(),filetype='pdf') as raw:
                    with fitz.open() as pdf:
                        page=pdf.new_page(width=468,height=510)
                        page.show_pdf_page(page.rect,raw,0);pdf.save(HEAD/(stem+'.pdf'))
    outputs={}
    for path in sorted(HEAD.iterdir()):
        if path.suffix not in ['.png','.pdf','.svg'] or 'full_1M_floor1000' not in path.name:
            continue
        alias=ROOT/'results'/path.name.replace('full_1M_floor1000','full_1M_topological_hashed')
        backup(alias);shutil.copy2(path,alias);outputs[str(alias.relative_to(ROOT))]=sha(alias)
    for src,dst in [('top_100_glp1_mimetics_full_1M_floor1000.csv','top_100_glp1_mimetics_full_1M_topological_hashed.csv'),
                    ('rerun.log','screening_full_1M_topological_hashed_external_run.log')]:
        alias=ROOT/'results'/dst;backup(alias);shutil.copy2(HEAD/src,alias)
        outputs[str(alias.relative_to(ROOT))]=sha(alias)
    stale_pid=ROOT/'results/screening_full_1M_topological_hashed_external_terminal.pid'
    if stale_pid.exists():
        assert stale_pid.resolve().is_relative_to((ROOT/'results').resolve())
        backup(stale_pid)
        stale_pid.unlink()
    note=('These diagnostics describe the freshly scored Stage-3 cohort and its top 20, not the final native top 20. '
          'The legacy pharmacophore_3d filename now contains the labeled PCA projection of the same receptor features. '
          'Native-ranked structures are in ACS_Omega_resubmission/figA4_top20_structures. All raster figures are rendered at 600 dpi; vector counterparts are supplied.\n')
    (HEAD/'FIGURES.md').write_text(note,encoding='utf-8')
    (ROOT/'results/FIGURES.md').write_text(note,encoding='utf-8')
    assert sha(source)==before
    (META/'headline_figure_refresh.json').write_text(json.dumps(dict(status='completed',source_sha256=before,
        source=str(source.relative_to(ROOT)),molecular_scores_changed=False,canonical_outputs=outputs),indent=2),encoding='utf-8')
    print('Regenerated three headline diagnostics at 600 dpi/vector resolution and refreshed all canonical figure/top-100 aliases.')


if __name__=='__main__':main()
