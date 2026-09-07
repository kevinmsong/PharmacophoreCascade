"""Graphical abstract for ACS Omega, in the journal's 3.25 x 1.75-inch TOC box.

The graphic answers the question the paper answers: the terminal native
peptide-contact score supplies the enrichment, on every system tested, while the
staged cascade upstream supplies the throughput. Configuration detail belongs in
the Methods, and every number here is read from the recorded outputs rather than
typed in, so the graphic cannot drift away from the tables it summarizes.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle as fs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'ACS_Omega_resubmission'
BENCH = ROOT / 'evidence/outputs'
HEADLINE = ROOT / 'results/screening_full_1M_topological_hashed_run_summary.json'

#: Benchmark directories, in the order the manuscript reports them. The en dash
#: matches the manuscript's spelling of the protein-protein pair.
SYSTEMS = {'glp1r': 'GLP-1R', 'ghsr': 'GHSR', 'ntsr1': 'NTSR1', 'mdm2': 'MDM2–p53'}
#: The two methods the graphic contrasts, terminal score against the conventional one.
ARMS = [('native_only', 'Native peptide-contact score'),
        ('standard_3d_pharmacophore', 'Conventional single-pass 3D')]


def measured():
    """Collect what the graphic states, from the recorded runs."""
    run = json.loads(HEADLINE.read_text(encoding='utf-8'))
    roc = {}
    for key in SYSTEMS:
        summary = pd.read_csv(BENCH / f'benchmark_{key}_full/benchmark_summary.csv',
                              index_col='method')['roc_auc']
        roc[key] = [float(summary[m]) for m, _ in ARMS]
    return dict(
        library=int(run['counts']['total_scanned']),
        final=int(run['native_rerank']['counts']['final_rows']),
        hours=run['timings']['total_pipeline_sec'] / 3600.0,
        roc=roc,
        wins=sum(v[0] > v[1] for v in roc.values()),
    )


def main():
    d = measured()
    fs.apply()
    green, grey = fs.OKABE_ITO['green'], '#9A9A9A'
    ink, mid = '#1A1A1A', '#4D4D4D'

    fig = plt.figure(figsize=(3.25, 1.75))

    # --- the question, and the scale it was answered at ---------------------
    fig.text(.5, .955, 'Where does the enrichment come from?', ha='center', va='top',
             fontsize=8.6, weight='bold', color=ink)
    fig.text(.5, .805, f'GLP-1R screen: {d["library"]:,} compounds  →  {d["final"]:,} ranked'
                       f'  in  {d["hours"]:.2f} h',
             ha='center', va='center', fontsize=5.9, color=mid)

    # --- retrieval on every system tested -----------------------------------
    ax = fig.add_axes([0.115, 0.325, 0.87, 0.415])
    labels = list(SYSTEMS.values())
    x = np.arange(len(labels))
    width = 0.37
    for j, (method, _) in enumerate(ARMS):
        vals = [d['roc'][k][j] for k in SYSTEMS]
        ax.bar(x + (j - 0.5) * width, vals, width=width,
               color=green if j == 0 else grey,
               hatch='' if j == 0 else '///',
               edgecolor='white', linewidth=0.45, zorder=2)
        for xi, v in zip(x + (j - 0.5) * width, vals):
            ax.text(xi, v + .02, f'{v:.2f}', ha='center', va='bottom',
                    fontsize=4.9, color=green if j == 0 else mid,
                    weight='bold' if j == 0 else 'normal')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.6)
    ax.set_ylim(0, 1.16)
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(axis='y', labelsize=5.2, length=2, pad=1)
    ax.tick_params(axis='x', length=0, pad=1.5)
    ax.set_ylabel('ROC-AUC', fontsize=5.8, labelpad=1.5)
    ax.grid(axis='x', visible=False)
    ax.grid(axis='y', linewidth=0.35, alpha=0.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=green if j == 0 else grey,
                             hatch='' if j == 0 else '///', edgecolor='white',
                             linewidth=0.45, label=lab)
               for j, (_, lab) in enumerate(ARMS)]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.5, 0.185),
               ncol=2, fontsize=5.3, frameon=False, handlelength=1.1,
               handletextpad=0.4, columnspacing=1.6)

    # --- the one-line reading -----------------------------------------------
    fig.text(.5, .055,
             f'The peptide-contact score wins on {d["wins"]} of {len(SYSTEMS)} interfaces; '
             f'the staged cascade supplies the throughput',
             ha='center', va='center', fontsize=5.6, color=ink)

    # Passing bbox_inches=None alone still consults the global "tight" setting.
    # Disable that setting explicitly to preserve the journal's exact box size.
    with plt.rc_context({'savefig.bbox': None, 'savefig.pad_inches': 0}):
        fig.savefig(OUT / 'toc_graphic.png', dpi=600, bbox_inches=None)
        fig.savefig(OUT / 'toc_graphic.pdf', bbox_inches=None)
    plt.close(fig)
    print('Graphical abstract: 1950 x 1050 pixels at 600 dpi; vector PDF 3.25 x 1.75 inches')


if __name__ == '__main__':
    main()
