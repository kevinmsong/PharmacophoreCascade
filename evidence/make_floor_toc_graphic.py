"""Graphical abstract for ACS Omega, in the journal's 3.25 x 1.75-inch TOC box.

The graphic answers the question the paper answers: the terminal native
peptide-contact score supplies the enrichment, and the staged cascade upstream
supplies tractability at library scale. Configuration detail (alert policy,
shortlist rule, capping behaviour) belongs in the Methods, not here.

Every number is read from the recorded outputs rather than typed in, so the
graphic cannot drift away from the tables it summarizes.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd

import figstyle as fs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'ACS_Omega_resubmission'
FLOOR = ROOT / 'evidence/outputs/absolute_floor'
BENCH = ROOT / 'evidence/outputs'
HEADLINE = ROOT / 'results/screening_full_1M_topological_hashed_run_summary.json'

#: Benchmark directories, in the order the manuscript reports them.
SYSTEMS = {'glp1r': 'GLP-1R', 'ghsr': 'GHSR', 'ntsr1': 'NTSR1', 'mdm2': 'MDM2-p53'}


def measured():
    """Collect the four quantities the graphic states, from the recorded runs."""
    run = json.loads(HEADLINE.read_text(encoding='utf-8'))
    counts, timings = run['counts'], run['timings']

    # Native-only beats the conventional single-pass 3D pharmacophore on N of 4.
    wins = 0
    for key in SYSTEMS:
        summary = pd.read_csv(BENCH / f'benchmark_{key}_full/benchmark_summary.csv',
                              index_col='method')['roc_auc']
        wins += summary['native_only'] > summary['standard_3d_pharmacophore']

    # Shortlist depth against final active retention, from the paired runs.
    policy = pd.read_csv(FLOOR / 'policy_comparison.csv')
    retention = {}
    for key in ['ghsr', 'ntsr1']:
        arm = policy[policy.system == key].set_index('policy')
        retention[SYSTEMS[key]] = (int(arm.loc['percentage', 'final_ranked_actives']),
                                   int(arm.loc['floor1000', 'final_ranked_actives']),
                                   int(arm.loc['percentage', 'n_actives']))
    return dict(
        library=int(counts['total_scanned']),
        final=int(run['native_rerank']['counts']['final_rows']),
        hours=timings['total_pipeline_sec'] / 3600.0,
        wins=int(wins),
        retention=retention,
    )


def main():
    d = measured()
    fs.apply()
    blue, green = fs.OKABE_ITO['blue'], fs.OKABE_ITO['green']
    ink, grey = '#1A1A1A', '#4D4D4D'

    fig, ax = plt.subplots(figsize=(3.25, 1.75))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')

    # --- the question, and the scale it was answered at ---------------------
    ax.text(.5, .965, 'Where does the enrichment come from?', ha='center', va='top',
            fontsize=8.8, weight='bold', color=ink)
    ax.text(.5, .800, f'GLP-1R  ·  {d["library"]:,} compounds  →  {d["final"]:,} ranked'
                      f'  ·  {d["hours"]:.2f} h',
            ha='center', va='center', fontsize=6.1, color=grey)

    # --- the two halves of the answer ---------------------------------------
    panels = [(.020, .430, '#EAF3F8', blue, 'Staged cascade', 'makes the scale tractable'),
              (.545, .435, '#E8F4EF', green, 'Native peptide-contact score',
               'supplies the enrichment')]
    for x, w, fill, edge, head, tail in panels:
        ax.add_patch(FancyBboxPatch((x, .478), w, .215,
                                    boxstyle='round,pad=.012,rounding_size=.03',
                                    facecolor=fill, edgecolor=edge, lw=.9))
        ax.text(x + w / 2, .630, head, ha='center', va='center',
                fontsize=6.8, weight='bold', color=edge)
        ax.text(x + w / 2, .535, tail, ha='center', va='center',
                fontsize=6.3, color=ink)

    ax.add_patch(FancyArrowPatch((.462, .5855), (.532, .5855), arrowstyle='-|>',
                                 mutation_scale=8, lw=1.0, color=ink))

    ax.text(.5, .394, f'beats a single-pass 3D pharmacophore on {d["wins"]} of 4 systems',
            ha='center', va='center', fontsize=6.5, weight='bold', color=ink)

    # --- and what the tractability costs ------------------------------------
    ax.add_patch(FancyBboxPatch((.020, .075), .960, .200,
                                boxstyle='round,pad=.010,rounding_size=.025',
                                facecolor='#FFFFFF', edgecolor='#C9C9C9', lw=.7))
    ax.text(.5, .224, 'Shortlist depth sets how many actives survive',
            ha='center', va='center', fontsize=6.3, color=ink)
    pairs = '        '.join(f'{name}  {a}/{n} → {b}/{n}'
                            for name, (a, b, n) in d['retention'].items())
    ax.text(.5, .124, pairs, ha='center', va='center',
            fontsize=6.5, weight='bold', color=blue)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    # Passing bbox_inches=None alone still consults the global "tight" setting.
    # Disable that setting explicitly to preserve the journal's exact box size.
    with plt.rc_context({'savefig.bbox': None, 'savefig.pad_inches': 0}):
        fig.savefig(OUT / 'toc_graphic.png', dpi=600, bbox_inches=None)
        fig.savefig(OUT / 'toc_graphic.pdf', bbox_inches=None)
    plt.close(fig)
    print('Graphical abstract: 1950 x 1050 pixels at 600 dpi; vector PDF 3.25 x 1.75 inches')


if __name__ == '__main__':
    main()
