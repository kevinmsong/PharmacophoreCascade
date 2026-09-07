"""Render completed native diagnostics at 600 dpi and as vector PDFs.

Apply presentation-only substitutions in memory to the established diagnostic
module. The executing screening engine and its pinned source files are untouched.
Numerical CSVs are regenerated in a temporary folder and compared before figures
are promoted into a completed run bundle.
"""
import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'GLP1_top_ligand_analysis/scripts/plot_top100_native_correlation.py'


def renderer():
    source = SOURCE.read_text(encoding='utf-8')
    substitutions = {
        'fig.savefig(output_path, dpi=320)':
            'fig.savefig(output_path, dpi=600)\n    fig.savefig(output_path.with_suffix(".pdf"))',
        'cmap="RdYlBu"':
            'cmap=__import__("matplotlib").colors.LinearSegmentedColormap.from_list("orange_white_blue", ["#E69F00", "#F7F7F7", "#0072B2"])',
        '"#3b6fb6"': '"#0072B2"', '"#b24a3d"': '"#D55E00"',
        '"#e4bdb4"': '"#F4D6B5"', '"#2c699a"': '"#0072B2"',
        '"#8c3c33"': '"#994400"',
    }
    for old, new in substitutions.items():
        if old not in source:
            raise ValueError(f'Diagnostic presentation source changed: {old}')
        source = source.replace(old, new)
    namespace = {'__name__': 'floor_diagnostic_renderer', '__file__': str(SOURCE)}
    exec(compile(source, str(SOURCE), 'exec'), namespace)
    from floor_diagnostic_figures import install
    install(namespace)
    return namespace['run_correlation_analysis']


def render_bundle(bundle, run):
    summary = bundle / 'analysis/ligand_best_native_mapping_summary.csv'
    reports = bundle / 'reports'
    n = len(pd.read_csv(summary))
    target=next((label for key,label in [('glp1r','GLP-1R'),('ghsr','GHSR'),('ntsr1','NTSR1'),('mdm2','MDM2-p53')]
                 if key in str(bundle).lower()),'GLP-1R')
    policy='5% + floor' if 'floor1000' in str(bundle) else '5% only' if 'percentage' in str(bundle) else 'benchmark'
    label=f'{target} | {policy} | {n:,} native-scored ligands'
    with tempfile.TemporaryDirectory(prefix='floor_diagnostics_') as folder:
        temp = Path(folder)
        result = run(summary, temp/'analysis', temp/'reports', 2000,
                     label)
        for name in ['native_correlation_metrics.csv', 'native_correlation_quartiles.csv', 'native_rank_comparison.csv']:
            assert_frame_equal(pd.read_csv(bundle/'analysis'/name), pd.read_csv(temp/'analysis'/name),
                               check_exact=False, rtol=1e-12, atol=1e-12)
        for path in (temp/'reports').glob('*'):
            if path.suffix=='.md':
                text=path.read_text(encoding='utf-8')
                text=text.replace('GLP-1','peptide-contact').replace('GLP1','peptide-contact')
                path.write_text(text,encoding='utf-8')
            (reports/path.name).write_bytes(path.read_bytes())
    manifest = {'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                'presentation_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'figure_module_sha256': hashlib.sha256((ROOT/'evidence/floor_diagnostic_figures.py').read_bytes()).hexdigest(),
                'raster_dpi': 600, 'vector_pdf': True,
                'palette': 'Okabe-Ito blue/vermillion; orange-white-blue signed rank shifts',
                'numerical_csv_comparison': 'unchanged to 1e-12', 'ligands': n}
    (reports/'figure_rendering_manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'Rendered and numerically verified: {bundle}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundles', nargs='+', type=Path)
    args = parser.parse_args()
    run = renderer()
    for bundle in args.bundles:
        render_bundle(bundle.resolve(), run)


if __name__ == '__main__':
    main()
