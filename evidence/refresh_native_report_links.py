"""Repair native diagnostic report links and distinguish coverage-rank diagnostics."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
META=ROOT/'evidence/outputs/alert_disabled_revision'


def main():
    queue=json.loads((META/'native_figure_queue.json').read_text())
    audited=[]
    for relative in queue['rendered_or_verified_bundles']:
        bundle=ROOT/relative
        report=bundle/'reports/native_correlation_analysis.md'
        manifest=bundle/'manifests/analysis_manifest.json'
        assert report.exists() and manifest.exists()
        text=report.read_text(encoding='utf-8')
        text=text.replace('See also [methodology](methodology.md) for the exact native 3D scoring procedure.',
            'See the [execution manifest](../manifests/analysis_manifest.json) for the native configuration, selection settings, and final sorting contract.')
        marker='## Scope of the diagnostic ranks'
        if marker not in text:
            text+='\n'+marker+'\n\nThese descriptive plots compare coverage-only rankings. They break equal coverage values by ligand ID and therefore differ from the full production sorting contract, which also uses fit geometry and other criteria. They do not establish active enrichment. The current retrospective enrichment statistics instead retain equal-status/equal-score ties and average early-retrieval metrics over within-tie orders.\n'
        report.write_text(text,encoding='utf-8')
        assert '(methodology.md)' not in text
        audited.append(relative)
    (META/'native_report_link_audit.json').write_text(json.dumps(dict(
        status='completed' if queue['status']=='completed' else 'partial',
        reports=audited, remaining_experiments=queue['pending_experiments']),indent=2))
    print(f'Repaired and clarified {len(audited)} native reports; figure queue is {queue["status"]}')


if __name__=='__main__':main()
