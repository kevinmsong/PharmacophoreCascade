"""Build final document artifacts only after the fresh scientific build succeeds.

Rendering is preparation for human/agent visual QA, not a claim of visual approval.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / 'ACS_Omega_resubmission'
META = ROOT / 'evidence/outputs/alert_disabled_revision'
DOCS = ['main', 'main_highlighted', 'supporting_information', 'response_to_reviewers']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wait', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    # MiKTeX's latexmk launcher requires Perl. Use the existing Git runtime when
    # the workstation's default PATH does not expose it.
    if shutil.which('perl') is None:
        git_perl=Path('C:/Program Files/Git/usr/bin/perl.exe')
        if git_perl.exists():
            os.environ['PATH']=str(git_perl.parent)+os.pathsep+os.environ.get('PATH','')
    (META/'document_tool_paths.json').write_text(json.dumps(
        {name:shutil.which(name) for name in ['latexmk','pdflatex','biber','perl']},indent=2),encoding='utf-8')
    while True:
        source = META / 'deliverable_build.json'
        build = json.loads(source.read_text()) if source.exists() else {}
        failed = [k for k, v in build.items() if v['status'] == 'failed']
        if failed:
            raise RuntimeError(f'Scientific build failed: {failed}')
        if build.get('figure_accessibility', {}).get('status') == 'completed':
            break
        if not args.wait:
            raise RuntimeError('Scientific build is incomplete')
        time.sleep(30)
    py = sys.executable
    phases = [
        ('metric_assertions', [py, 'evidence/verify_current_metrics.py'], ROOT),
        ('reviewer_coverage', [py, 'evidence/check_reviewer_coverage.py', '--dir', str(SUB)], ROOT),
        ('documentation', [py, 'evidence/write_alert_disabled_readmes.py'], ROOT),
    ]
    phases += [(f'latex_{name}', ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', name + '.tex'], SUB) for name in DOCS]
    phases += [
        ('word_build', [py, 'evidence/make_acs_docx.py'], ROOT),
        ('crossrefs', [py, 'evidence/check_response_crossrefs.py', '--dir', str(SUB)], ROOT),
        ('word_render', [py, 'evidence/render_word_documents.py', *[str(SUB / (name + '.docx')) for name in DOCS],
                         '--output-dir', str(ROOT / 'tmp/alert_disabled_revision/final_word')], ROOT),
    ]
    statepath = META / 'document_build.json'
    states = json.loads(statepath.read_text()) if statepath.exists() else {}
    if states and not args.resume:
        raise RuntimeError('Document build state exists; inspect and use --resume')
    for name, command, cwd in phases:
        if states.get(name, {}).get('status') == 'completed':
            continue
        print(f'Final document phase: {name}', flush=True)
        started = time.perf_counter()
        with open(META / f'document_{name}.log', 'w', encoding='utf-8') as log:
            result = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
        states[name] = dict(status='completed' if result.returncode == 0 else 'failed',
                            returncode=result.returncode, command=command,
                            elapsed_seconds=time.perf_counter()-started,
                            finished_utc=datetime.now(timezone.utc).isoformat())
        statepath.write_text(json.dumps(states, indent=2))
        if result.returncode:
            raise RuntimeError(f'{name} failed; inspect document_{name}.log')
    print('Documents rendered; EVERY page still requires visual QA before delivery.', flush=True)


if __name__ == '__main__':
    main()
