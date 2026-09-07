"""Package the visually approved current journal documents and figure sources."""
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / 'ACS_Omega_resubmission'
META = ROOT / 'evidence/outputs/alert_disabled_revision'


def main():
    qa = json.loads((META / 'final_visual_qa.json').read_text(encoding='utf-8'))
    assert qa['status'] == 'passed', 'Final page-by-page visual review is incomplete'
    for document in qa['documents']:
        path = ROOT / document['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == document['sha256'], f'Document changed after visual review: {path}'
    for name in ['final_embedded_image_audit.json', 'final_output_figure_resolution.json']:
        assert json.loads((META / name).read_text(encoding='utf-8'))['status'] == 'passed'
    archive = ROOT / 'ACS_Omega_resubmission.zip'
    if archive.exists():
        backup = ROOT / 'tmp/alert_disabled_revision/documentation_before/ACS_Omega_resubmission.zip'
        backup.parent.mkdir(parents=True, exist_ok=True)
        if not backup.exists():
            shutil.copy2(archive, backup)
    suffixes = {'.docx', '.pdf', '.tex', '.bib', '.csl', '.png', '.svg', '.txt', '.md'}
    files = sorted(p for p in SUB.iterdir() if p.is_file() and p.suffix.lower() in suffixes
                   and not p.name.startswith('_'))
    required = [SUB / (stem + ext) for stem in
                ['main', 'main_highlighted', 'supporting_information', 'response_to_reviewers']
                for ext in ['.docx', '.pdf', '.tex']]
    assert all(p in files for p in required)
    assert SUB / 'abstract.txt' in files
    rows = []
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
        for path in files:
            relative = path.relative_to(ROOT).as_posix()
            bundle.write(path, relative)
            rows.append(dict(path=relative, bytes=path.stat().st_size,
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    with zipfile.ZipFile(archive) as bundle:
        assert bundle.testzip() is None
    record = dict(status='completed', created_utc=datetime.now(timezone.utc).isoformat(),
                  archive=archive.name, bytes=archive.stat().st_size,
                  sha256=hashlib.sha256(archive.read_bytes()).hexdigest(), files=rows,
                  scope='Current journal documents, figure assets, and LaTeX sources. Scientific run data remain in the project output directories.')
    (META / 'submission_package.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(f'Packaged and verified {len(files)} current submission files: {archive}')


if __name__ == '__main__':
    main()
