"""Render final LaTeX PDFs for page-by-page layout review."""
import json
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / 'ACS_Omega_resubmission'
META = ROOT / 'evidence/outputs/alert_disabled_revision'


def main():
    state = json.loads((META / 'document_build.json').read_text(encoding='utf-8'))
    rows = []
    for name in ['main', 'main_highlighted', 'supporting_information', 'response_to_reviewers']:
        assert state['latex_' + name]['status'] == 'completed'
        folder = ROOT / 'tmp/alert_disabled_revision/final_latex' / name
        folder.mkdir(parents=True, exist_ok=True)
        with fitz.open(SUB / (name + '.pdf')) as document:
            for number, page in enumerate(document, 1):
                fitz.TOOLS.store_shrink(100)  # Avoid stale pattern resources across pages/documents.
                page.get_pixmap(dpi=110, alpha=False).save(folder / f'page-{number:03d}.png')
                words = page.get_text('words')
                outside = [w[:5] for w in words if w[0] < -1 or w[1] < -1 or
                           w[2] > page.rect.width + 1 or w[3] > page.rect.height + 1]
                rows.append(dict(document=name, page=number, words=len(words), text_outside_page=outside))
            print(f'{name}.pdf: {len(document)} pages rendered', flush=True)
    (META / 'pdf_render_inventory.json').write_text(json.dumps(dict(
        pages=rows, visual_review='Separate page-by-page review required'), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
