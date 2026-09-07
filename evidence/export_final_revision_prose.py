"""Export final abstract and reviewer-response prose after the Word build."""
import json
import shutil
from pathlib import Path

import pypandoc
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / 'ACS_Omega_resubmission'
META = ROOT / 'evidence/outputs/alert_disabled_revision'


def replace(path, text):
    backup = ROOT / 'tmp/alert_disabled_revision/documentation_before' / path.relative_to(ROOT)
    if path.exists() and not backup.exists():
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, backup)
    path.write_text(text, encoding='utf-8')


def main():
    state = json.loads((META / 'document_build.json').read_text(encoding='utf-8'))
    assert state['word_build']['status'] == 'completed'
    paragraphs = Document(SUB / 'main.docx').paragraphs
    start = next(i for i, p in enumerate(paragraphs) if p.text.strip() == 'Abstract') + 1
    content = []
    for paragraph in paragraphs[start:]:
        if paragraph.style.name.startswith('Heading'):
            break
        if paragraph.text.strip():
            content.append(paragraph.text.strip())
    abstract = '\n\n'.join(content) + '\n'
    assert '1,000' in abstract and '158/160' in abstract
    assert 150 < len(abstract.split()) < 400
    replace(SUB / 'abstract.txt', abstract)
    response = pypandoc.convert_file(str(SUB / '_response_docx.tex'), 'gfm', format='latex',
                                     extra_args=['--wrap=none'])
    assert len(response.split()) > 1000
    replace(ROOT / 'author_response.md', response)
    print(f'Exported {len(abstract.split())}-word abstract and current reviewer-response Markdown.')


if __name__ == '__main__':
    main()
