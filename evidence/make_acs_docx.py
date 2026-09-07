#!/usr/bin/env python3
"""Build editable Word versions of the ACS Omega resubmission documents.

Each .docx is generated from the same LaTeX source that produces the .pdf, so
there is one source of truth and the two cannot drift apart.

Pandoc reads plain ``article`` LaTeX but not the ACS-flavored constructs this
submission uses, so this script first rewrites each source into a document
pandoc understands, preserving all content:

* ``\\input{...}`` files are inlined, so table sources land in the Word file.
* The ``scheme`` float becomes a figure whose caption is prefixed "Scheme 1.",
  keeping ACS's separate Scheme numbering visible in Word.
* ``\\resizebox{..}{!}{ tabular }`` is unwrapped; pandoc's table reader cannot
  see through it and would drop the table.
* ``threeparttable`` wrappers and their ``tablenotes`` are lifted out, since
  pandoc drops the enclosing table otherwise.
* ``\\multirow`` cells are expanded and full-width ``\\multicolumn`` note rows
  are moved out of the table into a following paragraph.
* ``\\new{...}`` revision markup is either unwrapped (clean copy) or turned
  into ``\\textcolor{RevBlue}{...}`` (highlighted copy), matched with a
  brace-counting scanner rather than a regex because the argument nests and
  spans paragraphs.
* Figures point at the 600 dpi PNGs, which Word embeds reliably.
* Equations stay LaTeX; pandoc converts them to native Word (OMML) equations,
  which remain editable in Word's equation editor.
* Citations resolve through ``references.bib`` with ``acs.csl``, giving
  superscript numerals and an ACS-formatted reference list.

Usage:
    python evidence/make_acs_docx.py                # every document
    python evidence/make_acs_docx.py main response  # only these
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUB = ROOT / "ACS_Omega_resubmission"

#: Colors named in the LaTeX preamble that pandoc will not see.
COLOR_DEFS = r"""
\usepackage{xcolor}
\definecolor{RevBlue}{RGB}{0,70,190}
"""


# --------------------------------------------------------------------------
# Small LaTeX manipulation helpers
# --------------------------------------------------------------------------
def strip_comments(text: str) -> str:
    """Drop full-line LaTeX comments (keep escaped \\% inside content)."""
    return "\n".join(l for l in text.split("\n") if not l.lstrip().startswith("%"))


def read_group(s: str, i: int) -> tuple[str, int]:
    """Read a balanced ``{...}`` group starting at ``s[i] == '{'``.

    Returns the group's contents and the index just past its closing brace.
    """
    if s[i] != "{":
        raise ValueError(f"expected '{{' at {i}, found {s[i]!r}")
    depth, j = 0, i
    while j < len(s):
        if s[j] == "\\":            # skip an escaped character
            j += 2
            continue
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError("unbalanced braces")


def rewrite_macro(text: str, name: str, wrap) -> str:
    """Rewrite every ``\\name{...}`` call, honoring nested braces.

    ``wrap`` is a callable receiving the argument, or None to unwrap. A
    callable rather than a format string, because LaTeX replacements are full
    of literal braces that ``str.format`` would read as fields.
    """
    token, out, pos = "\\" + name + "{", [], 0
    while True:
        k = text.find(token, pos)
        if k == -1:
            out.append(text[pos:])
            return "".join(out)
        out.append(text[pos:k])
        arg, j = read_group(text, k + len(token) - 1)
        arg = rewrite_macro(arg, name, wrap)      # handle nesting
        out.append(arg if wrap is None else wrap(arg))
        pos = j


def inline_inputs(text: str, base: Path, depth: int = 0) -> str:
    """Replace ``\\input{file}`` with the file's contents, recursively."""
    if depth > 5:
        return text

    def repl(m: re.Match) -> str:
        name = m.group(1).strip()
        p = base / name
        if not p.suffix:
            p = p.with_suffix(".tex")
        if not p.exists():
            print(f"  warning: \\input{{{name}}} not found, dropped")
            return ""
        return inline_inputs(p.read_text(encoding="utf-8"), base, depth + 1)

    return re.sub(r"\\input\{([^}]+)\}", repl, text)


def convert_scheme(text: str) -> str:
    """Turn the ACS ``scheme`` float into a figure with a "Scheme n." caption."""
    n = 0

    def repl(m: re.Match) -> str:
        nonlocal n
        n += 1
        body = m.group(1)
        body = re.sub(r"\\caption\{", f"\\\\caption{{Scheme {n}. ", body, count=1)
        return "\\begin{figure}[htbp]" + body + "\\end{figure}"

    text, k = re.subn(r"\\begin\{scheme\}(?:\[[^\]]*\])?(.*?)\\end\{scheme\}",
                      repl, text, flags=re.S)
    print(f"  schemes converted: {k}")
    return text


def unwrap_resizebox(text: str) -> str:
    """Drop ``\\resizebox{w}{h}{...}`` wrappers, keeping their contents.

    Pandoc's table reader cannot see a tabular nested inside resizebox and
    silently drops the whole table.
    """
    token, out, pos, n = "\\resizebox", [], 0, 0
    while True:
        k = text.find(token, pos)
        if k == -1:
            out.append(text[pos:])
            break
        j = k + len(token)
        try:
            while j < len(text) and text[j] in " \n":
                j += 1
            _w, j = read_group(text, j)
            while j < len(text) and text[j] in " \n":
                j += 1
            _h, j = read_group(text, j)
            while j < len(text) and text[j] in " \n%":
                j += 1
            body, j = read_group(text, j)
        except (ValueError, IndexError):
            out.append(text[pos:k + len(token)])
            pos = k + len(token)
            continue
        out.append(text[pos:k])
        out.append(body)
        pos = j
        n += 1
    if n:
        print(f"  resizebox unwrapped: {n}")
    return "".join(out)


def flatten_threeparttable(text: str) -> str:
    """Remove threeparttable wrappers and lift tablenotes into paragraphs."""
    notes: list[str] = []

    def grab(m: re.Match) -> str:
        body = re.sub(r"\\item\s*", " ", m.group(1))
        notes.append(" ".join(body.split()))
        return ""

    text = re.sub(r"\\begin\{tablenotes\}(?:\[[^\]]*\])?(.*?)\\end\{tablenotes\}",
                  grab, text, flags=re.S)
    text = text.replace("\\begin{threeparttable}", "")
    text = text.replace("\\end{threeparttable}", "")
    if notes:
        # Append each lifted note after the table it belonged to.
        parts = text.split("\\end{table}")
        rebuilt = parts[0]
        for part in parts[1:]:
            note = f"\n\n\\noindent {notes.pop(0)}\n" if notes else ""
            rebuilt += "\\end{table}" + note + part
        text = rebuilt
        print(f"  tablenotes lifted: {len(notes) if notes else 'all'}")
    return text


def lift_multicolumn_notes(text: str) -> str:
    """Move full-width ``\\multicolumn`` note rows out of tables.

    Pandoc's table reader cannot parse a spanning note row and silently drops
    the whole table, so the note is emitted as an ordinary paragraph after it.
    """
    notes, out, pos = [], [], 0
    while True:
        k = text.find("\\multicolumn{", pos)
        if k == -1:
            out.append(text[pos:])
            break
        try:
            _span, j = read_group(text, k + len("\\multicolumn"))
            _spec, j = read_group(text, j)
            content, j = read_group(text, j)
        except (ValueError, IndexError):
            out.append(text[pos:k + 13])
            pos = k + 13
            continue

        tail = text[j:j + 40]
        if "\\\\" not in tail or len(content) < 60:   # a label, not a note row
            out.append(text[pos:j])
            pos = j
            continue

        head = re.sub(r"\\midrule\s*$", "", text[pos:k])
        out.append(head)
        notes.append(" ".join(content.split()))
        pos = j + tail.index("\\\\") + 2

    text = "".join(out)
    if notes:
        for env in ("\\end{longtable}", "\\end{tabular}"):
            if not notes:
                break
            parts = text.split(env)
            rebuilt = parts[0]
            for part in parts[1:]:
                note = f"\n\n\\noindent {notes.pop(0)}\n" if notes else ""
                rebuilt += env + note + part
            text = rebuilt
        print("  multicolumn notes lifted")
    return text


def expand_multirow(text: str) -> str:
    """Expand ``\\multirow{n}{*}{X}`` into a plain cell.

    Pandoc's table reader ignores multirow spans, which would otherwise drop
    the row labels.
    """
    return re.sub(r"\\multirow\{[^}]*\}\{[^}]*\}\{([^{}]*)\}", r"\1", text)


def expand_multicolumn(text: str) -> str:
    """Expand a surviving ``\\multicolumn{n}{spec}{X}`` into X plus n-1 blanks.

    Pandoc's table reader drops an entire table when it meets a column span in
    the body, so every table needs to reach it rectangular. Note rows have
    already been lifted out by ``lift_multicolumn_notes``; what is left here is
    a genuine spanning cell, and padding it with empty cells preserves both the
    text and the column count.
    """
    out, pos, n = [], 0, 0
    while True:
        k = text.find("\\multicolumn{", pos)
        if k == -1:
            out.append(text[pos:])
            break
        try:
            span, j = read_group(text, k + len("\\multicolumn"))
            _spec, j = read_group(text, j)
            content, j = read_group(text, j)
            width = int(span.strip())
        except (ValueError, IndexError):
            out.append(text[pos:k + 13])
            pos = k + 13
            continue
        out.append(text[pos:k])
        out.append(content + " &" * (width - 1))
        pos = j
        n += 1
    if n:
        print(f"  multicolumn spans expanded: {n}")
    return "".join(out)


def drop_partial_rules(text: str) -> str:
    """Remove ``\\cmidrule`` and friends, which pandoc has no equivalent for."""
    return re.sub(r"\\cmidrule(?:\([^)]*\))?\{[^}]*\}", "", text)


def figures_to_png(text: str) -> str:
    """Point every \\includegraphics at the PNG render of the same figure."""
    def repl(m: re.Match) -> str:
        name = Path(m.group(2)).with_suffix(".png").name
        return "\\includegraphics[width=\\linewidth]{%s}" % name
    return re.sub(r"\\includegraphics(?:\[([^\]]*)\])?\{([^}]+)\}", repl, text)


#: Sentinels bracketing revised text on its way through pandoc.
#:
#: Pandoc's docx writer discards \textcolor entirely (it survives as a style
#: attribute the writer has no mapping for), so the color has to be applied to
#: the finished .docx. These markers are plain ASCII, contiguous, and absent
#: from the manuscript, so they reach the Word runs intact and can be found
#: again afterwards.
MARK_START = "@@NEWSTART@@"
MARK_END = "@@NEWEND@@"
REV_RGB = "0046BE"


def apply_revision_markup(text: str, highlighted: bool) -> str:
    """Resolve \\new{...} the same way the corresponding LaTeX build does."""
    wrap = (lambda a: MARK_START + a + MARK_END) if highlighted else None
    return rewrite_macro(text, "new", wrap)


def _iter_paragraphs(doc):
    """Every paragraph in the document, including those inside tables."""
    yield from doc.paragraphs
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                yield from cell.paragraphs
                for inner in cell.tables:
                    for r in inner.rows:
                        for c in r.cells:
                            yield from c.paragraphs


def colorize_marked_runs(path: Path) -> int:
    """Color the text between the sentinels blue, then delete the sentinels.

    Runs are split rather than rewritten, so bold, italic, and superscript
    formatting inside a revised passage is preserved. The inside/outside state
    carries across runs and paragraphs, because a marked passage can span
    both.
    """
    import copy
    from docx import Document
    from docx.oxml.ns import qn

    doc = Document(str(path))
    inside = False
    colored = 0

    for para in _iter_paragraphs(doc):
        for run in list(para.runs):
            text = run.text
            if not text:
                continue
            if MARK_START not in text and MARK_END not in text:
                if inside and text.strip():
                    _set_color(run, REV_RGB, qn)
                    colored += len(text)
                continue

            # Split this run's text into (segment, is_revised) pieces.
            segments, buf = [], ""
            i = 0
            while i < len(text):
                if text.startswith(MARK_START, i):
                    segments.append((buf, inside)); buf = ""
                    inside = True
                    i += len(MARK_START)
                elif text.startswith(MARK_END, i):
                    segments.append((buf, inside)); buf = ""
                    inside = False
                    i += len(MARK_END)
                else:
                    buf += text[i]
                    i += 1
            segments.append((buf, inside))
            segments = [s for s in segments if s[0]]

            anchor = run._element
            for seg_text, seg_rev in segments:
                new_r = copy.deepcopy(run._element)
                for t in new_r.findall(qn("w:t")):
                    new_r.remove(t)
                t = new_r.makeelement(qn("w:t"), {})
                t.text = seg_text
                t.set(qn("xml:space"), "preserve")
                new_r.append(t)
                if seg_rev:
                    _set_color_element(new_r, REV_RGB, qn)
                    colored += len(seg_text)
                anchor.addnext(new_r)
                anchor = new_r
            run._element.getparent().remove(run._element)

    # Treat refreshed tables as complete revision units, including retained
    # headings. This matches the review-copy LaTeX table coloring.
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    for run in para.runs:
                        _set_color(run, REV_RGB, qn)
                        colored += len(run.text)
    for para in doc.paragraphs:
        if re.match(r'Table\s+\d+[.:]', para.text):
            for run in para.runs:
                _set_color(run, REV_RGB, qn)
    doc.save(str(path))
    return colored


def _set_color(run, rgb: str, qn) -> None:
    _set_color_element(run._element, rgb, qn)


def _set_color_element(r_element, rgb: str, qn) -> None:
    """Set <w:color> on a run element, creating <w:rPr> if needed."""
    rPr = r_element.find(qn("w:rPr"))
    if rPr is None:
        rPr = r_element.makeelement(qn("w:rPr"), {})
        r_element.insert(0, rPr)
    for old in rPr.findall(qn("w:color")):
        rPr.remove(old)
    color = rPr.makeelement(qn("w:color"), {qn("w:val"): rgb})
    rPr.insert(0, color)


# --------------------------------------------------------------------------
# Document builders
# --------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[11pt,letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{multirow}
\usepackage{url}
""" + COLOR_DEFS + r"""\newcommand{\angstrom}{\text{\AA}}
\newcommand{\figref}[1]{Figure~\ref{#1}}
\newcommand{\tabref}[1]{Table~\ref{#1}}
\newcommand{\schemeref}[1]{Scheme~\ref{#1}}
\newcommand{\eqnref}[1]{eq~\ref{#1}}
\begin{document}
"""


def resolve_document_numbers(body: str, aux_name: str) -> str:
    """Carry LaTeX's actual float/equation numbering into editable Word output."""
    aux = (SUB / aux_name).read_text(encoding="utf-8")
    numbers = dict(re.findall(r"\\newlabel\{([^}]+)\}\{\{([^}]+)\}", aux))

    def caption(match):
        block = match.group(0)
        label = re.search(r"\\label\{([^}]+)\}", block)
        cap = re.search(r"\\caption(?:\[[^\]]*\])?\{", block)
        if not label or not cap:
            return block
        key = label.group(1)
        if key not in numbers:
            raise ValueError(f"Unresolved caption {key}: build LaTeX first")
        kind = "Scheme" if key.startswith("sch:") else "Figure" if match.group(1)=="figure" else "Table"
        arg, end = read_group(block, cap.end()-1)
        if not re.match(r"(?:Scheme|Figure|Table)\s+\d", arg):
            arg = f"{kind} {numbers[key]}. " + arg
        return block[:cap.end()] + arg + block[end-1:]

    body = re.sub(r"\\begin\{(figure|table|longtable)\}.*?\\end\{\1\}", caption, body, flags=re.S)

    def equation(match):
        block = match.group(1)
        label = re.search(r"\\label\{([^}]+)\}", block)
        if not label: return match.group(0)
        number = numbers[label.group(1)]
        block = block[:label.start()] + block[label.end():]
        return r"\begin{equation}" + block + r"\qquad\text{(" + number + r")}\end{equation}"

    body = re.sub(r"\\begin\{equation\}(.*?)\\end\{equation\}", equation, body, flags=re.S)

    def reference(match):
        macro,key=match.groups()
        if key not in numbers:
            raise ValueError(f"Unresolved reference {key}: build LaTeX first")
        number=numbers[key]
        prefixes={"figref":"Figure~","tabref":"Table~","schemeref":"Scheme~","eqnref":"eq~"}
        return "("+number+")" if macro=="eqref" else prefixes.get(macro,"")+number

    return re.sub(r"\\(figref|tabref|schemeref|eqnref|eqref|ref)\{([^}]+)\}",reference,body)


def common_cleanup(body: str, highlighted: bool, aux_name: str = "main.aux") -> str:
    """Transformations every document needs, in dependency order."""
    body = inline_inputs(body, SUB)
    # Pandoc consumes numeric text after LaTeX array column modifiers.
    # Word column widths/alignment are assigned by the document formatter.
    body = re.sub(r'[<>]\{\\(?:raggedleft|raggedright|centering)\\arraybackslash\}', '', body)
    # A Word table repeats its own header; do not emit LaTeX continuation
    # headers as a second data row.
    body = re.sub(r'\\endfirsthead.*?\\endhead', r'\\endhead', body, flags=re.S)
    body = body.replace(
        r' & & \multicolumn{2}{c}{5\% only} & \multicolumn{2}{c}{5\% + floor} \\' + '\n' +
        r'System & Stage & Total & Actives & Total & Actives \\',
        r'System & Stage & 5\% total & 5\% actives & Floor total & Floor actives \\')
    body = re.sub(r"\\angstrom\b", "Å", body)
    body = re.sub(r"\{\\rm\s+([^{}]+)\}", lambda m: r"{\mathrm{" + m.group(1).strip() + "}}", body)
    body = apply_revision_markup(body, highlighted)
    body = convert_scheme(body)
    body = unwrap_resizebox(body)
    body = flatten_threeparttable(body)
    body = lift_multicolumn_notes(body)      # long note rows leave the table
    body = expand_multicolumn(body)          # remaining spans become padded cells
    body = expand_multirow(body)
    body = drop_partial_rules(body)
    body = figures_to_png(body)
    body = strip_comments(body)
    body = resolve_document_numbers(body, aux_name)
    # Environments and switches with no meaning in Word.
    for env in ("singlespace",):
        body = body.replace(f"\\begin{{{env}}}", "").replace(f"\\end{{{env}}}", "")
    for cmd in (r"\doublespacing", r"\FloatBarrier", r"\sloppy", r"\maketitle",
                r"\printbibliography", r"\newpage", r"\clearpage"):
        body = body.replace(cmd, "")
    return body


def build_manuscript(out_tex: Path, highlighted: bool) -> None:
    src = (SUB / "manuscript_body.tex").read_text(encoding="utf-8")
    pre = (SUB / "acs_preamble.tex").read_text(encoding="utf-8")

    title = " ".join(re.search(r"\\title\{(.*?)\}\s*\n", pre, re.S).group(1).split())
    if highlighted:
        title = MARK_START + title + MARK_END
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", src, re.S).group(1)
    abstract = apply_revision_markup(abstract.strip().replace("\\noindent", ""),
                                     highlighted)

    body = src.split("\\end{abstract}", 1)[1]
    body = common_cleanup(body, highlighted)

    banner = ""
    if highlighted:
        banner = (r"\textbf{Supporting Information for Review Only.} "
                  r"Text shown in " + MARK_START + "blue" + MARK_END +
                  r" is new or rewritten relative to the manuscript originally "
                  r"submitted to ACS Omega. Text in black is carried over "
                  r"substantially unchanged. Tables are marked blue as complete "
                  r"revised units, including retained headings. Figure graphics "
                  r"retain the publication palettes. The clean and highlighted "
                  r"copies contain the same scientific text and results."
                  "\n\n\\vspace{1em}\n\n")

    head = PREAMBLE + banner + r"""
\begin{center}
{\LARGE\textbf{%s}}

\vspace{1em}
Kevin Song, John Zhang, Lei Ye, Jianyi Zhang*

\textit{Department of Biomedical Engineering,
The University of Alabama at Birmingham}

*Email: jayzhang@uab.edu
\end{center}

\vspace{1em}

\subsection*{Abstract}
%s

\vspace{1em}
""" % (title, abstract)

    out_tex.write_text(head + body + "\n\\end{document}\n", encoding="utf-8")


def build_supplementary(out_tex: Path) -> None:
    src = (SUB / "supporting_information.tex").read_text(encoding="utf-8")
    body = src.split("\\begin{document}", 1)[1].split("\\end{document}", 1)[0]
    body = common_cleanup(body, highlighted=False, aux_name="supporting_information.aux")
    for cmd in (r"\setcounter{page}{1}", r"\setcounter{figure}{0}",
                r"\setcounter{table}{0}"):
        body = body.replace(cmd, "")
    out_tex.write_text(PREAMBLE + body + "\n\\end{document}\n", encoding="utf-8")


def build_response(out_tex: Path) -> None:
    src = (SUB / "response_to_reviewers.tex").read_text(encoding="utf-8")
    body = src.split("\\begin{document}", 1)[1].split("\\end{document}", 1)[0]
    body = strip_comments(body)

    # Reviewer quotes become italic block quotes; our own macros become plain
    # markup. Manuscript quotes stay roman so the two remain distinguishable in
    # Word, where the LaTeX shading and rule do not survive the conversion.
    body = body.replace("\\begin{reviewerquote}", "\\begin{quote}\\itshape")
    body = body.replace("\\end{reviewerquote}", "\\end{quote}")
    body = body.replace("\\begin{manuscriptquote}", "\\begin{quote}")
    body = body.replace("\\end{manuscriptquote}", "\\end{quote}")
    body = rewrite_macro(body, "comment", lambda a: "\\subsection*{" + a + "}")
    # Without this the quote runs straight on from the preceding paragraph and
    # the reader cannot tell where our reply ends and the manuscript begins.
    body = rewrite_macro(
        body, "quoted",
        lambda a: "\n\n\\textbf{Revised text} (" + a + "):\n")
    body = body.replace("\\response", "\n\n\\textbf{Response.} ")
    body = body.replace("\\changes", "\n\n\\textbf{Changes to the manuscript.} ")
    body = body.replace("\\hrule", "")
    for cmd in (r"\clearpage", r"\newpage"):
        body = body.replace(cmd, "")

    out_tex.write_text(PREAMBLE + body + "\n\\end{document}\n", encoding="utf-8")


# --------------------------------------------------------------------------
# Pandoc driver
# --------------------------------------------------------------------------
def run_pandoc(tex: Path, docx: Path, with_bib: bool) -> None:
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        import pypandoc
        pandoc = pypandoc.get_pandoc_path()

    cmd = [pandoc, tex.name, "-f", "latex", "-o", docx.name,
           "--resource-path", ".", "--wrap=preserve"]
    if with_bib:
        csl = "acs.csl" if (SUB / "acs.csl").exists() else "../IEEE_Access_submission/ieee.csl"
        cmd += ["--citeproc", "--bibliography=references.bib", f"--csl={csl}"]

    res = subprocess.run(cmd, cwd=tex.parent, capture_output=True,
                         text=True, encoding="utf-8", errors="replace")
    if res.returncode != 0:
        print(res.stdout or "")
        print(res.stderr or "", file=sys.stderr)
        raise SystemExit(f"pandoc failed for {tex.name}")
    for line in (res.stderr or "").strip().split("\n")[:8]:
        if line.strip():
            print(f"  pandoc: {line}")


#: name -> (builder, output file, resolve citations, apply revision color)
DOCS = {
    "main":          (lambda p: build_manuscript(p, highlighted=False),
                      "main.docx", True, False),
    "highlighted":   (lambda p: build_manuscript(p, highlighted=True),
                      "main_highlighted.docx", True, True),
    "supplementary": (build_supplementary, "supporting_information.docx", False, False),
    "response":      (build_response, "response_to_reviewers.docx", False, False),
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("docs", nargs="*", choices=list(DOCS),
                    help="documents to build (default: all)")
    args = ap.parse_args()

    for name in (args.docs or list(DOCS)):
        builder, out_name, with_bib, colorize = DOCS[name]
        print(f"{name}:")
        tex = SUB / f"_{name}_docx.tex"
        builder(tex)
        out = SUB / out_name
        run_pandoc(tex, out, with_bib)
        if colorize:
            n = colorize_marked_runs(out)
            print(f"  revised text colored: {n:,} characters")
        from format_acs_docx import format_document
        format_document(out)
        print(f"  wrote {out_name} ({out.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
