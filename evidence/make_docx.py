#!/usr/bin/env python3
"""Build an editable Word version of the manuscript from the LaTeX source.

The Word file is generated from ``main.tex`` rather than written separately, so
there is one source of truth and the two cannot drift apart.

Pandoc cannot read ``ieeeaccess.cls``, so this script first rewrites the source
into a plain ``article`` document that pandoc understands, preserving all
content:

* IEEE title-block macros (``\\history``, ``\\doi``, ``\\corresp``,
  ``\\authorrefmark``, ``\\address``, ``\\tfootnote``, ``\\markboth``) become an
  ordinary title block.
* The class's ``\\Figure[..](..)[..]{file}{caption}`` macro becomes a standard
  ``figure`` environment pointing at the 600 dpi PNG, so figures land inline in
  the Word file where Word can display and move them.
* ``\\multirow`` cells are expanded, since pandoc's table reader does not
  support them.
* Equations are left as LaTeX; pandoc converts them to native Word (OMML)
  equations, which stay editable in Word's equation editor.
* Citations resolve through ``references.bib`` with the IEEE CSL style, giving
  numbered citations and an IEEE-formatted reference list.

Usage:
    python evidence/make_docx.py                # main manuscript
    python evidence/make_docx.py --supplementary
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUB = ROOT / "IEEE_Access_submission"


def strip_comments(text: str) -> str:
    """Drop full-line LaTeX comments (keep escaped \\% in content)."""
    out = []
    for line in text.split("\n"):
        if line.lstrip().startswith("%"):
            continue
        out.append(line)
    return "\n".join(out)


def convert_figures(text: str) -> str:
    """Turn the class's \\Figure macro into a standard figure environment.

    Figures point at the PNG rather than the PDF: Word embeds raster images
    reliably, and the PNGs are the same 600 dpi renders used for the PDF.
    """
    pat = re.compile(
        r"\\Figure\[[^\]]*\]\([^)]*\)\[[^\]]*\]\{([^}]+)\}\s*\n?\s*\{(.*?)\}\s*(?=\n\n|\n\\)",
        re.S)

    def repl(m: re.Match) -> str:
        fname, cap = m.group(1), m.group(2).strip()
        png = Path(fname).with_suffix(".png").name
        return ("\\begin{figure}[htbp]\n\\centering\n"
                f"\\includegraphics[width=\\linewidth]{{{png}}}\n"
                f"\\caption{{{cap}}}\n\\end{{figure}}\n")

    text, n = pat.subn(repl, text)
    print(f"  figures converted: {n}")
    return text


def lift_multicolumn_notes(text: str) -> str:
    """Move full-width \\multicolumn note rows out of tables into paragraphs.

    Pandoc's table reader cannot parse a spanning note row inside a longtable
    and silently drops the whole table, so the note is lifted out and emitted
    as an ordinary paragraph immediately after it.
    """
    def read_group(s: str, i: int) -> tuple[str, int]:
        """Read a balanced {...} group starting at s[i] == '{'."""
        assert s[i] == "{"
        depth, j = 0, i
        while j < len(s):
            if s[j] == "{":
                depth += 1
            elif s[j] == "}":
                depth -= 1
                if depth == 0:
                    return s[i + 1:j], j + 1
            j += 1
        raise ValueError("unbalanced braces")

    notes: list[str] = []
    out, pos = [], 0
    while True:
        k = text.find("\\multicolumn{", pos)
        if k == -1:
            out.append(text[pos:])
            break
        try:
            _span, j = read_group(text, k + len("\\multicolumn"))
            _spec, j = read_group(text, j)
            content, j = read_group(text, j)
        except (ValueError, AssertionError):
            out.append(text[pos:k + 13])
            pos = k + 13
            continue

        tail = text[j:j + 40]
        if "\\\\" not in tail:                      # not a full note row
            out.append(text[pos:j])
            pos = j
            continue

        # Drop the row, plus the \midrule that introduced it, and keep the note.
        head = text[pos:k]
        head = re.sub(r"\\midrule\s*$", "", head)
        out.append(head)
        notes.append(" ".join(content.split()))
        pos = j + tail.index("\\\\") + 2

    text = "".join(out)

    if notes:
        parts = text.split("\\end{longtable}")
        rebuilt = parts[0]
        for part in parts[1:]:
            note = f"\n\n\\noindent {notes.pop(0)}\n" if notes else ""
            rebuilt += "\\end{longtable}" + note + part
        text = rebuilt
        print(f"  multicolumn notes lifted: {len(parts) - 1}")
    return text


def expand_multirow(text: str) -> str:
    """Expand \\multirow{n}{*}{X} into a plain cell.

    Pandoc's LaTeX table reader ignores multirow spans, which would otherwise
    silently drop the row labels.
    """
    text = re.sub(r"\\multirow\{\d+\}\{\*\}\{([^}]*)\}", r"\1", text)
    # Rows that previously relied on the span now need their label repeated.
    return text


def build_main(out_tex: Path) -> None:
    src = (SUB / "main.tex").read_text(encoding="utf-8")

    title = re.search(r"\\title\{(.*?)\}\s*\n\s*\n", src, re.S).group(1)
    title = " ".join(title.split())
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", src, re.S).group(1).strip()
    keywords = re.search(r"\\begin\{keywords\}(.*?)\\end\{keywords\}", src, re.S).group(1).strip()

    body = src.split("\\maketitle", 1)[1]
    body = body.split("\\bibliographystyle", 1)[0]

    # Biographies carry author placeholders; keep them as a plain section.
    body = re.sub(r"\\begin\{IEEEbiographynophoto\}\{([^}]*)\}",
                  r"\n\\subsection*{\1}", body)
    body = body.replace("\\end{IEEEbiographynophoto}", "")

    body = body.replace("\\PARstart{P}{eptide}", "Peptide")
    body = re.sub(r"\\PARstart\{(.)\}\{([^}]*)\}", r"\1\2", body)
    body = convert_figures(body)
    body = expand_multirow(body)
    body = strip_comments(body)

    preamble = r"""\documentclass[11pt,letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{url}
\newcommand{\angstrom}{\text{\AA}}
\newcommand{\headeretal}{et al.}
\title{%s}
\author{Kevin Song, John Zhang, Lei Ye, and Jianyi Zhang\\
Department of Biomedical Engineering,
The University of Alabama at Birmingham, Birmingham, AL 35294 USA\\
Corresponding author: Jianyi Zhang (jayzhang@uab.edu)}
\date{}
\begin{document}
\maketitle

\begin{abstract}
%s
\end{abstract}

\textbf{Index terms:} %s

\vspace{1em}

\noindent\textit{Funding: This work was supported in part by the National Heart, Lung,
and Blood Institute under Grants U01HL134764, P01HL160476, R01HL131017, and
R01HL149137.}

""" % (title, abstract, keywords)

    out_tex.write_text(preamble + body + "\n\\end{document}\n", encoding="utf-8")


def build_supplementary(out_tex: Path) -> None:
    src = (SUB / "supplementary.tex").read_text(encoding="utf-8")
    body = src.split("\\begin{document}", 1)[1].split("\\end{document}", 1)[0]

    # Inline the external attrition table so pandoc sees it.
    inc = SUB / "attrition_full_table.tex"
    if inc.exists():
        body = body.replace("\\input{attrition_full_table.tex}",
                            inc.read_text(encoding="utf-8"))

    body = lift_multicolumn_notes(body)
    body = re.sub(r"\\includegraphics\[([^\]]*)\]\{([^}]+)\}",
                  lambda m: "\\includegraphics[width=\\linewidth]{%s}"
                  % Path(m.group(2)).with_suffix(".png").name, body)
    body = strip_comments(body)

    preamble = r"""\documentclass[11pt,letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{longtable}
\newcommand{\angstrom}{\text{\AA}}
\begin{document}
"""
    out_tex.write_text(preamble + body + "\n\\end{document}\n", encoding="utf-8")


def run_pandoc(tex: Path, docx: Path, with_bib: bool) -> None:
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        import pypandoc
        pandoc = pypandoc.get_pandoc_path()

    cmd = [pandoc, str(tex.name), "-f", "latex", "-o", str(docx.name),
           "--resource-path", ".", "--wrap=preserve"]
    if with_bib:
        cmd += ["--citeproc", "--bibliography=references.bib", "--csl=ieee.csl"]
    res = subprocess.run(cmd, cwd=tex.parent, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        raise SystemExit(f"pandoc failed for {tex.name}")
    if res.stderr.strip():
        for line in res.stderr.strip().split("\n")[:10]:
            print(f"  pandoc: {line}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--supplementary", action="store_true")
    args = ap.parse_args()

    if args.supplementary:
        tex = SUB / "_supplementary_docx.tex"
        build_supplementary(tex)
        run_pandoc(tex, SUB / "supplementary.docx", with_bib=False)
        print(f"wrote {SUB / 'supplementary.docx'}")
    else:
        tex = SUB / "_main_docx.tex"
        build_main(tex)
        run_pandoc(tex, SUB / "main.docx", with_bib=True)
        print(f"wrote {SUB / 'main.docx'}")


if __name__ == "__main__":
    main()
