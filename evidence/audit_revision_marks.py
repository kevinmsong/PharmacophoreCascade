#!/usr/bin/env python3
"""Audit the change highlighting against the previously submitted manuscript.

The review copy marks revised passages with ``\\new{...}``. Those marks are only
trustworthy if they reflect what actually changed relative to the manuscript the
reviewers read, so this script compares the two sources sentence by sentence and
reports where the marking disagrees with the comparison:

* **under-marked** - a sentence with no close counterpart in the previous
  submission that is nevertheless left unmarked, so the review copy would show
  new text in black.
* **over-marked** - a sentence carried over nearly verbatim yet marked as new,
  which overstates the extent of the revision.

Matching is on normalized word sequences with a similarity ratio, because a
revision rewords far more often than it moves text intact; the threshold is
deliberately generous so that only near-verbatim carryover counts as unchanged.

Usage:
    python evidence/audit_revision_marks.py
    python evidence/audit_revision_marks.py --threshold 0.80 --show-all
"""
from __future__ import annotations

import argparse
import difflib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Similarity at or above which a sentence counts as carried over unchanged.
DEFAULT_THRESHOLD = 0.85
#: Sentences shorter than this (in words) are boilerplate; skip them.
MIN_WORDS = 6


def strip_latex(text: str) -> str:
    """Reduce LaTeX to readable prose, dropping floats, math, and markup."""
    text = "\n".join(l for l in text.split("\n") if not l.lstrip().startswith("%"))
    for env in ("figure", "table", "scheme", "equation", "aligned", "longtable",
                "tabular", "threeparttable", "abstract"):
        text = re.sub(r"\\begin\{" + env + r"\*?\}.*?\\end\{" + env + r"\*?\}",
                      " ", text, flags=re.S)
    text = re.sub(r"\\(cite|ref|label|tabref|figref|schemeref|eqnref|input|url)"
                  r"\{[^}]*\}", " ", text)
    text = re.sub(r"\$[^$]*\$", " NUM ", text)
    text = re.sub(r"\\(?:emph|textbf|textit|texttt|new|noindent)\s*", " ", text)
    text = re.sub(r"\\[a-zA-Z@]+\*?", " ", text)
    text = re.sub(r"[{}\\&~^_]", " ", text)
    return re.sub(r"\s+", " ", text)


def sentences(prose: str) -> list[str]:
    parts = re.split(r"(?<=[.:;])\s+(?=[A-Z(])", prose)
    return [" ".join(p.split()) for p in parts if len(p.split()) >= MIN_WORDS]


def key(sentence: str) -> list[str]:
    """Comparison form: lowercase words, punctuation and digits normalized."""
    s = sentence.lower()
    s = re.sub(r"[0-9][0-9,.]*", "#", s)
    s = re.sub(r"[^a-z#\s]", " ", s)
    return s.split()


def best_ratio(target: list[str], pool: list[list[str]],
               index: dict[str, set[int]]) -> float:
    """Highest similarity of ``target`` against any sentence in ``pool``.

    Candidates are narrowed by shared rare-ish words first; comparing every
    sentence against every other is quadratic and needlessly slow here.
    """
    candidates: set[int] = set()
    for w in set(target):
        if w in index and len(index[w]) < 400:
            candidates |= index[w]
    if not candidates:
        candidates = set(range(len(pool)))

    best = 0.0
    sm = difflib.SequenceMatcher(autojunk=False)
    sm.set_seq2(target)
    for i in candidates:
        sm.set_seq1(pool[i])
        if sm.real_quick_ratio() <= best or sm.quick_ratio() <= best:
            continue
        best = max(best, sm.ratio())
        if best >= 0.99:
            break
    return best


def marked_spans(body: str) -> list[tuple[int, int]]:
    """Character spans covered by \\new{...}, honoring nested braces."""
    spans, pos = [], 0
    token = "\\new{"
    while True:
        k = body.find(token, pos)
        if k == -1:
            return spans
        depth, j = 0, k + len(token) - 1
        while j < len(body):
            if body[j] == "\\":
                j += 2
                continue
            if body[j] == "{":
                depth += 1
            elif body[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        spans.append((k, j))
        pos = j + 1


#: Paragraphs carrying these are structural, not prose, and are left alone.
STRUCTURAL = ("\\begin{", "\\end{", "\\input{", "\\caption{", "\\label{",
              "\\section", "\\subsection", "\\includegraphics", "\\printbibliography",
              "\\maketitle", "\\newpage", "\\clearpage", "\\FloatBarrier",
              "\\doublespacing", "\\vspace", "\\smallskip", "\\centering")

#: Floats are marked through their captions, which carry their own \new{};
#: everything else that interrupts a paragraph is display math to write around.
FLOAT_ENVS = ("figure", "table", "scheme", "longtable", "tabular",
              "threeparttable", "abstract", "singlespace")

_ENV_RE = re.compile(r"\\begin\{([a-zA-Z*]+)\}.*?\\end\{\1\}", re.S)


def mark_prose_around_math(text: str, old_keys, index, threshold: float) -> tuple[str, int]:
    """Mark the changed prose runs of a paragraph that display math interrupts.

    Methods paragraphs wrap around numbered equations, so the paragraph cannot
    be marked as a unit without swallowing the equation. The equations are left
    untouched and the prose between them is marked run by run.
    """
    if any(f"\\begin{{{e}}}" in text for e in FLOAT_ENVS):
        return text, 0

    pieces, pos, changed = [], 0, 0
    for m in _ENV_RE.finditer(text):
        pieces.append(("prose", text[pos:m.start()]))
        pieces.append(("math", m.group(0)))
        pos = m.end()
    if not pieces:
        return text, 0
    pieces.append(("prose", text[pos:]))

    out = []
    for kind, seg in pieces:
        if kind == "math" or not seg.strip():
            out.append(seg)
            continue
        core = seg.strip()
        if core.startswith("\\new{") and core.endswith("}"):
            out.append(seg)
            continue
        sents = sentences(" ".join(strip_latex(core).split()))
        if not sents or all(best_ratio(key(s), old_keys, index) >= threshold
                            for s in sents):
            out.append(seg)
            continue
        lead = seg[:len(seg) - len(seg.lstrip())]
        tail = seg[len(seg.rstrip()):]
        out.append(f"{lead}\\new{{{core}}}{tail}")
        changed += 1
    return "".join(out), changed


def read_group(s: str, i: int) -> tuple[str, int]:
    """Read a balanced ``{...}`` group starting at ``s[i] == '{'``."""
    depth, j = 0, i
    while j < len(s):
        if s[j] == "\\":
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


def mark_captions(body: str, old_keys, index, threshold: float) -> tuple[str, int]:
    """Mark a whole caption whose text is not carried over from the old draft.

    Captions are marked as units rather than clause by clause: a caption
    describing a figure that has been redrawn is new in its entirety, and
    splitting it into black and blue fragments would suggest otherwise.
    """
    out, pos, n = [], 0, 0
    token = "\\caption{"
    while True:
        k = body.find(token, pos)
        if k == -1:
            out.append(body[pos:])
            return "".join(out), n
        arg, j = read_group(body, k + len(token) - 1)
        core = arg.strip()
        if core.startswith("\\new{") and core.endswith("}"):
            out.append(body[pos:j])
            pos = j
            continue
        sents = sentences(" ".join(strip_latex(core).split()))
        if sents and all(best_ratio(key(s), old_keys, index) >= threshold
                         for s in sents):
            out.append(body[pos:j])
            pos = j
            continue
        out.append(body[pos:k])
        out.append(token + "\\new{" + core + "}}")
        pos = j
        n += 1


def apply_paragraph_marks(path: Path, body: str, old_keys, index,
                          threshold: float) -> int:
    """Wrap each changed prose paragraph in ``\\new{}``.

    Marking whole paragraphs rather than individual sentences keeps the review
    copy readable: a paragraph reworded throughout should not be a checkerboard
    of black and blue. Nesting an existing inline ``\\new{}`` inside the wrap is
    harmless, since the macro just sets a color.
    """
    blocks = body.split("\n\n")
    changed = 0
    for i, block in enumerate(blocks):
        stripped = block.strip()
        if not stripped or stripped.startswith("%"):
            continue

        # A heading shares its block with the paragraph beneath it, and a
        # declaration such as \sloppy must stay outside the wrap: \new{} opens a
        # group, and a declaration inside it would expire before the paragraph
        # breaks. Split both off the front; the heading carries its own \new{}
        # where the title changed.
        head_lines: list[str] = []
        rest = stripped
        while rest.startswith(("\\section", "\\subsection", "\\sloppy",
                               "\\noindent\\textbf{Supporting")):
            line, _, rest = rest.partition("\n")
            head_lines.append(line)
            rest = rest.lstrip("\n")
        if not rest:
            continue

        if any(tok in rest for tok in STRUCTURAL):
            marked, n = mark_prose_around_math(rest, old_keys, index, threshold)
            if n:
                lead = block[:len(block) - len(block.lstrip())]
                head = ("\n".join(head_lines) + "\n") if head_lines else ""
                blocks[i] = f"{lead}{head}{marked}"
                changed += n
            continue
        if rest.startswith("\\new{") and rest.endswith("}"):
            continue                      # already wrapped whole
        sents = sentences(" ".join(strip_latex(rest).split()))
        if not sents:
            continue
        if all(best_ratio(key(s), old_keys, index) >= threshold for s in sents):
            continue                      # carried over unchanged

        lead = block[:len(block) - len(block.lstrip())]
        head = ("\n".join(head_lines) + "\n") if head_lines else ""
        blocks[i] = f"{lead}{head}\\new{{{rest}}}"
        changed += 1

    text = "\n\n".join(blocks)
    text, n_caps = mark_captions(text, old_keys, index, threshold)
    changed += n_caps
    if changed:
        path.write_text(text, encoding="utf-8")
    return changed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--new", default="ACS_Omega_resubmission",
                    help="folder holding manuscript_body.tex")
    ap.add_argument("--old", default="ACS_Omega_submission",
                    help="folder holding the previously submitted main.tex")
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    ap.add_argument("--show-all", action="store_true",
                    help="list every sentence with its similarity score")
    ap.add_argument("--apply", action="store_true",
                    help="wrap every changed prose paragraph in \\new{} in place")
    args = ap.parse_args()

    old_path = ROOT / args.old / "main.tex"
    new_path = ROOT / args.new / "manuscript_body.tex"
    for p in (old_path, new_path):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    old_sents = sentences(strip_latex(old_path.read_text(encoding="utf-8")))
    old_keys = [key(s) for s in old_sents]
    index: dict[str, set[int]] = {}
    for i, k in enumerate(old_keys):
        for w in set(k):
            index.setdefault(w, set()).add(i)

    body = new_path.read_text(encoding="utf-8")

    if args.apply:
        n = apply_paragraph_marks(new_path, body, old_keys, index, args.threshold)
        print(f"wrapped {n} changed prose paragraphs in \\new{{}}")
        body = new_path.read_text(encoding="utf-8")

    spans = marked_spans(body)

    def is_marked(pos: int) -> bool:
        return any(a <= pos <= b for a, b in spans)

    # Walk the new body sentence by sentence, keeping each one's offset so the
    # \new{} state at that point in the source can be recovered.
    prose_chunks: list[tuple[int, str]] = []
    offset, plain = 0, []
    for m in re.finditer(r"[^.;:]+[.;:]", body):
        prose_chunks.append((m.start(), m.group(0)))

    under, over, unchanged, new_ok = [], [], 0, 0
    for pos, chunk in prose_chunks:
        text = " ".join(strip_latex(chunk).split())
        if len(text.split()) < MIN_WORDS:
            continue
        k = key(text)
        ratio = best_ratio(k, old_keys, index)
        marked = is_marked(pos)
        carried = ratio >= args.threshold
        if args.show_all:
            flag = "M" if marked else " "
            print(f"  [{flag}] {ratio:.2f}  {text[:88]}")
        if carried and marked:
            over.append((ratio, text))
        elif carried and not marked:
            unchanged += 1
        elif not carried and marked:
            new_ok += 1
        else:
            under.append((ratio, text))

    print(f"\ncompared against {args.old}/main.tex "
          f"({len(old_sents)} sentences), threshold {args.threshold:.2f}\n")
    print(f"  marked new, and genuinely new     : {new_ok}")
    print(f"  unmarked, and genuinely carried   : {unchanged}")
    print(f"  UNDER-MARKED (new but unmarked)   : {len(under)}")
    print(f"  OVER-MARKED (carried but marked)  : {len(over)}")

    for label, rows in (("UNDER-MARKED", under), ("OVER-MARKED", over)):
        if not rows:
            continue
        print(f"\n{label}:")
        for ratio, text in sorted(rows, reverse=True)[:40]:
            print(f"  {ratio:.2f}  {text[:100]}")

    total = new_ok + unchanged + len(under) + len(over)
    if total:
        agree = 100 * (new_ok + unchanged) / total
        print(f"\nmarking agrees with the comparison on {agree:.1f}% of "
              f"{total} sentences")
    sys.exit(1 if under else 0)


if __name__ == "__main__":
    main()
