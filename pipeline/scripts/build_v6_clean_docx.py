"""
build_v6_clean_docx.py -- render the clean, submission-facing manuscript
(manuscript_AIIM_v6_CLEAN.md) into a .docx with tables and figures embedded
inline at the position they appear in the markdown.

Distinct from build_v6_docx.py, which reads the WORKING markdown (with its
internal notes and end-of-document Tables/Figures sections) and preserves the
v4 template's raw OOXML for formatting fidelity. This script reads the CLEAN
markdown instead, where notes are already stripped and tables/figures are
already inline; because of that it can use python-docx directly rather than
splicing XML, at the cost of the two documents having independently defined
(though matching) page/font settings rather than sharing one styles.xml.

Usage:
    python build_v6_clean_docx.py [--out PATH]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from docx.oxml.ns import qn



def _assert_figures_current(figdir):
    """Refuse to build if an embedded figure is older than what its script made.

    Nothing connects `30_manuscript/figures/`, where the figure scripts write,
    to `CURRENT_v5/figures_png/`, which this builder embeds: the copy is manual.
    Figure 6 was embedded five days stale on 2026-09-03, after gaining a panel,
    and the graphical abstract had been stale for five days once before
    (ADR-AIIM-160). G8 audits text, not images, so nothing else catches this.
    """
    src = figdir.parent / "figures"
    stale = []
    for png in sorted((figdir / "figures_png").glob("*.png")):
        s = src / png.name
        if s.is_file() and s.stat().st_mtime > png.stat().st_mtime + 1:
            stale.append(png.name)
    if stale:
        raise SystemExit(
            "REFUSED: newer versions of these figures exist in "
            "30_manuscript/figures/ than the copies this builder embeds:"
            + "".join(chr(10) + "  " + n for n in stale)
            + chr(10) + "Copy them into CURRENT_v5/figures_png/ and rebuild.")

def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SRC = ROOT / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6_CLEAN.md"
FIGDIR = ROOT / "30_manuscript/CURRENT_v5"
# No date in the name, deliberately. A date in a filename is a claim that
# nothing verifies: this file was called `..._2026-08-24.docx` while its
# content was rebuilt daily for twelve days, and the name was wrong for
# eleven of them. The build date is in the file's own mtime, the version
# is in the decision log, and both of those are checkable. An undated name
# also means the artifact has one path, so a stale copy cannot sit beside
# the current one wearing a different date.
DEFAULT_OUT = ROOT / "30_manuscript/CURRENT_v5/KDSS_A2AR_AIIM_v6_CLEAN.docx"

# Full-page figures (all four here) are wider than the 6.5in text column, so
# they are inserted at this width regardless of column and get a caption
# beneath rather than beside them -- matching how the source SVGs were laid
# out (190/193/180 mm, i.e. full journal-page width, not column width).
IMG_WIDTH_IN = 6.5


def add_runs(paragraph, text: str) -> None:
    """Inline markdown -> runs. Handles **bold**, *italic*, `code`, <sup>."""
    token = re.compile(r"(<sup>.+?</sup>|\*\*.+?\*\*|\*[^*\n]+?\*|`[^`\n]+?`)", re.S)
    pos = 0
    for m in token.finditer(text):
        if m.start() > pos:
            _plain_run(paragraph, text[pos:m.start()])
        s = m.group(0)
        if s.startswith("<sup>"):
            r = paragraph.add_run(s[5:-6])
            r.font.superscript = True
        elif s.startswith("**"):
            r = paragraph.add_run(s[2:-2])
            r.bold = True
        elif s.startswith("`"):
            r = paragraph.add_run(s[1:-1])
            r.font.name = "Courier New"
            r.font.size = Pt(10)
        else:
            r = paragraph.add_run(s[1:-1])
            r.italic = True
        pos = m.end()
    if pos < len(text):
        _plain_run(paragraph, text[pos:])


def _plain_run(paragraph, s: str) -> None:
    if not s:
        return
    s = re.sub(r"\s--\s", "\u2014", s)
    paragraph.add_run(s)


def add_table(doc: Document, rows: list[list[str]]) -> None:
    ncol = max(len(r) for r in rows)
    tbl = doc.add_table(rows=len(rows), cols=ncol)
    tbl.style = "Table Grid"
    for ri, row in enumerate(rows):
        for ci in range(ncol):
            cell_text = row[ci] if ci < len(row) else ""
            p = tbl.cell(ri, ci).paragraphs[0]
            if ri == 0 and not cell_text.startswith("**"):
                cell_text = f"**{cell_text}**"
            add_runs(p, cell_text)
            for run in p.runs:
                run.font.size = Pt(9)


def main() -> int:
    _assert_figures_current(FIGDIR)
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    md = SRC.read_text(encoding="utf-8")
    lines = md.split("\n")

    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    for m in ("top", "bottom", "left", "right"):
        setattr(section, f"{m}_margin", Inches(1))

    n_par = n_tbl = n_img = 0
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue

        # horizontal rule -> skip (section separators, not content)
        if line.strip() == "---":
            i += 1
            continue

        # image
        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", line.strip())
        if m:
            img_path = FIGDIR / m.group(2)
            if not img_path.is_file():
                raise SystemExit(f"missing image referenced in markdown: {img_path}")
            doc.add_picture(str(img_path), width=Inches(IMG_WIDTH_IN))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            n_img += 1
            i += 1
            continue

        # table
        if line.lstrip().startswith("|") and i + 1 < len(lines) and re.match(
            r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]
        ):
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                raw = lines[i].strip().strip("|")
                if not re.match(r"^[\s:|-]+$", raw):
                    rows.append([c.strip() for c in raw.split("|")])
                i += 1
            add_table(doc, rows)
            n_tbl += 1
            continue

        # bullet list item
        bm = re.match(r"^\s*[-*]\s+(.*)$", line)
        if bm:
            buf = [bm.group(1).strip()]
            i += 1
            while (i < len(lines) and lines[i].strip()
                   and not re.match(r"^\s*[-*]\s+", lines[i])
                   and not lines[i].lstrip().startswith(("|", "#", "!["))
                   and lines[i].startswith(("  ", "\t"))):
                buf.append(lines[i].strip())
                i += 1
            bp = doc.add_paragraph(style="List Bullet")
            add_runs(bp, " ".join(buf))
            n_par += 1
            continue
        # heading
        hm = re.match(r"^(#{1,4})\s+(.*)$", line)
        if hm:
            lvl = len(hm.group(1))
            head = hm.group(2).strip()
            style = "Heading 1" if lvl == 2 else "Heading 2"
            p = doc.add_paragraph(style=style)
            add_runs(p, head)
            i += 1
            continue

        # normal paragraph: gather until blank line, table start, heading, or image
        buf = [line]
        i += 1
        while (i < len(lines) and lines[i].strip()
               and not lines[i].lstrip().startswith(("|", "#", "!["))
               and lines[i].strip() != "---"):
            buf.append(lines[i].rstrip())
            i += 1
        text = " ".join(buf).strip()
        p = doc.add_paragraph()
        add_runs(p, text)
        n_par += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(out))
    print(f"paragraphs   : {n_par}")
    print(f"tables       : {n_tbl}")
    print(f"images       : {n_img}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
