"""
build_v6_docx.py -- render the v6 AIIM manuscript into a .docx that inherits the
v4 submission's formatting.

Why not reuse `build_v5_docx.py`: that script *splices* revised sections into
the v4 document at documented placement points, which works when the revision
patches an existing structure. v6 does not — it is a different structure (the
lean AIIM research-article form), so it replaces the body outright. What is
inherited is the formatting package, not the text: Normal = Times New Roman
12 pt, US Letter, 1 in margins, continuous line numbering, Heading1/Heading2,
TableGrid. Those live in styles.xml and in the body's trailing w:sectPr, both
of which are carried over untouched.

The v4 file is never modified; output goes to a new filename.

Two things are deliberately dropped from the rendered document:
  * the front-matter blockquote (draft status), which is scaffolding
  * the "Drafting notes (not manuscript text)" tail, which says so itself

`[PLACEHOLDER: ...]` markers are KEPT and rendered in bold. They are the
author's worklist, and a document that silently omitted them would look
finished when it is not.

Usage:
    python build_v6_docx.py [--out PATH]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_v5_docx import (  # noqa: E402
    W, WNS, blocks_to_elements, md_to_blocks,
)


def _root() -> Path:
    marker = "00_project_control/08_DECISION_LOG.md"
    for cand in (Path(r"D:\Agentic VS Workflow"), Path("/mnt/d/Agentic VS Workflow")):
        if (cand / marker).is_file():
            return cand
    raise SystemExit("project root not found")


ROOT = _root()
SRC = ROOT / "30_manuscript/CURRENT_v5/manuscript_AIIM_v6.md"
TEMPLATE = ROOT / ("30_manuscript/archive_v4/docx/"
                   "KDSS_A2AR_AIIM_Submission_v4_2026-08-11_reviewed.docx")
DEFAULT_OUT = ROOT / "30_manuscript/CURRENT_v5/KDSS_A2AR_AIIM_v6_2026-08-21_draft.docx"


def manuscript_markdown() -> str:
    """The prose only: front-matter blockquote and drafting-notes tail removed."""
    txt = SRC.read_text(encoding="utf-8")
    txt = txt.split("## Drafting notes", 1)[0]
    # Drop the leading '# Manuscript v6 ...' title line and the '>' status block;
    # the document's own title lives under '## Title'.
    lines = [l for l in txt.split("\n") if not l.startswith(">")]
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
    txt = "\n".join(lines)
    # Horizontal rules carry no meaning in the .docx; '---' would otherwise be
    # rendered as a literal paragraph of three hyphens.
    txt = re.sub(r"^---\s*$", "", txt, flags=re.M)
    # '##' is a top-level section in v6 ('## 2. Materials and methods') and
    # '###' a subsection ('### 2.1 ...'). md_to_blocks maps level 1 to Heading1
    # and everything deeper to Heading2, so promote one level.
    txt = re.sub(r"^####\s", "### ", txt, flags=re.M)
    txt = re.sub(r"^###\s", "## ", txt, flags=re.M)
    txt = re.sub(r"^##\s", "# ", txt, flags=re.M)
    return txt.strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--template", default=str(TEMPLATE))
    args = ap.parse_args()
    out = Path(args.out)

    tmp = Path(tempfile.mkdtemp(prefix="v6docx_"))
    try:
        with zipfile.ZipFile(args.template) as z:
            z.extractall(tmp)

        docxml = tmp / "word/document.xml"
        raw = docxml.read_text(encoding="utf-8")
        # Same namespace discipline as build_v5_docx: ElementTree only emits
        # declarations for prefixes it uses, but mc:Ignorable references
        # several that appear nowhere else. Register every prefix the file
        # declares and restore the original root tag verbatim afterwards, or
        # Word rejects the package.
        original_root_tag = re.search(r"<w:document\b[^>]*>", raw).group(0)
        for prefix, uri in re.findall(r'xmlns:([\w.-]+)="([^"]+)"', raw):
            try:
                ET.register_namespace(prefix, uri)
            except ValueError:
                pass

        tree = ET.parse(docxml)
        body = tree.getroot().find(WNS + "body")

        # The trailing w:sectPr carries page size, margins and line numbering.
        # It must survive and must remain last.
        sectpr = body.find(WNS + "sectPr")
        for child in list(body):
            body.remove(child)

        els = blocks_to_elements(md_to_blocks(manuscript_markdown()))
        for e in els:
            body.append(e)
        if sectpr is not None:
            body.append(sectpr)

        tree.write(docxml, encoding="utf-8", xml_declaration=True)
        fixed = docxml.read_text(encoding="utf-8")
        fixed = re.sub(r"<w:document\b[^>]*>", original_root_tag, fixed, count=1)
        docxml.write_text(fixed, encoding="utf-8")

        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists():
            out.unlink()
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
            # [Content_Types].xml must be the first entry in a valid package.
            ct = tmp / "[Content_Types].xml"
            z.write(ct, "[Content_Types].xml")
            for f in sorted(tmp.rglob("*")):
                if f.is_file() and f != ct:
                    z.write(f, str(f.relative_to(tmp)).replace("\\", "/"))

        n_par = len(body.findall(WNS + "p"))
        n_tbl = len(body.findall(WNS + "tbl"))
        n_ph = manuscript_markdown().count("PLACEHOLDER")
        print(f"paragraphs   : {n_par}")
        print(f"tables       : {n_tbl}")
        print(f"placeholders : {n_ph} (kept, in bold, as the author's worklist)")
        print(f"template     : {Path(args.template).name}")
        print(f"written -> {out}")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
