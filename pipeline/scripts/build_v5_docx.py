"""
build_v5_docx.py -- assemble the v5 manuscript draft sections into a .docx that
reuses the v4 submission's own formatting.

Approach: rather than generating a new document from scratch (which would lose
the v4 file's styles, tables, figures, header/footer, line numbering, and every
Methods/Results subsection the v5 rescue did not touch), this script unpacks the
v4 .docx, splices the v5 markdown sections into it at the placement points each
v5 file documents in its own "Placement note", and repacks to a NEW filename.
The v4 file is never modified.

Formatting is inherited, not re-specified: Normal = Times New Roman 12 pt,
US Letter, 1 in margins, continuous line numbering, Heading1/Heading2 built-in
styles, tables in the TableGrid style -- all already defined in the v4 package.

Usage:
    python build_v5_docx.py --unpacked <dir> --out <docx path>
"""
from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from kdss_paths import project_root, resolve, under

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
WNS = "{%s}" % W
ET.register_namespace("w", W)

ROOT = project_root()
MS = resolve("30_manuscript/CURRENT_v5/sections")
V4_TEMPLATE = resolve("30_manuscript/archive_v4/docx/KDSS_A2AR_AIIM_Submission_v4_2026-08-11_reviewed.docx")

TABLE_WIDTH = 9360


def esc(t: str) -> str:
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def md_content(path: Path) -> str:
    """Strip the drafting scaffolding: everything before the first '---' rule
    (status/placement/skill notes) and the 'Notes for the next drafting pass'
    tail. What remains is the manuscript prose itself."""
    txt = path.read_text(encoding="utf-8")
    parts = re.split(r"^---\s*$", txt, flags=re.M)
    # first chunk is front-matter; drop it. Drop any chunk that is the notes tail.
    body = []
    for chunk in parts[1:]:
        if re.search(r"##\s*Notes for the next drafting pass", chunk):
            chunk = re.split(r"##\s*Notes for the next drafting pass", chunk)[0]
            body.append(chunk)
            break
        body.append(chunk)
    return "\n".join(body).strip()


def runs_xml(text: str) -> str:
    """Inline markdown -> w:r runs. Handles **bold**, *italic*, `code`, <sup>.

    <sup> is matched BEFORE the italic alternative on purpose: an author
    block's corresponding-author marker is written <sup>*</sup>, and if the
    italic rule saw that asterisk first it would pair it with the next one in
    the paragraph and silently swallow both the marker and the text between
    them. That is exactly what happened to the first build of this author
    block, where <sup>*</sup> rendered as a stray '*' at the end of the line.
    """
    out = []
    token = re.compile(r"(<sup>.+?</sup>|\*\*.+?\*\*|\*[^*\n]+?\*|`[^`\n]+?`)", re.S)
    pos = 0
    for m in token.finditer(text):
        if m.start() > pos:
            out.append(("", text[pos:m.start()]))
        s = m.group(0)
        if s.startswith("<sup>"):
            out.append(("sup", s[5:-6]))
        elif s.startswith("**"):
            out.append(("b", s[2:-2]))
        elif s.startswith("`"):
            out.append(("c", s[1:-1]))
        else:
            out.append(("i", s[1:-1]))
        pos = m.end()
    if pos < len(text):
        out.append(("", text[pos:]))

    xml = []
    for kind, s in out:
        if not s:
            continue
        s = s.replace("\n", " ")
        if kind not in ("c", "sup"):
            # The v5 markdown uses ASCII '--' for an em dash; the v4 document
            # uses real typographic dashes, so match it.
            s = re.sub(r"\s--\s", "—", s)
            s = s.replace("...", "…")
        rpr = ""
        if kind == "b":
            rpr = "<w:rPr><w:b/></w:rPr>"
        elif kind == "i":
            rpr = "<w:rPr><w:i/></w:rPr>"
        elif kind == "c":
            rpr = '<w:rPr><w:rFonts w:ascii="Courier New" w:hAnsi="Courier New"/><w:sz w:val="20"/></w:rPr>'
        elif kind == "sup":
            rpr = '<w:rPr><w:vertAlign w:val="superscript"/></w:rPr>'
        xml.append(f'<w:r>{rpr}<w:t xml:space="preserve">{esc(s)}</w:t></w:r>')
    return "".join(xml) or '<w:r><w:t xml:space="preserve"></w:t></w:r>'


def para(text: str, style: str | None = None) -> str:
    ppr = f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>' if style else ""
    return f"<w:p>{ppr}{runs_xml(text)}</w:p>"


def table_xml(rows: list[list[str]]) -> str:
    ncol = max(len(r) for r in rows)
    w = TABLE_WIDTH // ncol
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for _ in range(ncol))
    body = []
    for ri, row in enumerate(rows):
        cells = []
        for ci in range(ncol):
            cell = row[ci] if ci < len(row) else ""
            bold_open = "**" if (ri == 0 and not cell.startswith("**")) else ""
            content = f"{bold_open}{cell}{'**' if bold_open else ''}" if cell else ""
            cells.append(
                f'<w:tc><w:tcPr><w:tcW w:type="dxa" w:w="{w}"/></w:tcPr>'
                f"{para(content)}</w:tc>"
            )
        body.append(f"<w:tr>{''.join(cells)}</w:tr>")
    return (
        "<w:tbl><w:tblPr><w:tblStyle w:val=\"TableGrid\"/>"
        '<w:tblW w:type="auto" w:w="0"/>'
        '<w:tblLook w:firstColumn="1" w:firstRow="1" w:lastColumn="0" '
        'w:lastRow="0" w:noHBand="0" w:noVBand="1" w:val="04A0"/></w:tblPr>'
        f"<w:tblGrid>{grid}</w:tblGrid>{''.join(body)}</w:tbl>"
    )


def md_to_blocks(md: str) -> list[str]:
    """Markdown -> list of OOXML block strings (w:p / w:tbl)."""
    blocks: list[str] = []
    lines = md.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
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
            blocks.append(table_xml(rows))
            continue
        # bullet list item -- ListBullet exists in the v4 template's styles.xml
        # and carries the bullet glyph and indent, so the "- " marker is
        # consumed rather than printed.
        bm = re.match(r"^\s*[-*]\s+(.*)$", line)
        if bm:
            buf = [bm.group(1).strip()]
            i += 1
            # a wrapped continuation line is indented and starts no new marker
            while (i < len(lines) and lines[i].strip()
                   and not re.match(r"^\s*[-*]\s+", lines[i])
                   and not lines[i].lstrip().startswith(("|", "#"))
                   and lines[i].startswith(("  ", "\t"))):
                buf.append(lines[i].strip())
                i += 1
            blocks.append(para(" ".join(buf), "ListBullet"))
            continue
        # headings
        m = re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            lvl = len(m.group(1))
            head = m.group(2).strip()
            if re.match(r"^Table\b", head):
                # Table titles are captions, not document headings -- emitting
                # them as Heading1 would put them in the manuscript's outline
                # between real sections.
                blocks.append(para(f"**{head}**"))
            elif head.lower() == "related work":
                blocks.append(para(head, "Heading1"))
            else:
                # '#' = top-level section, '##' and deeper = subsection. The v5
                # files use '##' for numbered subsections (2.x, 3.x, 3.y) and
                # '###' for sub-subsections (4.1-4.5).
                blocks.append(para(head, "Heading1" if lvl == 1 else "Heading2"))
            i += 1
            continue
        # italic caption line starting with *Table caption
        # normal paragraph: gather until blank line
        buf = [line]
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].lstrip().startswith(("|", "#")):
            buf.append(lines[i].rstrip())
            i += 1
        blocks.append(para(" ".join(buf).strip()))
    return blocks


def blocks_to_elements(blocks: list[str]) -> list[ET.Element]:
    els = []
    for b in blocks:
        wrapped = f'<w:root xmlns:w="{W}">{b}</w:root>'
        root = ET.fromstring(wrapped)
        els.extend(list(root))
    return els


def heading_index(kids, style_prefix, text_prefix):
    for i, k in enumerate(kids):
        if k.tag != WNS + "p":
            continue
        pr = k.find(WNS + "pPr")
        if pr is None:
            continue
        ps = pr.find(WNS + "pStyle")
        if ps is None or not ps.get(WNS + "val", "").startswith(style_prefix):
            continue
        t = "".join(e.text or "" for e in k.iter(WNS + "t")).strip()
        if t.lower().startswith(text_prefix.lower()):
            return i
    return None


def next_heading_after(kids, start):
    for i in range(start + 1, len(kids)):
        k = kids[i]
        if k.tag != WNS + "p":
            continue
        pr = k.find(WNS + "pPr")
        if pr is None:
            continue
        ps = pr.find(WNS + "pStyle")
        if ps is not None and ps.get(WNS + "val", "").startswith(("Heading1", "Heading2")):
            return i
    return len(kids)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unpacked", required=True)
    ap.add_argument("--report", default="")
    args = ap.parse_args()

    docxml = Path(args.unpacked) / "word/document.xml"
    # ElementTree only emits namespace declarations for prefixes it actually
    # uses, but the document root declares several (w14, w15, w16*, wp14...)
    # solely so that mc:Ignorable can reference them. Capture the original
    # opening tag verbatim and restore it after writing, or Word/validators
    # reject the file with "Namespace X in Ignorable but not declared".
    raw = docxml.read_text(encoding="utf-8")
    original_root_tag = re.search(r"<w:document\b[^>]*>", raw).group(0)
    # Register every prefix the original declares, so ElementTree round-trips
    # attributes like w14:paraId under their original prefix instead of
    # inventing ns0/ns1/... names that no declaration covers.
    # Scan the WHOLE file, not just the root tag: DrawingML prefixes (a:, pic:,
    # wp:) are declared inline on nested elements, and ET re-emits them at the
    # point of use only if the prefix is registered.
    for prefix, uri in re.findall(r'xmlns:([\w.-]+)="([^"]+)"', raw):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError:
            pass
    tree = ET.parse(docxml)
    body = tree.getroot().find(WNS + "body")
    kids = list(body)

    # ---- load v5 prose ----
    abs_intro = md_content(under("30_manuscript/CURRENT_v5/sections", "01_abstract_and_introduction_DRAFT_v5.md"))
    arch = md_content(under("30_manuscript/CURRENT_v5/sections", "02a_agentic_architecture_and_skills_DRAFT_v5.md"))
    related = md_content(under("30_manuscript/CURRENT_v5/sections", "02b_related_work_DRAFT_v5.md"))
    libexp = md_content(under("30_manuscript/CURRENT_v5/sections", "02c_library_expansion_methods_DRAFT_v5.md"))
    # resolve 02c's open placeholder (value established by ADR-AIIM-085/086)
    libexp = libexp.replace(
        "[AUTHOR ACTION / RESULT PENDING: the top N compounds by\nfast-pass score]",
        "the top 500 compounds by fast-pass score",
    ).replace(
        "[AUTHOR ACTION / RESULT PENDING: the top N compounds by fast-pass score]",
        "the top 500 compounds by fast-pass score",
    )
    sbdd = md_content(under("30_manuscript/CURRENT_v5/sections", "03a_sbdd_rescue_null_result_DRAFT_v5.md"))
    final = md_content(under("30_manuscript/CURRENT_v5/sections", "03b_final_prioritization_results_DRAFT_v5.md"))
    disc = md_content(under("30_manuscript/CURRENT_v5/sections", "04_discussion_and_conclusion_DRAFT_v5.md"))

    log = []

    def replace_section(head_style, head_text, new_md, keep_heading=True, label=""):
        idx = heading_index(kids, head_style, head_text)
        if idx is None:
            log.append(f"SKIP (heading not found): {label or head_text}")
            return
        end = next_heading_after(kids, idx)
        new_els = blocks_to_elements(md_to_blocks(new_md))
        if keep_heading:
            start = idx + 1
        else:
            start = idx
        removed = end - start
        for _ in range(removed):
            body.remove(kids[start])
            kids.pop(start)
        for off, el in enumerate(new_els):
            body.insert(start + off, el)
            kids.insert(start + off, el)
        log.append(f"REPLACED {label or head_text}: -{removed} +{len(new_els)} blocks")

    def insert_after_section(head_style, head_text, new_md, label=""):
        idx = heading_index(kids, head_style, head_text)
        if idx is None:
            log.append(f"SKIP (anchor not found): {label}")
            return
        end = next_heading_after(kids, idx)
        new_els = blocks_to_elements(md_to_blocks(new_md))
        for off, el in enumerate(new_els):
            body.insert(end + off, el)
            kids.insert(end + off, el)
        log.append(f"INSERTED {label} after '{head_text}': +{len(new_els)} blocks")

    # Abstract + Introduction come from one file; split on the '## 1. Introduction' heading
    m = re.split(r"^##\s*1\.\s*Introduction\s*$", abs_intro, flags=re.M)
    abstract_md = re.sub(r"^##\s*Abstract\s*$", "", m[0], flags=re.M).strip()
    abstract_md = re.sub(r"^\*\(.*?\)\*\s*$", "", abstract_md, flags=re.M | re.S).strip()
    intro_md = m[1].strip() if len(m) > 1 else ""

    replace_section("Heading1", "Abstract", abstract_md, label="Abstract")
    replace_section("Heading1", "1. Introduction", intro_md, label="1. Introduction")

    # 02a restates the '2.1 ...' heading its content belongs under; the docx
    # already has that heading and replace_section keeps it, so drop the
    # duplicate rather than emitting the same heading twice in a row.
    arch_md = re.sub(r"^#{1,3}\s*2\.1\b.*$", "", arch, flags=re.M)
    arch_md = re.sub(r"^#\s+.*$", "", arch_md, flags=re.M).strip()
    replace_section("Heading2", "2.1", arch_md, label="2.1 architecture+skills")

    insert_after_section("Heading1", "1. Introduction", related, label="Related Work")
    insert_after_section("Heading2", "2.7", libexp, label="2.x library expansion")
    insert_after_section("Heading2", "3.2", sbdd, label="3.x SBDD rescue")
    insert_after_section("Heading2", "3.9", final, label="3.y final prioritization")

    # Discussion / Conclusion / Limitations from 04_v5
    d_disc = re.split(r"^##\s*5\.\s*Conclusion\s*$", disc, flags=re.M)
    disc_body = re.sub(r"^##\s*4\.\s*Discussion\s*$", "", d_disc[0], flags=re.M).strip()
    rest = d_disc[1] if len(d_disc) > 1 else ""
    d_lim = re.split(r"^##\s*6\.\s*Limitations\s*$", rest, flags=re.M)
    concl_body = d_lim[0].strip()
    lim_body = d_lim[1].strip() if len(d_lim) > 1 else ""

    replace_section("Heading1", "4. Discussion", disc_body, label="4. Discussion")
    replace_section("Heading1", "5. Limitations", lim_body, label="5. Limitations")
    replace_section("Heading1", "6. Conclusion", concl_body, label="6. Conclusion")

    tree.write(docxml, encoding="UTF-8", xml_declaration=True)
    written = docxml.read_text(encoding="utf-8")
    # ElementTree hoists inline (locally scoped) namespace declarations such as
    # DrawingML's xmlns:a on <a:graphic>, but does not always re-emit them, which
    # leaves prefixes used-but-undeclared. Declaring every prefix the original
    # file used on the restored root element puts them all in scope document-wide;
    # a redundant declaration deeper in the tree is legal XML.
    root_tag = original_root_tag
    for prefix, uri in dict(re.findall(r'xmlns:([\w.-]+)="([^"]+)"', raw)).items():
        if f'xmlns:{prefix}=' not in root_tag:
            root_tag = root_tag[:-1] + f' xmlns:{prefix}="{uri}">'
    written = re.sub(r"<(?:\w+:)?document\b[^>]*>", root_tag, written, count=1)
    docxml.write_text(written, encoding="utf-8")

    print("\n".join(log))
    if args.report:
        Path(args.report).write_text("\n".join(log), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
