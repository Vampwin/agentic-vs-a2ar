"""Accepted-changes text of the v8 working documents, in the form G8 reads.
Manuscript -> markdown (headings, paragraphs, table rows). Supplement -> a .docx whose document.xml has
every tracked deletion removed, because G8's docx reader keeps w:delText."""
import re, sys, zipfile, shutil
src_main, src_supp, out = sys.argv[1], sys.argv[2], sys.argv[3]
def accepted_xml(path):
    t = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    t = re.sub(r"<w:del [^>]*/>", "", t)
    return re.sub(r"<w:del [^>]*>.*?</w:del>", "", t, flags=re.S)
# --- supplement docx
zin = zipfile.ZipFile(src_supp)
zout = zipfile.ZipFile(out + "/Supplementary_Information.docx", "w", zipfile.ZIP_DEFLATED)
for it in zin.infolist():
    data = zin.read(it.filename)
    if it.filename == "word/document.xml":
        data = accepted_xml(src_supp).encode("utf-8")
    zout.writestr(it, data)
zout.close()
# --- manuscript md
x = accepted_xml(src_main)
body = re.search(r"<w:body>(.*)</w:body>", x, re.S).group(1)
lines = []
def ptext(p): return "".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", p)).replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">").replace("&quot;", '"').replace("&apos;", "'")
for m in re.finditer(r"<w:tbl>.*?</w:tbl>|<w:p[ >].*?</w:p>", body, re.S):
    blk = m.group(0)
    if blk.startswith("<w:tbl>"):
        for tr in re.findall(r"<w:tr[ >].*?</w:tr>", blk, re.S):
            cells = [ptext(tc).strip() for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, re.S)]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
        continue
    t = ptext(blk).strip()
    if not t: continue
    st = re.search(r'<w:pStyle w:val="([^"]+)"', blk)
    st = st.group(1) if st else ""
    if st == "Title": lines.append("# " + t)
    elif st == "Heading1": lines.append("## " + t)
    elif st == "Heading2": lines.append("### " + t)
    elif st == "Heading3": lines.append("#### " + t)
    else: lines.append(t)
    lines.append("")
open(out + "/manuscript_AIIM_v8_accepted.md", "w", encoding="utf-8", newline="\n").write("\n".join(lines))
print("md lines", len(lines))
