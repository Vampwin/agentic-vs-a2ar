"""render_svg_to_png.py -- rasterise a figure SVG at the project's 600 dpi.

The figure generators emit SVG only, and the manuscript embeds PNG. Nothing in
this repository recorded how the conversion was done, and no converter is
installed: no cairosvg, no rsvg-convert, no Inkscape, no ImageMagick. Installing
one needs approval under this project's rules, so this uses a browser that is
already on the machine.

**The scale factor is not a preference.** The existing PNGs are 8.333x their
SVG viewBox — 4259 x 1717 against a 511 x 206 canvas — which is 600 dpi, since
one point is 1/72 inch and 600/72 = 8.333. Matching it is what keeps a rebuilt
figure interchangeable with the ones beside it in the same manuscript.

The SVG is wrapped in a minimal HTML page that pins the element to its viewBox
size in CSS pixels, because the figures declare their width in millimetres and a
browser lays that out at 96 dpi rather than 72, which would silently rescale
every rebuilt figure by 96/72 against the ones it sits next to.

Usage:
    python render_svg_to_png.py <in.svg> <out.png> [--dpi 600]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BROWSERS = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]
POINTS_PER_INCH = 72.0


def viewbox(svg_text: str) -> tuple[float, float]:
    m = re.search(r'viewBox="\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+'
                  r'([\d.eE+-]+)\s+([\d.eE+-]+)\s*"', svg_text)
    if not m:
        raise SystemExit("no viewBox in the SVG; cannot size the render")
    return float(m.group(3)), float(m.group(4))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("svg")
    ap.add_argument("png")
    ap.add_argument("--dpi", type=float, default=600.0)
    args = ap.parse_args()

    src = Path(args.svg)
    dst = Path(args.png)
    svg_text = src.read_text(encoding="utf-8")
    w_pt, h_pt = viewbox(svg_text)
    scale = args.dpi / POINTS_PER_INCH

    browser = next((b for b in BROWSERS if Path(b).is_file()), None)
    if browser is None:
        raise SystemExit("no Chrome or Edge found; cannot rasterise without "
                         "installing a converter, which needs approval")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # Pin the SVG to viewBox units in CSS pixels. Without this the browser
        # honours the millimetre width at 96 dpi and the render comes out 4/3
        # too large relative to every other figure in the manuscript.
        page = td / "page.html"
        page.write_text(
            "<!doctype html><meta charset='utf-8'>"
            "<style>html,body{margin:0;padding:0;background:#fff}"
            f"svg{{width:{w_pt}px;height:{h_pt}px;display:block}}</style>"
            + svg_text, encoding="utf-8")
        shot = td / "shot.png"
        cmd = [browser, "--headless", "--disable-gpu", "--hide-scrollbars",
               "--default-background-color=FFFFFFFF",
               f"--force-device-scale-factor={scale}",
               f"--window-size={int(round(w_pt))},{int(round(h_pt))}",
               f"--screenshot={shot}", page.as_uri()]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if not shot.is_file():
            raise SystemExit("browser produced no screenshot:\n"
                             + (r.stderr or r.stdout)[-600:])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(shot, dst)

    try:
        from PIL import Image
        with Image.open(dst) as im:
            got = im.size
        want = (int(round(w_pt * scale)), int(round(h_pt * scale)))
        print(f"{src.name} -> {dst.name}  {got[0]}x{got[1]}  "
              f"(expected {want[0]}x{want[1]} at {args.dpi:.0f} dpi)")
        if abs(got[0] - want[0]) > 2 or abs(got[1] - want[1]) > 2:
            print("  WARNING: size differs from the expected render by more "
                  "than rounding; do not ship this beside the other figures")
            return 1
    except ImportError:
        print(f"{src.name} -> {dst.name} (Pillow absent; size unverified)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
