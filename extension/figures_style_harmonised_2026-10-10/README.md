# figures_style_harmonised_2026-10-10 — ADR-AIIM-367

`figstyle.py` is the one style of the manuscript's figures (Arial; nothing below 7.0 pt; panel labels bold uppercase without
brackets, 4 pt above the top-left corner; axes 0.7, reference lines 0.8, threshold lines 1.2 dashed, main data 1.6; one red, one blue,
fixed compound colours; every figure saved exactly 190 mm wide at 600 dpi plus SVG). `make_harmonised.py` writes `fig3_screening.py`,
`fig4_ad.py` and `fig5_permutation.py` from the originals in `10_new_system/scripts` (frozen, not modified); Figure 6, Figure 7 and
Supplementary Figure S1 are made by the scripts of `figures_docking_md_2026-10-09/analysis` (their earlier versions kept as
`*.before_style_harmonisation`). Figures 1 and 2 (Arial SVG at 190 mm, 7 pt floor) and the graphical abstract were not redrawn.
Outputs of Figures 3 to 5 are in `out/`. The scripts use absolute `D:\Agentic VS Workflow` paths, as the other extension scripts do.


## ADR-AIIM-368 addendum

`fig3_redesign.py` is applied by `make_harmonised.py` to the generated `fig3_screening.py`: one legend above the three panels and a broken x axis for the y-randomisation panel. Re-run `python make_harmonised.py` and then `python fig3_screening.py --out out` to rebuild Figure 3.
