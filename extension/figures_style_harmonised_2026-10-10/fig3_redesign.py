"""ADR-AIIM-368: Figure 3 redesign, applied by make_harmonised.py to the generated fig3_screening.py.

One legend for the whole figure above the panels (panel C reuses panel A's colours), and the y-randomisation panel drawn on a
broken x axis so the permuted null (0.46-0.54) is readable next to the real model (0.911).
"""


def _between(src, a, b):
    i = src.index(a)
    j = src.index(b, i)
    return i, j


FIG = '''    fig = plt.figure(figsize=(7.48, 3.05))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.06, 1.0, 1.0], wspace=0.40,
                          left=0.062, right=0.988, top=0.80, bottom=0.17)

'''

PANEL_B = '''    # ---- B: y-randomisation (broken axis) ---------------------------------
    sub = gs[0, 1].subgridspec(1, 2, width_ratios=[6.0, 1.0], wspace=0.10)
    axB = fig.add_subplot(sub[0, 0]); axB2 = fig.add_subplot(sub[0, 1])
    perm = np.asarray(yr["permuted_aucs"], dtype=float)
    assert perm.min() > 0.44 and perm.max() < 0.56 and 0.90 < yr["real_cv_auc"] < 0.92
    axB.hist(perm, bins=np.arange(0.44, 0.5601, 0.005), color="#b0b8c4", edgecolor="white", linewidth=fs.LW["bar_edge"])
    axB.axvline(0.5, color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH)
    axB.set_xlim(0.44, 0.56); axB.set_xticks([0.45, 0.50, 0.55])
    axB.set_ylim(0, axB.get_ylim()[1] * 1.06)
    axB.set_ylabel("permutations", fontsize=fs.FS["label"])
    axB.tick_params(labelsize=fs.FS["tick"])
    axB2.axvline(yr["real_cv_auc"], color=fs.COL["alert"], lw=fs.LW["main"])
    axB2.set_xlim(0.90, 0.92); axB2.set_xticks([round(yr["real_cv_auc"], 3)]); axB2.set_ylim(axB.get_ylim())
    axB2.spines["left"].set_visible(False); axB2.set_yticks([])
    axB2.tick_params(labelsize=fs.FS["tick"])
    axB2.annotate("real\\nlabels", xy=(yr["real_cv_auc"], axB.get_ylim()[1] * 0.55), xytext=(-5, 0), textcoords="offset points",
                  ha="right", va="center", fontsize=fs.FS["tick"], color=fs.COL["alert"])
    for ax_, x_ in ((axB, 1.0), (axB2, 0.0)):
        ax_.plot([x_ - 0.04, x_ + 0.04], [-0.02, 0.02], transform=ax_.transAxes, color=fs.COL["axis"], lw=fs.LW["axis"], clip_on=False)
    axB.set_xlabel("cross-validated ROC-AUC", fontsize=fs.FS["label"])
    axB.xaxis.set_label_coords(0.58, -0.135)
    fs.panel_label(axB, "B")

'''

LEGEND = '''    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    hA, lA = axA.get_legend_handles_labels()
    hA = hA + [Patch(facecolor=METHODS[0][3], alpha=0.13, linewidth=0), Line2D([], [], color="#999999", lw=fs.LW["ref"], ls=fs.REF_DASH)]
    lA = lA + ["95% bootstrap band, Vina alone", "chance"]
    fig.legend(hA, lA, loc="upper center", ncol=4, frameon=False, fontsize=fs.FS["legend"], bbox_to_anchor=(0.5, 1.0),
               handlelength=1.6, columnspacing=1.5, handletextpad=0.5, labelspacing=0.35)
    out = Path(args.out)
    out.mkdir'''


def apply(s):
    i, j = _between(s, "    # Taller than the panels need", "    # ---- A: ROC")
    s = s[:i] + FIG + s[j:]
    i, j = _between(s, "    axA.legend(", "    fs.panel_label(axA")
    s = s[:i] + s[j:]
    i, j = _between(s, "    # ---- B: y-randomisation", "    # ---- C: early recognition")
    s = s[:i] + PANEL_B + s[j:]
    i, j = _between(s, "    axC.legend(", "    fs.panel_label(axC")
    s = s[:i] + s[j:]
    old = "    out = Path(args.out)\n    out.mkdir"
    assert old in s
    s = s.replace(old, LEGEND, 1)
    return s
